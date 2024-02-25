u"""
Export a Salome Mesh to OpenFOAM.

It handles all types of cells. Use
salomeToOpenFOAM.exportToFoam(Mesh_1)
to export. Optionally an output dir can be given as argument.

It's also possible to select a mesh in the object browser and
run the script via file->load script (ctrl-T).

Groups of volumes will be treated as cellZones. If they are
present they will be put in the file cellZones. In order to convert
to regions use the OpenFOAM tool
splitMeshRegions - cellZones

No sorting of faces is done so you'll have to run
renumberMesh -overwrite
In order to use the mesh.
"""
#Copyright 2013
#Author Nicolas Edh,
#Nicolas.Edh@gmail.com,
#or user "nsf" at cfd-online.com
#
#License
#
#    This script is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    salomeToOpenFOAM  is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with hexBlocker.  If not, see <http://www.gnu.org/licenses/>.
#
#    The license is included in the file LICENSE.
#

import sys
import salome
import SMESH
from salome.smesh import smeshBuilder
import os,time
# import salome_pluginsmanager
# from platform import system
try:
    from PyQt4 import QtGui,QtCore
    from PyQt4.QtGui import *
    from PyQt4.QtCore import *
except:
    from PyQt5.QtWidgets import QWidget, QMessageBox
    from PyQt5 import QtCore, QtGui
    import PyQt5.QtCore as QtCore
    from PyQt5.QtWidgets import *
    from PyQt5.QtCore import Qt

#different levels of verbosities, 0 all quiet,
#higher values means more information

debug=1

globalMeshName = ""

verify=True
"""verify face order, migt take longer time"""

#Note: to skip renumberMesh just sort owner
#while moving positions also move neighbour,faces, and bcfaces
#will probably have to first sort the internal faces then bc-faces within each bc

#obj=theStudy.FindObjectByName('name').GetObject()

class MeshBuffer(object):
    """
    Limits the calls to Salome by buffering the face and key details of volumes to speed up exporting
    """
    
    def __init__(self, mesh, v):
        i = 0
        faces = list()
        keys = list()
        fnodes = mesh.GetElemFaceNodes(v, i)

        while fnodes:                           #While not empty list
            faces.append(fnodes)                #Face list
            keys.append(tuple(sorted(fnodes)))  #Buffer key
            i += 1
            fnodes=mesh.GetElemFaceNodes(v, i)
        
        self.v = v         #The volume
        self.faces = faces #The sorted face list
        self.keys = keys
        self.fL = i        #The number of faces
    
    @staticmethod
    def Key(fnodes):
        """Takes the nodes and compresses them into a hashable key"""
        return tuple(sorted(fnodes))
    
    @staticmethod
    def ReverseKey(fnodes):
        """Takes the nodes and compresses them into a hashable key reversed for baffles"""
        if type(fnodes) is tuple:
            return tuple(reversed(fnodes))
        else:
            return tuple(sorted(fnodes, reverse=True)) 


def exportToFoam(mesh, dirname='polyMesh'):
    """
    Export a mesh to OpenFOAM.
    
    args: 
        +    mesh: The mesh
        + dirname: The mesh directory to write to
    
    The algorithm works as follows:
    [1] Loop through the boundaries and collect all faces in each group.
        Faces that don't have a group will be added to the group defaultPatches.
    
    [2] Loop through all cells (volumes) and each face in the cell. 
        If the face has been visited before it is added to the neighbour list
        If not, then check it is a boundary face. 
            If it is, add the cell to the end of owner.
            If not a boundary face and has not yet been visited add it to the list of internal faces. 
    
    To check if a face has been visited a dictionary is used. 
    The key is the sorted list of face nodes converted to a string.
    The value is the face id. Eg: facesSorted[key] = value
    """
    starttime=time.time()
    #try to open files
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    try:
        filePoints = open(dirname + '/points', 'w')
        fileFaces = open(dirname + '/faces', 'w')
        fileOwner = open(dirname + '/owner', 'w')
        fileNeighbour = open(dirname + '/neighbour', 'w')
        fileBoundary = open(dirname + '/boundary', 'w')
    except Exception:
        print('could not open files aborting')
        return

    #Get salome properties
    smesh = smeshBuilder.New()

    __debugPrint__('Number of nodes: %d\n' %(mesh.NbNodes()))
    volumes=mesh.GetElementsByType(SMESH.VOLUME)
    __debugPrint__('Number of cells: %d\n' %len(volumes))
    __debugPrint__('Counting number of faces:\n')

    #Filter faces
    filter = smesh.GetFilter(SMESH.EDGE, SMESH.FT_FreeFaces)
    extFaces = set(mesh.GetIdsFromFilter(filter))
    nrBCfaces = len(extFaces)
    nrExtFaces = len(extFaces)
    buffers=list()
    nrFaces = 0

    for v in volumes:
        b = MeshBuffer(mesh, v) 
        nrFaces += b.fL
        buffers.append(b)

    #all internal faces will be counted twice, external faces once
    #so:
    nrFaces = int((nrFaces + nrExtFaces) / 2)
    nrIntFaces = int(nrFaces - nrBCfaces)
    __debugPrint__('total number of faces: %d, internal: %d, external %d\n'  \
        %(nrFaces, nrIntFaces, nrExtFaces))

    __debugPrint__('Converting mesh to OpenFOAM\n')
    faces = [] #list of internal face nodes ((1 2 3 4 ... ))
    facesSorted = dict() #each list of nodes is sorted.
    bcFaces = [] #list of bc faces
    bcFacesSorted = dict()
    owner = [] #owner file, (of face id, volume id)
    neighbour = [] #neighbour file (of face id, volume id) only internal faces

    #Loop over all salome boundary elemets (faces) 
    # and store them inte the list bcFaces
    grpStartFace = [] # list of face ids where the BCs starts
    grpNrFaces = [] # list of number faces in each BC
    grpNames = [] #list of the group name.
    ofbcfid = 0   # bc face id in openfoam
    nrExtFacesInGroups = 0

    for gr in mesh.GetGroups():
        if gr.GetType() == SMESH.FACE:
            grpNames.append(gr.GetName())
            __debugPrint__('found group \"%s\" of type %s, %d\n' \
                           %(gr.GetName(), gr.GetType(), len(gr.GetIDs())), 2)
            grIds = gr.GetIDs()
            nr = len(grIds)
            if nr > 0:
                grpStartFace.append(nrIntFaces+ofbcfid)
                grpNrFaces.append(nr)

            #loop over faces in group
            for sfid in grIds:
                fnodes = mesh.GetElemNodes(sfid)
                key = MeshBuffer.Key(fnodes)
                if not key in bcFacesSorted:
                    bcFaces.append(fnodes)
                    bcFacesSorted[key] = ofbcfid
                    ofbcfid += 1
                else:
                    raise Exception(\
                        'Error the face, elemId %d, %s belongs to two ' % (sfid, fnodes)  +\
                            'or more groups. One is : %s'  % (gr.GetName()))

            #if the group is a baffle then the faces should be added twice
            if __isGroupBaffle__(mesh, gr, extFaces): #, grIds):
                nrBCfaces += nr
                nrFaces += nr
                nrIntFaces -= nr
                #since nrIntFaces is reduced all previously grpStartFaces are 
                #out of sync
                grpStartFace = [x - nr for x in grpStartFace]
                grpNrFaces[-1] = nr*2
                for sfid in gr.GetIDs():
                    fnodes = mesh.GetElemNodes(sfid)
                    key = MeshBuffer.ReverseKey(fnodes)
                    bcFaces.append(fnodes)
                    bcFacesSorted[key] = ofbcfid
                    ofbcfid += 1
            else:
                nrExtFacesInGroups += nr

    __debugPrint__('total number of faces: %d, internal: %d, external %d\n'  \
        %(nrFaces, nrIntFaces, nrExtFaces), 2)
    #Do the defined groups cover all BC-faces?
    if nrExtFacesInGroups < nrExtFaces:
        __debugPrint__('Warning, some elements don\'t have a group (BC). ' +\
                       'Adding to a new group called defaultPatches\n', 1)
        grpStartFace.append(nrIntFaces + ofbcfid)
        grpNrFaces.append(nrExtFaces - nrExtFacesInGroups)
        salomeIDs = []
        for face in extFaces:
            fnodes = mesh.GetElemNodes(face)
            key = MeshBuffer.Key(fnodes)
            try:
                bcFacesSorted[key]
            except KeyError:
                #if not in dict then add to default patches
                bcFaces.append(fnodes)
                bcFacesSorted[key] = ofbcfid
                salomeIDs.append(face)
                ofbcfid += 1
        newGrpName = 'defaultPatches'
        nri = 1
        while newGrpName in grpNames:
            newGrpName = "defaultPatches_%d" % nri
            nri += 1
        grpNames.append(newGrpName)
        #function might have different name
        try:
            defGroup = mesh.CreateGroup(SMESH.FACE, newGrpName)
        except AttributeError:
            defGroup = mesh.CreateEmptyGroup(SMESH.FACE, newGrpName)

        defGroup.Add(salomeIDs)
        smesh.SetName(defGroup, newGrpName)

        if salome.sg.hasDesktop():
            if sys.version_info.major < 3:
                salome.sg.updateObjBrowser(1)
            else:
                salome.sg.updateObjBrowser()

    #initialise the list faces vs owner/neighbour cells
    owner = [-1] * nrFaces
    neighbour = [-1] * nrIntFaces
    __debugPrint__('Finished processing boundary faces\n')
    __debugPrint__('bcFaces: %d\n' % len(bcFaces), 2)
    __debugPrint__(str(bcFaces) + '\n', 3)
    __debugPrint__('bcFacesSorted: %d\n' % len(bcFacesSorted), 2)
    __debugPrint__(str(bcFacesSorted) + '\n', 3)
    __debugPrint__('owner: %d\n' % len(owner), 2)
    __debugPrint__(str(owner) + '\n', 3)
    __debugPrint__('neighbour: %d\n' % len(neighbour), 2)
    __debugPrint__(str(neighbour) + '\n', 3)


    offid = 0
    ofvid = 0 #volume id in openfoam
    for b in buffers:
        
        nodes = mesh.GetElemNodes(b.v)
        if debug > 2: #Salome call only if verbose
            nodes = mesh.GetElemNodes(b.v)
            __debugPrint__('volume id: %d, num nodes %d, nodes:%s \n' %(b.v, len(nodes), nodes), 3)
        
        fi = 0 #Face index
        while fi < b.fL:
            fnodes = b.faces[fi]
            key = b.keys[fi]
            #Check if the node is already in list
            try:
                fidinof = facesSorted[key]
                #if faceSorted didn't throw an exception then the face is 
                #already in the dict. Its an internal face and should be added 
                # to the neighbour list
                #print('fidinof %d' % fidinof)
                #--ES: I have seen the following give a: list assignment index out of range
                neighbour[fidinof] = ofvid
                __debugPrint__('\tan owner already exist for %d, %s, cell %d\n' %(fi, fnodes, ofvid), 3)
            except KeyError:
                #the face is not in the list of internal faces
                #it might a new face or a BCface.
                try:
                    bcind = bcFacesSorted[key]
                    #if no exception was trown then it's a bc face
                    __debugPrint__('\t found bc face: %d, %s, cell %d\n' %(bcind, fnodes, ofvid), 3)
                    #if the face belongs to a baffle then it exits twice in owner
                    #check dont overwrite owner
                    if owner[nrIntFaces + bcind] == -1:
                        owner[nrIntFaces + bcind] = ofvid
                        bcFaces[bcind] = fnodes
                    else:
                        #build functions that looks for baffles in bclist. with bcind
                        key = MeshBuffer.ReverseKey(fnodes)
                        bcind = bcFacesSorted[key]
                        #make sure the faces has the correct orientation
                        bcFaces[bcind] = fnodes
                        owner[nrIntFaces + bcind] = ofvid
                except KeyError:
                    #the face is not in bc list either so it's a new internal face
                    __debugPrint__('\t a new face was found, %d, %s, cell %d\n' %(fi, fnodes, ofvid), 3)
                    if verify:
                        if not __verifyFaceOrder__(mesh, nodes, fnodes):
                            __debugPrint__('\t face has bad order, reversing order\n', 3)
                            fnodes.reverse()
                    faces.append(fnodes)
                    key = b.keys[fi]
                    facesSorted[key] = offid
                    owner[offid] = ofvid
                    offid += 1
                    if nrFaces > 50 and offid % (nrFaces/50) == 0:
                        if offid % ((nrFaces/50)*10) == 0:
                            __debugPrint__(':', 1)
                        else:
                            __debugPrint__('.', 1)
            fi += 1
            
        ofvid += 1
        # end for v in volumes

    nrCells = ofvid
    __debugPrint__('Finished processing volumes.\n')
    __debugPrint__('faces: %d\n' % len(faces), 2)
    __debugPrint__(str(faces) + '\n', 3)
    __debugPrint__('facesSorted: %d\n' % len(facesSorted), 2)
    __debugPrint__(str(facesSorted) + '\n', 3)
    __debugPrint__('owner: %d\n' %(len(owner)), 2)
    __debugPrint__(str(owner) + '\n', 3)
    __debugPrint__('neighbour: %d\n' %(len(neighbour)), 2)
    __debugPrint__(str(neighbour) + '\n', 3)


    #Convert to "upper triangular order"
    #owner is sorted, for each cell sort faces it's neighbour faces
    # i.e. change 
    # owner   neighbour      owner   neighbour
    #     0          15         0           3
    #     0           3  to     0          15
    #     0          17         0          17
    #     1           5         1           5
    # any changes made to neighbour are repeated to faces.
    __debugPrint__('Sorting faces in upper triangular order\n', 1)
    ownedfaces = 1
    quickrange = range if sys.version_info.major > 2 else xrange

    for faceId in quickrange(0, nrIntFaces):
        cellId = owner[faceId]
        nextCellId = owner[faceId + 1] #np since len(owner) > nrIntFaces
        if cellId == nextCellId:
            ownedfaces += 1
            continue
        
        if ownedfaces > 1:
            sId = faceId - ownedfaces + 1 #start ID
            eId = faceId #end ID
            inds = range(sId, eId + 1)

            if sys.version_info.major > 2:
                sorted(inds, key=neighbour.__getitem__)
            else:
                inds.sort(key = neighbour.__getitem__)

            neighbour[sId:eId + 1] = map(neighbour.__getitem__, inds)
            faces[sId:eId + 1] = map(faces.__getitem__, inds)

        ownedfaces = 1
    converttime = time.time() - starttime

    #WRITE points to file
    __debugPrint__('Writing the file points\n')
    __writeHeader__(filePoints, 'points')
    points = mesh.GetElementsByType(SMESH.NODE)
    nrPoints = len(points)
    filePoints.write('\n%d\n(\n' % nrPoints)
    for n, ni in enumerate(points):
        pos = mesh.GetNodeXYZ(ni)
        filePoints.write('\t(%.10g %.10g %.10g)\n' % (pos[0], pos[1], pos[2]))
    filePoints.write(')\n')
    filePoints.flush()
    filePoints.close()

    #WRITE faces to file
    __debugPrint__('Writing the file faces\n')
    __writeHeader__(fileFaces, 'faces')
    fileFaces.write('\n%d\n(\n' % nrFaces)
    for node in faces:
        fileFaces.write('\t%d(' % len(node))
        for p in node:
            #salome starts to count from one, OpenFOAM from zero
            fileFaces.write('%d ' % (p - 1))
        fileFaces.write(')\n')
    #internal nodes are done output bcnodes
    for node in bcFaces:
        fileFaces.write('\t%d(' % len(node))
        for p in node:
            #salome starts to count from one, OpenFOAM from zero
            fileFaces.write('%d ' % (p - 1))
        fileFaces.write(')\n')
    fileFaces.write(')\n')
    fileFaces.flush()
    fileFaces.close()

    #WRITE owner to file
    __debugPrint__('Writing the file owner\n')
    __writeHeader__(fileOwner, 'owner', nrPoints, nrCells, nrFaces, nrIntFaces)
    fileOwner.write('\n%d\n(\n' % len(owner))
    for cell in owner:
        fileOwner.write(' %d \n' % cell)
    fileOwner.write(')\n')
    fileOwner.flush()
    fileOwner.close()

    #WRITE neighbour
    __debugPrint__('Writing the file neighbour\n')
    __writeHeader__(fileNeighbour, 'neighbour', nrPoints, nrCells, nrFaces, nrIntFaces)
    fileNeighbour.write('\n%d\n(\n' %(len(neighbour)))
    for cell in neighbour:
        fileNeighbour.write(' %d\n' %(cell))
    fileNeighbour.write(')\n')
    fileNeighbour.flush()
    fileNeighbour.close()

    #WRITE boundary file
    __debugPrint__('Writing the file boundary\n')
    __writeHeader__(fileBoundary, 'boundary')
    fileBoundary.write('%d\n(\n' %len(grpStartFace))
    for ind, gname in enumerate(grpNames):
        fileBoundary.write('\t%s\n\t{\n' %gname)
        fileBoundary.write('\t\ttype\t\t')
        fileBoundary.write(str(bound[ind].currentText()) + ";\n")
        fileBoundary.write('\t\tnFaces\t\t%d;\n' %grpNrFaces[ind])
        fileBoundary.write('\t\tstartFace\t%d;\n' %grpStartFace[ind])
        fileBoundary.write('\t}\n')
    fileBoundary.write(')\n')
    fileBoundary.close()

    #WRITE cellZones
#Count number of cellZones
    nrCellZones = 0
    cellZonesName = list()

    for grp in mesh.GetGroups():
        if grp.GetType() == SMESH.VOLUME:
            nrCellZones += 1
            cellZonesName.append(grp.GetName())

    if nrCellZones > 0:
        try:
            fileCellZones = open(dirname + '/cellZones', 'w')
        except Exception:
            print('Could not open the file cellZones, other files are ok.')
        __debugPrint__('Writing file cellZones\n')
        #create a dictionary where salomeIDs are keys
        #and OF cell ids are values.
        scToOFc = dict([sa, of] for of, sa in enumerate(volumes))
        __writeHeader__(fileCellZones, 'cellZones')
        fileCellZones.write('\n%d(\n' %nrCellZones)

        for grp in mesh.GetGroups():
            if grp.GetType() == SMESH.VOLUME:
                fileCellZones.write(grp.GetName() + '\n{\n')
                fileCellZones.write('\ttype\tcellZone;\n')
                fileCellZones.write('\tcellLabels\tList<label>\n')
                cellSalomeIDs = grp.GetIDs()
                nrGrpCells = len(cellSalomeIDs)
                fileCellZones.write('%d\n(\n' %nrGrpCells)
                for csId in cellSalomeIDs:
                    ofID = scToOFc[csId]
                    fileCellZones.write('%d\n' %ofID)

                fileCellZones.write(');\n}\n')
        fileCellZones.write(')\n')
        fileCellZones.flush()
        fileCellZones.close()

    totaltime = time.time() - starttime
    __debugPrint__('Finished writing to %s \n' % dirname)
    __debugPrint__('Converted mesh in %.3fs\n' % (converttime), 1)
    __debugPrint__('Wrote mesh in %.3fs\n' % (totaltime - converttime), 1)
    __debugPrint__('Total time: %.3fs\n' % totaltime, 1)
                   

def __writeHeader__(file,fileType,nrPoints=0,nrCells=0,nrFaces=0,nrIntFaces=0):
    """Write a header for the files points, faces, owner, neighbour"""
    file.write("/*" + "-"*68 + "*\\\n" )
    file.write("|" + " "*70 + "|\n")
    file.write("|" + " "*4 + "File exported from Salome Platform" +\
                   " using SalomeToFoamExporter" +" "*5 +"|\n")
    file.write("|" + " "*70 + "|\n")
    file.write("\*" + "-"*68 + "*/\n")
    file.write("FoamFile\n{\n")
    file.write("\tversion\t\t2.0;\n")
    file.write("\tformat\t\tascii;\n")
    file.write("\tclass\t\t")
    if(fileType =="points"):
        file.write("vectorField;\n")
    elif(fileType =="faces"):
        file.write("faceList;\n")
    elif(fileType =="owner" or fileType=="neighbour"):
        file.write("labelList;\n")
        file.write("\tnote\t\t\"nPoints: %d nCells: %d nFaces: %d nInternalFaces: %d\";\n" \
                       %(nrPoints,nrCells,nrFaces,nrIntFaces))
    elif(fileType == "boundary"):
        file.write("polyBoundaryMesh;\n")
    elif(fileType=="cellZones"):
        file.write("regIOobject;\n")
    file.write("\tlocation\t\"constant/polyMesh\";\n")
    file.write("\tobject\t\t" + fileType +";\n")
    file.write("}\n\n")

def __debugPrint__(msg,level=1):
    """Print only if level >= debug """
    if(debug >= level ):
      print(msg)
      #message = QDialog()
      #message.setWindowTitle("DEBUG PRINT MESSAGES")
      #Box_L = QVBoxLayout(message)
      #text = QLabel(msg)
      #Box_L.addWidget(text)
      #message.show()

def __verifyFaceOrder__(mesh,vnodes,fnodes):
    """
    Verify if the face order is correct. I.e. pointing out of the cell
    calc vol center
    calc f center
    calc ftov=fcenter-vcenter
     calc fnormal=first to second cross first to last
    if ftov dot fnormal >0 reverse order
    """
    vc=__cog__(mesh,vnodes)
    fc=__cog__(mesh,fnodes)
    fcTovc=__diff__(vc,fc)
    fn=__calcNormal__(mesh,fnodes)
    if(__dotprod__(fn,fcTovc)>0.0):
        return False
    else:
        return True

def __cog__(mesh,nodes):
    """
    calculate the center of gravity.
    """
    c=[0.0,0.0,0.0]
    for n in nodes:
        pos=mesh.GetNodeXYZ(n)
        c[0]+=pos[0]
        c[1]+=pos[1]
        c[2]+=pos[2]
    c[0]/=len(nodes)
    c[1]/=len(nodes)
    c[2]/=len(nodes)
    return c

def __calcNormal__(mesh,nodes):
    """
    Calculate and return face normal.
    """
    p0=mesh.GetNodeXYZ(nodes[0])
    p1=mesh.GetNodeXYZ(nodes[1])
    pn=mesh.GetNodeXYZ(nodes[-1])
    u=__diff__(p1,p0)
    v=__diff__(pn,p0)
    return __crossprod__(u,v)

def __diff__(u,v):
    """
    u - v, in 3D
    """
    res=[0.0]*3
    res[0]=u[0]-v[0]
    res[1]=u[1]-v[1]
    res[2]=u[2]-v[2]
    return res

def __dotprod__(u,v):
    """
    3D scalar dot product
    """
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def __crossprod__(u,v):
    """
    3D cross product
    """
    res=[0.0]*3
    res[0]=u[1]*v[2]-u[2]*v[1]
    res[1]=u[2]*v[0]-u[0]*v[2]
    res[2]=u[0]*v[1]-u[1]*v[0]
    return res

def findMeshByName(name):
    meshes=list()
    smesh = smeshBuilder.New()
    selobjID=salome.myStudy.FindObject(name)

    if selobjID is None:
        QMessageBox.critical(None,'Error',"Mesh `{}` is not found".format(name),QMessageBox.Abort)
        return None

    selobj=selobjID.GetObject()
    mesh=smesh.Mesh(selobj)
    meshes.append(mesh)
    return meshes

def __isGroupBaffle__(mesh,group,extFaces):
    for sid in group.GetIDs():
        if not sid in extFaces:
            __debugPrint__("group %s is a baffle\n" %group.GetName(),1)
            return True
    return False

def run():
    global dialog
    """
    Main function. Export the selected mesh.
    Will try to find the selected mesh.
    """
    dialog.setEnabled(False)
    meshes=findMeshByName(globalMeshName)
    for mesh in meshes:
        if not mesh == None:
            mName=mesh.GetName()
            #outdir=os.getcwd()+"/"+mName+"/constant/polyMesh"
            outdir=str(le_direcOutput.text())+"/constant/polyMesh"
            exportToFoam(mesh,outdir)
            __debugPrint__("finished exporting",1)
            QMessageBox.information(None,'Information',"Finish: Mesh export in " + le_direcOutput.text())
            dialog.close()
            #QMessageBox.information(None,'Information',)

def hide():
    global dialog
    dialog.hide()

def meshFile():
    global le_direcOutput
    PageName = QFileDialog.getExistingDirectory(qApp.activeWindow(),'Select output directory ')
    le_direcOutput.setText(str(PageName))

def meshSelected():
    global globalMeshName
    globalMeshName = meshLineEdit.text()
    showMainDialog()

def showMainDialog():
    global dialog
    dialog = QDialog()
    dialog.resize(300,100)
    dialog.setWindowTitle("Salome to OpenFOAM")
    layout = QGridLayout(dialog)
    meshes = findMeshByName(globalMeshName)
    l_direcOutput   = QLabel("Output Directory:")
    layout.addWidget(l_direcOutput,1,0)
    pb_direcOutput = QPushButton()
    pb_direcOutput.setText("...")
    layout.addWidget(le_direcOutput,2,0)
    layout.addWidget(pb_direcOutput,2,1)
    
    for mesh in meshes:
        if not mesh == None:
            mName=mesh.GetName()
            l_selectMesh = QLabel("Selected Mesh")
            le_selectMesh = QLineEdit(mName)
            le_selectMesh.setEnabled(False)
            l_groups = QLabel("Groups Mesh:")
            l_boundary = QLabel("Boundary:")
            layout.addWidget(l_selectMesh,3,0)
            layout.addWidget(le_selectMesh,4,0)
            cb_verify  = QCheckBox("Verify face order")
            layout.addWidget(cb_verify,5,0)
            layout.addWidget(l_groups,6,0)
            layout.addWidget(l_boundary,6,1)
            for gr in mesh.GetGroups():
                l_groupMesh = QLabel(gr.GetName())
                layout.addWidget(l_groupMesh)
                cmb_bounds = QComboBox()
                cmb_bounds.addItems(["patch","wall","symmetry","empty","wedge","cyclic"])
                bound.append(cmb_bounds)
                layout.addWidget(cmb_bounds)
                    
    okbox = QDialogButtonBox(dialog)
    okbox.setOrientation(QtCore.Qt.Horizontal)
    okbox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
    layout.addWidget(okbox)
    okbox.accepted.connect(run)
    okbox.rejected.connect(hide)
    pb_direcOutput.clicked.connect(meshFile)
    QtCore.QMetaObject.connectSlotsByName(dialog)
    dialog.show()


# globals ...
le_direcOutput = QLineEdit()
bound  = [] #list of boundaries

# Mesh selector
dialog = QDialog()
dialog.resize(300,50)

dialog.setWindowTitle("Select mesh by name")
meshLayout = QGridLayout(dialog)

meshLabel  = QLabel("Mesh name:")

meshLineEdit = QLineEdit()

meshSelect = QPushButton()
meshSelect.setText("Select")
meshSelect.clicked.connect(meshSelected)

meshLayout.addWidget(meshLabel,1,0)
meshLayout.addWidget(meshLineEdit,1,1)
meshLayout.addWidget(meshSelect,2,1)

QtCore.QMetaObject.connectSlotsByName(dialog)

dialog.show()

