#!/usr/bin/python

#=======================================================================#
#  															   							    # 
#  Oliver Gindele, UCL, 2015									   				   				    #
#  Version 1.1																	                            #
#  Calculate polarization in each B-centred unit cell compared to a reference system (inp.cell)				    #
#=======================================================================#


import numpy as np
import sys,string, os, math
import argparse

version=1.1

#get command arguments
parser = argparse.ArgumentParser(description='Calculate polarization from .cell file')
parser.add_argument('--file', '-f', type=str,  help="Input .cell file name", dest="inputFile")
parser.add_argument('--defect', '-d', type=int,  help="Input contains defects. Enter number of atoms of ideal structure.", dest="defect")
args = parser.parse_args()


########################### F U N C T I O N S ###########################

def getTotalAtom(inputFile):   
	"""get number of atoms .cell file file and save in array"""
	atomList=['Pb', 'Ti','Ba','Zr','Mg','Nb', 'Sc', 'Fe','O', 'F', 'V']

        i=0                                                                    		 					#counter for position
        try:
        	inputFile=open(inputFile,'r')
        except: 
        	print "No file named ", inputFile
        	sys.exit("Exit script")

        line=inputFile.readline()  
	
        while line:
		if np.size(string.split(line))>2:	
			if string.split(line)[0] in atomList    :             
		                i+=1
		line=inputFile.readline()
        inputFile.close()

        return i
        
def getCellVector(inputFile):   
	"""get cell vectors from .cell file file and save in array"""
        cellVector=np.empty(shape=(3,3),dtype=object)
                                                       		 					
        inputFile=open(inputFile,'r')
        line=inputFile.readline()  
	
        while line:
		if np.size(string.split(line))>1:	
			if string.split(line)[1]=='LATTICE_CART' or string.split(line)[1]=='lattice_cart'   :
				line=inputFile.readline();line=inputFile.readline()
				cellVector[0,0]=float(string.split(line)[0]);cellVector[0,1]=float(string.split(line)[1]); cellVector[0,2]=float(string.split(line)[2]);   
				line=inputFile.readline()   
				cellVector[1,0]=float(string.split(line)[0]);cellVector[1,1]=float(string.split(line)[1]); cellVector[1,2]=float(string.split(line)[2]);   
				line=inputFile.readline()   
				cellVector[2,0]=float(string.split(line)[0]);cellVector[2,1]=float(string.split(line)[1]); cellVector[2,2]=float(string.split(line)[2]);   
				break
		line=inputFile.readline()  
        inputFile.close()
        return cellVector

def getPos(inputFile,totalAtom):   
	"""get positions from .cell file file and save in array"""
	atomList=['Pb', 'Ti', 'Ba', 'Zr','Mg','Nb', 'Sc', 'Fe','O', 'F','V']
        position=np.empty(shape=(totalAtom,4),dtype=object)

        i=0                                                                    		 					#counter for position
        inputFile=open(inputFile,'r')
        line=inputFile.readline()  
	
        while line:
		if np.size(string.split(line))>2:	
			if string.split(line)[0] in atomList    :             
			       	position[i,0]=str(string.split(line)[0])             				 #get atom name
	                        
		                position[i,1]=float(string.split(line)[1])            #get x pos
		                position[i,2]=float(string.split(line)[2])            #get y pos
		                position[i,3]=float(string.split(line)[3])            #get z pos
		                i+=1
		line=inputFile.readline()
        inputFile.close()
        return position
        
def calcLattice(PbList,latticeDirection):
	"""calculate local lattice parameters from 8 corner Pb atoms"""
	directions={'x':0, 'X':0, 'y':1, 'Y':1, 'z':2, 'Z':2}									#dictionary for directions
	topPb=np.empty(3); bottomPb=np.empty(3);
	
	if latticeDirection=='y' or latticeDirection=='Y':
		PbList=PbList[PbList[:,2].argsort()]																	#sort along Z direction
		bottomPb=np.average(np.vstack((PbList[0:2],PbList[4:6])),axis=0); 	topPb=np.average(np.vstack((PbList[2:4],PbList[6:8])),axis=0);	#manually build list for Y direction 
		#(top and bottom can be degenerate) 				

	else:
		PbList=PbList[PbList[:,directions[latticeDirection]].argsort()]											#sort along direction
		bottomPb=np.average(PbList[0:4],axis=0); 	topPb=np.average(PbList[4:8],axis=0);						#take 4 atoms each as top and bottom layer

	localCellVector=np.array([topPb[0]-bottomPb[0], topPb[1]-bottomPb[1], topPb[2]-bottomPb[2]])
	localCellParameter=np.sqrt((topPb[0]-bottomPb[0])**2+(topPb[1]-bottomPb[1])**2+(topPb[2]-bottomPb[2])**2)	#calculate cell parameter difference between top and bottom layer
	
	return localCellParameter, localCellVector     

def calcTransformationMatrix(Va,Vb):
	"""create transformation (rotation) matrix between two vectors.
	transformation matrix rotates vectot Va onto Vb"""
	
	Va=Va/np.linalg.norm(Va)		#first vector = 0 0 c
	Vb=Vb/np.linalg.norm(Vb)		#second vector = actual vector of unit cell

	GMatrix=[
	[np.dot(Va,Vb), np.linalg.norm(np.cross(Va,Vb))*-1, 0],
	[np.linalg.norm(np.cross(Va,Vb)), np.dot(Va,Vb) , 0],
	[0 , 0 ,1]]

	if np.linalg.norm(Vb-np.dot(np.dot(Va,Vb),Va))==0:				#if Va is along Vb this is zero: catch this warning
		basisMatrix=[[1,0,0],[0,1,0],[0,0,1]]
	else:
		basisMatrix=np.transpose([ Va, (Vb-np.dot(np.dot(Va,Vb),Va)) / np.linalg.norm(Vb-np.dot(np.dot(Va,Vb),Va)),  np.cross(Vb,Va) ])
	
	try:
		transformationMatrix=np.dot(basisMatrix,np.dot(GMatrix,np.linalg.inv(basisMatrix)))
	except np.linalg.linalg.LinAlgError:
		transformationMatrix=np.array([[1,0,0],[0,1,0],[0,0,1]])	#in case of cell aligned alon c already 
	
	return transformationMatrix
	
def printResult(result,name):
	"""print result array to file"""
		
	np.set_printoptions(threshold=np.nan)											#print whole array
	
	#print to file
	f = open(name+'.dat', 'w')
	print >>f, "#position" ,name
	for line in xrange(totalAtom/5):
		print >>f, str(result[line]).replace('[','').replace(']','')
	f.close()
     
###########################    M A I N    R O U T I N E     ###########################
inputFile=args.inputFile						#parse input arguments (file name)
defect=args.defect

#header
print "***** D I F F P O L *****"
print "\nVersion: ", version, "\n"

totalAtom=getTotalAtom(inputFile); totalAtomOriginal=totalAtom
print 'Number of atoms: ',totalAtom
if defect!=None:
	totalAtom=defect

cellVector=getCellVector(inputFile)
print '\nCell vectors:\n'
print str(cellVector).replace('[','').replace(']','') ,'\n'
print 'Cell volume per formula unit: ', np.dot(cellVector[0,:],np.cross(cellVector[1,:],cellVector[2,:]))/totalAtom*5

position=getPos(inputFile,totalAtom)

#Born Effective Charge tensors Cubic:
BqTi=np.array([[7.05221,0,0],  [0, 7.05221,0], [0,0, 7.05221]])
BqNb=np.array([[8.05221,0,0],  [0, 8.05221,0], [0,0, 8.05221]])
BqV=np.array([[8.05221,0,0],  [0, 8.05221,0], [0,0, 8.05221]])
BqFe=np.array([[6.05221,0,0],  [0, 6.05221,0], [0,0, 6.05221]])
BqSc=np.array([[6.05221,0,0],  [0, 6.05221,0], [0,0, 6.05221]])
BqFe=np.array([[6.05221,0,0],  [0, 6.05221,0], [0,0, 6.05221]])

BqPb=np.array([[ 3.96963,0,0], [0, 3.96963,0], [0,0, 3.96963]])
BqBa=np.array([[2.75936,0.,0.],[0.,2.75936,0.],[0.,0.,2.75936]])

BqO1=np.array([[-2.61443 ,0,0],[0,-2.61443 ,0],[0,0,-5.79298]]) # 0.5 0.5 0.0
BqO2=np.array([[ -5.79298,0,0],[0, -2.61443,0],[0,0,-2.61443]]) # 0.0 0.5 0.5
BqO3=np.array([[ -2.61443,0,0],[0, -5.79298,0],[0,0,-2.61443]]) # 0.5 0.0 0.5

############           multiplacte system to account for PBC               ############

extendedPosition=np.zeros(shape=(totalAtom*27,4),dtype=object)
k=0 #counter
for i in xrange(totalAtomOriginal):
	extendedPosition[k,:]=position[i,:]; k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]-1; extendedPosition[k,2]=position[i,2]-1;  extendedPosition[k,3]=position[i,3]-1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]; extendedPosition[k,2]=position[i,2]-1;  extendedPosition[k,3]=position[i,3]-1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]+1; extendedPosition[k,2]=position[i,2]-1;  extendedPosition[k,3]=position[i,3]-1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]-1; extendedPosition[k,2]=position[i,2];  extendedPosition[k,3]=position[i,3]-1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]; extendedPosition[k,2]=position[i,2];  extendedPosition[k,3]=position[i,3]-1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]+1; extendedPosition[k,2]=position[i,2];  extendedPosition[k,3]=position[i,3]-1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]-1; extendedPosition[k,2]=position[i,2]+1;  extendedPosition[k,3]=position[i,3]-1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]; extendedPosition[k,2]=position[i,2]+1;  extendedPosition[k,3]=position[i,3]-1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]+1; extendedPosition[k,2]=position[i,2]+1;  extendedPosition[k,3]=position[i,3]-1;  k=k+1
	
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]-1; extendedPosition[k,2]=position[i,2]-1;  extendedPosition[k,3]=position[i,3];  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]; extendedPosition[k,2]=position[i,2]-1;  extendedPosition[k,3]=position[i,3];  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]+1; extendedPosition[k,2]=position[i,2]-1;  extendedPosition[k,3]=position[i,3];  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]-1; extendedPosition[k,2]=position[i,2];  extendedPosition[k,3]=position[i,3];  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]+1; extendedPosition[k,2]=position[i,2];  extendedPosition[k,3]=position[i,3];  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]-1; extendedPosition[k,2]=position[i,2]+1;  extendedPosition[k,3]=position[i,3];  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]; extendedPosition[k,2]=position[i,2]+1;  extendedPosition[k,3]=position[i,3];  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]+1; extendedPosition[k,2]=position[i,2]+1;  extendedPosition[k,3]=position[i,3];  k=k+1
	
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]-1; extendedPosition[k,2]=position[i,2]-1;  extendedPosition[k,3]=position[i,3]+1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]; extendedPosition[k,2]=position[i,2]-1;  extendedPosition[k,3]=position[i,3]+1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]+1; extendedPosition[k,2]=position[i,2]-1;  extendedPosition[k,3]=position[i,3]+1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]-1; extendedPosition[k,2]=position[i,2];  extendedPosition[k,3]=position[i,3]+1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]; extendedPosition[k,2]=position[i,2];  extendedPosition[k,3]=position[i,3]+1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]+1; extendedPosition[k,2]=position[i,2];  extendedPosition[k,3]=position[i,3]+1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]-1; extendedPosition[k,2]=position[i,2]+1;  extendedPosition[k,3]=position[i,3]+1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]; extendedPosition[k,2]=position[i,2]+1;  extendedPosition[k,3]=position[i,3]+1;  k=k+1
	extendedPosition[k,0]=position[i,0]; extendedPosition[k,1]=position[i,1]+1; extendedPosition[k,2]=position[i,2]+1;  extendedPosition[k,3]=position[i,3]+1;  k=k+1
	
	
############            Get atoms                ############
BDist=np.zeros(shape=(totalAtom/5,1,7),dtype=float); PbDist=np.zeros(shape=(totalAtom/5,8,7),dtype=float); 
O1Dist=np.zeros(shape=(totalAtom/5,2,7),dtype=float);O2Dist=np.zeros(shape=(totalAtom/5,2,7),dtype=float); O3Dist=np.zeros(shape=(totalAtom/5,2,7),dtype=float);
polarization=np.zeros(shape=(totalAtom/5,7),dtype=float); volume=np.zeros(shape=(totalAtom/5,4),dtype=float)
polAngle=np.zeros(shape=(totalAtom/5,6),dtype=float)
incompleteCell=[]; 

m=0
for j in xrange(totalAtomOriginal):
	
	if position[j][0]=='Ti' or position[j][0]=='Zr' or position[j][0]=='Nb' or position[j][0]=='Fe' or position[j][0]=='Sc' or position[j][0]=='V':
	
		#initialize counters and lists
		PbCounter=0; OCounter=0
		PbList=np.zeros(shape=(8,3),dtype=float); OList=np.empty(shape=(6,3),dtype=float); OListAngles=np.zeros(shape=(6),dtype=float);TiList=np.zeros(shape=(1,3),dtype=float)
		unitCellParameter=np.zeros(3)
		unitCellVector=np.zeros(shape=(3,3),dtype=float)
		PbRefList=np.zeros(shape=(8,3),dtype=float);	O1RefList=np.zeros(shape=(2,3),dtype=float); O2RefList=np.zeros(shape=(2,3),dtype=float); O3RefList=np.zeros(shape=(2,3),dtype=float)
		polarizationTi=0; polarizationPb=0; polarizationO=0;
		
		BName=position[j][0]
		BList=np.dot(cellVector,position[j][1:4])
		#create Lists for Pb and O in unit cell in absolute coordinates
		for i in xrange(np.size(extendedPosition)/4):
			if extendedPosition[i][0]=='Pb' or extendedPosition[i][0]=='Ba':
				if np.sqrt(((extendedPosition[i][1]-position[j][1])*cellVector[0,0])**2+((extendedPosition[i][2]-position[j][2])*cellVector[1,1])**2+ ((extendedPosition[i][3]-position[j][3])*cellVector[2,2])**2)  <=4.0:
					PbList[PbCounter]=np.dot(cellVector,extendedPosition[i][1:4])
					PbCounter+=1
			
			if extendedPosition[i][0]=='O' or extendedPosition[i][0]=='F':
				if np.sqrt(((extendedPosition[i][1]-position[j][1])*cellVector[0,0])**2+((extendedPosition[i][2]-position[j][2])*cellVector[1,1])**2+ ((extendedPosition[i][3]-position[j][3])*cellVector[2,2])**2)  <=2.8:
					OList[OCounter]=np.dot(cellVector,extendedPosition[i][1:4])
					OCounter+=1

		#check for missing atoms				
		if  PbCounter <8:
			print "Missing Pb atom!"
			incompleteCell=np.append(incompleteCell,m)
		if  OCounter <6:
			print "Missing O atom!"
			incompleteCell=np.append(incompleteCell,m)

		#calculate reference cell
		BRefList=np.array([np.mean(PbList[:,0]),np.mean(PbList[:,1]),np.mean(PbList[:,2])])
		
		#get local unit cell parameters
		unitCellParameter[0],unitCellVector[0]=calcLattice(PbList,'x'); unitCellParameter[1],unitCellVector[1]=calcLattice(PbList,'y'); unitCellParameter[2],unitCellVector[2]=calcLattice(PbList,'z')

		#get primary cell
		primeAxis=unitCellVector[np.argsort(unitCellParameter)[2]]
		
		#calculate rotation matrix
		transformationMatrix=calcTransformationMatrix([0,0,unitCellParameter[2]], primeAxis)

		#find planar and axial oxygens
		for i in xrange(6):
			OListAngles[i]=np.dot((OList[i]-BRefList)/np.linalg.norm(OList[i]-BRefList), primeAxis/np.linalg.norm(primeAxis))	#calculate angle between O-B and prime cell Axis to distinguish bewteen oxygens
			
		O1List=OList[np.argsort(np.abs(OListAngles))][4:6]
		OpList=OList[np.argsort(np.abs(OListAngles))][0:4]
		O3List=OpList[np.argsort(OpList[:,0])][1:3]; O2List=OpList[np.argsort(OpList[:,1])][1:3]; 
		
		############            Build Reference Lists                ############
		
		#build Pb reference list
		PbRefList[0,:]=[-0.5,-0.5,-0.5]; PbRefList[1,:]=[+0.5,-0.5,-0.5]; PbRefList[2,:]=[-0.5,+0.5,-0.5]; PbRefList[3,:]=[+0.5,+0.5,-0.5];
		PbRefList[4,:]=[-0.5,-0.5,+0.5]; PbRefList[5,:]=[+0.5,-0.5,+0.5]; PbRefList[6,:]=[-0.5,+0.5,+0.5]; PbRefList[7,:]=[+0.5,+0.5,+0.5];
		PbRefList[:,0]=PbRefList[:,0]*unitCellParameter[0]; PbRefList[:,1]=PbRefList[:,1]*unitCellParameter[1]; PbRefList[:,2]=PbRefList[:,2]*unitCellParameter[2]; 
		
		for i in xrange(8):
			PbRefList[i]=np.dot(transformationMatrix,PbRefList[i])+BRefList			#rotate reference atoms and shift to B-cation
		PbList=np.sort(PbList,axis=0); PbRefList=np.sort(PbRefList,axis=0)			#sort actual PbList so it matches the reference PbList
		
		#build O1 (axial Oxygen) reference list
		O1RefList[0,:]=[0,0,-0.5]; O1RefList[1,:]=[0,0,0.5]; 
		O1RefList[:,0]=O1RefList[:,0]*unitCellParameter[0]; O1RefList[:,1]=O1RefList[:,1]*unitCellParameter[1]; O1RefList[:,2]=O1RefList[:,2]*unitCellParameter[2]; 
		
		for i in xrange(2):
			O1RefList[i]=np.dot(transformationMatrix,O1RefList[i])+BRefList			#rotate reference atoms and shift to B-cation
		O1List=np.sort(O1List,axis=0); O1RefList=np.sort(O1RefList,axis=0)			#sort actual PbList so it matches the reference PbList
			
		#build O2 (planar Oxygen) reference list
		O2RefList[0,:]=[-0.5,0,0]; O2RefList[1,:]=[+0.5,0,0]; 
		O2RefList[:,0]=O2RefList[:,0]*unitCellParameter[0]; O2RefList[:,1]=O2RefList[:,1]*unitCellParameter[1]; O2RefList[:,2]=O2RefList[:,2]*unitCellParameter[2]; 

		for i in xrange(2):
			O2RefList[i]=np.dot(transformationMatrix,O2RefList[i])+BRefList			#rotate reference atoms and shift to B-cation
		O2List=np.sort(O2List,axis=0); O3RefList=np.sort(O2RefList,axis=0)			#sort actual PbList so it matches the reference PbLis
		
		#build O3 (planar Oxygen) reference list
		O3RefList[0,:]=[0,-0.5,0]; O3RefList[1,:]=[0,+0.5,0]; 
		O3RefList[:,0]=O3RefList[:,0]*unitCellParameter[0]; O3RefList[:,1]=O3RefList[:,1]*unitCellParameter[1]; O3RefList[:,2]=O3RefList[:,2]*unitCellParameter[2]; 

		for i in xrange(2):
			O3RefList[i]=np.dot(transformationMatrix,O3RefList[i])+BRefList			#rotate reference atoms and shift to B-cation
		O3List=np.sort(O3List,axis=0); O3RefList=np.sort(O3RefList,axis=0)			#sort actual PbList so it matches the reference PbLis
		
		
		############            Calculate Displacements                ############
		#for B
		for i in xrange(1):
			try:
				BDist[m,i,0:3]=BList
			except IndexError:
				print "\nWarning: Missing atom!  \nAre there defects in the system? Try using -d (defect) option and specify number of atoms in defect free system"
				print "Exit"
				sys.exit()
				
			BDist[m,i,3:6]=BList-BRefList
			BDist[m,i,6]=np.dot(BDist[m,i,3:6],primeAxis/np.linalg.norm(primeAxis))		

		#for Pb
		for i in xrange(8):
			PbDist[m,i,0:3]=BList
			PbDist[m,i,3:6]=PbList[i]-PbRefList[i]
			PbDist[m,i,6]=np.dot(PbDist[m,i,3:6],primeAxis/np.linalg.norm(primeAxis))
	
		#for O1 
		for i in xrange(2):
			O1Dist[m,i,0:3]=BList
			O1Dist[m,i,3:6]=O1List[i]-O1RefList[i]
			O1Dist[m,i,6]=np.dot(O1Dist[m,i,3:6],primeAxis/np.linalg.norm(primeAxis))		
				
		#for O2
		for i in xrange(2):
			O2Dist[m,i,0:3]=BList
			O2Dist[m,i,3:6]=O2List[i]-O2RefList[i]
			O2Dist[m,i,6]=np.dot(O2Dist[m,i,3:6],primeAxis/np.linalg.norm(primeAxis))		
				
		#for O3
		for i in xrange(2):
			O3Dist[m,i,0:3]=BList
			O3Dist[m,i,3:6]=O3List[i]-O3RefList[i]
			O3Dist[m,i,6]=np.dot(O3Dist[m,i,3:6],primeAxis/np.linalg.norm(primeAxis))		
			
		############            Calculate Polarizaion                ############
		#B contribution
		if BName=='Ti':
			polarizationTi=np.dot(BqTi,BDist[m,0,3:6])				
		elif BName=='Nb':
			polarizationTi=np.dot(BqNb,BDist[m,0,3:6])				
		elif BName=='Sc':
			polarizationTi=np.dot(BqSc,BDist[m,0,3:6])
		elif BName=='Fe':
			polarizationTi=np.dot(BqFe,BDist[m,0,3:6])
		
		for i in xrange(8):					
			polarizationPb+=np.dot(BqPb,PbDist[m,i,3:6])/8		#Pb contribution
			#polarizationPb+=np.dot(BqBa,PbDist[m,i,3:6])/8		#Pb contribution
		for i in xrange(2):					
			polarizationO+=np.dot(BqO1,O1Dist[m,i,3:6])/2		#O1 contribution
		for i in xrange(2):					
			polarizationO+=np.dot(BqO2,O2Dist[m,i,3:6])/2		#O2 contribution
		for i in xrange(2):					
			polarizationO+=np.dot(BqO3,O3Dist[m,i,3:6])/2		#O3 contribution
			
		polarization[m,0:3]=BList
		polarization[m,3:6]=(polarizationTi+polarizationPb+polarizationO)*float(1.602177*10**(-19)*float(10**6))/(10**(-16))/np.abs((np.dot(unitCellVector[0,:],np.cross(unitCellVector[1,:],unitCellVector[2,:]))))
		polarization[m,6]=np.sqrt(polarization[m,3]**2+polarization[m,4]**2+polarization[m,5]**2)
		
		#calculate Angle
		polAngle[m,0:3]=BList
		for i in range(3,6):
			polAngle[m,i]=math.degrees(np.arccos(polarization[m,i]/np.linalg.norm(polarization[m,3:6])))
			
		#save cell volume
		volume[m,0:3]=BList
		volume[m,3]=np.abs((np.dot(unitCellVector[0,:],np.cross(unitCellVector[1,:],unitCellVector[2,:]))))
		
		m+=1											#B-cation counter

#set polarization on imcomplete cells to zero
if incompleteCell!=[]:
	for i in incompleteCell:
		polarization[i,3:7]=[0,0,0,0]

############            Print               ############	
#print options
np.set_printoptions(precision=6)
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=150)

print  "\nPolarization calculated for ", m, "B-cations" 
print  " Average polarization value: ", str(np.average(polarization[:,3:7],axis=0)[0:3]).replace('[','').replace(']',''), '\t', str(np.average(polarization[:,3:7],axis=0)[3]).replace('[','').replace(']','') 

############            Print   to file            ############
printResult(polarization[np.argsort(polarization[:,1])], 'polarization')
printResult(PbDist[np.argsort(PbDist[:,0,1])], 'Pb_dist')
printResult(BDist[np.argsort(PbDist[:,0,1])], 'Ti_dist')
printResult(O1Dist[np.argsort(PbDist[:,0,1])], 'O1_dist')
printResult(O2Dist[np.argsort(PbDist[:,0,1])], 'O2_dist')
printResult(O3Dist[np.argsort(PbDist[:,0,1])], 'O3_dist')			
printResult(polAngle[np.argsort(polAngle[:,1])], 'angle_polarization')
printResult(volume[np.argsort(volume[:,1])], 'volume')
     
