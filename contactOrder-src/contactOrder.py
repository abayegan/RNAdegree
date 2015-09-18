#! /usr/bin/python2.7

""" Amir Bayegan, May 2015
RNA contact order for the secondary and tertiary structures are computed. 
The secondary structure contact order is computed for a secondary structure obtainded from maximizing the number of base pairs predicted by RNAVIEW from pdb files. 
3D contact order is calculated after parsing a pdb file and selecting all the atoms in contact(<6 angstrom)

Input: a file containing paths of the pdbfiles. Each line of the file should include just one path.
Output: 1) a folder named RNAVIEW_out containing the output files of rnaview.
		 2) two folders named 3DcontactList and compatible3DcontactList. For each pdb file there is a contactList in each folder keeping information about
		 the nucleotides which are in contact.
		 3)files with extension _2DcontactOrder _3DcontactOrder and viewCompat3DcontactOrder containing the calculated contact order values

"""
			
 
 
import os,sys,glob,runRnaviewPeter, secStrFromRNAVIEW, misc, math, shutil

PRINTCONTACTLIST=1

purineAtoms=['N9','C8','N7','C5','C6','O6','N6','N1','C2','N2','N3','C4']
pyrimidineAtoms=['N1','C2','O2','N3','C4','N4','O4','C5','C6']

threeLetterAAcodes = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
 'SER', 'THR', 'TRP', 'TYR', 'VAL' ]

def isPurine(nuc):
	if (nuc=='A' or nuc=='G'):
		return 1;
	else:
		return 0;

def isPyrimidine(nuc):
	if (nuc=='U' or nuc=='C'):
		return 1;
	else:
		return 0;

def createDirectories():
	dirList = ["RNAVIEW_out","3DcontactList","compatible3DcontactList"]
	for f in dirList:
		if not os.path.exists(f):
			os.makedirs(f)
		else:
			r = glob.glob(f+"/*")
			for i in r:
				os.remove(i)
	return
	
def get2DcontactOrder(secStr):
	co=0.0
	bpList = misc.basePairList(secStr)
	for bp in bpList:
		co += abs(bp[1]-bp[0])
	if len(bpList) !=0:
		a_co = float(co)/len(bpList)
		r_co = co/(len(secStr)*len(bpList))
	else:
		a_co = float("inf")
		r_co = float("inf")
	return a_co,r_co
	
def getRecSpec(line):
	if len(line[17:20])==1:
		res = line[17:20].strip()
	elif line[17:20] in threeLetterAAcodes:
		res = -1
	else:
		res = line[19]
	return line[22:26].strip(),res ,line[12:16].strip(),line[21].strip(),line[30:38].strip(),line[38:46].strip(),line[46:55].strip()

def euclideanDist(c1,c2):
	return math.sqrt(math.pow(float(c2[2])-float(c1[2]),2) + math.pow(float(c2[1])-float(c1[1]),2) + math.pow(float(c2[0])-float(c1[0]),2))

def get3DpairwiseDistance(pdbName):
	d={};rnaDict={};ATOMseq={};numAtomsInCont=0;contactList=[];viewCompatibleContactList=[];prevPos=-1
	pdbId = pdbName.split('/')[-1].split(".")[0]
	pdbFile=open(pdbName,'r')
	line = pdbFile.readline();
	while line:
		if line[:6]=="ATOM  " or line[:6]=="HETATM":
			pos,nuc,atom,chainId,x,y,z = getRecSpec(line);
			if nuc==-1: #an amino acid is reached
				break
			#print pos,nuc,atom,chainId,x,y,z 
			if (pos,chainId) not in d.keys():
				d[(pos,chainId)]={}
			if (isPyrimidine(nuc) and atom in pyrimidineAtoms) or (isPurine(nuc) and atom in purineAtoms):
				d[(pos,chainId)][atom] = (nuc,x,y,z)
				#print pos,chainId,d[(pos,chainId)][atom]
				if prevPos!=pos:
					if chainId not in ATOMseq.keys():
						ATOMseq[chainId]=nuc
					else: 
						ATOMseq[chainId]+=nuc
				prevPos=pos
		elif line[:6]=="ENDMDL": #only the 1st model is condidered if more than one model exists
			break
		line = pdbFile.readline();
	#print ATOMseq
	rnaSeqDict = runRnaviewPeter.getRNAfromSEQRES(pdbName)
	#print rnaSeqDict
	
	key = d.keys()
	key.sort()
	#print key
	#find and report discrepencies between SEQRES and ATOM fields
	for c in ATOMseq.keys():
		rnaDict[c]=[rnaSeqDict[c],ATOMseq[c]]
		if len(ATOMseq[c])!=len(rnaSeqDict[c]):
			print "ATOM records in %s , chain %s have %d missing positions!" % (pdbName,c,len(rnaSeqDict[c])-len(ATOMseq[c]))

	#compute pairwise distance between all atoms of all bases
	for pos1 in key:
		for pos2 in key:
			if pos1[0] != pos2[0] and pos1[1]==pos2[1] : # check for different positions and the same chain id
				for atom1 in d[pos1].keys():
					for atom2 in d[pos2].keys():
						if euclideanDist(d[pos1][atom1][1:],d[pos2][atom2][1:]) < 6: # atoms are in contact if their distance is less than 6 Angstroms
							numAtomsInCont +=1
							break;  #if one atom is in contacts with at least one other heavy atom it is sufficients
				if numAtomsInCont>=1: #change here to put a threshold for the minimum number of atoms in contact between to base pairs
					contactList.append([pdbId,pos1[1],d[pos1][atom1][0],pos1[0],pos2[0],numAtomsInCont])
				numAtomsInCont=0
	if(PRINTCONTACTLIST==1):			
		outFile=open("3DcontactList/%s_3DContactList"%pdbId,'w')
		outFile.write("#Atom1\tAtom2\tPosition1\tPosition2\tNumber of atoms in contact\n")
		for c in contactList:
			outFile.write("%s\t%s\t%s\t%s\t%d\n" % (c[1],c[2],c[3],c[4],c[5]))
		outFile.close()
		pdbFile.close()
	
	#compute distance between atoms which appear in RNAVIEW	
	numAtomsInCont=0;viewList=[]
	for chain in ATOMseq.keys():
		viewName = pdbId+chain+".rnaview.out"
		if os.path.exists(viewName):
			viewFile=open(viewName,'r')
			for line in viewFile.readlines()[1:]:
				pos1 = line.split('\t')[0]
				pos2 = line.split('\t')[2].split()[0]
				if pos2!='0':
					viewList.append((pos1,chain))
					viewList.append((pos2,chain))
			viewFile.close()
	
	#posList = [(i[0]) for i in key]
	#print posList
	for pos1 in set(viewList):
		for pos2 in set(viewList):
			if pos1!=pos2:
				if pos1 in key and pos2 in key:
					for atom1 in d[pos1]:
						for atom2 in d[pos2]:
							if euclideanDist(d[pos1][atom1][1:],d[pos2][atom2][1:]) < 6: # atoms are in contact if their distance is less than 6 Angstroms
								numAtomsInCont +=1
								break;  #if one atom is in contacts with at least one other heavy atom it is sufficients
					if numAtomsInCont>=1: #change here to put a threshold for the minimum number of atoms in contact between to base pairs
						viewCompatibleContactList.append([pdbId,pos1[1],d[pos1][atom1][0],pos1[0],pos2[0],numAtomsInCont])
					numAtomsInCont=0
				else:
					print "position %s or %s of rnaview for %s not found in pdb" % (pos1[0],pos2[0],pdbId)
					return
	if(PRINTCONTACTLIST==1):		
		outFile=open("compatible3DcontactList/%s_viewCompat3DContactList"%(pdbId),'w')
		outFile.write("#Atom1\tAtom2\tPosition1\tPosition2\tNumber of atoms in contact\n")
		for c in viewCompatibleContactList:
			outFile.write("%s\t%s\t%s\t%s\t%d\n" % (c[1],c[2],c[3],c[4],c[5]))
		outFile.close()
	return contactList,rnaDict,viewCompatibleContactList

def get3DcontactOrder(contactList,rnaDict):
	#hainIdList= set([int(i[0]) for i in contactList[1]])
	d={};cnt=0;COlist=[]
	for l in contactList:
		diff = abs(float(l[3])-float(l[4]))
		if diff > 1: #adjacent nucleotides are not considered in the definition of contact order
			if (l[0],l[1]) in d: #(pdb,chainId) 
				d[(l[0],l[1])] += diff
			else:
				d[(l[0],l[1])] = diff
			cnt += 1
	for key in d:
		abs_co = d[key]/cnt
		rel_co = abs_co/len(rnaDict[key[1]][1]) #normalized by the sequence length obtained from ATOM
		COlist.append((key[0],key[1],abs_co,rel_co,rnaDict[key[1]][0],len(rnaDict[key[1]][1])))
	return COlist

#rnaview program is called to get the base pairing information from pdb files. 
#Then a secondary structure with maximal number of such base pairs is computed using Peter's code
#The output files from rnaview are saved at RNAVIEW_out folder.

def run2DContactOrder(pdbPath):
	allowedBP = ["WWc","WWt","SHc","SHt","WSc","WSt","SSt","SSc","HHt","HHc","WHc","WHt","PP","MM"]
	co2Ddict={};rnaDict={}
	shutil.copyfile(pdbPath,"RNAVIEW_out/"+pdbPath)
	os.chdir("RNAVIEW_out")
	runRnaviewPeter.RNAview(pdbPath,allowedBP)
	for viewOutFile in glob.glob("*.rnaview.out"):
		pdbId = viewOutFile.split('/')[-1].split(".")[0]
		rna,secStr = secStrFromRNAVIEW.runSecStr(viewOutFile)
		rnaDict[pdbId]=rna
		abs_co,rel_co = get2DcontactOrder(secStr)
		co2Ddict[pdbId] = [rna,secStr,abs_co,rel_co]
	s = "pdb\tAbs_Co\tRelative_Co\trna\tsecStr\tlength\n"
	key = co2Ddict.keys()
	key.sort()
	for k in key:
		s += "%s\t%f\t%f\t%s\t%s\t%d\n" %(k,co2Ddict[k][2],co2Ddict[k][3],co2Ddict[k][0],co2Ddict[k][1],len(co2Ddict[k][0]))
	os.chdir("..")
	out = open(pdbPath+"_2DcontactOrder",'w')
	out.write(s)
	out.close()

def run3DContactOrder(pdbPath):
	COlist =[];viewList=[]
	fileList = map(lambda x:x[:-1], open(pdbPath,"r").readlines())
	for f in fileList:
		if os.path.exists(f):
			contactList,rnaDict,viewCompatibleContactList = get3DpairwiseDistance(f)
			if contactList:              
				COlist += get3DcontactOrder(contactList,rnaDict)
			else:
				print "%s is excluded from 3D computation. Pdb file does not seem like RNA!" % f
				continue
			if viewCompatibleContactList:
				viewList += get3DcontactOrder(viewCompatibleContactList,rnaDict)
			else:
				print "No base pairs reported for %s by RNAVIEW. It is excluded from 2D and pseudo-3D!" % f
				continue
		else: print "file %s not found!" %f               
	s = "pdb\tAbs_Co\tRelative_Co\trna\tlength(ATOM)\n"
	for l in COlist:
		s += "%s%s\t%f\t%f\t%s\t%d\n" %(l[0],l[1],l[2],l[3],l[4],l[5])
	out = open(pdbPath+"_3DcontactOrder",'w')
	out.write(s)
	out.close()
	
	s = "pdb\tAbs_Co\tRelative_Co\trna\tlength(ATOM)\n"
	for l in viewList:
		s += "%s%s\t%f\t%f\t%s\t%d\n" %(l[0],l[1],l[2],l[3],l[4],l[5])
	out = open(pdbPath+"_viewCompat3DcontactOrder",'w')
	out.write(s)
	out.close()

if __name__ == "__main__":
	if len(sys.argv)<2:
		print "USAGE:%s pdbPath"%sys.argv[0]
		sys.exit(1)
	pdbPath=sys.argv[1]
	createDirectories()
	run2DContactOrder(pdbPath)
	run3DContactOrder(pdbPath)
	
	
