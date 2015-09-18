#! /usr/bin/python

# runRnaviewPeter.py

# Sat Dec  2 10:51:54 EST 2006

# F.Ferre
# P.Clote modified in 2 ways
# 1) ensure that triples of form (i,nucl,j) output where pos i and j base pair
#    previously the program output (i,nucl1,nucl2). Perhaps this original
#    output is necessary for DIAL.py but perhaps this is a mistake in DIAL
#    as well
# 2) retrieve RNA nucleotide sequence from SEQRES rather than ATOM records
#    this is important, since the ATOM records can have gaps. This happens
#    in PDB files 2hok and 2hoj for instance. So this is IMPORTANT as well
#    for DIAL

# run rnaview on a list of PDB files,
# parse the output and for each chain produces an
# output file in bpseq format.
# By default, consider only standard cis Watson-Crick base pairs.
# In verbose mode, output contains base-pair type following Leontis notation 
# and base pair configuration (cis/trans).
# NB: rnaview creates output files in the directory where input pdb files are.
# If variable DELETE_OUT is set to 1, these files are deleted.

import os,sys,string,subprocess,shlex

#RNAVIEW    = '/usr/local/bin/rnaview'
LIST       = ""
VERBOSE    = 1
DELETE_OUT = 0
BPconvert  = {("W/W","tran"):"WWt",
              ("W/W","cis") :"WWc",
              ("S/H","cis") :"SHc",
              ("H/S","cis") :"SHc",
              ("S/H","tran"):"SHt",
              ("H/S","tran"):"SHt",
              ("S/W","tran"):"WSt",
              ("W/S","tran"):"WSt",
              ("S/W","cis") :"WSc",
              ("W/S","cis") :"WSc",
              ("S/S","tran"):"SSt",
              ("S/S","cis") :"SSc",
              ("H/H","tran"):"HHt",
              ("H/H","cis") :"HHc",
              ("W/H","tran"):"WHt",
              ("H/W","tran"):"WHt",
              ("W/H","cis") :"WHc",
              ("H/W","cis") :"WHc",
              ("+/+","cis") :"PP",
              ("-/-","cis") :"MM"}

def usage():
  print '''Usage: %s 
  -a listOfPDBfiles
  -v [0/1] writes base-pair class and configuration in output files [default 0]
  -MM  [0/1] consider standard AU,AT Watson-Crick/Watson-Crick (-/- in rnaview) [default 1]
  -PP  [0/1] consider standard GC Watson-Crick/Watson-Crick (+/+ in rnaview) [default 1]
  -WWc [0/1] consider Watson-Crick/Watson-Crick cis [default 0]
  -WWt [0/1] consider Watson-Crick/Watson-Crick trans [default 0]
  -WSc [0/1] consider Watson-Crick/Sugar cis [default 0]
  -WSt [0/1] consider Watson-Crick/Sugar trans [default 0]
  -WHc [0/1] consider Watson-Crick/Hoogsten cis [default 0]
  -WHt [0/1] consider Watson-Crick/Hoogsten trans [default 0]
  -SHc [0/1] consider Sugar/Hoogsten cis [default 0]
  -SHt [0/1] consider Sugar/Hoogsten trans [default 0]
  -SSc [0/1] consider Sugar/Sugar cis [default 0]
  -SSt [0/1] consider Sugar/Sugar trans [default 0]
  -HHc [0/1] consider Hoogsten/Hoogsten cis [default 0]
  -HHt [0/1] consider Hoogsten/Hoogsten trans [default 0]
  '''% sys.argv[0]

def writeBpseqFile(BP,PDB,filename,chain):
  OutputList = [] #list of triples (i,nucl,j)
  if chain.strip() == "":
    chain = "_"
  c = 0
  pdbcode = filename.split("/")[-1].split(".")[0]
  O = open("%s%s.rnaview.out" % (pdbcode,chain),"w")
  O.write("#%s %s\n" % (pdbcode,chain))
  K = PDB.keys()
  K.sort()
  for k in K:
    tmpList = []  #temporary list to later convert to triple (i,nucl,j)
    if BP.has_key(k):
      O.write("%d\t%s" % (k,PDB[k]))
      tmpList = [k,BP[k][0][0]]
      for j in BP[k]:
        if VERBOSE:
          O.write("\t%s %s %s" % (j[2],j[3],j[4]))
#          print k,BP[k][0][0],j
        else:
          O.write("\t%s" % j[2])
          tmpList.append(j[2])
      O.write("\n")
      c += 1
#      del BP[k]
    else:
      O.write("%d\t%s\t0\n" % (k,PDB[k]))
      tmpList = [k,PDB[k],0]
    OutputList.append( tuple(tmpList) )
  O.close()
  #~ print PDB
  #~ print BP
  #~ assert c==len(BP)
  return OutputList

def parseOutput(filename,allowedBP):
  if not os.path.exists(filename):
    print "rnaview output not created"
    sys.exit(1)
  I = map(lambda x:x[:-1], open(filename,"r").readlines())
  D = {}
  start = 0
  for i in I:
    if i[:15] == "BEGIN_base-pair":
      start = 1
    elif i[:13] == "END_base-pair":
      start = 0
    else:
      if start:
        items = i.split()      
        if ("!" not in i) and ("stacked" not in i):
          chainA    = items[1][0].strip()  # chain ID
          if chainA == ":":
            chainA = "" #Clote 11/27/06
          posA      = int(items[2])  # sequence progressive id
          nucA      = items[3].split("-")[0].upper()  # nucleotide
          chainB    = items[5][0].strip()          
          if chainB == ":":
            chainB = "" #Clote 11/27/06
          posB      = int(items[4])
          nucB      = items[3].split("-")[1].upper()
          bpLeontis = items[6] 
          bpConform = items[7]
          if chainA == chainB and posA != posB: 
            if BPconvert.get((bpLeontis,bpConform),"NA") in allowedBP:
              D.setdefault(chainA,{})
              D[chainA].setdefault(posA,[])
              D[chainA].setdefault(posB,[])
              if (nucA,nucB,posB,bpLeontis,bpConform) not in D[chainA][posA]:
                D[chainA][posA].append((nucA,nucB,posB,bpLeontis,bpConform))
              if (nucB,nucA,posA,bpLeontis,bpConform) not in D[chainA][posB]:
                D[chainA][posB].append((nucB,nucA,posA,bpLeontis,bpConform))
#            else:
#              print bpLeontis,bpConform,BPconvert.get((bpLeontis,bpConform),"NA")
  return D


#parse SEQRES records for nucleotide sequence. POINT: Data obtained from SEQRES since ATOM may have gaps
def getRNAfromSEQRES(filename):
  f = open(filename)
  line = f.readline()
  RNAnuclSeqDict = {}  #RNAnuclSeqDict[chainID] = RNA sequence in that chain
  while line:
    if line[:6] == 'SEQRES':
      chainID = line[11]
      if chainID not in RNAnuclSeqDict.keys():
        RNAnuclSeqDict[chainID] = []
      words   = line[19:70].split()
      for w in words:
		  if len(w)==3:
			  w=w[2]
		  RNAnuclSeqDict[chainID]  += w 
    elif line[:4] in ['ATOM','HETA']:
      break #exit this while loop
    line = f.readline()
  rnaSeqDict = {}
  for chainID in RNAnuclSeqDict.keys():
    rnaSeqDict[chainID] = string.join(RNAnuclSeqDict[chainID],"")
    #print("chain,RNA:%s\t%s\n"%(chainID,rnaSeqDict[chainID]))
    #now rnaSeqDict[chainID] is a string of the RNA sequence
  f.close()
  return rnaSeqDict
  
def parsePDB(filename):
  rnaSeqDict = getRNAfromSEQRES(filename)
  file = open(filename)
  line = file.readline()
  OffsetDict = {}
  #now parse ATOM to get offset of first nucleotide
  while line:
    #parse ATOM records for nucleotide and position
    if line[:4] in ['ATOM','HETA']:
      chainID = line[21]
      if chainID not in OffsetDict.keys():
        resSeqNum = int(line[22:26].strip())
        OffsetDict[chainID]  = resSeqNum
    elif line[:6]=="ENDMDL": #only the 1st model is condidered if more than one model exists
			break
    line = file.readline()
    
  file.close()
  #for key in OffsetDict.keys():
	#print ("offset of first nucleotide:\t%s,%d\n"%(key,OffsetDict[key]))
  #Now return dictionary D[chainID][resSeqNum]=nucl
  D = {}
  for chainID in rnaSeqDict.keys():
    D[chainID] = {}
    rna        = rnaSeqDict[chainID]
    for i in range(len(rna)):
      D[chainID][i+OffsetDict[chainID]] = rna[i]
  return(D)
  
def runRnaview(filename):
  cmd = "rnaview  %s > /dev/null" % (filename)
  os.popen(cmd)
  #subprocess.check_call(cmd,shell=True)
  if os.path.exists("%s.out" % filename):
    outname = "%s.out" % filename
  elif os.path.exists("%s_nmr.pdb.out" % filename):
    outname = "%s_nmr.pdb.out" % filename
  return outname
  
  
def RNAview(LIST,allowedBP):
  if not os.path.exists(LIST):
    print "File %s not found" % LIST
    sys.exit(1) 
  L = map(lambda x:x[:-1], open(LIST,"r").readlines())
  for l in L:
    if os.path.exists(l):
      print l
      #print "running rnaview"
      outname = runRnaview(l)
      #print "parsing pdb file"
      PDB = parsePDB(l)
      #print "parsing rnaview output"
      BP = parseOutput(outname,allowedBP)
      #print "writing output"
      for k in BP.keys():  ## k is chain identifier
        writeBpseqFile(BP[k],PDB[k],l,k)
      if DELETE_OUT:
        os.remove(outname)
    else:
      print "File %s not found" % l
  return

#~ if __name__ == "__main__":
  #~ if len(sys.argv) < 3:
    #~ usage()
    #~ sys.exit(1)
#~ 
  #~ for i in range(len(sys.argv)):	
    #~ if sys.argv[i] == "-a":	
      #~ LIST = sys.argv[i+1]
    #~ elif sys.argv[i] == "-v":
      #~ VERBOSE = 1
    #~ elif sys.argv[i] == "-MM":
      #~ allowedBP.remove("MM") 
    #~ elif sys.argv[i] == "-PP":
      #~ allowedBP.remove("PP")
    #~ elif sys.argv[i] == "-WWc":
      #~ allowedBP.append("WWc")
    #~ elif sys.argv[i] == "-WWt":
      #~ allowedBP.append("WWt")
    #~ elif sys.argv[i] == "-WSc":
      #~ allowedBP.append("WSc")
    #~ elif sys.argv[i] == "-WSt":
      #~ allowedBP.append("WSt")
    #~ elif sys.argv[i] == "-WHc":
      #~ allowedBP.append("WHc")
    #~ elif sys.argv[i] == "-WHt":
      #~ allowedBP.append("WHt")
    #~ elif sys.argv[i] == "-SHc":
      #~ allowedBP.append("SHc")
    #~ elif sys.argv[i] == "-SHt":
      #~ allowedBP.append("SHt")
    #~ elif sys.argv[i] == "-SSc":
      #~ allowedBP.append("SSc")
    #~ elif sys.argv[i] == "-SSt":
      #~ allowedBP.append("SSt")
    #~ elif sys.argv[i] == "-HHc":
      #~ allowedBP.append("HHc")
    #~ elif sys.argv[i] == "-HHt":
      #~ allowedBP.append("HHt")
#~ 
  #~ if not LIST:
    #~ usage()
    #~ sys.exit(1)
  #~ RNAview(LIST,allowedBP,VERBOSE)
  #~ L = map(lambda x:x[:-1], open(LIST,"r").readlines())
  #~ for l in L:
    #~ if os.path.exists(l):
      #~ print l
      #~ print "running rnaview"
      #~ outname = runRnaview(l)
      #~ print "parsing pdb file",outname
      #~ PDB = parsePDB(l)
      #~ print "parsing rnaview output"
      #~ BP = parseOutput(outname)
      #~ print "writing output" ,BP
      #~ for k in BP.keys():  ## k is chain identifier
        #~ writeBpseqFile(BP[k],PDB[k],l,k)
        #~ print BP
        #~ a=PDB[k].keys()
        #~ a.sort()
        #~ print PDB[k]
      #~ if DELETE_OUT:
        #~ os.remove(outname)
    #~ else:
      #~ print "File %s not found" % l
