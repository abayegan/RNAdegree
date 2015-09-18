#! /usr/bin/env python

#A. Bayegan
#bruteForce computation of the uniform expected number of neighbors for MS1 and MS2

import sys,os,shlex,re,subprocess,math
from itertools import groupby

#find shift moves
def getNumShifts(s):
	bpList , baseList = getBpList(s)
	shiftList = []
	for l in bpList:
			x=l[0];y=l[1]
			for k in range(len(s)):
				if k not in baseList and k!=x and k!=y:
					if k < x and bp(k,x):
						if isCompatible(k,x,bpList):
							shiftList.append((k,x))
					if k > x and bp(x,k):
						if isCompatible(x,k,bpList):
							shiftList.append((x,k))
					if k < y and bp(k,y):
						if isCompatible(k,y,bpList):
							shiftList.append((k,y))
					if k > y and bp(y,k):
						if isCompatible(y,k,bpList):
							shiftList.append((y,k))
	#print "shiftList:",shiftList
	return shiftList 
			
#find compatible basepairs can be added to structure s 	
def getNumBpAdd(s):
	addList = []
	bpList , baseList = getBpList(s)
	#print bpList,baseList
	for i in range(len(s)):
		if i not in baseList:
			for j in range(i+4, len(s)):
				if j not in baseList and bp(i,j):
					if isCompatible(i,j,bpList):
						addList.append((i,j))
	#print "addList:",addList
	return addList
	
#get the total number of basepairs that can be deleted from a structure
def getNumBpDel(s):
	#print "delCnt=%s" %s.count('(')
	return s.count('(')

#check if structure L remains valid after adding basepair(i,j)	
def isCompatible(i,j,L):
	for l in L:
		x=l[0];y=l[1]
		if (i < x < j < y or
			x < i < y < j):
			return 0
	return 1				

def getBoltzman(e):
	return math.exp(-e/RT)
	
#compute partition function from the structure list containing energies
def PartitionFunction(strucList):
	z = 0
	for l in strucList[1:]:
		E = float(l.split()[1])
		#print "E:" , E
		z += math.exp(-E/RT)
	#print "z= " , z
	return z
					
#check if i and j could form a base pair				
def bp(i,j):
	if(j-i<4):
	    return 0;
	elif (SEQ[i]=='A' and SEQ[j]=='U')or(SEQ[i]=='U' and SEQ[j]=='A'):
	    return 1;
	elif (SEQ[i]=='C' and SEQ[j]=='G')or(SEQ[i]=='G' and SEQ[j]=='C'):
	    return 1;
	elif (SEQ[i]=='U' and SEQ[j]=='G')or(SEQ[i]=='G' and SEQ[j]=='U'):
	    return 1;
	else:
		return 0;

#given a structure get all the basepairs
def getBpList(s):
  BpList = []; q = []; BaseList = []
  for i in range(len(s)):
    if s[i]=="(":
      q.append(i)
    elif s[i]==")":
      p = q.pop()
      BpList.append((p,i))
      BaseList.append(p); BaseList.append(i);
  #print BpList[::-1], pBaseList[::-1]
  return BpList[::-1],BaseList[::-1]
				
def main():
	#get all possible structures
	NumNborNoShift = 0; NumNborWithShift = 0; expNoShift=0;expWithShift=0;expNoshift=0;totNumNborNoShift=0;totNumNborWithShift=0
	noShiftStrList = []; shiftStrList = []
	cmdEcho = "echo %s" % SEQ
	#cmdSubopt = "RNAsubopt -e 100 -d 0" #Turner 99
	cmdSubopt = "/home/amir/Vienna/ViennaRNA-2.1.7/Progs/RNAsubopt -e 100 -d 0"
	echo = subprocess.Popen(shlex.split(cmdEcho), stdout=subprocess.PIPE)
	subopt = subprocess.Popen(shlex.split(cmdSubopt),stdin=echo.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = subopt.communicate()
	if err != '':
		print repr(err)
		sys.exit(1)
	strucList = out.strip().split('\n')
	z = PartitionFunction(strucList)
	#for every possible structure get the number of possible bp additions, deletions and shifts
	#f = open(fName,'w')
	for st in strucList[1:]:
		p = getBoltzman(float(st.split()[1]))
		#print st
		#print"BF=%.6f" %p
		addingBpList = getNumBpAdd(st.split()[0])	
		remBpNum = getNumBpDel(st.split()[0])
		shiftBpList = getNumShifts(st.split()[0]);
		NumNborNoShift = len(addingBpList) + remBpNum
		NumNborWithShift = len(addingBpList)+ remBpNum + len(shiftBpList)
		expWithShift += p*NumNborWithShift
		expNoShift += p*NumNborNoShift
		#totExpNumNborNoShift += expNoShift
		#totExpNumNborWithShift += expWithShift
		totNumNborNoShift +=  NumNborNoShift
		totNumNborWithShift += NumNborWithShift
		#noShiftStrList.append(NumNborNoShift)
		#shiftStrList.append(NumNborWithShift)
		
		#print number of nbors of every structure
		print("%s\t%f"%(st,NumNborWithShift))
	print ("number of structures:%d"% (len(strucList)-1))
	print ("Q(shift-Turner): %f\tQ(shift-uniform):%d\nQ(no shift-Turner):%f\tQ(no shift-uniform):%d\npf:%f\n" %(expWithShift,totNumNborWithShift,expNoShift,totNumNborNoShift,z))
	print ("expNumNbor(shift):%f\texpNumNbor(no shift):%f" % (expWithShift/z,expNoShift/z))
	#print "adding base pair list:%s\nshifting base pairs list %s \n" % (ad,sh)
	#~ print "Num structures: %s "% (len(strucList)-1)
	#~ print "Total Num nbors without shift: %s " % totNumNborNoShift
	#~ print "Total Num nbors with shift: %s" % totNumNborWithShift
	#~ print "Unifrom Exp num nbors without shift: %lf" % (float(totNumNborNoShift)/float(len(strucList)-1))
	#~ print "Uniform Exp num nbors with shift:%lf" % (float(totNumNborWithShift)/float(len(strucList)-1))
	#~ print "Turner Exp num nbors without shift: %s" % expNoShift
	#~ print "Turner Exp num nbors with shift:%s" % expWithShift

	
	#get the frequency of add/del or add/del/shift for all structures
	#~ noShiftStrList.sort();shiftStrList.sort()
	#~ addFile = open('structureFreqNoShift.txt','w')
	#~ addFile.write("NumNbors without shift\tFrequency\n")
	#~ for k,g in groupby(addStrList):
		#~ addFile.write("%s\t%s\n" % (k,(len(list(g)))))
	#~ 
	#~ shiftFile = open('structureFreqWithShift.txt','w')
	#~ shiftFile.write("NumNbors with shift\t\tFrequency\n")
	#~ for k, g in groupby(shiftStrList):
		#~ shiftFile.write("%s\t%s\n" % (k,(len(list(g)))))
	#~ shiftFile.close();addFile.close()
	
		
if __name__ == '__main__':
	if len(sys.argv) != 2:
		print "Usage: %s \"Sequence\" " % sys.argv[0]
		sys.exit(1)
	global SEQ,fName
	SEQ = sys.argv[1]
	global RT
	RT = 0.00198717 * 310.15
	main()
