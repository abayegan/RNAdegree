#! /usr/bin/env python

# secStrFromRNAVIEW.py
# P.Clote

# Applies RNAVIEW by running runRnaviewPeter.py
# then applies my implementation of Yann Ponty's algorithm using a modified
# NussJac to determine that (planar) sec str which has maximal
# size within a general sec str


import sys,os,tempfile,string
from runRnaviewPeter import runRnaview
from misc import THRESHOLD
from misc import list2string,basePairList

RNAVIEW    = '/usr/src/RNAVIEW/bin/rnaview'
DEBUG = 0

def printList(L):
  for x in L: print L

 
def initialEnergyMatrix(E,n):
   #E is dictionary of form E[(i,j)] = energy of s[i],...,s[j]
   #This function initializes energy to be 0
   for i in range(n):
     for j in range(i,n):
       E[(j,i)]=-1
       E[(i,j)]=0
       #Note: E[(i,i)] = 0


def printMatrix(E,n):
  for i in range(n):
    for j in range(n):
      sys.stdout.write("%d\t" % E[(i,j)])
    sys.stdout.write("\n")




def computeMaxSecStrInGenSecStr(L,n,leftBracket='(',rightBracket=')'):
  #-------------------------------
  #  L is list of generalized base pairs (allow triples, nonnestedness)
  #  paren is parenthesis expression for maximal RNA secondary structure
  #-------------------------------

   
  #----- initialization  --#
  E = {} #empty dictionary
  initialEnergyMatrix(E,n)
  paren = [] #paren is mutable list
  for i in range(n): paren.append('.')
  computeEnergyMatrix(E,L,n)
  backtrack(0,n-1,E,paren,leftBracket,rightBracket)
  maxNumBasePairs = -E[(0,n-1)]
  if (DEBUG):
    print "Number of nucleotides = %d" % n
    printMatrix(E,n)
  return list2string(paren),maxNumBasePairs


def computeEnergyMatrix(E,L,n):
   #convert L=[(i,j),...] into D[(i-1,j-1)] = -1 or 0
   #WARNING: base pairs in L satisfy 1 <= i < j <= n
   #         base pairs in D satisfy 1 <= i < j <= n
   #         However, indices in loops range over 0,...,n-1
   #However, my NussJac original code had indices 0<=i<n
   #so an expression like "D[(i+1,j+1)] + E[(i+1,j-1)]"
   #really means the sum of NussJac energy of (i+1,j+1) base pair plus
   #the energy E of interval [i+2,j] 
   D = {}
   for i in range(1,n+1):
     for j in range(i+1,n+1):
       if (i,j) in L:
         D[(i,j)] = -1
       else:
         D[(i,j)] =  0
   for d in range(THRESHOLD+1,n):
     for i in range(n):
        #print "d=%d" % d
        j=i+d
        if j < n:
          min =0
          index=n
           #-------------------------------------
           # if index<n at end of for-loop, then this
           # means that index and j form a base pair,
           # and this is noted by E[j][i]=index.
           # if index=n at end of for-loop, then this
           # means that j is not base paired.
           #-------------------------------------

          if E[(i,j-1)]<min: 
             min = E[(i,j-1)]
             index = n
             #j not basepaired with some k such that i<k<j
          if D[(i+1,j+1)] + E[(i+1,j-1)] < min:
             min = D[(i+1,j+1)] + E[(i+1,j-1)]
             index=i

          for k in range(i+1,j-THRESHOLD):
             val = D[(k+1,j+1)] + E[(i,k-1)] + E[(k+1,j-1)]
             if val < min:
                min = val
                index=k
          E[(i,j)]=min
          if (index<n):
             E[(j,i)]=index
          else:
             E[(j,i)]=-1
  #endComputeEnergyMatrix



def backtrack(i,j,E,paren,leftBracket,rightBracket):
   if j-i>THRESHOLD:
     k = E[(j,i)]
     if k != -1:
       if DEBUG: print "Base pair k=%d,j=%d\n" % (k,j)
       paren[k] = leftBracket
       paren[j] = rightBracket
       if THRESHOLD <= (j-1)-(k+1):
         backtrack(k+1,j-1,E,paren,leftBracket,rightBracket)
       if THRESHOLD <= k-1-i: 
         backtrack(i,k-1,E,paren,leftBracket,rightBracket)
     else: #k==-1
       if DEBUG: print "k=-1"
       if THRESHOLD <= j-1-i:
          #printf "i=%d j=%d\n" % (i,j-1)
          backtrack(i,j-1,E,paren,leftBracket,rightBracket)
       else:
          return 0
  #endBacktrack

def readInputList(filename):
  L    = []
  file = open(filename)
  n = int(file.readline())  #first line contains value n
  line = file.readline()
  while line: 
    words = line.split()
    i     = int(words[0])
    j     = int(words[1])
    L.append( (i,j) )
    line = file.readline()
  file.close()
  return L,n


def runSecStr(filename):
  OutputList = []
  file    = open(filename)
  line    = file.readline()  #discard first line
  line    = file.readline()  
  while line:
    words   = line.split()
    i       = int(words[0]) 
    nucl       = words[1]
    j       = int(words[2]) 
    OutputList.append( (i,nucl,j) )
    line    = file.readline()  

  SecStrList = []; SeqList = []

  #process first triple to get offset
  offset= OutputList[0][0]-1
  #now process output list of triples
  for (i,nucl,j) in OutputList:
    i     = i-offset
    j     = j-offset
    SeqList.append(nucl)
    if i<j: SecStrList.append( (i,j) )
  rna = string.join(SeqList,"")
  n   = len(rna)
  secStr,numBP = computeMaxSecStrInGenSecStr(SecStrList,n)
  #print rna
  #print secStr
  return rna,secStr

#~ if __name__ == '__main__':
  #~ if len(sys.argv) != 2:
    #~ print "Usage: %s runRnaviewOutputFile" % sys.argv[0]
    #~ sys.exit(1)

