#! /usr/bin/env python

# bruteForceNumNborsWithShift.py
# P.Clote

# Computes the uniform expected number of neighbors for an RNA sequence
# a1...an (not a homopolymer). Method exhaustion, using RNAsubopt
# This corresponds to setting all energy terms to zero.

import sys,math,copy,os
from misc import basePair,basePairList

PRINT_OUTPUT                   = 1
PRINT                          = 0
VERBOSE_PRINT                  = 0
PRINTbasepairadditionsremovals = 0
THETA                          = 3

def basePairList2dotBracketNotation(rna,L):
  #WARNING: 1-indexed
  n = len(rna)
  secStrList = '.'*(n+1)
  secStrList = list(secStrList)
  for (i,j) in L:
    secStrList[i] = '('
    secStrList[j] = ')'
  secStr = "".join(secStrList[1:])
  return secStr

def isSecStr(BasePairList):
  numBasePairs = len(BasePairList)
  for k in range(numBasePairs):
    (i,j) = BasePairList[k]
    for l in range(k+1,numBasePairs):
      (x,y) = BasePairList[l]
      if (i,j)!=(x,y):
        if (i in [x,y] or j in [x,y]): return False
        if (i<x<j<y or x<i<y<j): return False
  return True

def numBasePairAdditionsRemovals(rna,secStr): #brute force
  #rna and secStr is 1-indexed
  if PRINTbasepairadditionsremovals: 
    print "Sec str: %s" % secStr
  n      = len(rna)
  rna0   = rna
  rna    = '$'+rna
  secStr = '$'+secStr
  SS     = basePairList(secStr)
  num    = len(SS); tempSS = copy.deepcopy(SS)
  for x in range(1,n-THETA):
    for y in range(x+THETA+1,n+1):
      if basePair(rna[x],rna[y]):
        if (x,y) not in tempSS:
          tempSS.append( (x,y) )
          if isSecStr(tempSS): 
            num += 1
            if PRINTbasepairadditionsremovals: 
              ss = basePairList2dotBracketNotation(rna0,tempSS)
              print '(%d,%d)\t%s' % (x,y,ss)
          tempSS.remove( (x,y) )
  return num


def numShifts(rna,secStr): #brute force computation of number shifts
  #rna and secStr is 1-indexed
  if VERBOSE_PRINT:
    print "Sec str: %s" % secStr
  n      = len(rna)
  rna0   = rna
  secStr0= secStr
  rna    = '$'+rna
  secStr = '$'+secStr
  SS   = basePairList(secStr)
  num  = 0; tempSS = copy.deepcopy(SS)
  for (i,j) in SS:
    tempSS.remove( (i,j) )
    if VERBOSE_PRINT:
      print "Remove (%d,%d) from %s" % (i,j,secStr[1:])
    for x in range(1,n+1):
      if abs(i-x)>THETA and x!=j and basePair(rna[i],rna[x]):
        if i<x:
          tempSS.append( (i,x) )
          if isSecStr(tempSS): 
            num += 1
            if VERBOSE_PRINT:
              ss = basePairList2dotBracketNotation(rna0,tempSS)
              print '(%d,%d)\t%s' % (i,x,ss)
          tempSS.remove( (i,x) )
        else: #i>x
          tempSS.append( (x,i) )
          if isSecStr(tempSS): 
            num += 1
            if VERBOSE_PRINT:
              ss = basePairList2dotBracketNotation(rna0,tempSS)
              print '(%d,%d)\t%s' % (x,i,ss)
          tempSS.remove( (x,i) )
      elif abs(j-x)>THETA and x!= i and basePair(rna[j],rna[x]):
        if j<x:
          tempSS.append( (j,x) )
          if isSecStr(tempSS): 
            num += 1
            if VERBOSE_PRINT:
              ss = basePairList2dotBracketNotation(rna0,tempSS)
              print '(%d,%d)\t%s' % (j,x,ss)
          tempSS.remove( (j,x) )
        else: #j>x
          tempSS.append( (x,j) )
          if isSecStr(tempSS): 
            num += 1
            if VERBOSE_PRINT:
              ss = basePairList2dotBracketNotation(rna0,tempSS)
              print '(%d,%d)\t%s' % (x,j,ss)
          tempSS.remove( (x,j) )
    tempSS.append( (i,j) ) #put back the base pair temporarily removed
  if PRINT: 
    print "%s has %s shifts" % (secStr[1:],num) #recall '$' prepended to secStr
  return num

def partFunInner(rna):  
  #Z[(i,j)] is number of sec str of rna[i...j], using array n x n
  Z   = {}
  n   = len(rna)
  rna = '$'+rna
  #initialize
  for i in range(1,n+1):
    for j in range(1,n+1):
      Z[(i,j)] = 1L
  #fill
  for d in range(THETA+1,n): #d in [4,n-1]
    for i in range(1,n-THETA): #i in [1,n-THETA-1]
      j = i+d
      if j>n: break
      sum = Z[(i,j-1)]
      if (j-i>THETA):
        if basePair(rna[i],rna[j]):
          sum += Z[(i+1,j-1)]
        for k in range(i+1,j-THETA):# k in [i+1..j-THETA-1]
          if basePair(rna[k],rna[j]):
            sum += Z[(i,k-1)]*Z[(k+1,j-1)]
        Z[(i,j)] = sum
  return Z

def numNborsNoShift(rna,Z):  
  #Q[(i,j)] = sum_x N(x), where x is sec str of rna[i:j+1] without shift
  #Z[(i,j)] = number of sec str of rna[i:j+1]
  Q = {}; 
  n   = len(rna)
  rna = '$'+rna #add dummy character
  #base case j-i in [0,1,2,3]
  for i in range(1,n+1):
    for d in range(0,min(THETA,n-i)+1):
      j = i+d
      Q[(i,j)]=0.0
  #inductive case j-i > 3, compute Q[(i,j)] all 1<=i<=j<=n
  for d in range(THETA+1,n+1):
    for i in range(1,n+1):
      j = i+d
      if (j>n): break
      q = Q[(i,j-1)] 
      if basePair(rna[i],rna[j]):
        q += 2*Z[(i+1,j-1)]+Q[(i+1,j-1)]
      for k in range(i+1,j-THETA):
        if basePair(rna[k],rna[j]):
          q += 2*Z[(i,k-1)]*Z[(k+1,j-1)]+Z[(i,k-1)]*Q[(k+1,j-1)]+Q[(i,k-1)]*Z[(k+1,j-1)]
      Q[(i,j)] = q
  return Q[(1,n)]


def main(rna):
  numberNborsWithShift         = 0
  numberNborsNoShiftBruteForce = 0
  cmd  = 'echo %s | RNAsubopt -e 100' % rna
  file = os.popen(cmd)
  line = file.readline() #discard first line
  line = file.readline() #line contains secStr \t Energy
  while line:
    secStr = line.split()[0]
    shifts = numShifts(rna,secStr)
    addRem = numBasePairAdditionsRemovals(rna,secStr)
    numberNborsWithShift += shifts
    numberNborsNoShiftBruteForce += addRem
    if PRINT_OUTPUT:
      print "%s\t%s\t%s\t%s" % (secStr,shifts,addRem,(shifts+addRem))
    line   = file.readline() 
  n        = len(rna) 
  Z        = partFunInner(rna)
  print "%s of length %s with total num str of %s" % (rna,n,Z[(1,n)])
  numberNborsNoShift   = numNborsNoShift(rna,Z)
  n                    = len(rna)
  numStr               = Z[(1,n)]
  expNumNborsWithoutShift = numberNborsNoShift/float(numStr)
  expNumNborsWithShift    = numberNborsWithShift/float(numStr)
  print "Brute force enumeration results "
  print "Num nbors without shift: ", numberNborsNoShift
  print "Num nbors without shift (brute force): ", numberNborsNoShiftBruteForce
  print "Num nbors with shift: ",    numberNborsWithShift
  print "Num nbors: ", ( numberNborsNoShift+ numberNborsWithShift )
  print "Num structures: ", numStr
#  print "Exp num nbors without shift: ", expNumNborsWithoutShift
#  print "Exp num nbors with shift: ",    expNumNborsWithShift


if __name__ == '__main__':
  if len(sys.argv)!= 2:
    print "Usage: %s RNA" % sys.argv[0]
    sys.exit(1)
  rna = sys.argv[1].upper()
  main(rna)


