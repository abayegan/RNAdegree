#! /usr/bin/env python

# bruteForceComputeFandG.py
# P.Clote, Nov 22, 2014

# Computes F(1,n,ch,x) = num str having x visible occurrences of ch in [1,n]
# and G(1,n,ch,x)      = num str having x visible occurrences of ch in [1,n-4]

import sys,math,copy,os
from misc import basePair,basePairList

PRINT         = 1
VERBOSE_PRINT = 0
THETA         = 3
NUCL          = ['A','C','G','U']

def numVisible(rna,secStr):
  #secStr in dot bracket notation
  assert( len(rna)==len(secStr) ),"rna and secStr have different sizes"
  n      = len(rna)
  rna    = '$'+rna
  secStr = '$'+secStr
  #----------------------------------------------------
  #Initializations
  a = 0; c = 0; g = 0; u = 0
  NumVis       = {}  #NumVis['A'] is num visible A's in secStr in [1,n]
  NumVisPrefix = {}  #NumVisPrefix['A'] is num visible A's in secStr in [1,n-4]
  for ch in NUCL: 
    NumVis[ch]=0 
    NumVisPrefix[ch]=0 
  L = []; stack = [] #L is set of base pairs in secStr, 1-indexed
  #----------------------------------------------------
  #Recursion
  for i in range(1,n+1):  #i in [1,n] is 1-index in secStr string
    if len(stack)==0 and secStr[i]!='(': 
      NumVis[rna[i]] += 1
      if i<=n-4: NumVisPrefix[rna[i]] += 1
    elif secStr[i]=="(":
      stack.append(i)
    elif secStr[i]==")":
      x = stack.pop() #return and pop end of stack
      L.append( (x,i) )
  return NumVis,NumVisPrefix

def checkFandG(rna,F,Z):
  #check that sum of F[A][x] taken over x is Z[1,n]
  n = len(rna)
  for ch in NUCL:
    sum = 0
    for x in F[ch].keys():
      sum += F[ch][x]
    pf = Z[(1,n)]
    assert (sum==pf), "sum_x F[%s][x]=%s != Z[1,%s]=%s" % (ch,sum,n,pf)
 

def bruteForceComputeFandG(rna):
  F = {}; G = {}; n=len(rna); totalNumStr = 0
  #F['A'][x] is number str having x visible occ of A in [1,n], etc.
  #G['A'][x] is number str having x visible occ of A in [1,n-4], etc.
  for ch in NUCL: 
    F[ch] = {}; G[ch] = {}  
  for ch in NUCL:
    for x in range(n+1): 
      F[ch][x] = 0; G[ch][x] = 0 
  cmd = "echo %s | RNAsubopt -e 100" % rna
  file = os.popen(cmd)
  line = file.readline() #discard first line
  line = file.readline()
  while line:
    totalNumStr        += 1
    secStr              = line.split()[0]
    NumVis,NumVisPrefix = numVisible(rna,secStr)
    for ch in NUCL: 
      F[ch][NumVis[ch]]       += 1
      G[ch][NumVisPrefix[ch]] += 1
    if VERBOSE_PRINT:
      #print "\n",rna
      print secStr,NumVis,NumVisPrefix
    line = file.readline()
  return F,G



     
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


def numShifts(rna,secStr): #brute force
  #rna and secStr is 1-indexed
  if VERBOSE_PRINT:
    print "Sec str: %s" % secStr
  n      = len(rna)
  rna0   = rna
  rna    = '$'+rna
  secStr = '$'+secStr
  SS   = basePairList(secStr)
  num  = 0; tempSS = copy.deepcopy(SS)
  for (i,j) in SS:
    tempSS.remove( (i,j) )
    if VERBOSE_PRINT:
      print "Remove (%d,%d) from %s" % (i,j,secStr[1:])
    for x in range(1,n+1):
      if abs(i-x)>THETA and x!=j:
        if i<x and basePair(rna[i],rna[x]):
          tempSS.append( (i,x) )
          if isSecStr(tempSS): 
            num += 1
            if VERBOSE_PRINT:
              ss = basePairList2dotBracketNotation(rna0,tempSS)
              print '(%d,%d)\t%s' % (i,x,ss)
          tempSS.remove( (i,x) )
        elif basePair(rna[x],rna[i]): #i>x
          tempSS.append( (x,i) )
          if isSecStr(tempSS): 
            num += 1
            if VERBOSE_PRINT:
              ss = basePairList2dotBracketNotation(rna0,tempSS)
              print '(%d,%d)\t%s' % (x,i,ss)
          tempSS.remove( (x,i) )
      elif abs(j-x)>THETA and x!= i:
        if j<x and basePair(rna[j],rna[x]):
          tempSS.append( (j,x) )
          if isSecStr(tempSS): 
            num += 1
            if VERBOSE_PRINT:
              ss = basePairList2dotBracketNotation(rna0,tempSS)
              print '(%d,%d)\t%s' % (j,x,ss)
          tempSS.remove( (j,x) )
        elif basePair(rna[x],rna[j]): #j>x
          tempSS.append( (x,j) )
          if isSecStr(tempSS): 
            num += 1
            if VERBOSE_PRINT:
              ss = basePairList2dotBracketNotation(rna0,tempSS)
              print '(%d,%d)\t%s' % (x,j,ss)
          tempSS.remove( (x,j) )
    tempSS.append( (i,j) ) #put back the base pair temporarily removed
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
  n   = len(rna)
  F,G = bruteForceComputeFandG(rna)
  Z   = partFunInner(rna)
  checkFandG(rna,F,Z)
  if PRINT:
    print "F(1,n,ch,x) -- unshown values are zero"
    for ch in NUCL:
      keys = F[ch].keys()
      keys.sort
      for x in keys:
        if F[ch][x]>0: print "F[%s][%s]=%s" % (ch,x,F[ch][x])


if __name__ == '__main__':
  if len(sys.argv)!= 2:
    print "Usage: %s RNA" % sys.argv[0]
    sys.exit(1)
  rna = sys.argv[1].upper()
  main(rna)


