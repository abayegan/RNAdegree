#! /usr/bin/env python

# totalNumNborsAddRemovalShiftUnifProb.py 
# P.Clote, Nov 23, 2014

# Computes the uniform expected number of neighbors for an RNA sequence
# a1...an (not a homopolymer).
# This corresponds to setting all energy terms to zero.

# WARNING: There is a constant in misc.py called HOMOPOLYMER
# If HOMOPOLYMER = 1 then we treat the homopolymer model, otherwise
# we treat real RNA. In any case, the energy is zero.

import sys,math,os
from misc import basePair
#from computeFwithDP import computeE,computeEprime
from computeFwithDP import computeF,computeH,computeEleft,computeEright
from computeFwithDP import computeFF,computeJ,computeEleft,computeEright
from computeFwithDP import computeERprime

PRINT           = 0
PRINT_RNASUBOPT = 0
PRINT_OUTPUT    = 0
THETA           = 3
NUCL            = ['A','C','G','U']

XNUCL = ['A','C','G','U','X'] #extended nucleotides
CHAR2INT = {'A':0,'C':1,'G':2,'U':3,'X':4}
INT2CHAR = {0:'A',1:'C',2:'G',3:'U',4:'X'}


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

def int2char(x):
  #x in [0,624]=[0,(5**4)-1] and decode(x) is [c1,c2,c3,c4]
  #where ci in NUCL+['X']
  L = []
  for i in range(4):
    L.append( INT2CHAR[x%5] )
    x = x/5
  L.reverse()
  return L    

def char2int(L):
  #L is [c1,c2,c3,c4] where ci in NUCL+['0']
  sum = 0
  for ch in L:
    sum = 5*sum+CHAR2INT[ch]
  return sum


def computeEold(rna,Z): 
  #function E(i,j) = sum_S sum_{(x,y)} I[ (x,y) external in S] 
  E   = {}
  n   = len(rna)
  rna = '$'+rna
  #initialize
  for i in range(1,n+1):
    for j in range(1,n+1):
      E[(i,j)] = 0L
  #fill
  for d in range(THETA+1,n): #d in [4,n-1]
    for i in range(1,n-THETA): #i in [1,n-THETA-1]
      j = i+d
      if j>n: break
      sum = E[(i,j-1)]
      if basePair(rna[i],rna[j]):
        sum += Z[(i+1,j-1)] 
        #Tricky: contribution comes from 3rd term in derivation of E
        #Contribution from 2nd term in derivation of E is ZERO
      for k in range(i+1,j-THETA):# k in [i+1..j-THETA-1]
        if basePair(rna[k],rna[j]):
          sum += (E[(i,k-1)]+Z[(i,k-1)])*Z[(k+1,j-1)]
      E[(i,j)] = sum
  return E


def computeQ(rna,Z):
  FF      = computeFF(rna,Z) 
  J       = computeJ(rna,Z,FF)
  EL      = computeEleft(rna,Z)
  ER      = computeEright(rna,Z)
  ERprime = computeERprime(rna,Z,ER)
  #--------------------------------------------
  #E[(i,j) = sum_S sum_{(x,y)} I[(x,y) external in S]
  # 
  #Ebis[(i,j) = sum_S sum_{(x,y)} I[(x,y) external in S, y<n]
  #
  #Eprime[(i,j) = sum_S sum_{(x,y)} 
  #        I[(x,y) external in S, y <= j-4, j unpaired]
  #EL[(i,j) = sum_S sum_{(x,y)} 
  #        I[(x,y) external in S, x basepairs with rna[n]]
  #ER[(i,j) = sum_S sum_{(x,y)} 
  #        I[(x,y) external in S, y basepairs with rna[n]]
  #ERprime[(i,j) = sum_S sum_{(x,y)} 
  #  I[(x,y) external in S, y basepairs with rna[n], y<=j-4, n unpaired in S]
  Q = {}
  n   = len(rna)
  rna = '$'+rna
  #-------------------------------------------
  #Initialization to zero
  for d in range(n):         #d in [0,n-1]
    for i in range(1,n+1-d): #i in [1,n-d]
      j = i+d
      Q[(i,j)]=0
  #-------------------------------------------
  #Recursions
  for d in range(THETA+1,n):         #d in [4,n-1]
    for i in range(1,n+1-d):         #i in [1,n-d]
      j = i+d
      #Case 1: first term Q(i,j-1)
      Q[(i,j)] = Q[(i,j-1)]
      if PRINT:
        print "\n\n",Q[(i,j-1)], "\tQ(%d,%d)" % (i,j-1)
      #Case 2: second term 2 * sum_k z(i,k-1)*z(k+1,j-1)
      sum = 0
      if basePair(rna[i],rna[j]):
        sum += Z[(i+1,j-1)]
      for k in range(i+1,j-3): #k in [i+1,j-4]
        if basePair(rna[k],rna[j]):
          sum += Z[(i,k-1)]*Z[(k+1,j-1)]
      Q[(i,j)] += 2*sum
      if PRINT:
        text = "\t2 sum_{k=%s}^{%s-4} z(%s,k-1)*z(k+1,%s-1)" 
        print (2*sum),text % (i,j,i,j)
      #Case 3: 2*EL[(i,j-1)]+2*ERprime[(i,j)]
      sum = 2*EL[(i,j-1,rna[j])]+2*ERprime[(i,j,rna[j])]
      Q[(i,j)] += sum
      if PRINT:
        text1 = "%s\t" % sum
        text2 = "2*EL[(%s,%s-1,%s)]+2*ERprime[(%s,%s,%s)]" 
        text2 = text2 % (i,j,rna[j],i,j,rna[j])
        print text1,text2
        text1 = " %s\t" % 2*EL[(i,j-1,rna[j])]
        text2 = "2*EL[(%s,%s-1,%s)]" % (i,j,rna[j])
        print text1,text2
        text1 = " %s\t" % 2*ERprime[(i,j,rna[j])]
        text2 = "2*ERprime[(%s,%s,%s)]" % (i,j,rna[j])
        print text1,text2
      #Case 4: fourth term is sum_{x=2}^{n-4} x*(x-1)*H(1,n,ch,x)
      sum = 0
      for x in range(2,j-i+1-3): #x in [2,j-i+1-4]
        sum += x*(x-1)*J[(i,j,rna[j],x)]
      Q[(i,j)] += sum 
      newsum = 0
      if PRINT:
        text = "\tch\tx\tx-1\tJ(%s,%s,ch,x)\tsum_ch sum_x x(x-1)J(%s,%s,ch,x)"
        print sum,"\tch\tx\tx-1\tJ(1,n,ch,x)\tsum_ch sum_x x(x-1)J(1,n,ch,x)"
        ch = rna[j]
        for x in range(2,j-i+1-3): #x in [2,j-i+1-4]
          print "\t%s\t%s\t%s\t%s = \t%s"%(ch,x,x-1,J[(1,n,ch,x)],x*(x-1)*J[(1,n,ch,x)])
      #Case 5: fifth term  
      #sum_{k=1}^{n-\theta-1} \left( z(k-1) \cdot Q(n-k-1) \right) +
      #\left( Q(k-1) \cdot z(n-k-1) \right)
      #Case 5a: (i,j) paired
      sum = 0
      if basePair(rna[i],rna[j]):
        sum += Q[(i+1,j-1)] 
      for k in range(i+1,j-3): #k in [i+1,j-4]
        if basePair(rna[k],rna[j]):
          sum += (Z[(i,k-1)]*Q[(k+1,j-1)])+(Q[(i,k-1)]*Z[(k+1,j-1)])
      Q[(i,j)] += sum
      if PRINT:
        text = "sum_k (Z[(%s,k-1)]*Q[(k+1,%s-1)])+(Q[(%s,k-1)]*Z[(k+1,%s-1)])"
        text = text % (i,j,i,j)
        print "%s\t%s" % (sum,text)
  return Q         
       

def printMatrix(rna,M):
  n    = len(rna)
  keys = M.keys()
  if len(keys) == n**2:
    m = n
  else:
    m = 2*n
  for i in range(1,n+1):
    for j in range(1,m+1):
      sys.stdout.write("%s " % M[(i,j)])
    sys.stdout.write("\n")

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
  n        = len(rna)
  Z        = partFunInner(rna)
#  print "%s of length %s with total num str of %s" % (rna,n,Z[(1,n)])
  EL       = computeEleft(rna,Z)
  ER       = computeEright(rna,Z)
  ERprime  = computeERprime(rna,Z,ER)
  FF       = computeFF(rna,Z)
  J        = computeJ(rna,Z,FF)
  Q  = computeQ(rna,Z)
  n  = len(rna)
  numberNborsNoShift   = numNborsNoShift(rna,Z)
  numberNborsWithShift = Q[(1,n)]
  numStr               = Z[(1,n)]
  expNumNborsWithoutShift = numberNborsNoShift/float(numStr)
  expNumNborsWithShift    = numberNborsWithShift/float(numStr)
  print "rna\t%s" % rna
  print "Z\t%s" % Z[(1,n)]
  print "Q\t%s" % Q[(1,n)]
  print numberNborsNoShift, numberNborsWithShift, numStr, expNumNborsWithoutShift, expNumNborsWithShift
  if PRINT_OUTPUT: 
    for i in range(1,n+1):
      for j in range(i,n+1):
        print "%s\t%s\t%s" % (i-1,j-1,Q[(i,j)])
    for c in NUCL:
      print("Char is %c" % c);
      print "\tEL\t%s" % EL[(1,n,c)]
      print "\tER\t%s" % ER[(1,n,c)]
      print "\tER1\t%s" % ERprime[(1,n,c)]
      for x in range(4):
        print "\tG(0,%s,%s,%s)=%s" % (n-1,c,x,J[(1,n,c,x)])
    print "Num nbors without shift: ", numberNborsNoShift
    print "Num nbors with shift: ",    (numberNborsWithShift-numberNborsNoShift)
    print "Num nbors: %s" % numberNborsWithShift
    print "Num structures: ", numStr
    print "Exp num nbors without shift: ", expNumNborsWithoutShift
    print "Exp num nbors with shift: ",    expNumNborsWithShift


if __name__ == '__main__':
  if len(sys.argv)!= 2:
    print "Usage: %s RNA" % sys.argv[0]
    sys.exit(1)
  rna = sys.argv[1].upper()
  main(rna)


