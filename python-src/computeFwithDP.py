#! /usr/bin/env python

# computeFwithDP.py
# P.Clote
# Nov 23, 2014. Very tricky algorithm and program. Finally correct!


#---------------------------------------------------------------------
# This program computes Q by dynamic programming, where for given rna
# Q(i,j) = sum_{s on [i,j]} N(s), where N(s) is number of base pair
# additions, removals and shifts that transform s to a valid sec str t.
# Additionally, the number of structures Z(i,j) on [i,j] is computed,
# which is the partition function when energy of a base pair equals 0.
# NOTE: Throughout, the min number THETA of unpaired bases in a hairpin
# is 3; sometimes the constant THETA is explicitly noted, and other times,
# for brevity, simply the number 3 is written. This often happens in
# limits for sums (bounds for for-loops).
# NOTE: In order to realize the HOMOPOLYMER model, simply set
# the function basePair(x,y)=1 for all x,y in NUCL=['A','C','G','U'].
# 
# Define the prefix of length n RNA sequence to be the region [1,n-4]
# and the suffix to be the region [n-3,n].
#
# Function Q(i,j) requires the following functions, where
# the expression sum_s means sum over s str on [i,j] below.
#
# 1) EL(i,j,ch) = sum_s sum_{(x,y)} 
#      I[(x,y) external in s, x basepairs with ch]
# 2) ER(i,j,ch) = sum_s sum_{(x,y)} 
#      I[(x,y) external in s, y basepairs with ch]
# 3) ERprime(i,j,ch) = sum_s sum_{(x,y)}
#      I[ s on [i,j], (x,y) external in s, y <= j-4, 
#         j unpaired in s and y basepairs with ch]
# 4) FF(i,j,ch,x),  = sum_s 
#      I[ s has EXACTLY x many visible occurrences of a 
#         a nucleotide that basepairs with ch]
# 5) J(i,j,ch,x) = sums_
#      I[s has EXACTLY x visible occurrences of a nucleotide
#        in [i,j-4] which can basepair with ch, and j unpaired in s]
#
# NOTE: In the C-program, we have the following correspondences with
# the Python program:
#      C-function		Python-function
#    	EL			EL
#	ER			ER
# 	ER1			ERprime
#	F			FF
#	G			J			
#
# The function Q is defined in totalNumNborsAddRemovalShiftUnifProb.py
# which calls the functions EL,ER,ERprime,FF,J in the current file.
# In the C-program, Q is defined in the file aux.c with the functions
# EL,ER,ER1,F,G.
#
# NOTE: Some additional precursor functions are defined in this file.
# These are E(i,j), Ebis(i,j), Eprime(i,j), F(i,j,ch,x),
# G(i,j,ch,x), H(i,j,ch,x) defined as follows.
#
# 1) E[(i,j) = sum_S sum_{(x,y)} I[(x,y) external in S]
# 2) Ebis[(i,j) = sum_S sum_{(x,y)} I[(x,y) external in S, y<n]
# 3) Eprime[(i,j) = sum_S sum_{(x,y)}
#        I[(x,y) external in S, y <= j-4, j unpaired]
# 4) F(i,j,ch,x) = sum_S
#        I[S has EXACTLY x visible occurrences of ch]
# 5) G(i,j,ch,x) = sum_S
#        I[S has EXACTLY x visible occurrences of ch lying in [i,j-4]]
# 6) H(i,j,ch,x) = sum_S
#        I[S has EXACTLY x visible occurrences of ch lying in [i,j-4],
#          and j is unpaired in S]
#---------------------------------------------------------------------


import sys,math,os
from misc import basePair

PRINT           = 0
PRINT_RNASUBOPT = 0
THETA           = 3
NUCL            = ['A','C','G','U']

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

def computeF(rna,Z):  
  #Compute the function F(i,j,ch,x) defined to be the
  #num str on [i,j] having x visible occurrences of nucleotide ch
  F   = {}
  n   = len(rna)
  rna = '$'+rna
  #initialization
  #----------------------------------------------
  #initialize all entries to 0
  for i in range(1,n+1):        #i in [1,n]
    for j in range(i,n+1):      #j in [i,n]
      for ch in NUCL:
        for x in range(n+1):    #x in [0,n]
          F[(i,j,ch,x)] = 0
  #----------------------------------------------
  #Base case: define F(i,j,ch,x) for i=j
  for i in range(1,n+1):        #i in [1,n]
    for ch in NUCL:
      if ch==rna[i]:
        F[(i,i,ch,1)] = 1
      else:
        F[(i,i,ch,0)] = 1
  #----------------------------------------------
  #Base case: define F(i,j,ch,x) for i<j<i+4
  for i in range(1,n+1):        #i in [1,n]
    for j in range(i+1,min(n,i+3)+1):    #j in [i+1,min(n,i+3)]
      for ch in NUCL:
        for x in range(j-i+2): #x in [0,j-i+1]
          if ch==rna[j]:
            if x>0:
              F[(i,j,ch,x)] = F[(i,j-1,ch,x-1)]
          else:
            F[(i,j,ch,x)] = F[(i,j-1,ch,x)]
  #----------------------------------------------
  #Inductive case: define F(i,j,ch,x) for j>=i+4
  for d in range(4,n):       #d in [4,n-1]
    for i in range(1,n-d+1): #i in [1,n-d]
      j = i+d                #note that [i,j] contains d+1 elements
      #Case 1: j unpaired in [i,j]
      for ch in NUCL:
        for x in range(d+2):              #x in [0,d+1]
          if ch==rna[j]:
            if x>0:
              F[(i,j,ch,x)] = F[(i,j-1,ch,x-1)]
          else:
            F[(i,j,ch,x)] = F[(i,j-1,ch,x)]
      #Case 2: (i,j) is base pair in S
      if basePair(rna[i],rna[j]):
        for ch in NUCL:
          F[(i,j,ch,0)] += Z[(i+1,j-1)] 
      #Case 3: (k,j) is base pair in S for k in [i+1,j-4]
      for k in range(i+1,j-3):       #k in [i+1,j-4]
        if basePair(rna[k],rna[j]):
          for ch in NUCL:
            for x in range(k-i+1):  #x in [0,k-i], where k-i num elem in [i,k-1]
              F[(i,j,ch,x)] += F[(i,k-1,ch,x)]*Z[(k+1,j-1)]
  return F

def computeFF(rna,Z):  
  #Compute the function FF(i,j,ch,x) defined to be the
  #num str on [i,j] having x visible occurrences of a nucleotide 
  #which can basepair with x. Note the difference with computeF(rna,Z)!!!
  FF  = {}
  n   = len(rna)
  rna = '$'+rna
  #initialization
  #----------------------------------------------
  #initialize all entries to 0
  for i in range(1,n+1):        #i in [1,n]
    for j in range(i,n+1):      #j in [i,n]
      for ch in NUCL:
        for x in range(n+1):    #x in [0,n]
          FF[(i,j,ch,x)] = 0
  #----------------------------------------------
  #Base case: define FF(i,j,ch,x) for i=j
  for i in range(1,n+1):        #i in [1,n]
    for ch in NUCL:
      if basePair(ch,rna[i]):
        FF[(i,i,ch,1)] = 1
      else:
        FF[(i,i,ch,0)] = 1
  #----------------------------------------------
  #Base case: define FF(i,j,ch,x) for i<j<i+4
  for i in range(1,n+1):        #i in [1,n]
    for j in range(i+1,min(n,i+3)+1):    #j in [i+1,min(n,i+3)]
      for ch in NUCL:
        for x in range(j-i+2): #x in [0,j-i+1]
          if basePair(ch,rna[j]):
            if x>0:
              FF[(i,j,ch,x)] = FF[(i,j-1,ch,x-1)]
          else:
            FF[(i,j,ch,x)]   = FF[(i,j-1,ch,x)]
  #----------------------------------------------
  #Inductive case: define FF(i,j,ch,x) for j>=i+4
  for d in range(4,n):       #d in [4,n-1]
    for i in range(1,n-d+1): #i in [1,n-d]
      j = i+d                #note that [i,j] contains d+1 elements
      #Case 1: j unpaired in [i,j]
      for ch in NUCL:
        for x in range(d+2):              #x in [0,d+1]
          if basePair(ch,rna[j]):
            if x>0:
              FF[(i,j,ch,x)] = FF[(i,j-1,ch,x-1)]
          else:
            FF[(i,j,ch,x)]   = FF[(i,j-1,ch,x)]
      #Case 2: (i,j) is base pair in S
      if basePair(rna[i],rna[j]):
        for ch in NUCL:
          FF[(i,j,ch,0)] += Z[(i+1,j-1)] 
      #Case 3: (k,j) is base pair in S for k in [i+1,j-4]
      for k in range(i+1,j-3):       #k in [i+1,j-4]
        if basePair(rna[k],rna[j]):
          for ch in NUCL:
            for x in range(k-i+1):  #x in [0,k-i], where k-i num elem in [i,k-1]
              FF[(i,j,ch,x)] += FF[(i,k-1,ch,x)]*Z[(k+1,j-1)]
  return FF

def computeG(rna,F,Z):
  G = {}; n = len(rna); rna = '$'+rna
  #-------------------------------------------
  #Initialization to zero
  for d in range(n):         #d in [0,n-1]
    for i in range(1,n+1-d): #i in [1,n-d]
      j = i+d
      for ch in NUCL:
        for x in range(n+1):
          G[(i,j,ch,x)] = 0
  #-------------------------------------------
  #Recursions
  for d in range(4,n):       #d in [4,n-1]
    for i in range(1,n+1-d): #i in [1,n-d]
      j = i+d                #note that num elements in [i,j] is d+1
      for ch in NUCL:
        for x in range(d+1-4+1):   #x in [0,d+1-4]
          #Case 1: positions j-3,j-2,j-1,j are unpaired in S
          G[(i,j,ch,x)] = F[(i,j-4,ch,x)]
          #Case 2: j-4+u is paired, but j-4+u+1,...,j are unpaired in S
          for u in range(1,5):     #u in [1,4], consider 4-u unpaired in suffix
            #Subcase A: (i,j-4+u) is base pair in S
            if x==0 and j-4+u-i>THETA and basePair(rna[i],rna[j-4+u]):
              G[(i,j,ch,x)] += Z[(i+1,j-4+u-1)]
            #Subcase B: (k,j-4+u) is base pair in S, some k in [i+1,j-4+u-4]
            for k in range(i+1,j-4+u-4+1):     #k in [i+1,(j-4+u)-4]
              if j-4+u-k>THETA and basePair(rna[k],rna[j-4+u]):
                G[(i,j,ch,x)] += F[(i,k-1,ch,x)]*Z[(k+1,j-4+u-1)]
  return G

def computeH(rna,Z,F):
  #H(i,j,ch,x) is num str S on [i,j] having x visible occ of ch in [i,j-4]
  #and in which j is unpaired in S.
  #WARNING: H(i,j,ch,0) defined to be zero if i <= j <= i+3
  H = {}; n = len(rna); rna = '$'+rna
  #-------------------------------------------
  #Initialization to zero
  for d in range(n):         #d in [0,n-1]
    for i in range(1,n+1-d): #i in [1,n-d]
      j = i+d
      for ch in NUCL:
        for x in range(n+1):
          H[(i,j,ch,x)] = 0
  #-------------------------------------------
  #Recursions
  for d in range(4,n):       #d in [4,n-1]
    for i in range(1,n+1-d): #i in [1,n-d]
      j = i+d                #note that num elements in [i,j] is d+1
      for ch in NUCL:
        for x in range(d+1-4+1):   #x in [0,d+1-4]
          #Case 1: positions j-3,j-2,j-1,j are unpaired in S
          H[(i,j,ch,x)] = F[(i,j-4,ch,x)]
          #Case 2: j-4+u is paired, but j-4+u+1,...,j are unpaired in S
          for u in range(1,4):     #u in [1,3], consider 4-u unpaired in suffix
            #Subcase A: (i,j-4+u) is base pair in S
            if x==0 and j-4+u-i>THETA and basePair(rna[i],rna[j-4+u]):
              H[(i,j,ch,x)] += Z[(i+1,j-4+u-1)]
            #Subcase B: (k,j-4+u) is base pair in S, some k in [i+1,j-4+u-4]
            for k in range(i+1,j-4+u-4+1):     #k in [i+1,(j-4+u)-4]
              if j-4+u-k>THETA and basePair(rna[k],rna[j-4+u]):
                H[(i,j,ch,x)] += F[(i,k-1,ch,x)]*Z[(k+1,j-4+u-1)]
  return H

def computeJ(rna,Z,FF):
  #J(i,j,ch,x) is num str S on [i,j] having x visible occ of a nucl [i,j-4]
  #which can basepair with rna[j] and in which j is unpaired in S.
  #WARNING: J(i,j,ch,0) defined to be zero if i <= j <= i+3
  J = {}; n = len(rna); rna = '$'+rna
  #-------------------------------------------
  #Initialization to zero
  for d in range(n):         #d in [0,n-1]
    for i in range(1,n+1-d): #i in [1,n-d]
      j = i+d
      for ch in NUCL:
        for x in range(n+1):
          J[(i,j,ch,x)] = 0
  #-------------------------------------------
  #Recursions
  for d in range(4,n):       #d in [4,n-1]
    for i in range(1,n+1-d): #i in [1,n-d]
      j = i+d                #note that num elements in [i,j] is d+1
      for ch in NUCL:
        for x in range(d+1-4+1):   #x in [0,d+1-4]
          #Case 1: positions j-3,j-2,j-1,j are unpaired in S
          J[(i,j,ch,x)] = FF[(i,j-4,ch,x)]
          #Case 2: j-4+u is paired, but j-4+u+1,...,j are unpaired in S
          for u in range(1,4):     #u in [1,3], consider 4-u unpaired in suffix
            #Subcase A: (i,j-4+u) is base pair in S
            if x==0 and j-4+u-i>THETA and basePair(rna[i],rna[j-4+u]):
              J[(i,j,ch,x)] += Z[(i+1,j-4+u-1)]
            #Subcase B: (k,j-4+u) is base pair in S, some k in [i+1,j-4+u-4]
            for k in range(i+1,j-4+u-4+1):     #k in [i+1,(j-4+u)-4]
              if j-4+u-k>THETA and basePair(rna[k],rna[j-4+u]):
                J[(i,j,ch,x)] += FF[(i,k-1,ch,x)]*Z[(k+1,j-4+u-1)]
  return J 


def computeE(rna,Z): 
  #function E(i,j) = sum_S sum_{(x,y)} I[ (x,y) external in S] 
  E   = {}
  n   = len(rna)
  rna = '$'+rna
  #initialize
  for i in range(1,n+1):
    for j in range(i,n+1):
      E[(i,j)] = 0
  #fill
  for d in range(THETA+1,n): #d in [4,n-1]
    for i in range(1,n-d+1): #i in [1,n-d]
      j = i+d
      sum = E[(i,j-1)]
      if basePair(rna[i],rna[j]):
        sum += Z[(i+1,j-1)]
      for k in range(i+1,j-THETA):# k in [i+1..j-THETA-1]
        if basePair(rna[k],rna[j]):
          sum += (E[(i,k-1)]+Z[(i,k-1)])*Z[(k+1,j-1)]
      E[(i,j)] = sum
  return E

def computeEleft(rna,Z): 
  #EL(i,j,ch) = sum_S sum_{(x,y)} 
  #   I[(x,y) external in S, x basepairs with ch]
  EL  = {}
  n   = len(rna)
  rna = '$'+rna
  #initialize
  for i in range(1,n+1):
    for j in range(i,n+1):
      for ch in NUCL:
        EL[(i,j,ch)] = 0
  #fill
  for ch in NUCL:             
    for d in range(THETA+1,n): #d in [4,n-1]
      for i in range(1,n-d+1): #i in [1,n-d]
        j = i+d
        sum = EL[(i,j-1,ch)]
        if basePair(rna[i],rna[j]) and basePair(rna[i],ch):
          sum += Z[(i+1,j-1)]
        for k in range(i+1,j-THETA):# k in [i+1..j-THETA-1]
          if basePair(rna[k],rna[j]):
            sum += EL[(i,k-1,ch)]*Z[(k+1,j-1)]
            if basePair(rna[k],ch):
              sum += Z[(i,k-1)]*Z[(k+1,j-1)]
#            sum += (EL[(i,k-1,ch)]+Z[(i,k-1)])*Z[(k+1,j-1)]
        EL[(i,j,ch)] = sum
  return EL

def computeEright(rna,Z): 
  #ER(i,j,ch) = sum_S sum_{(x,y)} 
  #   I[(x,y) external in S, y basepairs with ch]
  ER  = {}
  n   = len(rna)
  rna = '$'+rna
  #initialize
  for i in range(1,n+1):
    for j in range(i,n+1):
      for ch in NUCL:
        ER[(i,j,ch)] = 0
  #fill
  for ch in NUCL:             
    for d in range(THETA+1,n): #d in [4,n-1]
      for i in range(1,n-d+1): #i in [1,n-d]
        j = i+d
        sum = ER[(i,j-1,ch)]
        if basePair(rna[i],rna[j]) and basePair(rna[j],ch):
          sum += Z[(i+1,j-1)]
        for k in range(i+1,j-THETA):# k in [i+1..j-THETA-1]
          if basePair(rna[k],rna[j]):
            sum += ER[(i,k-1,ch)]*Z[(k+1,j-1)]
            if basePair(rna[j],ch):
              sum += Z[(i,k-1)]*Z[(k+1,j-1)]
#            sum += (ER[(i,k-1,ch)]+Z[(i,k-1)])*Z[(k+1,j-1)]
        ER[(i,j,ch)] = sum
  return ER


def computeEprime(rna,Z,E): 
  #Eprime(i,j) is sum over str on [i,j] sum over base pairs (x,y) in S
  #I[(x,y) external in S and y <= j-4, AND j unpaired in S]
  Eprime = {}; n = len(rna); rna = '$'+rna
  #-------------------------------------------
  #Initialization to zero
  for d in range(n):         #d in [0,n-1]
    for i in range(1,n+1-d): #i in [1,n-d]
      j = i+d
      Eprime[(i,j)] = 0
  #-------------------------------------------
  #Recursions
  for d in range(4,n):       #d in [4,n-1]
    for i in range(1,n+1-d): #i in [1,n-d]
      j = i+d                #note that num elements in [i,j] is d+1
      #Case 1: positions j-3,j-2,j-1,j are unpaired in S
      Eprime[(i,j)] = E[(i,j-4)]
      #Case 2: j-4+u is paired, but j-4+u+1,...,j are unpaired in S
      for u in range(1,4):     #u in [1,3], consider 4-u unpaired in suffix
        #Subcase A: (i,j-4+u) is base pair in S
        if j-4+u-i>THETA and basePair(rna[i],rna[j-4+u]):
          Eprime[(i,j)] += 0   #left here for clarity
        #Subcase B: (k,j-4+u) is base pair in S, some k in [i+1,j-4+u-4]
        for k in range(i+1,j-4+u-4+1):     #k in [i+1,(j-4+u)-4]
          if j-4+u-k>THETA and basePair(rna[k],rna[j-4+u]):
            Eprime[(i,j)] += E[(i,k-1)]*Z[(k+1,j-4+u-1)]
  return Eprime

def computeERprime(rna,Z,ER): 
  #ERprime(i,j,ch) is sum over str on [i,j] sum over base pairs (x,y) in S 
  #I[ S on [i,j], (x,y) external in S, y <= j-4, j unpaired in S and
  #   y basepairs with ch]
  ERprime = {}; n = len(rna); rna = '$'+rna
  #-------------------------------------------
  #Initialization to zero
  for d in range(n):         #d in [0,n-1]
    for i in range(1,n+1-d): #i in [1,n-d]
      for ch in NUCL:
        j = i+d
        ERprime[(i,j,ch)] = 0
  #-------------------------------------------
  #Recursions
  for ch in NUCL:
    for d in range(4,n):       #d in [4,n-1]
      for i in range(1,n+1-d): #i in [1,n-d]
        j = i+d                #note that num elements in [i,j] is d+1
        #Case 1: positions j-3,j-2,j-1,j are unpaired in S
        ERprime[(i,j,ch)] = ER[(i,j-4,ch)]
        #Case 2: j-4+u is paired, but j-4+u+1,...,j are unpaired in S
        for u in range(1,4):     #u in [1,3], consider 4-u unpaired in suffix
          #Subcase A: (i,j-4+u) is base pair in S
          if j-4+u-i>THETA and basePair(rna[i],rna[j-4+u]):
            ERprime[(i,j,ch)] += 0   #left here for clarity
          #Subcase B: (k,j-4+u) is base pair in S, some k in [i+1,j-4+u-4]
          for k in range(i+1,j-4+u-4+1):     #k in [i+1,(j-4+u)-4]
            ok = j-4+u-k>THETA and basePair(rna[k],rna[j-4+u])
            #WARNING: It is an ERROR to additionally require that
            #basePair(ch,rna[j-4+u]), since j-4+u too close to j and so
            #can't base pair with it.
            if ok:
              ERprime[(i,j,ch)] += ER[(i,k-1,ch)]*Z[(k+1,j-4+u-1)]
  return ERprime


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


def checkSumOfF(rna,F,Z):
  n = len(rna)
  L = []
  for ch in NUCL:
    sum = 0
    for x in range(n+1):
      sum += F[(1,n,ch,x)]
    text = "sum_x F(1,%s,%s,x) = %s != %s" % (n,ch,sum,Z[(1,n)])
    assert sum==Z[(1,n)],text

def checkSumOfG(rna,G,Z):
  n = len(rna)
  L = []
  for ch in NUCL:
    sum = 0
    for x in range(n+1):
      sum += G[(1,n,ch,x)]
    text = "sum_x G(1,%s,%s,x) = %s != %s" % (n,ch,sum,Z[(1,n)])
    assert sum==Z[(1,n)],text

def main(rna):
  n = len(rna)
  Z  = partFunInner(rna)
  F = computeF(rna,Z)
  G = computeG(rna,F,Z)
  if PRINT:
    for ch in NUCL:
      for x in range(n+1):
        if F[(1,n,ch,x)]!=0:
          print "F(1,%d,%s,%s)=%s" % (n,ch,x,F[(1,n,ch,x)])
  if PRINT_RNASUBOPT:
    cmd  = "echo %s | RNAsubopt -e 100" % rna
    file = os.popen(cmd)
    print file.read()
  if PRINT:
    for ch in NUCL:
      for x in range(n+1-4):
        if G[(1,n,ch,x)]>0:
          print "G(1,%s,%s,%s)=%s" % (n,ch,x,G[(1,n,ch,x)])
  checkSumOfF(rna,F,Z)
  checkSumOfG(rna,G,Z)


if __name__ == '__main__':
  if len(sys.argv)!= 2:
    print "Usage: %s RNA" % sys.argv[0]
    sys.exit(1)
  rna = sys.argv[1].upper()
  if len(rna)<5:
    print "Enter RNA of length at least 5"
  else:
    main(rna)


