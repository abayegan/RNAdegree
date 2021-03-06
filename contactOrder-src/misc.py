#! /usr/bin/env python

# misc.py
# P.Clote
# Sat Dec  2 06:32:48 EST 2006

# Miscellaneous
# Functions defined in this program
# error(s)
# firstProj(L)
# printList(L)
# printDict(D)
# guPair(x,y)
# auPair(x,y)
# cgPair(x,y)
# basePair(x,y)
# watsonCrick(x,y)
# a(i,j,s)
# getRNA(rna,rnaIsString)
# list2string(L)
# basePairList(secStr)
# numBasePairs(secStr)
# compatible(rna,secStr)
# mountainLevel(s)
# basePairDist(S,T)      where S,T are lists of bp
# basePairAccuracy(S,T)  where S,T are lists of bp
# minNumRemovedFromLocallyOpt(rna,s)


import sys,os,tempfile,string,copy

DEBUG   = 0

BINDIR    = '/usr/local/' #directory of RNA exec on miles
RNAfold   = BINDIR+'bin/RNAfold -d2 '
RNAeval   = BINDIR+'bin/RNAeval -d2 '
RNAsubopt = BINDIR+'bin/RNAsubopt -d2 '
RNAshapes = BINDIR+'bin/RNAshapes '
NORMALIZE = 0 #normalize energies by sequence length

RgasConstant = 0.0019872370936902486  #kcal per mol per K
TEMP      = 310.1  #temperature 
_RT_      = 0.0019872370936902486 * 310.1 

GUbond  = -1
AUbond  = -1
CGbond  = -1
INF     = sys.maxint
THRESHOLD = 3
  #for i and j to basepair, |i-j| >= THRESHOLD


def error(s):
  print "%s" % s
  sys.exit(1)

def firstProj( L ): #projects first coordinate of tuple, list or sequence
  return L[0]  #used together with function map
                                                                                
def printList(L):
  for x in L: print x

def printDict(D):
  keys = D.keys()
  keys.sort()
  for key in keys:
    print key,D[key]


def guPair(x,y):
   return ( x=='G' and y=='U' or x=='U'and y=='G' )
def auPair(x,y):
   return ( x=='A' and y=='U' or x=='U'and y=='A' )
def cgPair(x,y):
   return ( x=='G' and y=='C' or x=='C'and y=='G' )
def basePair(x,y):
   return auPair(x,y)  or  cgPair(x,y)  or  guPair(x,y) 
def watsonCrick(x,y):
   return auPair(x,y)  or  cgPair(x,y) 

def a(i,j,s):
  if i<len(s) and i<j<len(s): 
    if auPair(s[i],s[j]):
      return AUbond
    elif cgPair(s[i],s[j]):
      return CGbond
    elif guPair(s[i],s[j]):
      return GUbond
    else:
      return INF
  else:
    error("Indices out of bounds for RNA sequence")


def getRNA(rna,rnaIsString):
  if rnaIsString:
    return rna.upper()
  else:
    file = open(rna,"r")
    line = file.readline()
    if line[0]==">":
      line = file.readline()
    return line.strip().upper()
   


def list2string(L):
  s = string.join(L,'')
  return s

def basePairList(secStr):
  L = []; stack = []
  for i in range(len(secStr)):
    if secStr[i]=="(":
      stack.append(i)
    elif secStr[i]==")":
      x = stack.pop() #return and pop end of stack 
      L.append( (x,i) ) 
  return L
  
def numBasePairs(secStr):
  count = 0
  for ch in secStr: 
    if ch=='(': count += 1
  return count

def compatible(rna,secStr):
  #whether secStr is compatible with rna sequence
  #secStr is string
  assert len(rna)==len(secStr)
  S = basePairList(secStr)
  for (i,j) in S:
    if not basePair(rna[i],rna[j]): return 0
  return 1


def mountainLevel(s):
  #s is sec str as a string, ie ((...))
  D = {}
  n = 0; D[0] = []
  levels = []
  for i in range(len(s)):
    if s[i]=="(":
      n += 1
    elif s[i]==")":
      n -= 1
    if n not in D.keys():
      D[n] = [i]
      levels.append(n)
    else: 
      D[n].append(i)
      levels.append(n)
  return levels,D


def basePairDist(S,T): 
  #S,T are lists of base pairs, produced by basePairList
  #compute and return base pair distance
  num=0
  for bp in S:
    if bp not in T: num+=1
  for bp in T:
    if bp not in S: num+=1
  return num

def basePairAccuracy(S,T): 
  #S,T are lists of base pairs, produced by basePairList
  #S is native sec str, T is model sec str (eg from RNAbor)
  #compute and return fraction of bp of S in T
  lenS = len(S); lenT = len(T)
  num = 0
  for bp in S:
    if bp in T: num+=1
  truePos        = num
  sensitivity    = float(num)/lenS
  num = 0
  for bp in T:
    if bp in S: num+=1
  posPredVal     = float(num)/lenT
  return sensitivity,posPredVal,truePos

def minNumRemovedFromLocallyOpt(rna,s):
  #Computes minimum number of base pairs which can be added to 
  #obtain a locally optimal sec str
  #rna is RNA sequence, uppercase, s is secStr as string
  levels,D = mountainLevel(s)
  keys = D.keys()
  keys.sort()
  count    = 0; minCount = len(rna); maxCount = -1; found = 0
  for key in keys:
    m = len(D[key])
    for i in range(m-1):
      for j in range(i+1,m):
        x = D[key][i]; y = D[key][j]
        free = (s[x]=="." and s[y]==".")
        if abs(y-x)>THRESHOLD and free and basePair(rna[x],rna[y]):
          found = 1
          sNew = s[:x]+"("+s[x+1:y]+")"+s[y+1:]
          (a,b) = minNumRemovedFromLocallyOpt(rna,sNew)
          count = a+1
          if count < minCount:
            minCount  = count
            LminCount = copy.deepcopy(b)
            LminCount.append( (x,y) )
  if not found: #found=1 means that current structure is loc opt
    return (0,[])
  else:
    return (minCount,LminCount)




def main(rna,s):
  levels,D = mountainLevel(s)
  print rna
  print s
  print levels 
  printDict(D)
  count = 0
  for k in D.keys():
    count += len(D[k])
  print "Total count: %d" % count
  

if __name__ == '__main__':
  if len(sys.argv) < 3:
    text = """Usage: %s rna secStr
           1) rna is RNA string
           2) secStr is balanced paren secondary structure representation
           """ 
    text = text % sys.argv[0]
    print text
    sys.exit(1)
  else:
    rna    = sys.argv[1].upper() #RNA sequence as a string
    secStr = sys.argv[2]
    print minNumRemovedFromLocallyOpt(rna,secStr)
    main(rna,secStr)
 
