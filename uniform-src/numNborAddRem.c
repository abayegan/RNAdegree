/*---------------------------------------------
 * numNborsAddRem.c
 *
 * WARNING: Assign HOMOPOLYMER to 0 or 1
 * 
 * P.Clote
 * Program computes sum_s P(s)*N(s), where sum
 * ranges over secondary structures s of given RNA
 * sequence, P(s) is Boltzmann probability and N(s)
 * is number of secondary structures at base pair
 * distance 1 from s. This code works for uniform
 * distribution only.
 * WARNING: SHIFTS are NOT allowed. These are computed
 * in the program numNborsAddRemShift.c
 *
 * WARNING: HOMOPOLYMER below is set to either 1 or 0, where
 * value of 1 means that any two nucleotides can base-pair, as
 * long as there are 3 unpaired nucleotides between them.
 * A COMMON error is not to check the value of HOMOPOLYMER!!!!
 * -----------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <limits.h> //for INT_MAX
#include <stdlib.h>
#include <ctype.h>  //for toupper()
#include <string.h>
#define LEN 1001
#define pf_index(i,j,N) ((i)*(N)+(j)-(i)*(1+(i))/2)
#define RgasConstant 0.0019872370936902486  //kcal per mol per K
#define TEMP  310.1  //temperature
#define RT    (0.0019872370936902486 * 310.1)
#define DEBUG 0

typedef long double DBL_TYPE;



int run_ms1(char rna[LEN]) {
  DBL_TYPE pf,numNbors,numStr,nhbors1,nhbors2,nhbors3;
  DBL_TYPE *QN,*QNB,*ZN,*Nhbors,*SS,val;
  int seqLen,arraySizeQ;
  int seq_int[LEN];

  int i, j, k, d;
  int ij,kj,ijminus1,iplus1jminus1,ikminus1,kplus1jminus1; //pf_index indices
  int zeroNminus1;

  rna[LEN-1] = '\0';

  seqLen = strlen(rna);

  // Allocate extern matrices
  arraySizeQ   = seqLen * (seqLen + 1) / 2 + seqLen + 1;
  //Nhbors is sum_{s in ss[i,j],(i,j) in s} N(s), where N(s) is num nhbors of s
  //SS is sum_{s in ss[i,j],(i,j) in s} 1; i.e.num str in [i,j]
  Nhbors     = (DBL_TYPE*)calloc(arraySizeQ, sizeof(DBL_TYPE));
  SS         = (DBL_TYPE*)calloc(arraySizeQ, sizeof(DBL_TYPE));
  
  // Base case for SS is non-zero, while for Nhbors, value is 0
  for (i = 0; i < seqLen; i++) 
    for (j = i; j < i+4; j++)
      SS[pf_index(i,j,seqLen)] = 1;

  // Inductive case
  for (d = 4; d < seqLen; d++) {
    for (i = 0; i < seqLen; i++) {
      j = i+d;
      if (j>=seqLen) break;
      ij              = pf_index(i,j,seqLen);
      ijminus1        = pf_index(i,j-1,seqLen);
      iplus1jminus1   = pf_index(i+1,j-1,seqLen);
      //case 0: j unpaired in [i,j]
      Nhbors[ij] = Nhbors[ijminus1];
      SS[ij]     = SS[ijminus1];

      //case 1: (i,j) forms a base pair
      if (BP(i,j,rna)){
        Nhbors[ij] += 2*SS[iplus1jminus1]+Nhbors[iplus1jminus1];
        SS[ij]     += SS[iplus1jminus1];
        }
      //case 2: k>i and (k,j) forms a base pair
      for (k = i+1; k<=j-4; k++) 
        if (BP(k,j,rna)){
          ikminus1      = pf_index(i,k-1,seqLen); 
          kplus1jminus1 = pf_index(k+1,j-1,seqLen); 
          kj            = pf_index(k,j,seqLen);
          SS[ij]  += SS[ikminus1]*SS[kplus1jminus1];
          //compute terms for Nhbors[ij]
          nhbors1       = 2*SS[ikminus1]*SS[kplus1jminus1];
          nhbors2       = Nhbors[ikminus1]*SS[kplus1jminus1];
          nhbors3       = SS[ikminus1]*Nhbors[kplus1jminus1];
          Nhbors[ij]   += nhbors1+nhbors2+nhbors3;
          }
    }
  }

  //  output
  zeroNminus1 = pf_index(0,seqLen-1,seqLen);
//  printf("rna\tnhbors (unif)\tnum str\tnhbors/numStr\n");
//  printf("%s\t%Lf\t%Lf\t",rna,Nhbors[zeroNminus1],SS[zeroNminus1]);
//  printf("%Lf\n",(Nhbors[zeroNminus1]/SS[zeroNminus1]));
  numNbors  = Nhbors[zeroNminus1];
  numStr    = SS[zeroNminus1];
  val = (Nhbors[zeroNminus1]/SS[zeroNminus1]);
  printf("%Lf\t%Lf\n",val,(val/seqLen));
 // printf("%.Lf\t%.Lf\t%Lf\t%d\t%Lf\n",numNbors,numStr,val,seqLen,(val/seqLen));
//  printf("%s\t%Lf\n",rna, (Nhbors[zeroNminus1]/SS[zeroNminus1]), (Nhbors[zeroNminus1]/SS[zeroNminus1]/(float) seqLen));
//  printf("%Lf\n",(Nhbors[zeroNminus1]/SS[zeroNminus1]), (Nhbors[zeroNminus1]/SS[zeroNminus1]/(float) seqLen));
 
  free(Nhbors);
  free(SS);
  Nhbors = NULL; SS = NULL;
  return(0);
}


