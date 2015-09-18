/*---------------------------------------------
 * expNumNhborsTurner.c
 *
 * P.Clote
 * Program computes sum_s P(s)*N(s), where sum
 * ranges over secondary structures s of given RNA
 * sequence, P(s) is Boltzmann probability and N(s)
 * is number of secondary structures at base pair
 * distance 1 from s. This code works for Turner energy.
 * -----------------------------------------*/

#include <stdio.h>
#include<assert.h>
#include <math.h>
#include <limits.h> //for INT_MAX
#include <stdlib.h>
#include <ctype.h>  //for toupper()
#include <string.h>
#define LEN 300
#define pf_index(i,j,N) ((i)*(N)+(j)-(i)*(1+(i))/2)
#define RgasConstant 0.0019872370936902486  //kcal per mol per K
#define TEMP  310.15  //temperature
//#define RT    (0.0019872370936902486 * 310.1)
#define RT   (0.00198717 * 310.15)
#define arraySizeQ   ((N)*((N)+1)/2+(N)+1)
#define ENERGY_IS_ZERO 0



#include"fold.h"
#include"energy_const.h"
#include"fold_vars.h"
#include"pair_mat.h"
#include"convert_Vienna.h"
#include"params.h"
#include"myConst.h"
#include"misc.h"
#include"RNAhairpin.h"
#include"McCaskill.h"
#include "expNumNhborsTurner.h"

typedef long double DBL_TYPE;

/*---------------------------------------------------------
  Function prototypes for following:
    encode_seq() from Vienna RNA Package 1.8.5
    MaxNumHairpin()
    HairpinPartition()
    sample(),sampleB(),sampleM(),sampleM1()
---------------------------------------------------------*/
PRIVATE short *encode_seq(const char *seq);

int main(int argc, char *argv[]){
  char sequence[MAXSIZE],filename[FILENAMELENGTH];
  char paren[MAXSIZE];
  FILE *inputfile;
  int i,j,k,h;
  int ok=0; 
  //ok used for checking that -s rnaSeq in command line arguments
  double totalpar; //total partition function
  double totalNumStr; //total number of structures
  double expNumNhbors, expNumNhborsNorm;

  /*---------------------------------------------------
  Check command line parameters
  ---------------------------------------------------*/
  if (argc<2){ 
    printf("Usage: %s rna\n",argv[0]);
    exit(1);
    }
  sequence[0]='@';
  strncpy(sequence+1,argv[1],strlen(argv[1])); 
  sequence[strlen(argv[1])+1]='\0'; //termination character for string
 
  /*---------------------------------------------------
  Set up computation
  ---------------------------------------------------*/

  CheckSequence(sequence);
  S0=encode_seq(sequence+1);
  seqlen=strlen(sequence)-1;
  Initialize_Params();
  make_pair_matrix();//needed for pair matching
  kT = (temperature+K0)*GASCONST/1000.0;
  //printf("%.15f\n",kT);
  if (ENERGY_IS_ZERO){
    ML_base=0;
    ML_close=0;
    }
  else{
    ML_base=(double)P->MLbase/100;
    ML_close=(double)P->MLclosing/100;
    }
  
  //allocate space for partition function 
  double **Z, **ZB, **ZM1, **ZM;
  Z   = Allocate2DMatrix( seqlen+1,seqlen+1);
  ZB  = Allocate2DMatrix( seqlen+1,seqlen+1);
  ZM1 = Allocate2DMatrix( seqlen+1,seqlen+1);
  ZM  = Allocate2DMatrix( seqlen+1,seqlen+1);
  for(i=0;i<=seqlen;i++) //necessary since Yang didn't use calloc
    for (j=0;j<=seqlen;j++) {
	  ZB[i][j]=0; Z[i][j]=0; ZM1[i][j]=0; ZM[i][j]=0;
	}

  //allocate space for Q partition function 
  double **Q, **QB, **QM1, **QM;
  Q   = Allocate2DMatrix( seqlen+1,seqlen+1);
  QB  = Allocate2DMatrix( seqlen+1,seqlen+1);
  QM1 = Allocate2DMatrix( seqlen+1,seqlen+1);
  QM  = Allocate2DMatrix( seqlen+1,seqlen+1);
  for(i=0;i<=seqlen;i++) //necessary since Yang didn't use calloc
    for (j=0;j<=seqlen;j++) {
	  QB[i][j]=0; Q[i][j]=0; QM1[i][j]=0; QM[i][j]=0;
	}

  //compute partition function 
  mcCaskill(sequence, Z, ZB, ZM1, ZM);
  mcCaskillQ(sequence,Q,QB,QM1,QM, Z, ZB, ZM1, ZM);
  expNumNhbors = Q[1][seqlen]/Z[1][seqlen];
  expNumNhborsNorm = expNumNhbors/seqlen;
//  printf("%s\t%lf\t%d\t%lf\n",sequence+1,expNumNhbors,seqLen,expNumNhborsNorm);
  printf("%lf\t%d\t%lf\n",expNumNhbors,seqlen,expNumNhborsNorm);
#if 0
  totalpar    = Z[1][seqlen];
  totalNumStr = Z[seqlen][1];
  printf("The RNA sequence is %s:\n",sequence+1);
  printf("Partition function:%lf\n",totalpar);
  printf("Num structures:%lf\n",totalNumStr);
  printf("Q[1][n]:%lf\n",Q[1][seqlen]);
//  printf("Q[n][1]:%lf\n",Q[seqlen][1]);
#endif
  return 0;

#if 0
  PartitionFunction *pf = mcCaskill(sequence);
  totalpar    = pf->Z[1][seqlen];
  totalNumStr = pf->Z[seqlen][1];
  printf("The RNA sequence is %s:\n",sequence+1);
#endif
}

/*Below modified from fold.c*/
//Further modified by P.Clote on 12 May 2014
PRIVATE short *encode_seq(const char *seq) {
  unsigned int k,l;
  short *S0_out;
  int i;
  l = strlen(seq);
  if ( (S0_out = (short *) calloc(1, sizeof(short)*(l+2) )) == NULL) {
    printf("Out of memory!\n");
    exit(1);
    }
  S0_out[0]= (short) l;
  for (k=1; k<=l; k++) { /* make numerical encoding of seq */
    S0_out[k]= (short) encode_char(toupper(seq[k-1]));
  }
  return S0_out;
}

