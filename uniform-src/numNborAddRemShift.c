/*---------------------------------------------
 * numNborAddRemShift.c
 *
 * P.Clote
 * Code to compute Q(i,j) = sum_s N(s), where N(s) is
 * number of base pair additions,removals and shifts to
 * transform s into a str t. This code requires aux.c
 * which latter contains functions EL, ER, ER1, F, G, Q.
 * 
 * WARNING: HOMOPOLYMER in aux.c is set to either 1 or 0, where
 * value of 1 means that any two nucleotides can base-pair, as
 * long as there are 3 unpaired nucleotides between them.
 * A COMMON error is not to check the value of HOMOPOLYMER!!!!
 * Check the aux.c file.
 * -----------------------------------------*/


#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <limits.h> //for INT_MAX
#include <stdlib.h>
#include <ctype.h>  //for toupper()
#include <string.h>
#define LEN 300
#define i2(i,j,N) ((i)*(N)+(j))
#define i3(i,j,c,N) ((i)*(N)+(j)+(c)*(N)*(N))
#define i4(i,j,c,x,N) ((i)*(N)+(j)+(c)*(N)*(N)+(x)*4*(N)*(N))
#define RgasConstant 0.0019872370936902486  //kcal per mol per K
#define TEMP  310.1  //temperature
#define RT    (0.0019872370936902486 * 310.1)
#define DEBUG 0
#define MAXSIZE 1100         /*upper bound on num char in RNA sequence */
#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b)) 
#define int2c(x) ((x)==0 ? 'A' : ((x)==1 ? 'C' : ((x)==2 ? 'G' : 'U')))
                                                                               
typedef long double DBL_TYPE;

/*----------------------------------------------------------
   Function prototypes
----------------------------------------------------------*/
//~ int basePair(int i,int j,char rna[MAXSIZE]);
//~ void computeZ(DBL_TYPE *Z, char *rna);
//~ void computeEL(DBL_TYPE *EL, DBL_TYPE *Z, char *rna);
//~ void computeER(DBL_TYPE *ER, DBL_TYPE *Z, char *rna);
//~ void computeER1(DBL_TYPE *ER1, DBL_TYPE *ER, DBL_TYPE *Z, char *rna);
//~ void computeF(DBL_TYPE *F, DBL_TYPE *Z, char *rna);
//~ void computeG(DBL_TYPE *G, DBL_TYPE *F, DBL_TYPE *Z, char *rna);
//~ void computeQ(DBL_TYPE *Q, DBL_TYPE *Z, DBL_TYPE *EL, DBL_TYPE *ER1, \
              //~ DBL_TYPE *G, char *rna);


int run_ms2(char rna[LEN]) {
  DBL_TYPE pf,nbors1,nbors2,nbors3,*Q,*Z,*EL,*ER,*ER1,*F,*G;
  DBL_TYPE numNbors,numStr,val;

  int N,Nz,Nc,Nx;
    //N=len(rna), Nz=arraySizeZ, Nc=4*arraySizeZ, Nx=(N+1)*c*arraySizeZ
  //char rna[LEN];
  int seq_int[LEN];

  int i, j, k, d;
  int c,x;

  rna[LEN-1] = '\0';

  N = strlen(rna);
  // Convert to upper case nucleotides A,C,G,U
  for (i=0;i<N;i++)
    rna[i] = toupper(rna[i]);

  // Allocate extern matrices
  Nz = N*N;
  Nc = 4*Nz; 
  Nx = (N+1)*Nc; 
  Z         = (DBL_TYPE*)calloc(Nz, sizeof(DBL_TYPE));
  //Z = num str in [i,j]
  EL        = (DBL_TYPE*)calloc(Nc, sizeof(DBL_TYPE));
  //EL(i,j,c)= sum_s [ num (x,y) ext in s, BP(x,c) ]
  ER        = (DBL_TYPE*)calloc(Nc, sizeof(DBL_TYPE));
  //ER(i,j,c)= sum_s [ num (x,y) ext in s, BP(y,c) ]
  ER1       = (DBL_TYPE*)calloc(Nc, sizeof(DBL_TYPE));
  //ER1(i,j,c)= sum_s [num (x,y) ext in s, BP(y,c), y<j-4, j unpaired in s ]
  F         = (DBL_TYPE*)calloc(Nx, sizeof(DBL_TYPE));
  //F(i,j,c,x)= sum_s [ s has exactly x visible occ of nt that pairs with c]
  G         = (DBL_TYPE*)calloc(Nx, sizeof(DBL_TYPE));
  //G(i,j,c,x)= sum_s [ s has exactly x visible occ of nt in [1,j-4] that 
  //                    can pair with c, and j unpaired in s]
  Q         = (DBL_TYPE*)calloc(Nz, sizeof(DBL_TYPE));
  //Q(i,j)= sum_s N(s), where N(s) is the number of base pair additions,
  //removals and shifts.

  /*--------------------------------------------------------------
  Main Recursions
  --------------------------------------------------------------*/
  
  computeZ(Z,rna);
  computeEL(EL,Z,rna);
  computeER(ER,Z,rna);
  computeER1(ER1,ER,Z,rna);
  computeF(F,Z,rna);
  computeG(G,F,Z,rna);
  computeQ(Q,Z,EL,ER1,G,rna);
  numNbors = Q[i2(0,N-1,N)];
  numStr   = Z[i2(0,N-1,N)];
  val      = numNbors/numStr;
  printf("%Lf\t%Lf\n",val,(val/N));
  //printf("%.Lf\t%.Lf\t%Lf\t%d\t%Lf\n",numNbors,numStr,val,N,(val/N));
#if 0
  printf("rna\t%s\n",rna);
  printf("Z\t%Lf\n",Z[i2(0,N-1,N)]);
  printf("Q\t%Lf\n",Q[i2(0,N-1,N)]);
  for (i=0;i<N;i++)
    for (j=i;j<N;j++)
      printf("%d\t%d\t%.Lf\n",i,j,Q[i2(i,j,N)]);
#endif
#if 0
  for (c=0;c<4;c++){
    printf("Char is %c\n",int2c(c));
    printf("EL\t%Lf\n",EL[i3(0,N-1,c,N)]);
    printf("ER\t%Lf\n",ER[i3(0,N-1,c,N)]);
    printf("ER1\t%Lf\n",ER1[i3(0,N-1,c,N)]);
    for (x=0;x<=3;x++){
      printf("\tG(0,%d,%c,%d)=%.Lf\n",N-1,int2c(c),x,G[i4(0,N-1,c,x,N)]);
    }
  }
#endif
  free(Q); free(Z); free(EL); free(ER); free(ER1); free(F); free(G);
  Q = NULL; Z = NULL; EL = NULL; ER = NULL; ER1 = NULL; F = NULL; G = NULL;
  return(0);
}


