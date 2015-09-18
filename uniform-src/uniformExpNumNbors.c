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
//~ #define i2(i,j,N) ((i)*(N)+(j))
//~ #define i3(i,j,c,N) ((i)*(N)+(j)+(c)*(N)*(N))
//~ #define i4(i,j,c,x,N) ((i)*(N)+(j)+(c)*(N)*(N)+(x)*4*(N)*(N))
//~ #define RgasConstant 0.0019872370936902486  //kcal per mol per K
//~ #define TEMP  310.1  //temperature
//~ #define RT    (0.0019872370936902486 * 310.1)
//~ #define DEBUG 0
#define MAXSIZE 1100         /*upper bound on num char in RNA sequence */
//~ #define min(a, b)  (((a) < (b)) ? (a) : (b))
//~ #define max(a, b)  (((a) > (b)) ? (a) : (b)) 
//~ #define int2c(x) ((x)==0 ? 'A' : ((x)==1 ? 'C' : ((x)==2 ? 'G' : 'U')))
                                                                               
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


int main(int argc, char *argv[]) {
  char rna[LEN];
  int seq_int[LEN];

  int i, j, k, d;
  int c,x;

  if (argc != 3) {
    fprintf(stderr,"Usage: %s RNA -ms1/ms2 \n",argv[0]);
    exit(1);
  }
  if (!strcmp(argv[2],"-ms1"))
	run_ms1(argv[1]);
  else if(!strcmp(argv[2],"-ms2"))
	run_ms2(argv[1]);
return 0;
}
  
  
