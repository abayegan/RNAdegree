/*
/*
 * expNumNborTurner.c
 *
 *  Created on: Apr 29, 2015
 *      Author: A. Bayegan, P. Clote
 */

// C headers
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <limits.h> //for INT_MAX
#include <stdlib.h>
#include <ctype.h>  //for toupper()
#include <string.h>

// ViennaRNA headers
#include "fold.h"
#include "fold_vars.h"
#include "HP.h"
#include "read_epars.h"


//Program headers
#include "convert_Vienna.h"
#include "myConst.h"
#include "misc.h"
#include "ms2_Turner.h"

#define LEN 300
#define pf_index(i,j,N) ((i)*(N)+(j)-(i)*(1+(i))/2)
#define RgasConstant 0.0019872370936902486  //kcal per mol per K
//#define TEMP 310.15
#define RT   (0.00198717 * 310.15)
#define arraySizeQ   ((N)*((N)+1)/2+(N)+1)
#define ENERGY_IS_ZERO 0
extern paramT *P;

int ARC[4][MAXSIZE][MAXSIZE];
PRIVATE short *encode_seq(const char *seq);

int main(int argc, char *argv[]){
	DBL_TYPE pf,nbors1,nbors2,nbors3, *EL, *ER, *ER1, *G, *F;
	DBL_TYPE numNbors,numStr,val;
	char sequence[MAXSIZE],filename[FILENAMELENGTH];
	char paren[MAXSIZE];
	FILE *inputfile;
	int i,j,k,h,c;
	int N,Nz,Nc,Nx;
	int ok=0;
	//ok used for checking that -s rnaSeq in command line arguments
	double totalpar; //total partition function
	double totalNumStr; //total number of structures
	double expNumNhbors, expNumNhborsNorm;
	int t99Flag = 0, VERBOSE = 0;
	double TEMP=37.;
	const char turner99_eng[] = "rna_turner1999.par";
	/*---------------------------------------------------
	  Check command line parameters
	  ---------------------------------------------------*/
	if (argc<2){
		dispUsage(argv[0]);
	}
	if (argc>2)
		for(i=2;i<argc;i++)
			if (argv[i][0] == '-')
				switch (argv[i][1]){
					case 't':
						if(argc>i+1 && argv[i+1][0]!='-'){
							char *e;
							TEMP = strtod(argv[i+1], &e);
							if(e==argv[i+1]){  //error in conversion
								printf("\nerror in the input temperature");
								dispUsage(argv[0]);
							}
						}
						else
						{
							printf("\nerror in the input temperature");
							dispUsage(argv[0]);
						}
						break;
					case 'e':
						if(argc>i+1 && !strcmp(argv[i+1],"99"))
							t99Flag = 1;
						else if (argc>i+1 && !strcmp(argv[i+1],"04"))
							t99Flag = 0;
						else
						{
							printf("\nerror in the input Turner flag");
							dispUsage(argv[0]);
						}
						break;
					case 'v':
						VERBOSE=1;
						break;
					case 'h':
						displayReadMe();
					default : 
						dispUsage(argv[0]);					
				}
	if(!strcmp(argv[1],"-h"))
		displayReadMe();
	
	sequence[0]='@';
	strncpy(sequence+1,argv[1],strlen(argv[1]));
	sequence[strlen(argv[1])+1]='\0'; //termination character for string
	  //printf("seq:%s\n",sequence);


	/*---------------------------------------------------
	  Set up computation
	  ---------------------------------------------------*/

	CheckSequence(sequence);
	temperature = TEMP;
	if (t99Flag != 0)
		read_parameter_file(turner99_eng); //use turner99 energy model
	S0=encode_seq(sequence+1);
	seqlen=strlen(sequence)-1;
	Initialize_Params();
	make_pair_matrix();//needed for pair matching
	kT = (temperature+K0)*GASCONST/1000.0;
	
	update_fold_params();
	IL_initialize(seqlen);
	HP_init(seqlen);
	
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
		for (j=0;j<=seqlen;j++)  {
			ZB[i][j]=0; Z[i][j]=0; ZM1[i][j]=0; ZM[i][j]=0;
		}

	double **Q, **QB, **QM1, **QM;
	Q   = Allocate2DMatrix( seqlen+1,seqlen+1);
	QB  = Allocate2DMatrix( seqlen+1,seqlen+1);
	QM1 = Allocate2DMatrix( seqlen+1,seqlen+1);
	QM  = Allocate2DMatrix( seqlen+1,seqlen+1);
	for(i=0;i<=seqlen;i++) //necessary since Yang didn't use calloc
		for (j=0;j<=seqlen;j++) {
			QB[i][j]=0; Q[i][j]=0; QM1[i][j]=0; QM[i][j]=0;
		}
	N = strlen(sequence);
	Nz = (N)*(N);
	Nc = 4*Nz; 	
	Nx = (N+1)*Nc; 
	//DBL_TYPE ***EL;
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

	computeARC(ARC,sequence);
	mcCaskill(sequence, Z, ZB, ZM1, ZM);
	computeEL(EL,Z,ZB,sequence);
	computeER(ER,Z,ZB,sequence);
	computeER1(ER1,ER,Z,ZB,sequence);
	computeF(F,Z,ZB,sequence);
	computeG(G,F,Z,ZB,sequence);
	mcCaskillQ(Q,QB,QM1,QM,Z,ZB,ZM1,ZM,EL,ER1,G,F,sequence);
	expNumNhbors = Q[1][seqlen]/Z[1][seqlen];
	expNumNhborsNorm = expNumNhbors/seqlen;
	if(VERBOSE)
		printf("RNA Sequence:%s\nExpected degree:%lf\tNormalized expected degree:%lf\nPartition function:%lf\tTotal number of structures:%lf\n",
				sequence+1,expNumNhbors,expNumNhborsNorm,Z[1][seqlen],Z[seqlen][1]);
	else
		printf("%lf\t%lf\n",expNumNhbors,expNumNhborsNorm);
#if 0
	totalpar    = Z[1][seqlen];
	totalNumStr = Z[seqlen][1];
	printf("The RNA sequence is %s:\n",sequence+1);
	printf("Partition function:%lf\n",totalpar);
	printf("Num structures:%lf\n",totalNumStr);
	printf("Q[1][%d]:%lf\tQB[1][n]:%lf\tQM[1][n]:%lf\tQM1[1][n]:%lf\n",seqlen,Q[1][seqlen],QB[1][seqlen],QM[1][seqlen],QM1[1][seqlen]);
	printf("Z[1][%d]:%lf\tZB[1][n]:%lf\tZM[1][n]:%lf\tZM1[1][n]:%lf\n",seqlen,Z[1][seqlen],ZB[1][seqlen],ZM[1][seqlen],ZM1[1][seqlen]);
	printf("Q[%d][1]:%lf\n",seqlen,Q[seqlen][1]);
#endif
#if 0
	printf("arc1A[1][%d]=%d\n",seqlen,arc1A(1,seqlen,sequence));
	printf("arc1B[1][%d]=%d\n",seqlen,arc1B(1,seqlen,sequence));
	printf("arc1C[1][%d]=%d\n",seqlen,arc1C(1,seqlen,sequence));
	printf("arc2A[1][%d]=%d\n",seqlen,arc2A(1,seqlen,sequence));
	printf("arc2B1[1][%d]=%d\n",seqlen,arc2B1(1,seqlen,sequence));
	printf("arc2B2[1][%d]=%d\n",seqlen,arc2B2(1,seqlen,sequence));
	printf("arc2C1[1][%d]=%d\n",seqlen,arc2C1(1,seqlen,sequence));
	printf("arc2C2[1][%d]=%d\n",seqlen,arc2C2(1,seqlen,sequence));



#endif
#if 0
printf("\nZ\ti\tj------\n");
  for (i=1;i<seqlen+1;i++)
    for (j=i;j<seqlen+1;j++)
      printf("%d\t%d\t%.5f\n",i,j,Z[i][j]);
printf("\nZB\ti\tj------\n");
  for (i=1;i<seqlen+1;i++)
    for (j=i;j<seqlen+1;j++)
      printf("%d\t%d\t%lf\n",i,j,ZB[i][j]);
printf("\nZM\ti\tj------\n");
  for (i=1;i<seqlen+1;i++)
    for (j=i;j<seqlen+1;j++)
      printf("%d\t%d\t%lf\n",i,j,ZM[i][j]);
printf("\nZM1\ti\tj------\n");
  for (i=1;i<seqlen+1;i++)
    for (j=i;j<seqlen+1;j++)
      printf("%d\t%d\t%lf\n",i,j,ZM1[i][j]);
#endif

#if 0
printf("\nQ\ti\tj------\n");
  for (i=1;i<seqlen+1;i++)
    for (j=i;j<seqlen+1;j++)
      printf("%d\t%d\t%.5f\n",i,j,Q[i][j]);
printf("\nQB\ti\tj------\n");
  for (i=1;i<seqlen+1;i++)
    for (j=i;j<seqlen+1;j++)
      printf("%d\t%d\t%lf\n",i,j,QB[i][j]);
printf("\nQM\ti\tj------\n");
  for (i=1;i<seqlen+1;i++)
    for (j=i;j<seqlen+1;j++)
      printf("%d\t%d\t%lf\n",i,j,QM[i][j]);
printf("\nQM1\ti\tj------\n");
  for (i=1;i<seqlen+1;i++)
    for (j=i;j<seqlen+1;j++)
      printf("%d\t%d\t%lf\n",i,j,QM1[i][j]);
#endif

#if 0
//~ ;
//~ int d;
//~ for (c=0;c<4;c++)
		//~ for (d = 4; d < N; d++)
			//~ for (i = 0; i < N-d; i++) {
				//~ j=i+d;
				//~ printf("EL**(%d,%d,%d)\t%Lf\n",i,j,c,EL[i3(i,j,c,N)]);
			//~ }
  int x;
  for (c=0;c<4;c++){
    printf("Char is %c\n",int2c(c));
    printf("EL(1,%d)\t%Lf\n",seqlen,EL[i3(1,seqlen,c,N)]);
    printf("ER\t%Lf\n",ER[i3(1,seqlen,c,N)]);
    printf("ER1\t%Lf\n",ER1[i3(1,seqlen,c,N)]);
    for (i=1;i<seqlen+1;i++)
		for (j=i;j<seqlen+1;j++)
			printf("ER[%d][%d][%c]\t%Lf\n",i,j,int2c(c),ER[i3(i,j,c,N)]);
			//~ for (x=0;x<=3;x++){
				//~ printf("\tG(%d,%d,%c,%d)=%.Lf\n",i,j,int2c(c),x,G[i4(i,j,c,x,N)]);
    //~ }
  }
#endif
	free(Q); free(Z); free(EL); free(ER); free(ER1); free(F); free(G);
	return 0;
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

