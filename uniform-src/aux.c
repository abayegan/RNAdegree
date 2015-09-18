/*---------------------------------------------
 * aux.c
 *
 * P.Clote
 * Code for EL,ER,ER1,F,G,Q to compute num bp add,rem,shifts.
 * WARNING: HOMOPOLYMER below is set to either 1 or 0, where
 * value of 1 means that any two nucleotides can base-pair, as
 * long as there are 3 unpaired nucleotides between them.
 * A COMMON error is not to check the value of HOMOPOLYMER!!!!
 * -----------------------------------------*/


/*---------------------------------------------------------------------
 This program computes Q by dynamic programming, where for given rna
 Q(i,j) = sum_{s on [i,j]} N(s), where N(s) is number of base pair
 additions, removals and shifts that transform s to a valid sec str t.
 Additionally, the number of structures Z(i,j) on [i,j] is computed,
 which is the partition function when energy of a base pair equals 0.
 NOTE: Throughout, the min number THETA of unpaired bases in a hairpin
 is 3; sometimes the constant THETA is explicitly noted, and other times,
 for brevity, simply the number 3 is written. This often happens in
 limits for sums (bounds for for-loops).
 NOTE: In order to realize the HOMOPOLYMER model, simply set
 the function basePair(x,y)=1 for all x,y in NUCL=['A','C','G','U'].
 
 Define the prefix of length n RNA sequence to be the region [1,n-4]
 and the suffix to be the region [n-3,n].

 Function Q(i,j) requires the following functions, where
 the expression sum_s means sum over s str on [i,j] below.

 1) EL(i,j,ch) = sum_s sum_{(x,y)} 
      I[(x,y) external in s, x basepairs with ch]
 2) ER(i,j,ch) = sum_s sum_{(x,y)} 
      I[(x,y) external in s, y basepairs with ch]
 3) ERprime(i,j,ch) = sum_s sum_{(x,y)}
      I[ s on [i,j], (x,y) external in s, y <= j-4, 
         j unpaired in s and y basepairs with ch]
 4) FF(i,j,ch,x),  = sum_s 
      I[ s has EXACTLY x many visible occurrences of a 
         a nucleotide that basepairs with ch]
 5) J(i,j,ch,x) = sums_
      I[s has EXACTLY x visible occurrences of a nucleotide
        in [i,j-4] which can basepair with ch, and j unpaired in s]

 NOTE: In the C-program, we have the following correspondences with
 the Python program:
      C-function		Python-function
    	EL			EL
	ER			ER
 	ER1			ERprime
	F			FF
	G			J			

 The function Q is defined in totalNumNborsAddRemovalShiftUnifProb.py
 which calls the functions EL,ER,ERprime,FF,J in the current file.
 In the C-program, Q is defined in the file aux.c with the functions
 EL,ER,ER1,F,G.

 NOTE: Some additional precursor functions are defined in this file.
 These are E(i,j), Ebis(i,j), Eprime(i,j), F(i,j,ch,x),
 G(i,j,ch,x), H(i,j,ch,x) defined as follows.

---------------------------------------------------------------------*/


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
#define c2int(c) ((c)=='A' ? 0 : ((c)=='C' ? 1 : ((c)=='G' ? 2 : 3)))

#define HOMOPOLYMER 0                                                                               
typedef long double DBL_TYPE;

int bp(int i,int j,char rna[MAXSIZE]){
  char x,y;
  int wc; //watson crick pair
  if (HOMOPOLYMER)
    return 1;
  x = toupper(rna[i]);
  y = toupper(rna[j]);
  wc =((x=='A'&&y=='U')||(x=='U'&&y=='A')||(x=='C'&&y=='G')||(x=='G'&&y=='C'));
  if (wc || (x=='G'&&y=='U')|| (x=='U'&&y=='G'))
    return 1;
  else
    return 0;
  }

int BP(int i,int j,char rna[LEN]){
  char x,y;
  int wc; //watson crick pair
  if (HOMOPOLYMER) 
    return 1;
  //nonhomopolymer case
  x = toupper(rna[i]);
  y = toupper(rna[j]);
  wc =((x=='A'&&y=='U')||(x=='U'&&y=='A')||(x=='C'&&y=='G')||(x=='G'&&y=='C'));
  if (wc || (x=='G'&&y=='U')|| (x=='U'&&y=='G'))
    return 1;
  else
    return 0;
  }

int basePair(char x, int c){
  //c codes a nt
  //0->A 1->C 2->G 3->U
  int wc; //watson crick pair
  //assert( ((0<=c)&&(c<4)) );
  if (HOMOPOLYMER)
    return 1;
  x = toupper(x);
  wc =((x=='A'&&c==3)||(x=='U'&&c==0)||(x=='C'&&c==2)||(x=='G'&&c==1));
  if (wc || (x=='G'&&c==3)|| (x=='U'&&c==2))
    return 1;
  else
    return 0;
  }


void computeZ(DBL_TYPE *Z, char *rna){
  /*---------------------------------------------------
  Z(i,j) = num str on [i,j]
  ---------------------------------------------------*/
  int i,j,d,k,N=strlen(rna);

  // Initialize all values to 0
  for (d = 0; d < N; d++) 
    for (i = 0; i < N-d; i++){
      j = i+d;
      Z[i2(i,j,N)] = 0;
      }

  // Base Case for Z is non-zero
  for (d = 0; d < 4; d++) 
    for (i = 0; i < N-d; i++){
      j = i+d;
      Z[i2(i,j,N)] = 1;
      }

  // Inductive Case
  for (d = 4; d < N; d++) {
    for (i = 0; i < N-d; i++) {
      j = i+d;
      //Case 0: j unpaired in [i,j]
      Z[i2(i,j,N)] = Z[i2(i,j-1,N)];
      //Case 1: (i,j) forms a base pair
      if (bp(i,j,rna)){
        Z[i2(i,j,N)] += Z[i2(i+1,j-1,N)];
        }
      //Case 2: k>i and (k,j) forms a base pair
      for (k = i+1; k<=j-4; k++) 
        if (bp(k,j,rna)){
          Z[i2(i,j,N)]  += Z[i2(i,k-1,N)]*Z[i2(k+1,j-1,N)];
          }

    }
  }
}


void computeEL(DBL_TYPE *EL, DBL_TYPE *Z, char *rna){
  /*---------------------------------------------------
  EL(i,j,c) = sum_S I[ (x,y) external bp in S, x base-pairs with c ]
  ---------------------------------------------------*/
  int c,i,j,d,k,N=strlen(rna);
  int Nz = N*N;
  int Nc = 4*Nz; 
  int Nx = (N+1)*Nc; 
  DBL_TYPE sum;

  // Initialize all values to 0
  for (c=0;c<4;c++)
    for (d = 0; d < N; d++) 
      for (i = 0; i < N-d; i++) {
        j = i+d;
        EL[i3(i,j,c,N)] = 0;
        }

  // Inductive Case
  for (c=0;c<4;c++)
    for (d = 4; d < N; d++) 
      for (i = 0; i < N-d; i++) {
        j = i+d;
        //Case 1: j unpaired in [i,j]
        sum = EL[i3(i,j-1,c,N)];
        //Case 2: (i,j) forms a base pair
        if (bp(i,j,rna) && basePair(rna[i],c))
          sum += Z[i2(i+1,j-1,N)];
        //Case 3: k>i and (k,j) forms a base pair
        for (k=i+1;k<=j-4;k++)
          if (bp(k,j,rna)){
            sum += EL[i3(i,k-1,c,N)]*Z[i2(k+1,j-1,N)];
            if (basePair(rna[k],c))
              sum += Z[i2(i,k-1,N)]*Z[i2(k+1,j-1,N)];
            }
        EL[i3(i,j,c,N)] = sum;
  }
}



void computeER(DBL_TYPE *ER, DBL_TYPE *Z, char *rna){
  /*---------------------------------------------------
  ER(i,j,c) = sum_S I[ (x,y) external bp in S, y base-pairs with c ]
  ---------------------------------------------------*/
  int c,i,j,d,k,N=strlen(rna);
  DBL_TYPE sum;
  int Nz = N*N;
  int Nc = 4*Nz; 
  int Nx = (N+1)*Nc; 

  // Initialize all values to 0
  for (c=0;c<4;c++)
    for (d = 0; d < N; d++) 
      for (i = 0; i < N-d; i++) {
        j = i+d;
        ER[i3(i,j,c,N)] = 0;
        }
  //---------------------------------------------------
  //Recursions 
  for (c=0;c<4;c++)
    for (d=4;d<N;d++)
      for (i=0;i<N-d;i++){
        j = i+d;
        sum = ER[i3(i,j-1,c,N)];
        if (bp(i,j,rna) && basePair(rna[j],c))
          sum += Z[i2(i+1,j-1,N)];
        for (k=i+1;k<=j-4;k++)
          if (bp(k,j,rna)){
            sum += ER[i3(i,k-1,c,N)]*Z[i2(k+1,j-1,N)];
            if (basePair(rna[j],c))
              sum += Z[i2(i,k-1,N)]*Z[i2(k+1,j-1,N)];
            }
        ER[i3(i,j,c,N)] = sum;
        }
  }



void computeER1(DBL_TYPE *ER1, DBL_TYPE *ER, DBL_TYPE *Z, char *rna){
  int c,i,j,d,k,u,N=strlen(rna);
  int Nz = N*N;
  int Nc = 4*Nz; 
  int Nx = (N+1)*Nc; 
  DBL_TYPE sum;

  // Initialize all values to 0
  for (c=0;c<4;c++)
    for (d = 0; d < N; d++) 
      for (i = 0; i < N-d; i++) {
        j = i+d;
        ER1[i3(i,j,c,N)] = 0;
        }
  //-------------------------------------------
  //Recursions
  for (c=0;c<4;c++)
    for (d=4;d<N;d++)
      for (i=0;i<N-d;i++){
        j = i+d;
        //Case 1: positions j-3,j-2,j-1,j are unpaired in S. Note rhs is ER!!
        sum = ER[i3(i,j-4,c,N)];
        //Case 2: j-4+u is paired, but j-4+u+1,...,j are unpaired in S
        for (u=1;u<=3;u++){
          //Subcase A: (i,j-4+u) is base pair in S
          if ((j-4+u-i>3) && bp(i,j-4+u,rna))
            sum += 0;   //left here for clarity
          //Subcase B: (k,j-4+u) is base pair in S, some k in [i+1,j-4+u-4]
          for (k=i+1;k<=j-4+u-4;k++)
            if ((j-4+u-k>3) && bp(k,j-4+u,rna))
              sum += ER[i3(i,k-1,c,N)]*Z[i2(k+1,j-4+u-1,N)];
          }
        ER1[i3(i,j,c,N)] = sum;
        }
  }



void computeF(DBL_TYPE *F, DBL_TYPE *Z, char *rna){
  /*---------------------------------------------------
  F(i,j,c,x) = sum_s I[ s has exactly x visible nt that can base-pair
               with c ]
  ---------------------------------------------------*/

  int c,i,j,d,k,x,N=strlen(rna);
  int Nz = N*N;
  int Nc = 4*Nz; 
  int Nx = (N+1)*Nc; 

  //Initialize all values to 0; calloc should do this, 
  //but this DIDN'T happen on my Mac Powerbook with /usr/local/bin/gcc-4.9 !

  for (d=0;d<N;d++)
    for (i=0;i<N-d;i++){
      j = i+d;
      for (c=0;c<4;c++)
        for (x=0;x<=N;x++)
          F[i4(i,j,c,x,N)] = 0;
      }

  //----------------------------------------------
  //Base case: define F(i,j,c,x) for i=j
  for (i=0;i<N;i++)
    for (c=0;c<4;c++)
      if (basePair(rna[i],c))
        F[i4(i,i,c,1,N)] = 1;
      else
        F[i4(i,i,c,0,N)] = 1;
  //----------------------------------------------
  //Base case: define F(i,j,c,x) for i<j<i+4
  for (i=0;i<N;i++)
    for (j=i+1;j<=min(N-1,i+3);j++)
      for (c=0;c<4;c++)
        for (x=0;x<=j-i+1;x++)
          if (basePair(rna[j],c)){
            if (x>0)
              F[i4(i,j,c,x,N)] = F[i4(i,j-1,c,x-1,N)];
              }
          else
            F[i4(i,j,c,x,N)] = F[i4(i,j-1,c,x,N)];
  //----------------------------------------------
  //Inductive case: define F(i,j,c,x) for j>=i+4
  for (d=4;d<N;d++)
    for (i=0;i<N-d;i++){
      j = i+d;                //note that [i,j] contains d+1 elements
      //Case 1: j unpaired in [i,j]
      for (c=0;c<4;c++)
        for (x=0;x<=j-i+1;x++)
          if (basePair(rna[j],c)){
            if (x>0)
              F[i4(i,j,c,x,N)] = F[i4(i,j-1,c,x-1,N)];
              }
          else
            F[i4(i,j,c,x,N)]   = F[i4(i,j-1,c,x,N)];
      //Case 2: (i,j) is base pair in S
      if (bp(i,j,rna))
        for (c=0;c<4;c++)
          F[i4(i,j,c,0,N)] += Z[i2(i+1,j-1,N)]; 
      //Case 3: (k,j) is base pair in S for k in [i+1,j-4]
      for (k=i+1;k<=j-4;k++)
        if (bp(k,j,rna))
          for (c=0;c<4;c++)
            for (x=0;x<=k-i;x++) //x in [0,k-i], where k-i num elem in [i,k-1]
              F[i4(i,j,c,x,N)]+=F[i4(i,k-1,c,x,N)]*Z[i2(k+1,j-1,N)];
      }
  }


void computeG(DBL_TYPE *G, DBL_TYPE *F, DBL_TYPE *Z, char *rna){
  /*---------------------------------------------------
  G(i,j,c,x) = sum_s I[ s has exactly x visible nt in [i,j-4] that 
               can base-pair with c, and j unpaired in s ]
  WARNING: G(i,j,c,0) defined to be zero if i <= j <= i+3
  ---------------------------------------------------*/

  int c,i,j,d,u,k,x,N=strlen(rna);
  int Nz = N*N;
  int Nc = 4*Nz; 
  int Nx = (N+1)*Nc; 

  // Initialization to 0
  for (d=0;d<N;d++)
    for (i=0;i<N-d;i++){
      j = i+d;
      for (c=0;c<4;c++)
        for (x=0;x<=j-i+1;x++)
          G[i4(i,j,c,x,N)] = 0;
      }

  //-------------------------------------------
  //Recursions
  for (d=4;d<N;d++)
    for (i=0;i<N-d;i++){
      j = i+d;                //note that num elements in [i,j] is d+1
      for (c=0;c<4;c++)
        for (x=0;x<=j-i+1-4;x++){
          //Case 1: positions j-3,j-2,j-1,j are unpaired in S
          G[i4(i,j,c,x,N)] = F[i4(i,j-4,c,x,N)];
          //Case 2: j-4+u is paired, but j-4+u+1,...,j are unpaired in S
          for (u=1;u<=3;u++){    //u in [1,3], consider 4-u unpaired in suffix
            //Subcase A: (i,j-4+u) is base pair in S
            if (x==0 && j-4+u-i>3 && bp(i,j-4+u,rna))
              G[i4(i,j,c,x,N)] += Z[i2(i+1,j-4+u-1,N)];
            //Subcase B: (k,j-4+u) is base pair in S, some k in [i+1,j-4+u-4]
            for (k=i+1;k<=j-4+u-4;k++) //k in [i+1,(j-4+u)-4]
              if (j-4+u-k>3 && bp(k,j-4+u,rna))
                G[i4(i,j,c,x,N)]+=F[i4(i,k-1,c,x,N)] * \
                                        Z[i2(k+1,j-4+u-1,N)];
            }
          }
      }
  }


void computeQ(DBL_TYPE *Q, DBL_TYPE *Z, DBL_TYPE *EL, DBL_TYPE *ER1, \
              DBL_TYPE *G, char *rna){
  /*---------------------------------------------------
  Q(i,j) = sum_s N(s), where N(s) is total number of base pair additions,
  removals, and shifts in s
  ---------------------------------------------------*/
  int i,j,d,k,x,N=strlen(rna);
  DBL_TYPE sum;
  int Nz = N*N;
  int Nc = 4*Nz; 
  int Nx = (N+1)*Nc; 
 
  //Initialize all values to 0
  for (d=0;d<N;d++)
    for (i=0;i<N-d;i++){
      j = i+d;
      Q[i2(i,j,N)] = 0;
      }

  //Recursions -- there are 5 cases
  for (d=4;d<N;d++)
    for (i=0;i<N-d;i++){
      j = i+d;

      //Case 1: first term Q(i,j-1)
      Q[i2(i,j,N)] = Q[i2(i,j-1,N)];
      //Case 2: second term 2 * sum_k z(i,k-1)*z(k+1,j-1)
      sum = 0;
      if (bp(i,j,rna))
        sum += Z[i2(i+1,j-1,N)];
      for (k=i+1;k<=j-4;k++)
        if (bp(k,j,rna))
          sum += Z[i2(i,k-1,N)]*Z[i2(k+1,j-1,N)];
      Q[i2(i,j,N)] += 2*sum;
      //Case 3: 2*EL[(i,j-1)]+2*ERprime[(i,j)]
      sum=2*EL[i3(i,j-1,c2int(rna[j]),N)]+2*ER1[i3(i,j,c2int(rna[j]),N)];
      Q[i2(i,j,N)] += sum;
      //Case 4: fourth term is sum_{x=2}^{n-4} x*(x-1)*G(1,n,ch,x)
      sum = 0;
      for (x=2;x<=j-i+1-4;x++)
        sum += x*(x-1)*G[i4(i,j,c2int(rna[j]),x,N)];
      Q[i2(i,j,N)] += sum;
      //Case 5: fifth term
      //sum_{k=1}^{n-4} z(k-1)*Q(n-k-1)+ Q(k-1)*z(n-k-1)
      //Case 5a: (i,j) paired
      sum = 0;
      if (bp(i,j,rna))
        sum += Q[i2(i+1,j-1,N)];
      for (k=i+1;k<=j-4;k++)
        if (bp(k,j,rna))
          sum+=Z[i2(i,k-1,N)]*Q[i2(k+1,j-1,N)]+Q[i2(i,k-1,N)]*Z[i2(k+1,j-1,N)];
      Q[i2(i,j,N)] += sum;
      }
  }

