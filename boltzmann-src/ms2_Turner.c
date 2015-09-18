/*
 * aux_Turner.c
 *
 *  Created on: Apr 29, 2015
 *      Author: A. Bayegan, P.Clote
 */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include"myConst.h"
#include"fold.h"
#include"energy_const.h"
#include"fold_vars.h"
#include"pair_mat.h"
#include"convert_Vienna.h"
#include"params.h"
#include<limits.h>
#include"misc.h"

extern int ARC[4][MAXSIZE][MAXSIZE];

//counts the number of neighbors obtained from adding (x,y)
int arc1A_ms2(int i, int j, char rna[MAXSIZE]){
	int x,sum=0;
	if(i>j) return 0;
	for(x=i;x<=j-4;x++)
		sum += ARC[ch2int(rna[x])][x+4][j];
	//~ printf("arc1A_ms2(1):%d\n",sum);
	//~ int y;
	//~ sum=0;
	//~ for (x=i;x<=j-4;x++)
		//~ for(y=x+4;y<=j;y++)
			//~ if(BP(x,y,rna))
			//~ {
				//~ printf("adding base pair(%d,%d)\n",x,y);
				//~ sum+=1;
			//~ }
	//~ printf("arc1A_ms2(2):%d\n",sum);
	return sum;
}

//counts the number of neighbors obtained from shifting (i,j) to (i,y)
int arc1B_ms2(int i, int j, char rna[MAXSIZE]){
	if(i>j) return 0;
	return ARC[ch2int(rna[i])][i+4][j-1];
}

//counts the number of neighbors obtained from shifting (i,j) to (x,j)
int arc1C_ms2(int i, int j,char rna[MAXSIZE]){
	if(i>j) return 0;
	return ARC[ch2int(rna[j])][i+1][j-4];
}

//arc2 counts the number of neighbors for internal loops bound by (i,j) and (l,r) where i < l < r < j
//counts the number of neighbors obtained from adding the base pair (x,y) to the internal loop
int arc2A_ms2(int i, int j, int l, int r, char rna[MAXSIZE]){
	int x,sum=0;
	if (!((i<l) && (l<r) && (r<j))) return 0;
	for(x=i+1;x<=l-1;x++)
		sum += ARC[ch2int(rna[x])][r+1][j-1];
	return sum;
}

//counts the number of neighbors obtained from shifting (i,j) to (i,y) in the internal loop
int arc2B1_ms2(int i, int j, int l, int r, char rna[MAXSIZE]){
	if (!((i<l) && (l<r) && (r<j))) return 0;
	return ARC[ch2int(rna[i])][r+1][j-1] + ARC[ch2int(rna[i])][i+4][l-1];
}

//counts the number of neighbors obtained from shifting (l,r) to (l,y) in the internal loop
int arc2B2_ms2(int i, int j, int l, int r, char rna[MAXSIZE]){
	if (!((i<l) && (l<r) && (r<j))) return 0;
	return ARC[ch2int(rna[l])][r+1][j-1] + ARC[ch2int(rna[l])][i+1][l-4];
}

int arc2B_ms2(int i, int j, int l, int r, char rna[MAXSIZE]){
	if (!((i<l) && (l<r) && (r<j))) return 0;
	return arc2B1_ms2(i,j,l,r,rna) + arc2B2_ms2(i,j,l,r,rna);
}

//counts the number of neighbors obtained from shifting (i,j) to (x,j) in the internal loop
int arc2C1_ms2(int i, int j, int l, int r, char rna[MAXSIZE]){
	if (!((i<l) && (l<r) && (r<j))) return 0;
	return ARC[ch2int(rna[j])][i+1][l-1] + ARC[ch2int(rna[j])][r+1][j-4];
}

//counts the number of neighbors obtained from shifting (l,r) to (x,r) in the internal loop
int arc2C2_ms2(int i, int j, int l, int r, char rna[MAXSIZE]){
	if (!((i<l) && (l<r) && (r<j))) return 0;
	return ARC[ch2int(rna[r])][i+1][l-1] + ARC[ch2int(rna[r])][r+4][j-1];
}

int arc2C_ms2(int i, int j, int l, int r, char rna[MAXSIZE]){
	if (!((i<l) && (l<r) && (r<j))) return 0;
	return arc2C1_ms2(i,j,l,r,rna)+arc2C2_ms2(i,j,l,r,rna);
}

int arc2_ms2(int i, int j, int l, int r, char rna[MAXSIZE]){
	if (!((i<l) && (l<r) && (r<j))) return 0;
	return arc2A_ms2(i,j,l,r,rna) + arc2B_ms2(i,j,l,r,rna) + arc2C_ms2(i,j,l,r,rna);
}

//counts the number of neighbors obtained from adding or shifting a base pair in the internal loop
int arc3_ms2(int i, int j, int l, int r, char rna[MAXSIZE]){
	if (!((i<l) && (l<r) && (r<j))) return 0;
	return arc1A_ms2(i+1,l-1,rna)+arc1A_ms2(r+1,j-1,rna)+arc2_ms2(i,j,l,r,rna);
}

//counts the number of neighbors obtained from shifting (i,j) to (i,y) where i<j<y<k
int arc4_ms2(int i, int j, int k, char rna[MAXSIZE]){
	if (!((i<j) && (j<k))) return 0;
	return ARC[ch2int(rna[i])][j+1][k];
}

//counts the number of neighbors obtained from shifting (i,j) to (j,y) where i<j<y<k
int arc5_ms2(int i, int j, int k, char rna[MAXSIZE]){
	if (!((i<j) && (j<k))) return 0;
	return ARC[ch2int(rna[j])][j+4][k];
}

void McGetZB_ms2(int i,int j, char sequence[MAXSIZE], double **ZM1, double **ZM, double **ZB)
{
	int l,r;
  ZB[i][j]+=exp(-HP_Energy(i,j,S0,sequence+1)/kT);
  ZB[j][i]+=1.0;
  for(l=i+1;l<min(i+30,j-5)+1;++l)
    {
      for(r=max(l+4,j-(30-(l-i)));r<j;++r)
	{if(BP(l,r,sequence))
	    {
	      ZB[i][j]+=ZB[l][r]*exp(-IL_Energy(i,j,l,r,S0)/kT);
	      ZB[j][i]+=ZB[r][l];
	    }
	}
    }
  for(r=i+6;r<j-4;++r)
    {
      ZB[i][j]+=exp(-(ML_close+MLbasepairAndAUpenalty(j,i,S0))/kT)*ZM[i+1][r-1]*ZM1[r][j-1];
      ZB[j][i]+=ZM[r-1][i+1]*ZM1[j-1][r];
    }
}

void McGetZM1_ms2(int i, int j, char sequence[MAXSIZE], double **ZM1, double **ZB){
  int r;
  for(r=i+4;r<j+1;++r)
    if (BP(i,r,sequence)) {
      ZM1[i][j]+=ZB[i][r]*exp(-(ML_base*(j-r)+MLbasepairAndAUpenalty(i,r,S0))/kT);
      ZM1[j][i]+=ZB[r][i];
    }
}

void McGetZM_ms2(int i, int j, char sequence[MAXSIZE],double **ZM1, double **ZM)
{
int r;
  for(r=i;r<j-3;++r)
    {ZM[i][j]+=ZM1[r][j]*exp(-ML_base*(r-i)/kT);
     ZM[j][i]+=ZM1[j][r];
    }
     
  for(r=i+5;r<j-3;++r)
    {
      ZM[i][j]+=ZM[i][r-1]*ZM1[r][j];
      ZM[j][i]+=ZM[r-1][i]*ZM1[j][r];
    }
}

void McGetZ_ms2(int i, int j, char sequence[MAXSIZE],double **Z,double **ZB)
{
	int r;
  if(j-i<4)
    {
      Z[i][j]=exp(0);
      Z[j][i]=1;
    }
  else
    {
      Z[i][j]+=Z[i][j-1];
      Z[j][i]+=Z[j-1][i];
      for(r=i;r<j-3;++r)
	{ 
	  if(BP(r,j,sequence))
	    {
	      if (r==i)
		{
		  Z[i][j]+=ZB[r][j]*exp(-AU_Penalty(i,j,S0)/kT);
		  Z[j][i]+=ZB[j][r];
		}
	      else
		{
		  Z[i][j]+=Z[i][r-1]*ZB[r][j]*exp(-AU_Penalty(r,j,S0)/kT);
		  Z[j][i]+=Z[r-1][i]*ZB[j][r];
		}
	    }
	}
    }
}

void mcCaskill_ms2(char sequence[MAXSIZE], double **Z, double **ZB, double **ZM1, double **ZM) {
  int i,j,d;
  for(d=4;d<seqlen;++d)
    for(i=1;i<=seqlen-d;++i){
	j=i+d;
	  if(BP(i,j,sequence))
	    McGetZB_ms2(i,j,sequence,ZM1,ZM,ZB);
	  McGetZM1_ms2(i,j,sequence,ZM1,ZB);
	  McGetZM_ms2(i,j,sequence,ZM1,ZM);
	}
  for(d=0;d<seqlen;++d)
      for(i=1;i<=seqlen-d;++i){
	j=i+d;
	McGetZ_ms2(i,j,sequence,Z,ZB);
	}
  return;
}

void McGetQB_ms2(int i,int j, char sequence[MAXSIZE], double **QM1, double **QM, double **QB, double **ZM1, double **ZM, double **ZB){
	int l,r;

	//A[i][j]
	QB[i][j]+=exp(-HP_Energy(i,j,S0,sequence+1)/kT)*((1+arc1A_ms2(i+1,j-1,sequence)) + arc1B_ms2(i,j,sequence) + arc1C_ms2(i,j,sequence));
	QB[j][i]+=1.0;

	//B[i][j]
	for(l=i+1;l<min(i+30,j-5)+1;++l) {
		for(r=max(l+4,j-(30-(l-i)));r<j;++r){
			if(BP(l,r,sequence)) {
				QB[i][j]+=exp(-IL_Energy(i,j,l,r,S0)/kT) * \
						(QB[l][r]+ZB[l][r]*(1+arc3_ms2(i,j,l,r,sequence)));
				QB[j][i]+=QB[r][l];
			}
		}
	}

	//C[i][j]
	for(r=i+6;r<j-4;++r) {
		QB[i][j] += exp(-(ML_close+MLbasepairAndAUpenalty(j,i,S0))/kT) * \
				(ZM[i+1][r-1]*ZM1[r][j-1] + QM[i+1][r-1]*ZM1[r][j-1] + \
						ZM[i+1][r-1]*QM1[r][j-1]);
		QB[j][i]+=QM[r-1][i+1]*QM1[j-1][r];
	}
}

void McGetQM1_ms2(int i, int j, char sequence[MAXSIZE], double **QM1, double **QB, double **ZM1, double **ZB){
	int r;
	for(r=i+4;r<j+1;++r)
		if (BP(i,r,sequence)) {
			QM1[i][j] += (QB[i][r] + ZB[i][r] * arc1A_ms2(r+1,j,sequence) + arc4_ms2(i,r,j,sequence) + arc5_ms2(i,r,j,sequence) ) \
					*exp(-(ML_base*(j-r)+MLbasepairAndAUpenalty(i,r,S0))/kT);
			QM1[j][i]+=QB[r][i];
		}
}

void McGetQM_ms2(int i, int j, char sequence[MAXSIZE], double **QM1, double **QM, double **ZM1, double **ZM){
	int r;
	double MLintern =   0.04;

	//Case 1: Only one component
	for(r=i;r<j-3;++r){
		QM[i][j] += (QM1[r][j] + ZM1[r][j] * (arc1A_ms2(i,r-1,sequence) + arc1C_ms2(i-1,r,sequence))) * \
				exp(-MLintern/kT)*exp(-ML_base*(r-i)/kT); //first clause
	}
	//Case 2: Two or more components
	for(r=i+5;r<j-3;++r) {
		QM[i][j] += (QM[i][r-1]*ZM1[r][j] + ZM[i][r-1]*QM1[r][j]) * \
				exp(-MLintern/kT);
		QM[j][i]+=QM[r-1][i]*QM1[j][r];
	}
}

void McGetQ_ms2(int i, int j, double **Q, double **Z, double **QB, double **ZB, DBL_TYPE *EL, DBL_TYPE *ER1, \
		DBL_TYPE *G, DBL_TYPE *F, char *rna){
	/*---------------------------------------------------
  Q(i,j) = sum_s N(s), where N(s) is total number of base pair additions,
  removals, and shifts in s
  ---------------------------------------------------*/
	int k,x,N=strlen(rna);
	DBL_TYPE sum;
	double c1_a=0,c1_b=0,c1_c=0,c1_d=0,c2_a=0,c2_b=0,c2_c=0,c2_d=0;
	//~ DBL_TYPE t1=0,t2=0,t3=0,t4=0,t5=0,t6=0;
	//int Nz = N*N;
	//int Nc = 4*Nz;
	//int Nx = (N+1)*Nc;

	if(j-i<4)
	{
		Q[j][i]=1;
		Q[i][j]=0; 
	}
	else
	{
		sum=0;
		//1(a) j unpaired, t obtained from ms2 in [i,j-1]
		c1_a = Q[i][j-1];
		//1(b):j unpaired, t obtained by adding (k,j)
		//2(a):(k,j) exists, t obtained from ms2 in [i,k-1]
		//2(b):(k,j) exists, t obtained from ms2 in [k,j]
		//2(c):(k,j) exists, t obtained by shifting (k,j) to (k',j)
		for (k=i;k<=j-4;k++){
			if(k==i){
				if(BP(i,j,rna)){
					c1_b += Z[i+1][j-1];
					c2_a += 0;
					c2_b += QB[k][j]*exp(-AU_Penalty(i,j,S0)/kT);
				}
			}
			else{
				if(BP(k,j,rna)){
					c1_b += Z[i][k-1]*Z[k+1][j-1];
					c2_a += Q[i][k-1] * ZB[k][j]*exp(-AU_Penalty(i,j,S0)/kT);
					c2_b += Z[i][k-1] * QB[k][j]*exp(-AU_Penalty(i,j,S0)/kT);
					for(x=1;x<=k-i;x++){
						c2_c += x*F[i4(i,k-1,ch2int(rna[j]),x,N)]*ZB[k][j]*exp(-AU_Penalty(i,j,S0)/kT);
						c2_d += x*G[i4(i,k,ch2int(rna[k]),x,N)]*ZB[k][j]*exp(-AU_Penalty(i,j,S0)/kT);
						//printf("G[%d][%d][%c][%d]=%lf\n",i,k,rna[k],x,G[i4(i,k,ch2int(rna[k]),x,N)]);
					}
				}
			}
		}
		//1(c): j unpaired, t obtained by shifting (x,y) to (x,j)
		c1_c = EL[i3(i,j-1,ch2int(rna[j]),N)];
		//1(d): j unpaired, t obtained by shifting (x,y) to (y,j)
		c1_d = ER1[i3(i,j,ch2int(rna[j]),N)];
		Q[i][j]+= (c1_a + c1_b + c1_c + c1_d + c2_a + c2_b + c2_c + c2_d);
		//printf("Q[%d,%d]:%lf\t1_a:%lf\t1_b:%lf\t1_c:%lf\t1_d:%lf\t2_a:%lf\t2_b:%lf\t2_c:%lf\n2_d:%lf\n",i,j,Q[i][j], c1_a, c1_b, c1_c , c1_d , c2_a , c2_b , c2_c , c2_d);		
	}
	return;
}

void mcCaskillQ_ms2(double **Q, double **QB, double **QM1, double **QM, \
		double **Z, double **ZB, double **ZM1, double **ZM,DBL_TYPE *EL, DBL_TYPE *ER1, DBL_TYPE *G,DBL_TYPE *F,char sequence[MAXSIZE]) {
	int i,j,d;
	for(d=4;d<seqlen;++d)
		for(i=1;i<=seqlen-d;++i){
			j=i+d;
			if(BP(i,j,sequence))
				McGetQB_ms2(i,j,sequence,QM1,QM,QB,ZM1,ZM,ZB);
			McGetQM1_ms2(i,j,sequence,QM1,QB,ZM1,ZB);
			McGetQM_ms2(i,j,sequence,QM1,QM,ZM1,ZM);
		}
	for(d=0;d<seqlen;++d)
		for(i=1;i<=seqlen-d;++i){
			j=i+d;
			McGetQ_ms2(i,j,Q,Z,QB,ZB,EL,ER1,G,F,sequence);
		}
	return;
}

void computeEL(DBL_TYPE *EL, double **Z, double **ZB, char *rna){
	/*---------------------------------------------------
  EL(i,j,c) = sum_S I[ (x,y) external bp in S, x base-pairs with c ]
  ---------------------------------------------------*/
	int c,i,j,d,k;
	//int Nz = N*N;
	//int Nc = 4*Nz;
	//int Nx = (N+1)*Nc;
	int N=strlen(rna);
	DBL_TYPE sum;

	// Initialize all values to 0
	for (c=0;c<4;c++)
		for (d = 0; d < N; d++)
			for (i = 0; i < N-d; i++) {
				j = i+d;
				EL[i3(i,j,c,N)] = 0;
				//printf("El[%d][%d][%d] = %lf\n",i,j,c,EL[i3(i,j,c,N+1)]);
			}

	// Inductive Case
	for (c=0;c<4;c++)
		for (d = 4; d < N; d++)
			for (i = 1; i < N-d; i++) {
				j = i+d;
				//Case 1: j unpaired in [i,j]
				sum = EL[i3(i,j-1,c,N)];
				//Case 2: (i,j) forms a base pair
				if (BP(i,j,rna) && basePair(rna[i],c)){
					sum += ZB[i][j]*exp(-AU_Penalty(i,j,S0)/kT); 
					//printf("EL case2: ZB[%d][%d] = %lf\n",i,j,ZB[i][j]);
				}
				//Case 3: k>i and (k,j) forms a base pair
				for (k=i+1;k<=j-4;k++)
					if (BP(k,j,rna)){
						//sum += EL[i3(i,k-1,c,N)]*Z[k+1][j-1]; //Amir
						sum += EL[i3(i,k-1,c,N)]*ZB[k][j]*exp(-AU_Penalty(i,j,S0)/kT);
						//printf("EL case3-1: ZB[%d][%d] = %lf\n",k,j,ZB[k][j]);
						if (basePair(rna[k],c))
							sum += Z[i][k-1]*ZB[k][j]*exp(-AU_Penalty(i,j,S0)/kT);
							//sum += Z[i][k-1]*Z[k+1][j-1]; //Amir
							//printf("EL case3-2: ZB[%d][%d] = %lf\n",k,j,ZB[k][j]);
					}
				EL[i3(i,j,c,N)] = sum;
				//printf("EL(%d,%d,%d)=%lf\n",i,j,c,EL[i3(i,j,c,N)]);
			}
	
}

void computeER(DBL_TYPE *ER, double **Z, double **ZB, char *rna){
	/*---------------------------------------------------
  ER(i,j,c) = sum_S I[ (x,y) external bp in S, y base-pairs with c ]
  ---------------------------------------------------*/
	int c,i,j,d,k,N=strlen(rna);
	DBL_TYPE sum;
	//int Nz = N*N;
	//int Nc = 4*Nz;
	//int Nx = (N+1)*Nc;

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
			for (i=1;i<N-d;i++){
				j = i+d;
				sum = ER[i3(i,j-1,c,N)];
				if (BP(i,j,rna) && basePair(rna[j],c))
					sum += ZB[i][j]*exp(-AU_Penalty(i,j,S0)/kT);
				for (k=i+1;k<=j-4;k++)
					if (BP(k,j,rna)){
						sum += ER[i3(i,k-1,c,N)]*ZB[k][j]*exp(-AU_Penalty(i,j,S0)/kT);
						//sum += ER[i3(i,k-1,c,N)]*Z[k+1][j-1]; //Amir
						if (basePair(rna[j],c))
							sum += Z[i][k-1]*ZB[k][j]*exp(-AU_Penalty(i,j,S0)/kT);
							//sum += Z[i][k-1]*Z[k+1][j-1]; //Amir
							
					}
				ER[i3(i,j,c,N)] = sum;
				//printf("ER(%d,%d,%s)=%lf\n",i,j,c,ER[i3(i,j,c,N)]);
			}
			
}

void computeER1(DBL_TYPE *ER1, DBL_TYPE *ER, double **Z, double **ZB, char *rna){
	int c,i,j,d,k,u,N=strlen(rna);
	//int Nz = N*N;
	//int Nc = 4*Nz;
	//int Nx = (N+1)*Nc;
	DBL_TYPE sum;

	// Initialize all values to 0
	for (c=0;c<4;c++)
		for (d = 0; d < N; d++)
			for (i = 0; i<N-d; i++) {
				j = i+d;
				ER1[i3(i,j,c,N)] = 0;
			}
	//-------------------------------------------
	//Recursions
	for (c=0;c<4;c++)
		for (d=4;d<N;d++)
			for (i=1;i<N-d;i++){
				j = i+d;
				//Case 1: positions j-3,j-2,j-1,j are unpaired in S. Note rhs is ER!!
				sum = ER[i3(i,j-4,c,N)];
				//Case 2: j-4+u is paired, but j-4+u+1,...,j are unpaired in S
				for (u=1;u<=3;u++){
					//Subcase A: (i,j-4+u) is base pair in S
					if ((j-4+u-i>3) && BP(i,j-4+u,rna))
						sum += 0;   //left here for clarity
					//Subcase B: (k,j-4+u) is base pair in S, some k in [i+1,j-4+u-4]
					for (k=i+1;k<=j-4+u-4;k++)
						if ((j-4+u-k>3) && BP(k,j-4+u,rna))
							sum += ER[i3(i,k-1,c,N)]*ZB[k][j-4+u]*exp(-AU_Penalty(i,j,S0)/kT);
				}
				ER1[i3(i,j,c,N)] = sum;
				//printf("ER1(%d,%d,%s)=%lf\n",i,j,c,ER1[i3(i,j,c,N)]);
			}
}



void computeF(DBL_TYPE *F, double **Z, double **ZB,char *rna){
	/*---------------------------------------------------
  F(i,j,c,x) = sum_s I[ s has exactly x visible nt that can base-pair
               with c ]
  ---------------------------------------------------*/

	int c,i,j,d,k,x,N=strlen(rna);
	//int Nz = N*N;
	//int Nc = 4*Nz;
	//int Nx = (N+1)*Nc;
	
	//Initialize all values to 0; calloc should do this,
	for (d=0;d<N;d++)
		for (i=0;i<N-d;i++){
			j = i+d;
			for (c=0;c<4;c++)
				for (x=0;x<=N;x++)
					F[i4(i,j,c,x,N)] = 0;
		}

	//----------------------------------------------
	//Base case: define F(i,j,c,x) for i=j
	for (i=1;i<N;i++)
		for (c=0;c<4;c++)
			if (basePair(rna[i],c))
				F[i4(i,i,c,1,N)] = 1;
			else
				F[i4(i,i,c,0,N)] = 1;
	//----------------------------------------------
	//Base case: define F(i,j,c,x) for i<j<i+4
	for (i=1;i<N;i++)
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
		for (i=1;i<N-d;i++){
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
			if (BP(i,j,rna))
				for (c=0;c<4;c++)
					F[i4(i,j,c,0,N)] += ZB[i][j]*exp(-AU_Penalty(i,j,S0)/kT);
			//Case 3: (k,j) is base pair in S for k in [i+1,j-4]
			for (k=i+1;k<=j-4;k++)
				if (BP(k,j,rna))
					for (c=0;c<4;c++)
						for (x=0;x<=k-i;x++) //x in [0,k-i], where k-i num elem in [i,k-1]
							F[i4(i,j,c,x,N)]+=F[i4(i,k-1,c,x,N)]*ZB[k][j]*exp(-AU_Penalty(i,j,S0)/kT);
		}
}


void computeG(DBL_TYPE *G, DBL_TYPE *F, double **Z, double **ZB, char *rna){
	/*---------------------------------------------------
  G(i,j,c,x) = sum_s I[ s has exactly x visible nt in [i,j-4] that
               can base-pair with c, and j unpaired in s ]
  WARNING: G(i,j,c,0) defined to be zero if i <= j <= i+3
  ---------------------------------------------------*/

	int c,i,j,d,u,k,x,N=strlen(rna);
	//int Nz = N*N;
	//int Nc = 4*Nz;
	//int Nx = (N+1)*Nc;

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
		for (i=1;i<N-d;i++){
			j = i+d;                //note that num elements in [i,j] is d+1
			for (c=0;c<4;c++)
				for (x=0;x<=j-i+1-4;x++){
					//Case 1: positions j-3,j-2,j-1,j are unpaired in S
					G[i4(i,j,c,x,N)] = F[i4(i,j-4,c,x,N)];
					//Case 2: j-4+u is paired, but j-4+u+1,...,j are unpaired in S
					for (u=1;u<=3;u++){    //u in [1,3], consider 4-u unpaired in suffix
						//Subcase A: (i,j-4+u) is base pair in S
						if (x==0 && j-4+u-i>3 && BP(i,j-4+u,rna))
							G[i4(i,j,c,x,N)] += ZB[i][j-4+u]*exp(-AU_Penalty(i,j,S0)/kT);
						//Subcase B: (k,j-4+u) is base pair in S, some k in [i+1,j-4+u-4]
						for (k=i+1;k<=j-4+u-4;k++) //k in [i+1,(j-4+u)-4]
							if (j-4+u-k>3 && BP(k,j-4+u,rna))
								G[i4(i,j,c,x,N)]+=F[i4(i,k-1,c,x,N)] * \
								ZB[k][j-4+u]*exp(-AU_Penalty(i,j,S0)/kT);
					}
				//printf("\tG(%d,%d,%c,%d)=%.Lf\n",i,j,int2c(c),x,G[i4(i,j,c,x,N)]);
				}
		}
}
