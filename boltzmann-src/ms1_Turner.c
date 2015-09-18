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
//#include"RNAhairpin.h"
#include"misc.h"

/* P. Clote, Aug 2014 */


/*----------------------------------------
Following from Vienna 1.5 beta, except multiplied by 1/100.
WARNING: These are global variables.
----------------------------------------*/
//PUBLIC double ML_BASE37 = 0;
//PUBLIC double ML_closing37 = 3.40;
//PUBLIC double ML_intern37 =  0.40;
//PUBLIC double MLintern =  0.40;


/*----------------------------------------
WARNING: Arrays in this file (Turner energy) are 1-indexed,
while those for the uniform and Nussinov case are 0-indexed.
----------------------------------------*/

extern int ARC[4][MAXSIZE][MAXSIZE];

typedef struct PartitionFunction {
   double **Z;
   double **ZB;
   double **ZM1;
   double **ZM;
} PartitionFunction;

//arc1,arc2,arc3 assume that sequence is 1-indexed
//S0 and sequence are 1-indexed; however, sequence is 0-indexed in fold.c
int arc1(int i, int j, char sequence[MAXSIZE]){
  int x,y,sum=0;
  if (i>j) return 0;
  for (x=i;x<=j-4;x++)
	sum += ARC[ch2int(sequence[x])][x+4][j];
    //for (y=x+4;y<=j;y++)
      //if (BP(x,y,sequence)) sum++;
  return sum;
  }

int arc2(int i, int j, int l, int r, char sequence[MAXSIZE]){
  int x,y,sum=0;
  if (i>j || l>r) return 0;
  for (x=i+1;x<=l-1;x++)
	sum += ARC[ch2int(sequence[x])][r+1][j-1];
    //for (y=r+1;y<=j-1;y++)
      //if (BP(x,y,sequence)) sum++;
  return sum;
  }

int arc3(int i, int j, int l, int r, char sequence[MAXSIZE]){
  return (arc1(i+1,l-1,sequence)+arc1(r+1,j-1,sequence)+arc2(i,j,l,r,sequence));
  }


void mcCaskill(char sequence[MAXSIZE], double **Z, double **ZB, double **ZM1, double **ZM) {
  int i,j,d;
  for(d=4;d<seqlen;++d)
    for(i=1;i<=seqlen-d;++i){
	j=i+d;
	  if(BP(i,j,sequence))
	    McGetZB(i,j,sequence,ZM1,ZM,ZB);
	  McGetZM1(i,j,sequence,ZM1,ZB);
	  McGetZM(i,j,sequence,ZM1,ZM);
	}
  for(d=0;d<seqlen;++d)
      for(i=1;i<=seqlen-d;++i){
	j=i+d;
	McGetZ(i,j,sequence,Z,ZB);
	}
  return;
}

void mcCaskillQ(char sequence[MAXSIZE], double **Q, double **QB, double **QM1, double **QM, double **Z, double **ZB, double **ZM1, double **ZM) {
  int i,j,d;
  for(d=4;d<seqlen;++d)
    for(i=1;i<=seqlen-d;++i){
	j=i+d;
	  if(BP(i,j,sequence))
	    McGetQB(i,j,sequence,QM1,QM,QB,ZM1,ZM,ZB); 
	  McGetQM1(i,j,sequence,QM1,QB,ZM1,ZB);
	  McGetQM(i,j,sequence,QM1,QM,ZM1,ZM);
	}
  for(d=0;d<seqlen;++d)
     for(i=1;i<=seqlen-d;++i){
	j=i+d;
	McGetQ(i,j,sequence,Q,QB,Z,ZB);
	}
  return;
}

#if 0
PartitionFunction *mcCaskill(char sequence[MAXSIZE]) {
  printf("hi\n");
  static PartitionFunction *pf;
  int i,j,d;
  double ** McZ;
  double ** McZB;
  double ** McZM;
  double ** McZM1;
  McZ=Allocate2DMatrix( seqlen+1,seqlen+1);
  McZB=Allocate2DMatrix( seqlen+1,seqlen+1);
  McZM=Allocate2DMatrix( seqlen+1,seqlen+1);
  McZM1=Allocate2DMatrix( seqlen+1,seqlen+1);
  for(i=0;i<seqlen+1;++i)
    {for (j=0;j<seqlen+1;++j)
	{
	  McZB[i][j]=0;
	  McZ[i][j]=0;
	  McZM[i][j]=0;
	  McZM1[i][j]=0;
	}
     }
  for(d=4;d<seqlen;++d)
    for(i=1;i<=seqlen-d;++i){
	j=i+d;
	  if(BP(i,j,sequence))
	    McGetZB(i,j,sequence,McZM1,McZM,McZB);
	  McGetZM1(i,j,sequence,McZM1,McZB);
	  McGetZM(i,j,sequence,McZM1,McZM);
	}
  for(d=0;d<seqlen;++d)
      for(i=1;i<=seqlen-d;++i){
	j=i+d;
	McGetZ(i,j,sequence,McZ,McZB);
	}
  pf->Z   = McZ; pf->ZB  = McZB; pf->ZM1 = McZM1; pf->ZM  = McZM;
  return pf;
}
#endif


double** runMcCaskill(char sequence[MAXSIZE]) {
  int i,j,d;
  double ** McZ;
  double ** McZB;
  double ** McZM;
  double ** McZM1;
  McZ=Allocate2DMatrix( seqlen+1,seqlen+1);
  McZB=Allocate2DMatrix( seqlen+1,seqlen+1);
  McZM=Allocate2DMatrix( seqlen+1,seqlen+1);
  McZM1=Allocate2DMatrix( seqlen+1,seqlen+1);
  for(i=0;i<seqlen+1;++i)
    {for (j=0;j<seqlen+1;++j)
	{
	  McZB[i][j]=0;
	  McZ[i][j]=0;
	  McZM[i][j]=0;
	  McZM1[i][j]=0;
	}
     }
  for(d=4;d<seqlen;++d)
    {for(i=1;i<=seqlen-d;++i)
	{j=i+d;
	  if(BP(i,j,sequence))
	    McGetZB(i,j,sequence,McZM1,McZM,McZB);
	  McGetZM1(i,j,sequence,McZM1,McZB);
	  McGetZM(i,j,sequence,McZM1,McZM);
	}
    }
  for(d=0;d<seqlen;++d)
    {
      for(i=1;i<=seqlen-d;++i)
	{
	  j=i+d;
	  McGetZ(i,j,sequence,McZ,McZB);
	}
    }
  return McZ;
}

int McGetZB(int i,int j, char sequence[MAXSIZE], double **ZM1, double **ZM, double **ZB)
{ int l,r;
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

int McGetQB(int i,int j, char sequence[MAXSIZE], double **QM1, double **QM, double **QB, double **ZM1, double **ZM, double **ZB){
  int l,r;
  QB[i][j]+=exp(-HP_Energy(i,j,S0,sequence+1)/kT)*(1+arc1(i+1,j-1,sequence));
  QB[j][i]+=1.0;
  for(l=i+1;l<min(i+30,j-5)+1;++l) {
      for(r=max(l+4,j-(30-(l-i)));r<j;++r){
	if(BP(l,r,sequence)) {
	      QB[i][j]+=exp(-IL_Energy(i,j,l,r,S0)/kT) * \
                        (QB[l][r]+ZB[l][r]*(1+arc3(i,j,l,r,sequence)));
	      QB[j][i]+=QB[r][l];
	    }
	}
    }
  for(r=i+6;r<j-4;++r) {
      QB[i][j]+=exp(-(ML_close+MLbasepairAndAUpenalty(j,i,S0))/kT) * \
                (ZM[i+1][r-1]*ZM1[r][j-1] + QM[i+1][r-1]*ZM1[r][j-1] + \
                ZM[i+1][r-1]*QM1[r][j-1]);
      QB[j][i]+=QM[r-1][i+1]*QM1[j-1][r];
    }
}

int McGetZM1(int i, int j, char sequence[MAXSIZE], double **ZM1, double **ZB){
  int r;
  for(r=i+4;r<j+1;++r)
    if (BP(i,r,sequence)) {
      ZM1[i][j]+=ZB[i][r]*exp(-(ML_base*(j-r)+MLbasepairAndAUpenalty(i,r,S0))/kT);
      ZM1[j][i]+=ZB[r][i];
    }
}

int McGetQM1(int i, int j, char sequence[MAXSIZE], double **QM1, double **QB, double **ZM1, double **ZB){
  int r;
  for(r=i+4;r<j+1;++r)
    if (BP(i,r,sequence)) {
      QM1[i][j]+=(QB[i][r]+ZB[i][r]*arc1(r+1,j,sequence))*exp(-(ML_base*(j-r)+MLbasepairAndAUpenalty(i,r,S0))/kT);
      QM1[j][i]+=QB[r][i];
    }
}

int McGetZM(int i, int j, char sequence[MAXSIZE],double **ZM1, double **ZM)
{ int r;
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

int McGetQM(int i, int j, char sequence[MAXSIZE], double **QM1, double **QM, double **ZM1, double **ZM){
  int r;
  double MLintern =  0.04; 

  //Case 1: Only one component
  for(r=i;r<j-3;++r){
    QM[i][j] += (QM1[r][j] + ZM1[r][j]*arc1(i,r-1,sequence)) * \
                exp(-MLintern/kT)*exp(-ML_base*(r-i)/kT); //first clause
    }
     
  //Case 2: Two or more components
  for(r=i+5;r<j-3;++r) {
    QM[i][j] += (QM[i][r-1]*ZM1[r][j] + ZM[i][r-1]*QM1[r][j]) * \
                exp(-MLintern/kT);
      QM[j][i]+=QM[r-1][i]*QM1[j][r];
    }
}

int McGetZ(int i, int j, char sequence[MAXSIZE],double **Z,double **ZB)
{ int r;
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


int McGetQ(int i, int j, char sequence[MAXSIZE], double **Q,double **QB, double **Z,double **ZB){
  int r;
  double val;
  if(j-i<4) {
      Q[j][i]=1;
      Q[i][j]=0; 
      //very tricky bug in previous version where these two statements
      //were reversed!!
    }
  else {
      Q[i][j]+=Q[i][j-1];
      Q[j][i]+=Q[j-1][i];
      for(r=i;r<j-3;++r) { 
	  if(BP(r,j,sequence)) {
	      if (r==i) {
                  val  = Z[i+1][j-1];
                  val += QB[i][j]*exp(-AU_Penalty(i,j,S0)/kT);
                  Q[i][j]+=val;
		  Q[j][i]+=QB[j][r];
		}
	      else {
                  val  = Q[i][r-1]*ZB[r][j]+Z[i][r-1]*QB[r][j];
                  val  = val*exp(-AU_Penalty(r,j,S0)/kT);
		  val += Z[i][r-1]*Z[r+1][j-1];
		  Q[i][j]+=val;
		  Q[j][i]+=Q[r-1][i]*QB[j][r];
		}
	    }
	}
    }

}


