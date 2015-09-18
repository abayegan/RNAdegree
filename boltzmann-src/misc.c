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
#include"misc.h"
#include<limits.h>

#ifdef _WIN32
	#include <windows>
#elif defined(__APPLE__) || defined(__linux)  || defined(__unix)  || defined(__posix) 
	#include <unistd.h>
	#include <limits.h>
#endif


extern int ARC[4][MAXSIZE][MAXSIZE];


void computeARC( int ARC[4][MAXSIZE][MAXSIZE], char rna[MAXSIZE]){
	//indices start at 1 up to seqLen
	//ARC['G',i,j]=number of k where rna[k] in ['C','U']
	int i,j,d,alpha;
	int seqLen=strlen(rna); //rna starts with @, so indices 1...,seqLen-1

	//Set all entries to zero -- unnecessary if use calloc
	for (alpha=0; alpha<4; alpha++)
		for (i=1; i<seqLen; i++)
			for (j=1; j<seqLen; j++)
				ARC[alpha][i][j]=0;

	//Base case: i=j
	for (i=1;i<seqLen;i++){
		if (rna[i]=='A')
			ARC[3][i][i]+=1; //ARC['U',i,i] += 1
		else if (rna[i]=='C')
			ARC[2][i][i]+=1; //ARC['G',i,i] += 1
		else if (rna[i]=='G'){
			ARC[1][i][i]+=1; //ARC['C',i,i] += 1
			ARC[3][i][i]+=1; //ARC['U',i,i] += 1
		}
		else if (rna[i]=='U'){
			ARC[0][i][i]+=1; //ARC['A',i,i] += 1
			ARC[2][i][i]+=1; //ARC['G',i,i] += 1
		}
	}
	//Inductive case: i<j
	for (d=1;d<seqLen;d++)
		for (i=1;i<seqLen-1;i++){
			j = i+d;
			if (j>=seqLen) break;
			for (alpha=0;alpha<4;alpha++)
				ARC[alpha][i][j]+=ARC[alpha][i][j-1];
			if (rna[j]=='A')
				ARC[3][i][j]+=1; //ARC['U',i,j]+=1
			else if (rna[j]=='C')
				ARC[2][i][j]+=1; //ARC['G',i,j]+=1
			else if (rna[j]=='G'){
				ARC[1][i][j]+=1; //ARC['C',i,j]+=1
				ARC[3][i][j]+=1; //ARC['U',i,i]+=1
			}
			else if (rna[j]=='U'){
				ARC[0][i][j]+=1; //ARC['A',i,j]+=1
				ARC[2][i][j]+=1; //ARC['G',i,j]+=1
			}
		}
}

int ch2int(char ch){
	ch = toupper(ch);
	if (ch=='A')
		return 0;
	if (ch=='C')
		return 1;
	if (ch=='G')
		return 2;
	if (ch=='U')
		return 3;
	printf("Error: ch2int(%c)\n",ch);
}

int basePair(char x, int c){
  //c codes a nt
  //0->A 1->C 2->G 3->U
  int wc; //watson crick pair
  //assert( ((0<=c)&&(c<4)) );
  x = toupper(x);
  wc =((x=='A'&&c==3)||(x=='U'&&c==0)||(x=='C'&&c==2)||(x=='G'&&c==1));
  if (wc || (x=='G'&&c==3)|| (x=='U'&&c==2))
    return 1;
  else
    return 0;
  }

/*The following function checks if the ith and jth positions in the sequence could base pair*/
int BP(int i,int j,char sequence[MAXSIZE])  //safe
{if(j-i<4)
    return 0;
  else if ((sequence[i]=='A'&& sequence[j]=='U')||(sequence[i]=='U'&& sequence[j]=='A'))
    return 1;
  else if ((sequence[i]=='C'&& sequence[j]=='G')||(sequence[i]=='G'&& sequence[j]=='C'))
    return 1;
  else if ((sequence[i]=='U'&& sequence[j]=='G')||(sequence[i]=='G'&& sequence[j]=='U'))
    return 1;
  else
    return 0;
}

int ElementInList(char subshape[SHAPELENGTH],char shapeList[SHAPELENGTH][SHAPELENGTH],int listlen) //check if a subshape is in list, if yes,return the index of the subshape,if not return -1
{int i;
 for(i=0;i<listlen;++i)
    {if(strcmp(shapeList[i],subshape)==0)
        return i;
    }
 return -1;
}

int NumInList(int a, int L[SHAPELENGTH],int listlen) //check if a number is in this list,if not return -1
{int i;
 for(i=0;i<listlen;++i)
    {if(a==L[i])
        return i;
    }
 return -1;
}

/*This function returns the largest value in an array.*/
double MaxInArray(double array[],int arraylen)   //safe
{ double largest=array[0];
  int i;
  for (i=0;i<arraylen;++i)
    {if (array[i]>largest)
        largest=array[i];
    }
  return largest;
}

/*The following function check if the input string is a valid RNA sequence*/
int CheckSequence(char sequence[MAXSIZE]){     //safe
  int i;
  for (i=1;i<strlen(sequence);++i)
    {
     if (toupper(sequence[i])!='A' && toupper(sequence[i])!='U' && toupper(sequence[i])!='C' && toupper(sequence[i])!='G'&& toupper(sequence[i])!='T') //check if there are invalid characters
        {
         printf("The input string should only contain A,U,T,C,G!\n");
         exit(1);
        }
     else if (toupper(sequence[i])=='T') //change T to U
        sequence[i]='U';
     else
        sequence[i]=toupper(sequence[i]); // change lower case to upper case
    }
  return 0;
}

double ***Allocate3DMatrix(int a, int b, int c)
// allocate a matrix of size a*b*c
{ int i;
  int j;
  double ***Matrix;
  Matrix=(double ***) malloc((a) * sizeof (double **));
  if(Matrix == NULL)
        {printf("out of memory\n");
         exit(1);}
  for (i=0; i<a; i++){
    Matrix[i]=(double **) malloc((b)*sizeof(double *));
    if(Matrix[i] == NULL)
        {printf("out of memory\n");
         exit(1);}
  }
  for(i=0;i<a;i++)
    {for(j=0;j<b;j++)
        {
        Matrix[i][j]=(double *) malloc((c)*sizeof(double));
        if(Matrix[i][j] == NULL)
        {printf("out of memory\n");
         exit(1);}
        }
    }
  return Matrix;
}

int Free3DMatrix(double ***Matrix, int a, int b, int c)
// Free a matrix of size a*b*c
{ int i;
  int j;
  for(i=0;i<a;i++)
	{
	for(j=0;j<b;j++)
	free(Matrix[i][j]);
	}
  for(i=0;i<a;i++)
	free(Matrix[i]);
  free(Matrix);
  return 0;
}

double **Allocate2DMatrix(int a, int b)
// allocate a matrix of size a*b
{ int i;
  double **Matrix;
  Matrix=(double **) malloc((a) * sizeof (double *));
  if(Matrix == NULL)
        {printf("out of memory\n");
         exit(1);}
  for (i=0; i<a; i++){
    Matrix[i]=(double *) malloc((b)*sizeof(double));
    if(Matrix[i] == NULL)
        {printf("out of memory\n");
         exit(1);}
  }
  return Matrix;
}

int Free2DMatrix(double **Matrix,int a, int b)
// free a matrix of size a*b
{ int i;
  for (i=0; i<a; i++){
    free(Matrix[i]);
  }
  free(Matrix);
  return 0;
}

void displayReadMe()
{
	int c;
	FILE *readMe;
	readMe = fopen("README.txt", "r");
	if (readMe) {
		while ((c = getc(readMe)) != EOF)
			putchar(c);
		fclose(readMe);
	}
	else
		printf ("No README.txt file to display\n");
	exit(1);
}


void dispUsage(char * name)
{
	printf("\nUSAGE: %s rnaSequence [-ms1|-ms2] [-t temperature] [-e turner99|turner04|andronescu07] [-v] [-h]\n", name);
	exit(1);
}

char* getExecPath(char* argv0){
	char* exec_path=(char*) malloc (sizeof(char)*PATH_MAX);
	exec_path[0]='.';
	exec_path[1]='\0';
	int path_length=1;
	FILE* file;
	#if defined(_WIN32)
	// GetModuleFileName
	//_pgmptr
	#elif defined(__APPLE__) || defined(__linux)  || defined(__unix)  || defined(__posix) 
	char buff[PATH_MAX];
	int bufsize = PATH_MAX-1;
	if(file = fopen("/proc/self/exe", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/self/exe", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}		
	}
	else if(file = fopen("/proc/curproc/file", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/curproc/file", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}		
	}
	else if(file = fopen("/proc/self/path/a.out", "r")){		
		fclose(file);
		ssize_t len = readlink("/proc/self/path/a.out", buff, bufsize);
		if (len != -1) {
			buff[len] = '\0';
			strcpy(exec_path, buff);
			path_length=len;
		}		
	}
	else{
		exec_path= argv0;
	}
	char* slash_pointer= strrchr(exec_path,'/');
	int slash_position = (int)(slash_pointer - exec_path);
	if(slash_position != path_length){			
		exec_path[slash_position]='\0';
	}
	else{
		exec_path[0]='.';
		exec_path[1]='\0';
	}

	#endif
	

	return exec_path;
}





