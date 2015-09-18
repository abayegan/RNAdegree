#define int2c(x) ((x)==0 ? 'A' : ((x)==1 ? 'C' : ((x)==2 ? 'G' : 'U')))
#define i2(i,j,N) ((i)*(N)+(j))
#define i3(i,j,c,N) ((i)*(N)+(j)+(c)*(N)*(N))
#define i4(i,j,c,x,N) ((i)*(N)+(j)+(c)*(N)*(N)+(x)*4*(N)*(N))
typedef long double DBL_TYPE;
void computeARC( int ARC[4][MAXSIZE][MAXSIZE], char rna[MAXSIZE]);
int ElementInList(char subshape[SHAPELENGTH],char shapeList[SHAPELENGTH][SHAPELENGTH],int listlen);
int NumInList(int a, int L[SHAPELENGTH],int listlen);
double MaxInArray(double array[],int arraylen);
int BP(int i,int j,char sequence[MAXSIZE]);
double ***Allocate3DMatrix(int a, int b, int c);
double **Allocate2DMatrix(int a, int b);
int CheckSequence(char sequence[MAXSIZE]);
int Free3DMatrix(double ***Matrix, int a, int b, int c);
int Free2DMatrix(double **Matrix,int a, int b);
int ch2int(char ch);
int basePair(char x,int c);
void displayReadMe();
void dispUsage(char * name);
char* getExecPath(char* );
