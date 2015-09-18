/*
 * aux_Turner.h
 *
 *  Created on: May 10, 2015
 *      Author: Amir Bayegan
 */
#ifndef AUX_TURNER_H_
#define AUX_TURNER_H_
void McGetZB_ms2(int ,int , char sequence[MAXSIZE], double **, double **, double **);
void McGetZM1_ms2(int i, int j, char sequence[MAXSIZE], double **ZM1, double **ZB);
void McGetZM_ms2(int i, int j, char sequence[MAXSIZE],double **ZM1, double **ZM);
void McGetZ_ms2(int i, int j, char sequence[MAXSIZE],double **Z,double **ZB);
void McGetQB_ms2(int i,int j, char sequence[MAXSIZE], double **QM1, double **QM, double **QB, double **ZM1, double **ZM, double **ZB);
void McGetQM1_ms2(int i, int j, char sequence[MAXSIZE], double **QM1, double **QB, double **ZM1, double **ZB);
void McGetQM_ms2(int i, int j, char sequence[MAXSIZE], double **QM1, double **QM, double **ZM1, double **ZM);
void mcCaskillQ_ms2(double **Q, double **QB, double **QM1, double **QM, \
		double **Z, double **ZB, double **ZM1, double **ZM,DBL_TYPE *EL, DBL_TYPE *ER1, DBL_TYPE *G, DBL_TYPE *F, char sequence[MAXSIZE]);
void McGetQ_ms2(int i, int j, double **Q, double **Z, double **QB, double **ZB, DBL_TYPE *EL, DBL_TYPE *ER1, \
		DBL_TYPE *G, DBL_TYPE *F, char *rna);

int arc1A_ms2(int i, int j, char rna[MAXSIZE]);
int arc1B_ms2(int i, int j, char rna[MAXSIZE]);
int arc1C_ms2(int i, int j, char rna[MAXSIZE]);
int arc2A_ms2(int i, int j, int l, int r, char rna[MAXSIZE]);
int arc2B1_ms2(int i, int j, int l, int r, char rna[MAXSIZE]);
int arc2B2_ms2(int i, int j, int l, int r, char rna[MAXSIZE]);
int arc2B_ms2(int i, int j, int l, int r, char rna[MAXSIZE]);
int arc2C1_ms2(int i, int j, int l, int r, char rna[MAXSIZE]);
int arc2C2_ms2(int i, int j, int l, int r, char rna[MAXSIZE]);
int arc2C_ms2(int i, int j, int l, int r, char rna[MAXSIZE]);
int arc2_ms2(int i, int j, int l, int r, char rna[MAXSIZE]);
int arc3_ms2(int i, int j, int l, int r, char rna[MAXSIZE]);
int arc4_ms2(int i, int j, int k, char rna[MAXSIZE]);
int arc5_ms2(int i, int j, int k, char rna[MAXSIZE]);


#endif /* AUX_TURNER_H_ */
