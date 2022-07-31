#ifndef __PASE_MATVEC_H__
#define __PASE_MATVEC_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdbool.h>
#include <float.h>

#define CLK_TCK 1000000

typedef struct pase_Matrix_struct
{
    int num_aux_vec;
    /* N_H阶的并行矩阵 */
    void *A_H;
    /* num_aux_vec个N_H阶的并行向量 */
    void **aux_Hh;
    /* num_aux_vec*num_aux_vec的数组 */
    double *aux_hh;
    void *factorization; //存分解后的矩阵
    void **AH_inv_aux_Hh;  //存A_H^-1*aux_Hh
    double *aux_hh_inv;
    bool if_sym;
    bool is_diag;
} pase_Matrix;
typedef struct pase_Matrix_struct *PASE_Matrix;

typedef struct pase_Vector_struct
{
    int num_aux_vec;
    /* N_H阶的并行向量 */
    void *b_H;
    /* num_aux_vec的数组 */
    double *aux_h;
    double *aux_h_tmp;

} pase_Vector;
typedef struct pase_Vector_struct *PASE_Vector;

typedef struct pase_MultiVector_struct
{
    int num_aux_vec;
    int num_vec;
    /* num_vec个N_H阶的并行向量 */
    void **b_H;
    /* num_aux_vec*num_vec的数组 */
    double *aux_h;
    double *aux_h_tmp;
} pase_MultiVector;
typedef struct pase_MultiVector_struct *PASE_MultiVector;

void Mat_Output(void *A, char *name);
void MultiVec_Output(void **MV, int vec_num, char *name);
void PASE_Mat_Output(PASE_Matrix A, char *name);
void PASE_MultiVec_Output(PASE_MultiVector MV, int vec_num, char *name);


#endif