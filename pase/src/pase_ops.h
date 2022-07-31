#ifndef _PASE_OPS_H_
#define _PASE_OPS_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "ops.h"
#include "pase_matvec.h"

#if OPS_USE_MPI
#include <mpi.h>
#endif

#if OPS_USE_SLEPC
#include "app_slepc.h"
#endif

typedef struct PASE_OPS_
{
    void (*Printf)(const char *fmt, ...);
    double (*GetWtime)(void);
    int (*GetOptionFromCommandLine)(
        const char *name, char type, void *data,
        int argc, char *argv[], struct OPS_ *ops);

    /* ------------- mat -------------- */
    void (*MatView)(void *mat, struct PASE_OPS_ *pase_ops);
    // Y = alpha X + beta Y
    void (*MatAxpby)(double alpha, void *matX, double beta, void *matY, struct PASE_OPS_ *pase_ops);

    /* -------------- vec -------------- */
    void (*VecCreateByMat)(void **des_vec, void *src_mat, struct PASE_OPS_ *pase_ops);
    void (*VecCreateByVec)(void **des_vec, void *src_vec, struct PASE_OPS_ *ops);
    void (*VecDestroy)(void **des_vec, struct PASE_OPS_ *pase_ops);
    void (*VecView)(void *x, struct PASE_OPS_ *pase_ops);
    // inner_prod = x'y
    void (*VecInnerProd)(void *x, void *y, double *inner_prod, struct PASE_OPS_ *pase_ops);
    // inner_prod = (x->b_H)'(y->b_H) for each proc
    void (*VecLocalInnerProd)(void *x, void *y, double *inner_prod, struct PASE_OPS_ *pase_ops);
    void (*VecSetRandomValue)(void *x, struct PASE_OPS_ *pase_ops);
    // y = alpha x + beta y
    void (*VecAxpby)(double alpha, void *x, double beta, void *y, struct PASE_OPS_ *pase_ops);
    // y = mat * x
    void (*MatDotVec)(void *mat, void *x, void *y, struct PASE_OPS_ *pase_ops);
    // y = mat' * x
    void (*MatTransDotVec)(void *mat, void *x, void *y, struct PASE_OPS_ *pase_ops);

    /* -------------- multivec -------------- */
    void (*MultiVecCreateByMat)(void ***multi_vec, int num_vec, void *src_mat, struct PASE_OPS_ *pase_ops);
    void (*MultiVecCreateByVec)(void ***multi_vec, int num_vec, void *src_vec, struct PASE_OPS_ *pase_ops);
    void (*MultiVecCreateByMultiVec)(void ***multi_vec, int num_vec, void **src_mv, struct PASE_OPS_ *pase_ops);
    void (*MultiVecDestroy)(void ***multi_vec, int num_vec, struct PASE_OPS_ *pase_ops);
    void (*GetVecFromMultiVec)(void **multi_vec, int col, void **vec, struct PASE_OPS_ *pase_ops);
    // *vec = multi_vec[col]
    void (*RestoreVecForMultiVec)(void **multi_vec, int col, void **vec, struct PASE_OPS_ *pase_ops);
    void (*MultiVecView)(void **x, int start, int end, struct PASE_OPS_ *pase_ops);
    // inner_prod = (x->b_H)'(y->b_H) for each proc
    void (*MultiVecLocalInnerProd)(char nsdIP, void **x, void **y, int is_vec, int *start, int *end,
                                   double *inner_prod, int ldIP, struct PASE_OPS_ *pase_ops);
    // inner_prod = x'y
    void (*MultiVecInnerProd)(char nsdIP, void **x, void **y, int is_vec, int *start, int *end,
                              double *inner_prod, int ldIP, struct PASE_OPS_ *pase_ops);
    void (*MultiVecSetRandomValue)(void **multi_vec, int start, int end, struct PASE_OPS_ *pase_ops);
    void (*MultiVecAxpby)(double alpha, void **x, double beta, void **y,
                          int *start, int *end, struct PASE_OPS_ *pase_ops);
    /* y = x coef + y diag(beta) , coef & beta : sequantial*/
    void (*MultiVecLinearComb)(void **x, void **y, int is_vec, int *start, int *end,
                               double *coef, int ldc, double *beta, int incb, struct PASE_OPS_ *pase_ops);

    // mat * x = y
    void (*MatDotMultiVec)(void *mat, void **x, void **y,
                           int *start, int *end, struct PASE_OPS_ *pase_ops);
    // mat^T * x = y
    void (*MatTransDotMultiVec)(void *mat, void **x, void **y,
                                int *start, int *end, struct PASE_OPS_ *pase_ops);

    /* -------------- other -------------- */
    // qAp = Qt A P
    void (*MultiVecQtAP)(char ntsA, char ntsdQAP,
                         void **mvQ, void *matA, void **mvP, int is_vec,
                         int *start, int *end, double *qAp, int ldQAP,
                         void **mv_ws, struct PASE_OPS_ *pase_ops);

    /* -------------- Dense matrix vector ops -------------- */
    struct OPS_ *lapack_ops;
    void (*DenseMatQtAP)(char ntluA, char nsdC,
                         int nrowsA, int ncolsA,
                         int nrowsC, int ncolsC,
                         double alpha, double *matQ, int ldQ,
                         double *matA, int ldA,
                         double *matP, int ldP,
                         double beta, double *matC, int ldC,
                         double *dbl_ws);
    void (*DenseMatOrth)(double *mat, int nrows, int ldm,
                         int start, int *end, double orth_zero_tol,
                         double *dbl_ws, int length, int *int_ws);

    /*TO DO*/
    /* linear solver */
    void (*LinearSolver)(void *mat, void *b, void *x,
                         struct OPS_ *ops);
    void *linear_solver_workspace;
    void (*MultiLinearSolver)(void *mat, void **b, void **x,
                              int *start, int *end, struct OPS_ *ops);
    void *multi_linear_solver_workspace;

    // orthonormal
    void (*MultiVecOrth)(void **x, int start_x, int *end_x,
                         void *B, struct PASE_OPS_ *pase_ops);
    void *orth_workspace;

    /* multi grid */
    /* get multigrid operator for num_levels = 4
     * P0     P1       P2
     * A0     A1       A2        A3
     * B0  P0'B0P0  P1'B1P1   P2'B2P2
     * A0 is the original matrix */
    void (*MultiGridCreate)(void ***A_array, void ***B_array, void ***P_array,
                            int *num_levels, void *A, void *B, struct OPS_ *ops);
    /* free A1 A2 A3 B1 B2 B3 P0 P1 P2
     * A0 and B0 are just pointers */
    void (*MultiGridDestroy)(void ***A_array, void ***B_array, void ***P_array,
                             int *num_levels, struct OPS_ *ops);
    void (*VecFromItoJ)(void **P_array, int level_i, int level_j,
                        void *vec_i, void *vec_j, void **vec_ws, struct OPS_ *ops);
    void (*MultiVecFromItoJ)(void **P_array, int level_i, int level_j,
                             void **multi_vec_i, void **multi_vec_j, int *startIJ, int *endIJ,
                             void ***multi_vec_ws, struct OPS_ *ops);

    /* eigen solver */
    void (*EigenSolver)(void *A, void *B, double *eval, void **evec,
                        int nevGiven, int *nevConv, struct OPS_ *ops);
    void *eigen_solver_workspace;

    /* for pas */
    OPS *gcge_ops;
} PASE_OPS;

void PASE_OPS_Create(PASE_OPS **pase_ops, OPS *ops);
void PASE_OPS_Destroy(PASE_OPS **pase_ops);

void PASE_MatrixCreate(PASE_Matrix *pase_matrix,
                       int num_aux_vec, void *A_H, PASE_OPS *pase_ops);
void PASE_MatrixDestroy(PASE_Matrix *matrix, PASE_OPS *pase_ops);

void PASE_DefaultMatView(void *mat, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMatAxpby(double alpha, void *matX, double beta, void *matY, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMatDotMultiVec(void *mat, void **x, void **y,
                                int *start, int *end, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMatTransDotMultiVec(void *mat, void **x, void **y,
                                     int *start, int *end, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMatDotVec(void *mat, void *x, void *y, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMatTransDotVec(void *mat, void *x, void *y, struct PASE_OPS_ *pase_ops);

void PASE_DefaultVecCreateByMat(void **des_vec, void *src_mat, struct PASE_OPS_ *pase_ops);
void PASE_DefaultVecCreateByVec(void **des_vec, void *src_vec, struct PASE_OPS_ *pase_ops);
void PASE_DefaultVecDestroy(void **des_vec, struct PASE_OPS_ *pase_ops);
void PASE_DefaultVecView(void *x, struct PASE_OPS_ *pase_ops);
void PASE_DefaultVecAxpby(double alpha, void *x, double beta, void *y, struct PASE_OPS_ *pase_ops);
void PASE_DefaultVecInnerProd(void *x, void *y, double *inner_prod, struct PASE_OPS_ *pase_ops);
void PASE_DefaultVecLocalInnerProd(void *x, void *y, double *inner_prod, struct PASE_OPS_ *pase_ops);
void PASE_DefaultVecSetRandomValue(void *x, struct PASE_OPS_ *pase_ops);

void PASE_DefaultMultiVecCreateByMat(void ***multi_vec, int num_vec, void *src_mat, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMultiVecCreateByVec(void ***multi_vec, int num_vec, void *src_vec, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMultiVecCreateByMultiVec(void ***multi_vec, int num_vec, void **src_mv, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMultiVecDestroy(void ***multi_vec, int num_vec, struct PASE_OPS_ *pase_ops);
void PASE_DefaultGetVecFromMultiVec(void **multi_vec, int col, void **vec, struct PASE_OPS_ *pase_ops);
void PASE_DefaultRestoreVecForMultiVec(void **multi_vec, int col, void **vec, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMultiVecView(void **x, int start, int end, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMultiVecLocalInnerProd(char nsdIP, void **x, void **y, int is_vec, int *start, int *end,
                                        double *inner_prod, int ldIP, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMultiVecInnerProd(char nsdIP, void **x, void **y, int is_vec, int *start, int *end,
                                   double *inner_prod, int ldIP, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMultiVecSetRandomValue(void **multi_vec, int start, int end, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMultiVecAxpby(double alpha, void **x, double beta, void **y,
                               int *start, int *end, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMultiVecLinearComb(void **x, void **y, int is_vec, int *start, int *end,
                                    double *coef, int ldc, double *beta, int incb, struct PASE_OPS_ *pase_ops);

void PASE_DefaultMultiVecQtAP(char ntsA, char ntsdQAP,
                              void **mvQ, void *matA, void **mvP, int is_vec,
                              int *startQP, int *endQP, double *qAp, int ldQAP,
                              void **mv_ws, struct PASE_OPS_ *pase_ops);
void PASE_DefaultMultiVecOrth(void **x, int start_x, int *end_x,
                              void *B, struct PASE_OPS_ *pase_ops);

void PASE_BLASMultiVecAxpby(double alpha, void **x, double beta, void **y,
                            int *start, int *end, struct PASE_OPS_ *pase_ops);

#endif