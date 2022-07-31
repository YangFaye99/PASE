#include "pase_ops.h"
#include <stdbool.h>

#define AUX_TIME_INFO 0
int one = 1;

#if OPS_USE_MPI
int SIZE_B, SIZE_E, LDA;
void user_fn_submatrix_sum_pase_ops(double *in, double *inout, int *len, MPI_Datatype *dptr)
{
    int i, j;
    double *b, *a;
    double one = 1.0;
    int inc = 1;
    for (i = 0; i < *len; ++i)
    {
        for (j = 0; j < SIZE_B; ++j)
        {
            b = inout + j * LDA;
            a = in + j * LDA;
            daxpy(&SIZE_E, &one, a, &inc, b, &inc);
        }
    }
}
#endif

void PASE_OPS_Create(PASE_OPS **pase_ops, OPS *ops)
{
    (*pase_ops) = (PASE_OPS *)malloc(sizeof(PASE_OPS));
    (*pase_ops)->gcge_ops = ops;

    (*pase_ops)->VecCreateByMat = PASE_DefaultVecCreateByMat;
    (*pase_ops)->VecCreateByVec = PASE_DefaultVecCreateByVec;
    (*pase_ops)->VecDestroy = PASE_DefaultVecDestroy;
    (*pase_ops)->VecView = PASE_DefaultVecView;
    (*pase_ops)->VecAxpby = PASE_DefaultVecAxpby;
    (*pase_ops)->VecInnerProd = PASE_DefaultVecInnerProd;
    (*pase_ops)->VecLocalInnerProd = PASE_DefaultVecLocalInnerProd;
    (*pase_ops)->VecSetRandomValue = PASE_DefaultVecSetRandomValue;

    (*pase_ops)->MultiVecCreateByMat = PASE_DefaultMultiVecCreateByMat;
    (*pase_ops)->MultiVecCreateByVec = PASE_DefaultMultiVecCreateByVec;
    (*pase_ops)->MultiVecCreateByMultiVec = PASE_DefaultMultiVecCreateByMultiVec;
    (*pase_ops)->MultiVecDestroy = PASE_DefaultMultiVecDestroy;
    (*pase_ops)->GetVecFromMultiVec = PASE_DefaultGetVecFromMultiVec;
    (*pase_ops)->RestoreVecForMultiVec = PASE_DefaultRestoreVecForMultiVec;
    (*pase_ops)->MultiVecView = PASE_DefaultMultiVecView;
    (*pase_ops)->MultiVecSetRandomValue = PASE_DefaultMultiVecSetRandomValue;
    (*pase_ops)->MultiVecLocalInnerProd = PASE_DefaultMultiVecLocalInnerProd;
    (*pase_ops)->MultiVecInnerProd = PASE_DefaultMultiVecInnerProd;
    // (*pase_ops)->MultiVecAxpby = PASE_DefaultMultiVecAxpby;
    (*pase_ops)->MultiVecAxpby = PASE_BLASMultiVecAxpby;
    (*pase_ops)->MultiVecLinearComb = PASE_DefaultMultiVecLinearComb;

    (*pase_ops)->MatView = PASE_DefaultMatView;
    (*pase_ops)->MatAxpby = PASE_DefaultMatAxpby;
    (*pase_ops)->MatDotVec = PASE_DefaultMatDotVec;
    (*pase_ops)->MatTransDotVec = PASE_DefaultMatTransDotVec;
    (*pase_ops)->MatDotMultiVec = PASE_DefaultMatDotMultiVec;
    (*pase_ops)->MatTransDotMultiVec = PASE_DefaultMatTransDotMultiVec;

    (*pase_ops)->lapack_ops = ops->lapack_ops;
    (*pase_ops)->MultiVecQtAP = PASE_DefaultMultiVecQtAP;
    (*pase_ops)->MultiVecOrth = PASE_DefaultMultiVecOrth;
    return;
}

void PASE_OPS_Destroy(PASE_OPS **pase_ops)
{
    free(*pase_ops);
    *pase_ops = NULL;
    return;
}

void PASE_MatrixCreate(PASE_Matrix *pase_matrix,
                       int num_aux_vec, void *A_H, PASE_OPS *pase_ops)
{
    *pase_matrix = (PASE_Matrix)malloc(sizeof(pase_Matrix));
    (*pase_matrix)->num_aux_vec = num_aux_vec;
    (*pase_matrix)->A_H = A_H;
    (*pase_matrix)->is_diag = 0;
    (*pase_matrix)->if_sym = 1;
    pase_ops->gcge_ops->MultiVecCreateByMat(&((*pase_matrix)->aux_Hh), num_aux_vec, A_H, pase_ops->gcge_ops);
    pase_ops->gcge_ops->MultiVecCreateByMat(&((*pase_matrix)->AH_inv_aux_Hh), num_aux_vec, A_H, pase_ops->gcge_ops);
    (*pase_matrix)->aux_hh = (double *)calloc(num_aux_vec * num_aux_vec, sizeof(double));
    (*pase_matrix)->aux_hh_inv = NULL;

    int idx;
    for (idx = 0; idx < num_aux_vec; idx++)
    {
        (*pase_matrix)->aux_hh[num_aux_vec * idx + idx] = 1.0;
    }

    return;
}

void PASE_MatrixDestroy(PASE_Matrix *matrix, PASE_OPS *pase_ops)
{
    int num_aux_vec = (*matrix)->num_aux_vec;
    pase_ops->gcge_ops->MultiVecDestroy(&((*matrix)->aux_Hh), num_aux_vec, pase_ops->gcge_ops);
    free((*matrix)->aux_hh);
    (*matrix)->aux_hh = NULL;
    free(*matrix);
    *matrix = NULL;
    return;
}

void PASE_DefaultMatView(void *mat, struct PASE_OPS_ *pase_ops)
{
    PASE_Matrix pase_mat = (PASE_Matrix)mat;
    OPS *ops = pase_ops->gcge_ops;
    int c_idx, r_idx;

    ops->Printf("PASE_Matrix View:\n");
    ops->Printf("(1) A_H:\n");
    ops->MatView(pase_mat->A_H, ops);
    ops->Printf("\n(2) aux_Hh:(%d)\n", pase_mat->num_aux_vec);
    if (pase_mat->is_diag)
    {
        ops->Printf("NULL(diag)\n");
    }
    else
    {
        ops->MultiVecView(pase_mat->aux_Hh, 0, pase_mat->num_aux_vec - 1, ops);
    }
    ops->Printf("\n(3) aux_hh:\n");
    for (r_idx = 0; r_idx < pase_mat->num_aux_vec; r_idx++)
    {
        for (c_idx = 0; c_idx < pase_mat->num_aux_vec; c_idx++)
        {
            ops->Printf("%f\t", pase_mat->aux_hh[c_idx * pase_mat->num_aux_vec + r_idx]);
        }
        ops->Printf("\n");
    }

    return;
}

void PASE_DefaultMatAxpby(double alpha, void *matX, double beta, void *matY, struct PASE_OPS_ *pase_ops)
{
    // pase_ops->gcge_ops->Printf("PASE_DefaultMatAxpby\n");
    void *Y_A_H = ((PASE_Matrix)matY)->A_H;
    void **Y_aux_Hh = ((PASE_Matrix)matY)->aux_Hh;
    double *Y_aux_hh = ((PASE_Matrix)matY)->aux_hh;
    int num_aux_vec = ((PASE_Matrix)matY)->num_aux_vec, inc = 1;
    int start[2] = {0, 0};
    int end[2] = {num_aux_vec, num_aux_vec};

    int length = num_aux_vec * num_aux_vec;
    if (beta == 0.0)
    {
        memset(Y_aux_hh, 0, length * sizeof(double));
    }
    else if (beta != 1.0)
    {
        dscal(&length, &beta, Y_aux_hh, &inc);
    }

    if (matX != NULL)
    {
        void *X_A_H = ((PASE_Matrix)matX)->A_H;
        void **X_aux_Hh = ((PASE_Matrix)matX)->aux_Hh;
        double *X_aux_hh = ((PASE_Matrix)matX)->aux_hh;
        assert(num_aux_vec == ((PASE_Matrix)matX)->num_aux_vec);
        // A_H
        pase_ops->gcge_ops->MatAxpby(alpha, X_A_H, beta, Y_A_H, pase_ops->gcge_ops);
        // aux_Hh
        pase_ops->gcge_ops->MultiVecAxpby(alpha, X_aux_Hh, beta, Y_aux_Hh, start, end, pase_ops->gcge_ops);
        // aux_hh
        daxpy(&length, &alpha, X_aux_hh, &inc, Y_aux_hh, &inc);
    }
    else
    {
        pase_ops->gcge_ops->MatAxpby(alpha, NULL, beta, Y_A_H, pase_ops->gcge_ops);
        pase_ops->gcge_ops->MultiVecAxpby(alpha, NULL, beta, Y_aux_Hh, start, end, pase_ops->gcge_ops);
    }
    return;
}

void PASE_DefaultMatDotVec(void *mat, void *x, void *y, struct PASE_OPS_ *pase_ops)
{
    pase_ops->gcge_ops->Printf(" TO DO PASE_DefaultMatDotVec ");
    return;
}

void PASE_DefaultMatTransDotVec(void *mat, void *x, void *y, struct PASE_OPS_ *pase_ops)
{
    pase_ops->gcge_ops->Printf(" TO DO PASE_DefaultMatTransDotVec ");
    return;
}

static void MatDotMultiVecDiag(void *mat, void **x, void **y,
                               int *start, int *end, struct PASE_OPS_ *pase_ops)
{
    void *A_H = ((PASE_Matrix)mat)->A_H;
    double *aux_hh = ((PASE_Matrix)mat)->aux_hh;
    int A_num_aux_vec = ((PASE_Matrix)mat)->num_aux_vec;

    void **x_b_H = ((PASE_MultiVector)x)->b_H;
    double *x_aux_h = ((PASE_MultiVector)x)->aux_h;
    int x_num_aux_vec = ((PASE_MultiVector)x)->num_aux_vec;

    void **y_b_H = ((PASE_MultiVector)y)->b_H;
    double *y_aux_h = ((PASE_MultiVector)y)->aux_h;
    int y_num_aux_vec = ((PASE_MultiVector)y)->num_aux_vec;

    assert(A_num_aux_vec == x_num_aux_vec);
    assert(end[0] - start[0] == end[1] - start[1]);

    int ncols = end[1] - start[1];
    double one = 1.0, zero = 0.0;
    char charN = 'N', charL = 'L', charT = 'T';
#if AUX_TIME_INFO
    extern double time_pase_matmv_x_GCGE_diag;
    time_pase_matmv_x_GCGE_diag -= DefaultGetWtime();
#endif
    pase_ops->gcge_ops->MatDotMultiVec(A_H, x_b_H, y_b_H, start, end, pase_ops->gcge_ops);
#if AUX_TIME_INFO
    time_pase_matmv_x_GCGE_diag += DefaultGetWtime();
    extern double time_pase_matmv_mu_blas_diag;
    time_pase_matmv_mu_blas_diag -= DefaultGetWtime();
#endif
    if (((PASE_Matrix)mat)->if_sym == 0)
    {

        dgemm(&charN, &charN, &A_num_aux_vec, &ncols, &A_num_aux_vec,
              &one, aux_hh, &A_num_aux_vec, x_aux_h + start[0] * x_num_aux_vec, &x_num_aux_vec,
              &zero, y_aux_h + y_num_aux_vec * start[1], &y_num_aux_vec);
    }
    else
    {
        dsymm(&charL, &charL, &A_num_aux_vec, &ncols,
              &one, aux_hh, &A_num_aux_vec, x_aux_h + start[0] * x_num_aux_vec, &x_num_aux_vec,
              &zero, y_aux_h + y_num_aux_vec * start[1], &y_num_aux_vec);
    }
#if AUX_TIME_INFO
    time_pase_matmv_mu_blas_diag += DefaultGetWtime();
#endif

    return;
}

static void MatDotMultiVecAll(void *mat, void **x, void **y,
                              int *start, int *end, struct PASE_OPS_ *pase_ops)
{
    OPS *gcge_ops = pase_ops->gcge_ops;
    void *A_H = ((PASE_Matrix)mat)->A_H;
    void **aux_Hh = ((PASE_Matrix)mat)->aux_Hh;
    double *aux_hh = ((PASE_Matrix)mat)->aux_hh;
    int A_num_aux_vec = ((PASE_Matrix)mat)->num_aux_vec;

    void **x_b_H = ((PASE_MultiVector)x)->b_H;
    double *x_aux_h = ((PASE_MultiVector)x)->aux_h;
    int x_num_aux_vec = ((PASE_MultiVector)x)->num_aux_vec;
    void **y_b_H = ((PASE_MultiVector)y)->b_H;
    double *y_aux_h = ((PASE_MultiVector)y)->aux_h;
    int y_num_aux_vec = ((PASE_MultiVector)y)->num_aux_vec;
    double *y_aux_h_tmp = ((PASE_MultiVector)y)->aux_h_tmp;
    if (y_aux_h_tmp == NULL)
    {
        y_aux_h_tmp = (double *)calloc(y_num_aux_vec * ((PASE_MultiVector)y)->num_vec, sizeof(double));
    }

    int mv_s[2];
    int mv_e[2];
    double one = 1.0;
    double zero = 0.0;
    int i = 0;

#if AUX_TIME_INFO
    extern double time_pase_matmv_x_GCGE;
    time_pase_matmv_x_GCGE -= DefaultGetWtime();
#endif
    // compute r_b_H = A_H * x_b_H
    gcge_ops->MatDotMultiVec(A_H, x_b_H, y_b_H, start, end, gcge_ops);
#if AUX_TIME_INFO
    time_pase_matmv_x_GCGE += DefaultGetWtime();
#endif

#if AUX_TIME_INFO
    extern double time_pase_matmv_mu_inner_l;
    time_pase_matmv_mu_inner_l -= DefaultGetWtime();
#endif
    // begin compute r->aux_h = aux_Hh^T * x->b_H
#if OPS_USE_MPI
    MPI_Request request;
    MPI_Status status;
    MPI_Datatype SUBMATRIX;
    MPI_Op SUBMATRIX_SUM;
    SIZE_B = end[1] - start[1];
    SIZE_E = A_num_aux_vec;
    LDA = y_num_aux_vec;
#endif
    mv_s[0] = 0;
    mv_e[0] = A_num_aux_vec;
    mv_s[1] = start[0];
    mv_e[1] = end[0];
    gcge_ops->MultiVecLocalInnerProd('N', aux_Hh, x_b_H, 0, mv_s, mv_e,
                                     y_aux_h + start[1] * y_num_aux_vec, y_num_aux_vec, gcge_ops);
#if OPS_USE_MPI
    MPI_Type_vector(SIZE_B, SIZE_E, LDA, MPI_DOUBLE, &SUBMATRIX);
    MPI_Type_commit(&SUBMATRIX);
    MPI_Op_create((MPI_User_function *)user_fn_submatrix_sum_pase_ops, 1, &SUBMATRIX_SUM);
    MPI_Iallreduce(MPI_IN_PLACE, y_aux_h + start[1] * y_num_aux_vec, 1, SUBMATRIX, SUBMATRIX_SUM, MPI_COMM_WORLD, &request);
#endif
#if AUX_TIME_INFO
    time_pase_matmv_mu_inner_l += DefaultGetWtime();
#endif

#if AUX_TIME_INFO
    extern double time_pase_matmv_x_blas;
    time_pase_matmv_x_blas -= DefaultGetWtime();
#endif
    // compute r_b_H += aux_Hh * aux_h
    mv_s[0] = 0;
    mv_e[0] = A_num_aux_vec;
    mv_s[1] = start[1];
    mv_e[1] = end[1];
    int ncols = end[0] - start[0];
#if OPS_USE_SLEPC
    const PetscScalar *A_mv_tmp, *y_b_H_data;
    int tmp_nrows, tmp_ncols, y_nrows, y_ncols;
    char charN = 'N', charL = 'L', charT = 'T';
    BVGetArrayRead((BV)aux_Hh, &A_mv_tmp);
    BVGetArrayRead((BV)y_b_H, &y_b_H_data);
    BVGetSizes((BV)aux_Hh, &tmp_nrows, NULL, &tmp_ncols);
    BVGetSizes((BV)y_b_H, &y_nrows, NULL, &y_ncols);
    if (tmp_nrows != 0)
    {
        dgemm(&charN, &charN, &tmp_nrows, &ncols, &x_num_aux_vec, &one,
              (double *)A_mv_tmp, &tmp_nrows, x_aux_h + start[0] * x_num_aux_vec, &x_num_aux_vec,
              &one, (double *)y_b_H_data + start[1] * y_nrows, &y_nrows);
    }
    BVRestoreArrayRead((BV)aux_Hh, &A_mv_tmp);
    BVRestoreArrayRead((BV)y_b_H, &y_b_H_data);
#elif
    gcge_ops->Printf("you can only use slepc for pase now.\n");
    exit();
#endif
#if AUX_TIME_INFO
    time_pase_matmv_x_blas += DefaultGetWtime();
#endif

#if AUX_TIME_INFO
    extern double time_pase_matmv_mu_blas;
    time_pase_matmv_mu_blas -= DefaultGetWtime();
#endif
    // compute y_aux_h_tmp = aux_hh^T * x->aux_h
    if (((PASE_Matrix)mat)->if_sym == 0)
    {
        dgemm(&charN, &charN, &A_num_aux_vec, &ncols, &A_num_aux_vec,
              &one, aux_hh, &A_num_aux_vec, x_aux_h + start[0] * x_num_aux_vec, &x_num_aux_vec,
              &zero, y_aux_h_tmp, &y_num_aux_vec);
    }
    else
    {
        dsymm(&charL, &charL, &A_num_aux_vec, &ncols,
              &one, aux_hh, &A_num_aux_vec, x_aux_h + start[0] * x_num_aux_vec, &x_num_aux_vec,
              &zero, y_aux_h_tmp, &y_num_aux_vec);
    }
#if AUX_TIME_INFO
    time_pase_matmv_mu_blas += DefaultGetWtime();
#endif

#if AUX_TIME_INFO
    extern double time_pase_matmv_mu_inner_c;
    time_pase_matmv_mu_inner_c -= DefaultGetWtime();
#endif
#if OPS_USE_MPI
    MPI_Wait(&request, &status);
    MPI_Op_free(&SUBMATRIX_SUM);
    MPI_Type_free(&SUBMATRIX);
#endif
#if AUX_TIME_INFO
    time_pase_matmv_mu_inner_c += DefaultGetWtime();
#endif

    // compute r->aux_h += y_aux_h_tmp
    int num = y_num_aux_vec * ncols;
    int inc = 1;
    daxpy(&num, &one, y_aux_h_tmp, &inc, y_aux_h + start[1] * y_num_aux_vec, &inc);
    return;
}

void PASE_DefaultMatDotMultiVec(void *mat, void **x, void **y,
                                int *start, int *end, struct PASE_OPS_ *pase_ops)
{
#if AUX_TIME_INFO
    extern double time_pase_matmv;
    extern int number_pase_matmv;
    extern int number_pase_matmv_computeW;
    extern int colnum_pase_matmv;
    extern int colnum_pase_matmv_computeW;
    extern bool if_compute_w;
    if (if_compute_w)
    {
        number_pase_matmv_computeW++;
        colnum_pase_matmv_computeW += end[0] - start[0];
    }
    number_pase_matmv++;
    colnum_pase_matmv += end[0] - start[0];
    time_pase_matmv -= DefaultGetWtime();
#endif
    if (((PASE_Matrix)mat)->is_diag == 0)
    {
#if AUX_TIME_INFO
        extern double time_pase_matmv_all;
        extern int number_pase_matmv_all;
        extern int colnum_pase_matmv_all;
        number_pase_matmv_all++;
        colnum_pase_matmv_all += end[0] - start[0];
        time_pase_matmv_all -= DefaultGetWtime();
#endif
        MatDotMultiVecAll(mat, x, y, start, end, pase_ops);
#if AUX_TIME_INFO
        time_pase_matmv_all += DefaultGetWtime();
#endif
    }
    else
    {
#if AUX_TIME_INFO
        extern double time_pase_matmv_diag;
        extern int number_pase_matmv_diag;
        extern int colnum_pase_matmv_diag;
        number_pase_matmv_diag++;
        colnum_pase_matmv_diag += end[0] - start[0];
        time_pase_matmv_diag -= DefaultGetWtime();
#endif
        MatDotMultiVecDiag(mat, x, y, start, end, pase_ops);
#if AUX_TIME_INFO
        time_pase_matmv_diag += DefaultGetWtime();
#endif
    }
#if AUX_TIME_INFO
    time_pase_matmv += DefaultGetWtime();
#endif
    return;
}

static void MatTransDotMultiVecDiag(void *mat, void **x, void **y,
                                    int *start, int *end, struct PASE_OPS_ *pase_ops)
{
    void *A_H = ((PASE_Matrix)mat)->A_H;
    double *aux_hh = ((PASE_Matrix)mat)->aux_hh;
    int A_num_aux_vec = ((PASE_Matrix)mat)->num_aux_vec;

    void **x_b_H = ((PASE_MultiVector)x)->b_H;
    double *x_aux_h = ((PASE_MultiVector)x)->aux_h;
    int x_num_aux_vec = ((PASE_MultiVector)x)->num_aux_vec;

    void **y_b_H = ((PASE_MultiVector)y)->b_H;
    double *y_aux_h = ((PASE_MultiVector)y)->aux_h;
    int y_num_aux_vec = ((PASE_MultiVector)y)->num_aux_vec;

    assert(A_num_aux_vec == x_num_aux_vec);
    assert(end[0] - start[0] == end[1] - start[1]);

    int ncols = end[1] - start[1];
    double one = 1.0, zero = 0.0;
    char charT = 'T', charL = 'L', charN = 'N';
    pase_ops->gcge_ops->MatTransDotMultiVec(A_H, x_b_H, y_b_H, start, end, pase_ops->gcge_ops);
    if (((PASE_Matrix)mat)->if_sym == 0)
    {
        dgemm(&charT, &charN, &A_num_aux_vec, &ncols, &A_num_aux_vec,
              &one, aux_hh, &A_num_aux_vec, x_aux_h + start[0] * x_num_aux_vec, &x_num_aux_vec,
              &zero, y_aux_h + y_num_aux_vec * start[1], &y_num_aux_vec);
    }
    else
    {
        dsymm(&charL, &charL, &A_num_aux_vec, &ncols,
              &one, aux_hh, &A_num_aux_vec, x_aux_h + start[0] * x_num_aux_vec, &x_num_aux_vec,
              &zero, y_aux_h + y_num_aux_vec * start[1], &y_num_aux_vec);
    }
    return;
}

static void MatTransDotMultiVecAll(void *mat, void **x, void **y,
                                   int *start, int *end, struct PASE_OPS_ *pase_ops)
{
    OPS *gcge_ops = pase_ops->gcge_ops;
    void *A_H = ((PASE_Matrix)mat)->A_H;
    void **aux_Hh = ((PASE_Matrix)mat)->aux_Hh;
    double *aux_hh = ((PASE_Matrix)mat)->aux_hh;
    int A_num_aux_vec = ((PASE_Matrix)mat)->num_aux_vec;

    void **x_b_H = ((PASE_MultiVector)x)->b_H;
    double *x_aux_h = ((PASE_MultiVector)x)->aux_h;
    int x_num_aux_vec = ((PASE_MultiVector)x)->num_aux_vec;
    void **y_b_H = ((PASE_MultiVector)y)->b_H;
    double *y_aux_h = ((PASE_MultiVector)y)->aux_h;
    int y_num_aux_vec = ((PASE_MultiVector)y)->num_aux_vec;
    double *y_aux_h_tmp = ((PASE_MultiVector)y)->aux_h_tmp;

    int mv_s[2];
    int mv_e[2];
    double one = 1.0;
    double zero = 0.0;
    int i = 0;

    // compute r_b_H = A_H * x_b_H
    gcge_ops->MatTransDotMultiVec(A_H, x_b_H, y_b_H, start, end, gcge_ops);

    // begin compute r->aux_h = aux_Hh^T * x->b_H
#if OPS_USE_MPI
    MPI_Request request;
    MPI_Status status;
    MPI_Datatype SUBMATRIX;
    MPI_Op SUBMATRIX_SUM;
    SIZE_B = end[1] - start[1];
    SIZE_E = A_num_aux_vec;
    LDA = y_num_aux_vec;
#endif
    mv_s[0] = 0;
    mv_e[0] = A_num_aux_vec;
    mv_s[1] = start[0];
    mv_e[1] = end[0];
    gcge_ops->MultiVecLocalInnerProd('N', aux_Hh, x_b_H, 0, mv_s, mv_e,
                                     y_aux_h + start[1] * y_num_aux_vec, y_num_aux_vec, gcge_ops);

#if OPS_USE_MPI
    MPI_Type_vector(SIZE_B, SIZE_E, LDA, MPI_DOUBLE, &SUBMATRIX);
    MPI_Type_commit(&SUBMATRIX);
    MPI_Op_create((MPI_User_function *)user_fn_submatrix_sum_pase_ops, 1, &SUBMATRIX_SUM);
    MPI_Iallreduce(MPI_IN_PLACE, y_aux_h + start[1] * y_num_aux_vec, 1, SUBMATRIX, SUBMATRIX_SUM, MPI_COMM_WORLD, &request);
#endif

    // compute r_b_H += aux_Hh * aux_h
    mv_s[0] = 0;
    mv_e[0] = A_num_aux_vec;
    mv_s[1] = start[1];
    mv_e[1] = end[1];
    int ncols = end[0] - start[0];
#if OPS_USE_SLEPC
    const PetscScalar *A_mv_tmp, *y_b_H_data;
    int tmp_nrows, tmp_ncols;
    char charN = 'N', charL = 'L', charT = 'T';
    BVGetArrayRead((BV)aux_Hh, &A_mv_tmp);
    BVGetArrayRead((BV)y_b_H, &y_b_H_data);
    BVGetSizes((BV)aux_Hh, &tmp_nrows, NULL, &tmp_ncols);
    dgemm(&charN, &charN, &tmp_nrows, &ncols, &x_num_aux_vec, &one,
          (double *)A_mv_tmp, &tmp_nrows, x_aux_h + start[0] * x_num_aux_vec, &x_num_aux_vec,
          &one, (double *)y_b_H_data + start[1] * y_num_aux_vec, &y_num_aux_vec);
    BVRestoreArrayRead((BV)aux_Hh, &A_mv_tmp);
    BVRestoreArrayRead((BV)y_b_H, &y_b_H_data);
#elif
    gcge_ops->Printf("you can only use slepc for pase now.\n");
    assert(0 == 1);
#endif

    // compute y_aux_h_tmp = aux_hh^T * x->aux_h
    if (((PASE_Matrix)mat)->if_sym == 0)
    {
        dgemm(&charT, &charN, &A_num_aux_vec, &ncols, &A_num_aux_vec,
              &one, aux_hh, &A_num_aux_vec, x_aux_h + start[0] * x_num_aux_vec, &x_num_aux_vec,
              &zero, y_aux_h_tmp, &y_num_aux_vec);
    }
    else
    {
        dsymm(&charL, &charL, &A_num_aux_vec, &ncols,
              &one, aux_hh, &A_num_aux_vec, x_aux_h + start[0] * x_num_aux_vec, &x_num_aux_vec,
              &zero, y_aux_h_tmp, &y_num_aux_vec);
    }
#if OPS_USE_MPI
    MPI_Wait(&request, &status);
    MPI_Op_free(&SUBMATRIX_SUM);
    MPI_Type_free(&SUBMATRIX);
#endif
    // compute r->aux_h += y_aux_h_tmp
    int num = y_num_aux_vec * ncols;
    int inc = 1;
    daxpy(&num, &one, y_aux_h_tmp, &inc, y_aux_h + start[1] * y_num_aux_vec, &inc);
    return;
}

void PASE_DefaultMatTransDotMultiVec(void *mat, void **x, void **y,
                                     int *start, int *end, struct PASE_OPS_ *pase_ops)
{
    if (((PASE_Matrix)mat)->is_diag == 0)
    {
        MatTransDotMultiVecAll(mat, x, y, start, end, pase_ops);
    }
    else
    {
        MatTransDotMultiVecDiag(mat, x, y, start, end, pase_ops);
    }
    return;
}

void PASE_DefaultVecAxpby(double alpha, void *x, double beta, void *y, struct PASE_OPS_ *pase_ops)
{
    PASE_Vector pase_y = (PASE_Vector)y;
    OPS *ops = pase_ops->gcge_ops;
    int inc = 1;
    int num_aux = pase_y->num_aux_vec;

    double *destin = pase_y->aux_h;
    if (beta == 0.0)
    {
        memset(destin, 0, num_aux * sizeof(double));
    }
    else if (beta != 1.0)
    {
        dscal(&num_aux, &beta, destin, &inc);
    }

    if (x != NULL)
    {
        PASE_Vector pase_x = (PASE_Vector)x;
        double *source = pase_x->aux_h;
        daxpy(&num_aux, &alpha, source, &inc, destin, &inc);
        pase_ops->gcge_ops->VecAxpby(alpha, pase_x->b_H, beta, pase_y->b_H, pase_ops->gcge_ops);
    }
    else
    {
        pase_ops->gcge_ops->VecAxpby(alpha, NULL, beta, pase_y->b_H, pase_ops->gcge_ops);
    }

    return;
}

// optimizable
void PASE_DefaultVecInnerProd(void *x, void *y, double *inner_prod, struct PASE_OPS_ *pase_ops)
{
    void *x_b_H = ((PASE_Vector)x)->b_H;
    double *x_aux_h = ((PASE_Vector)x)->aux_h;
    void *y_b_H = ((PASE_Vector)y)->b_H;
    double *y_aux_h = ((PASE_Vector)y)->aux_h;
    int num_aux_vec = ((PASE_Vector)y)->num_aux_vec;
    int inc = 1;
    pase_ops->gcge_ops->VecInnerProd(x_b_H, y_b_H, inner_prod, pase_ops->gcge_ops);
    *inner_prod += ddot(&num_aux_vec, x_aux_h, &inc, y_aux_h, &inc);
    return;
}

void PASE_DefaultVecLocalInnerProd(void *x, void *y, double *inner_prod, struct PASE_OPS_ *pase_ops)
{
#if AUX_TIME_INFO
    extern double time_pase_innerprod_local;
    time_pase_innerprod_local -= DefaultGetWtime();
#endif
    void *x_b_H = ((PASE_Vector)x)->b_H;
    void *y_b_H = ((PASE_Vector)y)->b_H;
    pase_ops->gcge_ops->VecLocalInnerProd(x_b_H, y_b_H, inner_prod, pase_ops->gcge_ops);
    return;
#if AUX_TIME_INFO
    extern double time_pase_innerprod_local;
    time_pase_innerprod_local += DefaultGetWtime();
#endif
}

void PASE_DefaultVecSetRandomValue(void *x, struct PASE_OPS_ *pase_ops)
{
    void *b_H = ((PASE_Vector)x)->b_H;
    double *aux_h = ((PASE_Vector)x)->aux_h;
    int num_aux_vec = ((PASE_Vector)x)->num_aux_vec;
    pase_ops->gcge_ops->VecSetRandomValue(b_H, pase_ops->gcge_ops);
    int i = 0;
    srand(time(NULL));
    for (i = 0; i < num_aux_vec; i++)
    {
        aux_h[i] = ((double)rand()) / ((double)RAND_MAX + 1);
    }
    return;
}

void PASE_DefaultVecDestroy(void **des_vec, struct PASE_OPS_ *pase_ops)
{
    pase_ops->gcge_ops->VecDestroy(&(((PASE_Vector)(*des_vec))->b_H), pase_ops->gcge_ops);
    free(((PASE_Vector)(*des_vec))->aux_h);
    free(((PASE_Vector)(*des_vec))->aux_h_tmp);
    ((PASE_Vector)(*des_vec))->aux_h = NULL;
    ((PASE_Vector)(*des_vec))->aux_h_tmp = NULL;
    free((PASE_Vector)(*des_vec));
    *des_vec = NULL;
    return;
}

void PASE_DefaultVecView(void *x, struct PASE_OPS_ *pase_ops)
{
    PASE_Vector pase_x = (PASE_Vector)x;
    OPS *ops = pase_ops->gcge_ops;
    int idx;

    ops->Printf("PASE_Vector View:\n");
    ops->Printf("(1) b_H:\n");
    ops->VecView(pase_x->b_H, ops);
    ops->Printf("\n(2) aux_h:(%d)\n", pase_x->num_aux_vec);
    for (idx = 0; idx < pase_x->num_aux_vec; idx++)
    {
        ops->Printf("%f\t", pase_x->aux_h[idx]);
        ops->Printf("\n");
    }
    return;
}

void PASE_DefaultVecCreateByMat(void **des_vec, void *src_mat, struct PASE_OPS_ *pase_ops)
{
    PASE_Vector v = (PASE_Vector)malloc(sizeof(pase_Vector));
    void *A_H = ((PASE_Matrix)src_mat)->A_H;
    int num_aux_vec = ((PASE_Matrix)src_mat)->num_aux_vec;

    v->num_aux_vec = num_aux_vec;
    pase_ops->gcge_ops->VecCreateByMat(&(v->b_H), A_H, pase_ops->gcge_ops);
    v->aux_h = (double *)calloc(num_aux_vec, sizeof(double));
    v->aux_h_tmp = (double *)calloc(num_aux_vec, sizeof(double));
    *des_vec = (void *)v;
    return;
}

void PASE_DefaultVecCreateByVec(void **des_vec, void *src_vec, struct PASE_OPS_ *pase_ops)
{
    PASE_Vector vec = (PASE_Vector)malloc(sizeof(pase_Vector));
    void *b_H = ((PASE_Vector)src_vec)->b_H;
    double *aux_h = ((PASE_Vector)src_vec)->aux_h;
    int num_aux_vec = ((PASE_Vector)src_vec)->num_aux_vec;

    vec->num_aux_vec = num_aux_vec;
    pase_ops->gcge_ops->VecCreateByVec(&(vec->b_H), b_H, pase_ops->gcge_ops);
    vec->aux_h = (double *)calloc(num_aux_vec, sizeof(double));
    vec->aux_h_tmp = (double *)calloc(num_aux_vec, sizeof(double));
    *des_vec = (void *)vec;
    return;
}

void PASE_DefaultMultiVecCreateByMat(void ***multi_vec, int num_vec, void *src_mat, struct PASE_OPS_ *pase_ops)
{
    void *A_H = ((PASE_Matrix)src_mat)->A_H;
    int num_aux_vec = ((PASE_Matrix)src_mat)->num_aux_vec;

    PASE_MultiVector vecs = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    vecs->num_vec = num_vec;
    vecs->num_aux_vec = num_aux_vec;
    pase_ops->gcge_ops->MultiVecCreateByMat(&(vecs->b_H), num_vec, A_H, pase_ops->gcge_ops);
    vecs->aux_h_tmp = (double *)calloc(num_vec * num_aux_vec, sizeof(double));
    vecs->aux_h = (double *)calloc(num_vec * num_aux_vec, sizeof(double));
    *multi_vec = (void **)vecs;
    return;
}

void PASE_DefaultMultiVecCreateByVec(void ***multi_vec, int num_vec, void *src_vec, struct PASE_OPS_ *pase_ops)
{
    void *b_H = ((PASE_Vector)src_vec)->b_H;
    int num_aux_vec = ((PASE_Vector)src_vec)->num_aux_vec;
    PASE_MultiVector vecs = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    vecs->num_vec = num_vec;
    vecs->num_aux_vec = num_aux_vec;
    pase_ops->gcge_ops->MultiVecCreateByVec(&(vecs->b_H), num_vec, b_H, pase_ops->gcge_ops);
    vecs->aux_h = (double *)calloc(num_vec * num_aux_vec, sizeof(double));
    vecs->aux_h_tmp = (double *)calloc(num_vec * num_aux_vec, sizeof(double));
    *multi_vec = (void **)vecs;
    return;
}

void PASE_DefaultMultiVecCreateByMultiVec(void ***multi_vec, int num_vec, void **src_mv, struct PASE_OPS_ *pase_ops)
{
    void **b_H = ((PASE_MultiVector)src_mv)->b_H;
    int num_aux_vec = ((PASE_MultiVector)src_mv)->num_aux_vec;
    PASE_MultiVector vecs = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    vecs->num_vec = num_vec;
    vecs->num_aux_vec = num_aux_vec;
    pase_ops->gcge_ops->MultiVecCreateByMultiVec(&(vecs->b_H), num_vec, b_H, pase_ops->gcge_ops);
    vecs->aux_h = (double *)calloc(num_vec * num_aux_vec, sizeof(double));
    vecs->aux_h_tmp = (double *)calloc(num_vec * num_aux_vec, sizeof(double));
    *multi_vec = (void **)vecs;
    return;
}

void PASE_DefaultMultiVecDestroy(void ***multi_vec, int num_vec, struct PASE_OPS_ *pase_ops)
{
    pase_ops->gcge_ops->MultiVecDestroy(&(((PASE_MultiVector)(*multi_vec))->b_H),
                                        ((PASE_MultiVector)(*multi_vec))->num_vec, pase_ops->gcge_ops);
    free(((PASE_MultiVector)(*multi_vec))->aux_h);
    free(((PASE_MultiVector)(*multi_vec))->aux_h_tmp);
    ((PASE_MultiVector)(*multi_vec))->aux_h = NULL;
    ((PASE_MultiVector)(*multi_vec))->aux_h_tmp = NULL;
    free((PASE_MultiVector)(*multi_vec));
    *multi_vec = NULL;
    return;
}

void PASE_DefaultGetVecFromMultiVec(void **multi_vec, int col, void **vec, struct PASE_OPS_ *pase_ops)
{
    void **mv_b_H = ((PASE_MultiVector)multi_vec)->b_H;
    double *mv_aux_h = ((PASE_MultiVector)multi_vec)->aux_h;
    double *mv_aux_h_tmp = ((PASE_MultiVector)multi_vec)->aux_h_tmp;
    int num_aux_vec = ((PASE_MultiVector)multi_vec)->num_aux_vec;
    PASE_Vector pase_x = (PASE_Vector)malloc(sizeof(pase_Vector));
    pase_x->num_aux_vec = num_aux_vec;
    pase_ops->gcge_ops->GetVecFromMultiVec(mv_b_H, col, &(pase_x->b_H), pase_ops->gcge_ops);
    pase_x->aux_h = mv_aux_h + col * num_aux_vec;
    pase_x->aux_h_tmp = mv_aux_h_tmp + col * num_aux_vec;
    *vec = pase_x;
    return;
}

void PASE_DefaultRestoreVecForMultiVec(void **multi_vec, int col, void **vec, struct PASE_OPS_ *pase_ops)
{
    void **mv_b_H = ((PASE_MultiVector)multi_vec)->b_H;
    double *mv_aux_h = ((PASE_MultiVector)multi_vec)->aux_h;
    int num_aux_vec = ((PASE_MultiVector)multi_vec)->num_aux_vec;
    pase_ops->gcge_ops->RestoreVecForMultiVec(mv_b_H, col, &(((PASE_Vector)(*vec))->b_H), pase_ops->gcge_ops);
    free((PASE_Vector)(*vec));
    *vec = NULL;
    return;
}

void PASE_DefaultMultiVecView(void **x, int start, int end, struct PASE_OPS_ *pase_ops)
{
    PASE_MultiVector pase_mv = (PASE_MultiVector)x;
    OPS *ops = pase_ops->gcge_ops;
    int idx, idx_vec;

    ops->Printf("PASE_MultiVector View ( %d - %d ):\n", start, end);
    // ops->Printf("(1) b_H:\n");
    // ops->MultiVecView(pase_mv->b_H, start, end, ops);
    ops->Printf("(2) aux_h:\n");
    for (idx_vec = start; idx_vec < end; idx_vec++)
    {
        ops->Printf("No.%d: ", idx_vec);
        for (idx = 0; idx < pase_mv->num_aux_vec; idx++)
        {
            ops->Printf(" %f", pase_mv->aux_h[idx_vec * pase_mv->num_aux_vec + idx]);
        }
        ops->Printf("\n");
    }
    // int myrank, nprocs, i;
    // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    // MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    // for (i = 0; i < nprocs; i++)
    // {
    //     if (myrank == i)
    //     {
    //         printf("here is rank %d\n", myrank);
    //         for (idx_vec = start; idx_vec < end; idx_vec++)
    //         {
    //             printf("No.%d: ", idx_vec);
    //             for (idx = 0; idx < pase_mv->num_aux_vec; idx++)
    //             {
    //                 printf(" %f", pase_mv->aux_h[idx_vec * pase_mv->num_aux_vec + idx]);
    //             }
    //             printf("\n");
    //         }
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    return;
}

void PASE_DefaultMultiVecLocalInnerProd(char nsdIP, void **x, void **y, int is_vec, int *start, int *end,
                                        double *inner_prod, int ldIP, struct PASE_OPS_ *pase_ops)
{
    void **x_b_H = ((PASE_MultiVector)x)->b_H;
    void **y_b_H = ((PASE_MultiVector)y)->b_H;
    pase_ops->gcge_ops->MultiVecLocalInnerProd(nsdIP, x_b_H, y_b_H, is_vec, start, end,
                                               inner_prod, ldIP, pase_ops->gcge_ops);

    int num_aux_vec = ((PASE_MultiVector)y)->num_aux_vec;
    int nrows = end[0] - start[0];
    int ncols = end[1] - start[1];
    double alpha = 1.0;
    double beta = 1.0;
    double *x_aux_h = ((PASE_MultiVector)x)->aux_h;
    double *y_aux_h = ((PASE_MultiVector)y)->aux_h;
    double *a_tmp = (double *)calloc(ldIP * (end[1] - start[1]), sizeof(double));
    LDA = ldIP;
#if OPS_USE_MPI
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0)
#endif
    {
        pase_ops->lapack_ops->DenseMatQtAP('L', nsdIP, num_aux_vec, num_aux_vec, nrows, ncols, alpha,
                                           x_aux_h + start[0] * num_aux_vec, num_aux_vec,
                                           NULL, num_aux_vec,
                                           y_aux_h + start[1] * num_aux_vec, num_aux_vec,
                                           beta, a_tmp, LDA, NULL);
    }
    int inc = 1;
    int mid = ncols * ldIP;
    daxpy(&mid, &alpha, a_tmp, &inc, inner_prod, &inc);
    free(a_tmp);
    a_tmp = NULL;
    return;
}

void PASE_DefaultMultiVecInnerProd(char nsdIP, void **x, void **y, int is_vec, int *start, int *end,
                                   double *inner_prod, int ldIP, struct PASE_OPS_ *pase_ops)
{
#if AUX_TIME_INFO
    extern double time_pase_innerprod_full;
    time_pase_innerprod_full -= DefaultGetWtime();
#endif
    void **x_b_H = ((PASE_MultiVector)x)->b_H;
    double *x_aux_h = ((PASE_MultiVector)x)->aux_h;
    double *a_tmp = (double *)calloc(ldIP * (end[1] - start[1]), sizeof(double));
    void **y_b_H;
    double *y_aux_h;
    int num_aux_vec;
    if (is_vec == 0)
    {
        y_b_H = ((PASE_MultiVector)y)->b_H;
        y_aux_h = ((PASE_MultiVector)y)->aux_h;
        num_aux_vec = ((PASE_MultiVector)y)->num_aux_vec;
    }
    else
    {
        y_b_H = &(((PASE_Vector)(y[0]))->b_H);
        y_aux_h = ((PASE_Vector)(y[0]))->aux_h;
        num_aux_vec = ((PASE_Vector)(y[0]))->num_aux_vec;
    }

#if OPS_USE_MPI
    MPI_Request request;
    MPI_Status status;
#endif
    pase_ops->gcge_ops->MultiVecLocalInnerProd(nsdIP, x_b_H, y_b_H, is_vec, start, end,
                                               inner_prod, ldIP, pase_ops->gcge_ops);
#if OPS_USE_MPI

    SIZE_B = end[1] - start[1];
    if (nsdIP == 'D')
    {
        SIZE_E = ldIP;
    }
    else
    {
        SIZE_E = end[0] - start[0];
    }
    LDA = ldIP;

    MPI_Datatype SUBMATRIX;
    MPI_Type_vector(SIZE_B, SIZE_E, LDA, MPI_DOUBLE, &SUBMATRIX);
    MPI_Type_commit(&SUBMATRIX);

    MPI_Op SUBMATRIX_SUM;
    MPI_Op_create((MPI_User_function *)user_fn_submatrix_sum_pase_ops, 1, &SUBMATRIX_SUM);
    MPI_Iallreduce(MPI_IN_PLACE, inner_prod, 1, SUBMATRIX, SUBMATRIX_SUM, MPI_COMM_WORLD, &request);
#endif
    // y_aux_h_tmp = x_aux_h * y_aux_h
    int nrows = end[0] - start[0];
    int ncols = end[1] - start[1];
    double alpha = 1.0;
    double beta = 1.0;
    // TO DO nrows==ncols==1的情况
    pase_ops->lapack_ops->DenseMatQtAP('L', nsdIP, num_aux_vec, num_aux_vec, nrows, ncols, alpha,
                                       x_aux_h + start[0] * num_aux_vec, num_aux_vec,
                                       NULL, num_aux_vec,
                                       y_aux_h + start[1] * num_aux_vec, num_aux_vec,
                                       beta, a_tmp, LDA, NULL);
#if OPS_USE_MPI
    MPI_Wait(&request, &status);
    MPI_Op_free(&SUBMATRIX_SUM);
    MPI_Type_free(&SUBMATRIX);
#endif
    //计算 a += y_aux_h_tmp
    int inc = 1;
    int mid = ncols * ldIP;
    daxpy(&mid, &alpha, a_tmp, &inc, inner_prod, &inc);
    free(a_tmp);
    a_tmp = NULL;
#if AUX_TIME_INFO
    time_pase_innerprod_full += DefaultGetWtime();
#endif
    return;
}

void PASE_DefaultMultiVecSetRandomValue(void **multi_vec, int start, int end, struct PASE_OPS_ *pase_ops)
{
    pase_ops->gcge_ops->MultiVecSetRandomValue(((PASE_MultiVector)multi_vec)->b_H,
                                               start, end, pase_ops->gcge_ops);
    int i = 0;
    srand((unsigned)time(NULL));
    int num_aux_vec = ((PASE_MultiVector)multi_vec)->num_aux_vec;
    int num_vec = ((PASE_MultiVector)multi_vec)->num_vec;
    double *aux_h = ((PASE_MultiVector)multi_vec)->aux_h;

    for (i = start * num_aux_vec; i < end * num_aux_vec; i++)
    {
        aux_h[i] = ((double)rand()) / ((double)RAND_MAX + 1);
    }
    return;
}

void PASE_BLASMultiVecAxpby(double alpha, void **x, double beta, void **y,
                            int *start, int *end, struct PASE_OPS_ *pase_ops)
{
#if AUX_TIME_INFO
    extern double time_pase_axpby;
    time_pase_axpby -= DefaultGetWtime();
#endif
    assert(end[0] - start[0] == end[1] - start[1]);
    void **y_b_H = ((PASE_MultiVector)y)->b_H;
    double *y_aux_h = ((PASE_MultiVector)y)->aux_h;
    int num_aux_vec = ((PASE_MultiVector)y)->num_aux_vec;
    int length = num_aux_vec * (end[0] - start[0]);
    int inc = 1;

    double *destin = y_aux_h + start[1] * num_aux_vec;
    if (beta == 0.0)
    {
        memset(destin, 0, length * sizeof(double));
    }
    else if (beta != 1.0)
    {
        dscal(&length, &beta, destin, &inc);
    }

    if (x != NULL)
    {
        void **x_b_H = ((PASE_MultiVector)x)->b_H;
        double *x_aux_h = ((PASE_MultiVector)x)->aux_h;
        double *source = x_aux_h + start[0] * num_aux_vec;
        daxpy(&length, &alpha, source, &inc, destin, &inc);
#if OPS_USE_SLEPC
        double *xarray;
        double *yarray;
        BVGetArray((BV)x_b_H, &xarray);
        BVGetArray((BV)y_b_H, &yarray);
        int localrow;
        BVGetSizes((BV)x_b_H, &localrow, NULL, NULL);
        int num = localrow * (end[0] - start[0]);
        dscal(&num, &beta, yarray + start[1] * localrow, &one);
        daxpy(&num, &alpha, xarray + start[0] * localrow, &one, yarray + start[1] * localrow, &one);
        BVRestoreArray((BV)x_b_H, &xarray);
        BVRestoreArray((BV)y_b_H, &yarray);
#endif
    }
    else
    {
        pase_ops->gcge_ops->MultiVecAxpby(alpha, NULL, beta, y_b_H, start, end, pase_ops->gcge_ops);
    }
#if AUX_TIME_INFO
    time_pase_axpby += DefaultGetWtime();
#endif
    return;
}

void PASE_DefaultMultiVecAxpby(double alpha, void **x, double beta, void **y,
                               int *start, int *end, struct PASE_OPS_ *pase_ops)
{
    assert(end[0] - start[0] == end[1] - start[1]);
    void **y_b_H = ((PASE_MultiVector)y)->b_H;
    double *y_aux_h = ((PASE_MultiVector)y)->aux_h;
    int num_aux_vec = ((PASE_MultiVector)y)->num_aux_vec;
    int length = num_aux_vec * (end[0] - start[0]);
    int inc = 1;

    double *destin = y_aux_h + start[1] * num_aux_vec;
    if (beta == 0.0)
    {
        memset(destin, 0, length * sizeof(double));
    }
    else if (beta != 1.0)
    {
        dscal(&length, &beta, destin, &inc);
    }

    if (x != NULL)
    {
        void **x_b_H = ((PASE_MultiVector)x)->b_H;
        double *x_aux_h = ((PASE_MultiVector)x)->aux_h;
        double *source = x_aux_h + start[0] * num_aux_vec;
        daxpy(&length, &alpha, source, &inc, destin, &inc);
        pase_ops->gcge_ops->MultiVecAxpby(alpha, x_b_H, beta, y_b_H, start, end, pase_ops->gcge_ops);
    }
    else
    {
        pase_ops->gcge_ops->MultiVecAxpby(alpha, NULL, beta, y_b_H, start, end, pase_ops->gcge_ops);
    }
    return;
}

void PASE_DefaultMultiVecLinearComb(void **x, void **y, int is_vec, int *start, int *end,
                                    double *coef, int ldc, double *beta, int incb, struct PASE_OPS_ *pase_ops)
{
    assert(is_vec == 0);
    void **x_b_H;
    double *x_aux_h;
    void **y_b_H = ((PASE_MultiVector)y)->b_H;
    double *y_aux_h = ((PASE_MultiVector)y)->aux_h;
    int num_aux_vec = ((PASE_MultiVector)y)->num_aux_vec;
    LAPACKVEC y_vec;
    y_vec.nrows = num_aux_vec;
    y_vec.ncols = ((PASE_MultiVector)y)->num_vec;
    y_vec.ldd = num_aux_vec;
    y_vec.data = y_aux_h;
    if (x != NULL)
    {
        x_b_H = ((PASE_MultiVector)x)->b_H;
        x_aux_h = ((PASE_MultiVector)x)->aux_h;
        LAPACKVEC x_vec;
        x_vec.nrows = num_aux_vec;
        x_vec.ncols = ((PASE_MultiVector)x)->num_vec;
        x_vec.ldd = num_aux_vec;
        x_vec.data = x_aux_h;
        pase_ops->lapack_ops->MultiVecLinearComb((void **)&x_vec, (void **)&y_vec, is_vec,
                                                 start, end, coef, ldc, beta, incb, pase_ops->lapack_ops);
        pase_ops->gcge_ops->MultiVecLinearComb(x_b_H, y_b_H, is_vec, start, end,
                                               coef, ldc, beta, incb, pase_ops->gcge_ops);
    }
    else
    {
        pase_ops->lapack_ops->MultiVecLinearComb(NULL, (void **)&y_vec, is_vec,
                                                 start, end, coef, ldc, beta, incb, pase_ops->lapack_ops);
        pase_ops->gcge_ops->MultiVecLinearComb(NULL, y_b_H, is_vec, start, end,
                                               coef, ldc, beta, incb, pase_ops->gcge_ops);
    }
    //通过gcge_ops计算b_H部分的线性组合
    return;
}

static void MultiVecQtAPDiag(char ntsA, char ntsdQAP,
                             void **mvQ, void *matA, void **mvP, int is_vec,
                             int *startQP, int *endQP, double *qAp, int ldQAP,
                             void **mv_ws, struct PASE_OPS_ *pase_ops)
{
#if AUX_TIME_INFO
    extern double time_pase_QTAP_matmv;
    extern int number_pase_QTAP_diagmv;
    number_pase_QTAP_diagmv++;
    time_pase_QTAP_matmv -= DefaultGetWtime();
#endif
    int start[2], end[2], nrows, ncols;
    start[0] = startQP[1];
    end[0] = endQP[1];
    start[1] = 0;
    end[1] = endQP[1] - startQP[1];
#if MAT3
    if (ntsA = 'S' &&ntsdQAP = 'N' &&)
    {
    }
#endif
    if (ntsA == 'N' || ntsA == 'S')
    {
        pase_ops->MatDotMultiVec(matA, mvP, mv_ws, start, end, pase_ops);
    }
    else if (ntsA == 'T')
    {
        pase_ops->MatTransDotMultiVec(matA, mvP, mv_ws, start, end, pase_ops);
    }
#if AUX_TIME_INFO
    time_pase_QTAP_matmv += DefaultGetWtime();
    extern double time_pase_QTAP_inner;
    time_pase_QTAP_inner -= DefaultGetWtime();
#endif
    if (ntsdQAP == 'T')
    {
        start[0] = 0;
        end[0] = endQP[1] - startQP[1];
        start[1] = startQP[0];
        end[1] = endQP[0];
        pase_ops->MultiVecInnerProd('N', mv_ws, mvQ, is_vec, start, end, qAp, ldQAP, pase_ops);
    }
    else
    {
        start[0] = startQP[0];
        end[0] = endQP[0];
        start[1] = 0;
        end[1] = endQP[1] - startQP[1];
        pase_ops->MultiVecInnerProd(ntsdQAP, mvQ, mv_ws, is_vec, start, end, qAp, ldQAP, pase_ops);
    }
#if AUX_TIME_INFO
    time_pase_QTAP_inner += DefaultGetWtime();
#endif
    return;
}

static void MultiVecQtAPAll(char ntsA, char ntsdQAP,
                            void **mvQ, void *matA, void **mvP, int is_vec,
                            int *startQP, int *endQP, double *qAp, int ldQAP,
                            void **mv_ws, struct PASE_OPS_ *pase_ops)
{
#if AUX_TIME_INFO
    extern double time_pase_QTAP_matmv;
    extern int number_pase_QTAP_allmv;
    number_pase_QTAP_allmv++;
    time_pase_QTAP_matmv -= DefaultGetWtime();
#endif
    int start[2], end[2], nrows, ncols;
    start[0] = startQP[1];
    end[0] = endQP[1];
    start[1] = 0;
    end[1] = endQP[1] - startQP[1];
    if (ntsA == 'N' || ntsA == 'S')
    {
        pase_ops->MatDotMultiVec(matA, mvP, mv_ws, start, end, pase_ops);
    }
    else if (ntsA == 'T')
    {
        pase_ops->MatTransDotMultiVec(matA, mvP, mv_ws, start, end, pase_ops);
    }
#if AUX_TIME_INFO
    time_pase_QTAP_matmv += DefaultGetWtime();
    extern double time_pase_QTAP_inner;
    time_pase_QTAP_inner -= DefaultGetWtime();
#endif
    if (ntsdQAP == 'T')
    {
        start[0] = 0;
        end[0] = endQP[1] - startQP[1];
        start[1] = startQP[0];
        end[1] = endQP[0];
        pase_ops->MultiVecInnerProd('N', mv_ws, mvQ, is_vec, start, end, qAp, ldQAP, pase_ops);
    }
    else
    {
        start[0] = startQP[0];
        end[0] = endQP[0];
        start[1] = 0;
        end[1] = endQP[1] - startQP[1];
        pase_ops->MultiVecInnerProd(ntsdQAP, mvQ, mv_ws, is_vec, start, end, qAp, ldQAP, pase_ops);
    }
#if AUX_TIME_INFO
    time_pase_QTAP_inner += DefaultGetWtime();
#endif
    return;
}

void PASE_DefaultMultiVecQtAP(char ntsA, char ntsdQAP,
                              void **mvQ, void *matA, void **mvP, int is_vec,
                              int *startQP, int *endQP, double *qAp, int ldQAP,
                              void **mv_ws, struct PASE_OPS_ *pase_ops)
{
    assert(is_vec == 0);

    int start[2], end[2], nrows, ncols;
    nrows = endQP[0] - startQP[0];
    ncols = endQP[1] - startQP[1];
    if (nrows <= 0 || ncols <= 0)
        exit(0);
    if (matA == NULL)
    {
        if (ntsdQAP == 'T')
        {
            start[0] = startQP[1];
            end[0] = endQP[1];
            start[1] = startQP[0];
            end[1] = endQP[0];
            pase_ops->MultiVecInnerProd('N', mvP, mvQ, is_vec, start, end, qAp, ldQAP, pase_ops);
        }
        else
        {
            pase_ops->MultiVecInnerProd(ntsdQAP, mvQ, mvP, is_vec, startQP, endQP, qAp, ldQAP, pase_ops);
        }
        return;
    }
#if AUX_TIME_INFO
    extern double time_pase_QTAP;
    extern int number_pase_QTAP;
    number_pase_QTAP++;
    time_pase_QTAP -= DefaultGetWtime();
#endif
    if (((PASE_Matrix)matA)->is_diag == 1)
    {
        MultiVecQtAPDiag(ntsA, ntsdQAP, mvQ, matA, mvP, is_vec,
                         startQP, endQP, qAp, ldQAP, mv_ws, pase_ops);
    }
    else
    {
        MultiVecQtAPAll(ntsA, ntsdQAP, mvQ, matA, mvP, is_vec,
                        startQP, endQP, qAp, ldQAP, mv_ws, pase_ops);
    }
#if AUX_TIME_INFO
    time_pase_QTAP += DefaultGetWtime();
#endif
    return;
}

void PASE_DefaultMultiVecOrth(void **x, int start_x, int *end_x,
                              void *B, struct PASE_OPS_ *pase_ops)
{
    return;
}
