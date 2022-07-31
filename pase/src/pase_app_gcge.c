#include "pase_app_gcge.h"

int pase_app_gcge_flag = 0;

void PASE_GCGE_MatView(void *mat, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMatView(mat, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MatAxpby(double alpha, void *matX, double beta, void *matY, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMatAxpby(alpha, matX, beta, matY, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MatDotMultiVec(void *mat, void **x, void **y,
                              int *start, int *end, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMatDotMultiVec(mat, x, y, start, end, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MatTransDotMultiVec(void *mat, void **x, void **y,
                                   int *start, int *end, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMatTransDotMultiVec(mat, x, y, start, end, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MatDotVec(void *mat, void *x, void *y, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMatDotVec(mat, x, y, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MatTransDotVec(void *mat, void *x, void *y, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMatTransDotVec(mat, x, y, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MultiVecCreateByMat(void ***multi_vec, int num_vec, void *src_mat, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMultiVecCreateByMat(multi_vec, num_vec, src_mat, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_VecCreateByMat(void **des_vec, void *src_mat, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultVecCreateByMat(des_vec, src_mat, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_VecCreateByVec(void **des_vec, void *src_vec, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultVecCreateByVec(des_vec, src_vec, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_VecDestroy(void **des_vec, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultVecDestroy(des_vec, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_VecView(void *x, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultVecView(x, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_VecAxpby(double alpha, void *x, double beta, void *y, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultVecAxpby(alpha, x, beta, y, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_VecInnerProd(void *x, void *y, double *inner_prod, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultVecInnerProd(x, y, inner_prod, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_VecLocalInnerProd(void *x, void *y, double *inner_prod, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultVecLocalInnerProd(x, y, inner_prod, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_VecSetRandomValue(void *x, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultVecSetRandomValue(x, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MultiVecCreateByVec(void ***multi_vec, int num_vec, void *src_vec, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMultiVecCreateByVec(multi_vec, num_vec, src_vec, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MultiVecCreateByMultiVec(void ***multi_vec, int num_vec, void **src_mv, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMultiVecCreateByMultiVec(multi_vec, num_vec, src_mv, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MultiVecDestroy(void ***multi_vec, int num_vec, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMultiVecDestroy(multi_vec, num_vec, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_GetVecFromMultiVec(void **multi_vec, int col, void **vec, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultGetVecFromMultiVec(multi_vec, col, vec, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_RestoreVecForMultiVec(void **multi_vec, int col, void **vec, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultRestoreVecForMultiVec(multi_vec, col, vec, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MultiVecView(void **x, int start, int end, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMultiVecView(x, start, end, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MultiVecLocalInnerProd(char nsdIP, void **x, void **y, int is_vec, int *start, int *end,
                                      double *inner_prod, int ldIP, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMultiVecLocalInnerProd(nsdIP, x, y, is_vec, start, end,
                                       inner_prod, ldIP, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MultiVecInnerProd(char nsdIP, void **x, void **y, int is_vec, int *start, int *end,
                                 double *inner_prod, int ldIP, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMultiVecInnerProd(nsdIP, x, y, is_vec, start, end,
                                  inner_prod, ldIP, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MultiVecSetRandomValue(void **multi_vec, int start, int end, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMultiVecSetRandomValue(multi_vec, start, end, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MultiVecAxpby(double alpha, void **x, double beta, void **y,
                             int *start, int *end, struct OPS_ *pase_ops_gcge)
{
    PASE_BLASMultiVecAxpby(alpha, x, beta, y, start, end, (PASE_OPS *)(pase_ops_gcge->pase_ops));
    // PASE_DefaultMultiVecAxpby(alpha, x, beta, y, start, end, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MultiVecLinearComb(void **x, void **y, int is_vec, int *start, int *end,
                                  double *coef, int ldc, double *beta, int incb, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMultiVecLinearComb(x, y, is_vec, start, end,
                                   coef, ldc, beta, incb, (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void PASE_GCGE_MultiVecQtAP(char ntsA, char ntsdQAP,
                            void **mvQ, void *matA, void **mvP, int is_vec,
                            int *startQP, int *endQP, double *qAp, int ldQAP,
                            void **mv_ws, struct OPS_ *pase_ops_gcge)
{
    PASE_DefaultMultiVecQtAP(ntsA, ntsdQAP, mvQ, matA, mvP,
                             is_vec, startQP, endQP, qAp, ldQAP, mv_ws,
                             (PASE_OPS *)(pase_ops_gcge->pase_ops));
}

void GCGE_PASE_SetOps(OPS *pase_ops_gcge, PASE_OPS *pase_ops)
{
    pase_ops_gcge->pase_ops = (void *)pase_ops;
    pase_ops_gcge->MatView = PASE_GCGE_MatView;
    pase_ops_gcge->MultiVecCreateByMat = PASE_GCGE_MultiVecCreateByMat;
    pase_ops_gcge->MatAxpby = PASE_GCGE_MatAxpby;
    pase_ops_gcge->MatDotMultiVec = PASE_GCGE_MatDotMultiVec;
    pase_ops_gcge->MatTransDotMultiVec = PASE_GCGE_MatTransDotMultiVec;
    pase_ops_gcge->MatDotVec = PASE_GCGE_MatDotVec;
    pase_ops_gcge->MatTransDotVec = PASE_GCGE_MatTransDotVec;
    pase_ops_gcge->VecCreateByMat = PASE_GCGE_VecCreateByMat;
    pase_ops_gcge->VecCreateByVec = PASE_GCGE_VecCreateByVec;
    pase_ops_gcge->VecDestroy = PASE_GCGE_VecDestroy;
    pase_ops_gcge->VecView = PASE_GCGE_VecView;
    pase_ops_gcge->VecAxpby = PASE_GCGE_VecAxpby;
    pase_ops_gcge->VecInnerProd = PASE_GCGE_VecInnerProd;
    pase_ops_gcge->VecLocalInnerProd = PASE_GCGE_VecLocalInnerProd;
    pase_ops_gcge->VecSetRandomValue = PASE_GCGE_VecSetRandomValue;
    pase_ops_gcge->MultiVecCreateByVec = PASE_GCGE_MultiVecCreateByVec;
    pase_ops_gcge->MultiVecCreateByMultiVec = PASE_GCGE_MultiVecCreateByMultiVec;
    pase_ops_gcge->MultiVecDestroy = PASE_GCGE_MultiVecDestroy;
    pase_ops_gcge->GetVecFromMultiVec = PASE_GCGE_GetVecFromMultiVec;
    pase_ops_gcge->RestoreVecForMultiVec = PASE_GCGE_RestoreVecForMultiVec;
    pase_ops_gcge->MultiVecView = PASE_GCGE_MultiVecView;
    pase_ops_gcge->MultiVecLocalInnerProd = PASE_GCGE_MultiVecLocalInnerProd;
    pase_ops_gcge->MultiVecInnerProd = PASE_GCGE_MultiVecInnerProd;
    pase_ops_gcge->MultiVecSetRandomValue = PASE_GCGE_MultiVecSetRandomValue;
    pase_ops_gcge->MultiVecAxpby = PASE_GCGE_MultiVecAxpby;
    pase_ops_gcge->MultiVecLinearComb = PASE_GCGE_MultiVecLinearComb;
    pase_ops_gcge->MultiVecQtAP = PASE_GCGE_MultiVecQtAP;

    return;
}

void PASE_Diag_BlockPCG(void *mat, void **mv_b, void **mv_x,
                        int *start_bx, int *end_bx, struct OPS_ *pase_ops_to_gcge)
{
    PASE_OPS *pase_ops = pase_ops_to_gcge->pase_ops;
    OPS *gcge_ops = pase_ops->gcge_ops;

    PASE_Matrix paseA = (PASE_Matrix)mat;
    void *A = paseA->A_H;
    void **A_inv_a = paseA->AH_inv_aux_Hh;
    double *alpha_inv = paseA->aux_hh_inv;
    PASE_MultiVector paseb = (PASE_MultiVector)mv_b;
    void *b = paseb->b_H;
    double *beta = paseb->aux_h;
    double *beta_tmp = paseb->aux_h_tmp;
    PASE_MultiVector pasex = (PASE_MultiVector)mv_x;
    void *x = pasex->b_H;
    double *gamma = pasex->aux_h;

    int size = paseA->num_aux_vec;
    int length = end_bx[0] - start_bx[0];
    int num = size * length;
    double p_one = 1.0, m_one = -1.0, zero = 0.0;
    int i_one = 1;
    int mv_s[2] = {0, start_bx[0]};
    int mv_e[2] = {paseA->num_aux_vec, end_bx[0]};
    char charN = 'N';
    gcge_ops->MultiVecInnerProd('N', A_inv_a, b, 0, mv_s, mv_e, beta_tmp, paseA->num_aux_vec, gcge_ops);
    daxpy(&num, &m_one, beta_tmp, &i_one, beta + (start_bx[0] * paseb->num_aux_vec), &i_one);
    dgemm(&charN, &charN, &size, &length, &size, &p_one, alpha_inv, &size,
          beta + (start_bx[0] * size), &size, &zero, gamma + (start_bx[1] * size), &size);

    PASE_BlockPCGSolver *pase_bpcg = (PASE_BlockPCGSolver *)(pase_ops_to_gcge->multi_linear_solver_workspace);
    PASE_MultiVector pase_mv_ws[3];
    pase_mv_ws[0] = (PASE_MultiVector)(pase_bpcg->mv_ws[0]);
    pase_mv_ws[1] = (PASE_MultiVector)(pase_bpcg->mv_ws[1]);
    pase_mv_ws[2] = (PASE_MultiVector)(pase_bpcg->mv_ws[2]);
    void **mv_ws[3];
    mv_ws[0] = pase_mv_ws[0]->b_H;
    mv_ws[1] = pase_mv_ws[1]->b_H;
    mv_ws[2] = pase_mv_ws[2]->b_H;
    MultiLinearSolverSetup_BlockPCG(pase_bpcg->max_iter, pase_bpcg->rate, pase_bpcg->tol, pase_bpcg->tol_type,
                                    mv_ws, pase_bpcg->dbl_ws, pase_bpcg->int_ws, pase_bpcg->pc,
                                    pase_bpcg->MatDotMultiVec, gcge_ops);
    pase_app_gcge_flag = 1;
    gcge_ops->MultiLinearSolver(A, b, x, start_bx, end_bx, gcge_ops);
    pase_app_gcge_flag = 0;

    const double *A_inv_a_tmp, *x_tmp;
    int nrows, ncols;
    BVGetArrayRead((BV)A_inv_a, &A_inv_a_tmp);
    BVGetArrayRead((BV)x, &x_tmp);
    BVGetSizes((BV)A_inv_a, &nrows, NULL, &ncols);
    dgemm(&charN, &charN, &nrows, &length, &ncols, &m_one,
          (double *)A_inv_a_tmp, &nrows, gamma + (start_bx[1] * size), &size,
          &p_one, (double *)x_tmp + (start_bx[1] * nrows), &nrows);
    BVRestoreArrayRead((BV)A_inv_a, &A_inv_a_tmp);
    BVRestoreArrayRead((BV)x, &x_tmp);

}

void PASE_Diag_BlockPCG_Setup(int max_iter, double rate, double tol,
                              const char *tol_type, void **mv_ws[3], double *dbl_ws, int *int_ws,
                              void *pc, void (*MatDotMultiVec)(void **x, void **y, int *, int *, void **z, int s, struct OPS_ *),
                              struct OPS_ *pase_ops_to_gcge)
{
    static PASE_BlockPCGSolver pase_bpcg_static = {
        .max_iter = 50, .rate = 1e-2, .tol = 1e-12, .tol_type = "abs", .mv_ws = {}, .pc = NULL, .MatDotMultiVec = NULL};
    pase_bpcg_static.max_iter = max_iter;
    pase_bpcg_static.rate = rate;
    pase_bpcg_static.tol = tol;
    strcpy(pase_bpcg_static.tol_type, tol_type);
    pase_bpcg_static.mv_ws[0] = mv_ws[0];
    pase_bpcg_static.mv_ws[1] = mv_ws[1];
    pase_bpcg_static.mv_ws[2] = mv_ws[2];
    pase_bpcg_static.dbl_ws = dbl_ws;
    pase_bpcg_static.int_ws = int_ws;
    pase_bpcg_static.MatDotMultiVec = MatDotMultiVec;

    pase_bpcg_static.niter = 0;
    pase_bpcg_static.residual = -1.0;

    pase_ops_to_gcge->multi_linear_solver_workspace = (void *)(&pase_bpcg_static);
    pase_ops_to_gcge->MultiLinearSolver = PASE_Diag_BlockPCG;
    return;
}
