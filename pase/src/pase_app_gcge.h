#ifndef __PASE_APP_GCGE_H__
#define __PASE_APP_GCGE_H__

#include "ops.h"
#include "ops_lin_sol.h"
#include "pase_ops.h"

typedef BlockPCGSolver PASE_BlockPCGSolver;

void GCGE_PASE_SetOps(OPS *ops, PASE_OPS *pase_ops);

void PASE_Diag_BlockPCG_Setup(int max_iter, double rate, double tol,
                              const char *tol_type, void **mv_ws[3], double *dbl_ws, int *int_ws,
                              void *pc, void (*MatDotMultiVec)(void **x, void **y, int *, int *, void **z, int s, struct OPS_ *),
                              struct OPS_ *ops);

#endif