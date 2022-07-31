#include "pase.h"

typedef struct AugLinSolver_
{
    OPS *gcge_ops;
    PASE_OPS *pase_ops;
    char amg_type[16];
    char aux_solver_type[16];

    int num_levels;
    int aux_level;
    void **A_array;
    void **P_array;
    void ***b_array;
    void ***x_array;
    void ***r_array;
    PASE_Matrix aux_A;
    PASE_MultiVector aux_b;
    PASE_MultiVector aux_x;
    double tol;

    double residual;
    int iterations;
    int pre_smooth_iter_num;
    int post_smooth_iter_num;

    void ***cg_res;
    void ***cg_p;
    void ***cg_w;
    double *double_tmp;
    int *int_tmp;
    PASE_MultiVector aux_gcg_mv[3];
} AugLinSolver;

void AugLinSolver_Create(AugLinSolver **solver, void *A, void **b, int num_levels, double residual);
void AugLinSolve(AugLinSolver *solver, void **x, char type[64]);
void AugLinSolver_Free(AugLinSolver **solver);

void Petsc_KSPCG_Solver(AugLinSolver *solver, void **x);
void GCGE_PCG_Solver(AugLinSolver *solver, void **x);

void Aug_Solver(AugLinSolver *solver, void **x, char type[64]);