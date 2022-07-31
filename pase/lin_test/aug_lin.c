#include "aug_lin.h"

#define PRINT_INFO 1

#define MAX_ITS 2e4
#define SMT_ITS 30
#define GCGEPCGLEVEL(x) 2

char amg_type[32] = "petsc_gamg";
char aux_solver_type[32] = "pcg";

static void amg_petsc_gamg(AugLinSolver *solver, void *A, void **b, int *num_levels)
{
    PetscErrorCode ierr;
    Mat petsc_A = (Mat)A;
    Mat *A_array = NULL, *P_array = NULL;
    Mat *Aarr = NULL, *Parr = NULL;

    PC pc;
    PCCreate(PETSC_COMM_WORLD, &pc);
    PCSetOperators(pc, petsc_A, petsc_A);
    PCSetType(pc, PCGAMG);
    PCGAMGSetNlevels(pc, *num_levels);
    PCGAMGSetType(pc, PCGAMGAGG);
    PCGAMGSetRepartition(pc, PETSC_TRUE);
    PCSetUp(pc);
    PCGetCoarseOperators(pc, num_levels, &Aarr);
    PCGetInterpolations(pc, num_levels, &Parr);
#if PRINT_INFO
    solver->gcge_ops->Printf("[petsc_gamg] multigrid 层数: %d levels\n", *num_levels);
#endif

    A_array = (Mat *)malloc(*num_levels * sizeof(Mat));
    P_array = (Mat *)malloc((*num_levels - 1) * sizeof(Mat));
    A_array[0] = petsc_A;
    int idx_level, m, n;
    MatGetSize(A_array[0], &m, &n);
#if PRINT_INFO
    solver->gcge_ops->Printf("·            size of level %d: %d x %d\n", 0, m, n);
#endif
    for (idx_level = 1; idx_level < *num_levels; idx_level++)
    {
        P_array[idx_level - 1] = Parr[(*num_levels) - idx_level - 1];
#if PRINT_INFO
        MatGetSize(P_array[idx_level - 1], &m, &n);
        solver->gcge_ops->Printf("·             |  interpolation %d: %d x %d\n", idx_level - 1, m, n);
#endif
        A_array[idx_level] = Aarr[(*num_levels) - idx_level - 1];
        MatGetSize(A_array[idx_level], &m, &n);
#if PRINT_INFO
        solver->gcge_ops->Printf("·            size of level %d: %d x %d\n", idx_level, m, n);
#endif
    }
    PetscFree(Aarr);
    PetscFree(Parr);
    PCDestroy(&pc);
    solver->A_array = (void **)A_array;
    solver->P_array = (void **)P_array;

    BV *b_array = (BV *)malloc(*num_levels * sizeof(BV));
    b_array[0] = (BV)b;
    for (idx_level = 1; idx_level < *num_levels; idx_level++)
    {
        solver->gcge_ops->MultiVecCreateByMat((void ***)(&(b_array[idx_level])), 1, A_array[idx_level], solver->gcge_ops);
        BVSetActiveColumns(b_array[idx_level - 1], 0, 1);
        BVSetActiveColumns(b_array[idx_level], 0, 1);
        BVMatMultTranspose(b_array[idx_level - 1], P_array[idx_level - 1], b_array[idx_level]);
    }
    solver->b_array = (void ***)b_array;

    BV *x = (BV *)malloc(*num_levels * sizeof(BV));
    BV *cg_res = (BV *)malloc(*num_levels * sizeof(BV));
    BV *cg_p = (BV *)malloc(*num_levels * sizeof(BV));
    BV *cg_w = (BV *)malloc(*num_levels * sizeof(BV));
    OPS *gcge_ops = solver->gcge_ops;
    for (idx_level = 0; idx_level < *num_levels; idx_level++)
    {
        gcge_ops->MultiVecCreateByMat((void ***)(&(x[idx_level])), 1,
                                      (void **)A_array[idx_level], gcge_ops);
        gcge_ops->MultiVecCreateByMat((void ***)(&(cg_res[idx_level])), 1,
                                      (void **)A_array[idx_level], gcge_ops);
        gcge_ops->MultiVecCreateByMat((void ***)(&(cg_p[idx_level])), 1,
                                      (void **)A_array[idx_level], gcge_ops);
        gcge_ops->MultiVecCreateByMat((void ***)(&(cg_w[idx_level])), 1,
                                      (void **)A_array[idx_level], gcge_ops);
    }
    solver->x_array = (void ***)x;
    solver->cg_res = (void ***)cg_res;
    solver->cg_p = (void ***)cg_p;
    solver->cg_w = (void ***)cg_w;

    return;
}

void AugLinSolve(AugLinSolver *solver, void **x, char type[64])
{
    if (strcmp(type, "petsc_kspcg") == 0)
    {
        Petsc_KSPCG_Solver(solver, x);
        return;
    }
    if (strcmp(type, "gcge_pcg") == 0)
    {
        GCGE_PCG_Solver(solver, x);
        return;
    }
    if (strcmp(type, "aug_gcge_pcg") == 0)
    {
        Aug_Solver(solver, x, "gcge_pcg");
        return;
    }
}

void AugLinSolver_Create(AugLinSolver **solver, void *A, void **b, int num_levels, double tol)
{
    (*solver) = (AugLinSolver *)malloc(sizeof(AugLinSolver));

    OPS *gcge_ops;
    OPS_Create(&gcge_ops);
    OPS_SLEPC_Set(gcge_ops);
    OPS_Setup(gcge_ops);
    (*solver)->gcge_ops = gcge_ops;

    (*solver)->tol = tol;
    (*solver)->iterations = 0;
    (*solver)->double_tmp = (double *)malloc(128 * sizeof(double));
    (*solver)->int_tmp = (int *)malloc(128 * sizeof(int));

    if (num_levels == 1)
    {
        (*solver)->num_levels = 1;
        (*solver)->aux_level = -1;
        (*solver)->A_array = (void **)malloc(sizeof(void *));
        (*solver)->A_array[0] = A;
        (*solver)->P_array = NULL;
        (*solver)->b_array = (void ***)malloc(sizeof(void **));
        (*solver)->b_array[0] = b;
        (*solver)->x_array = (void ***)malloc(sizeof(void **));
        (*solver)->r_array = NULL;
        (*solver)->aux_A = NULL;
        (*solver)->aux_b = NULL;
        (*solver)->aux_x = NULL;
        (*solver)->cg_res = (void ***)malloc(sizeof(void **));
        (*solver)->cg_p = (void ***)malloc(sizeof(void **));
        (*solver)->cg_w = (void ***)malloc(sizeof(void **));
        gcge_ops->MultiVecCreateByMat(&((*solver)->cg_res[0]), 1, A, gcge_ops);
        gcge_ops->MultiVecCreateByMat(&((*solver)->cg_p[0]), 1, A, gcge_ops);
        gcge_ops->MultiVecCreateByMat(&((*solver)->cg_w[0]), 1, A, gcge_ops);
        return;
    }
    else
    {
        strcpy((*solver)->amg_type, amg_type);
        strcpy((*solver)->aux_solver_type, aux_solver_type);
        if (strcmp((*solver)->amg_type, "petsc_gamg") == 0)
        {
            amg_petsc_gamg((*solver), A, b, &num_levels);
        }
        (*solver)->num_levels = num_levels;
    }
}

void AugLinSolver_Free(AugLinSolver **solver)
{
    free((*solver)->A_array);
    free((*solver)->b_array);
    free((*solver)->double_tmp);
    free((*solver)->int_tmp);
    free(*solver);
    *solver = NULL;
}

void Petsc_KSPCG_Solver(AugLinSolver *solver, void **x)
{
    Vec Vec_b, Vec_x;
    BVCreateVec((BV)(solver->b_array[0]), &Vec_b);
    BVCreateVec((BV)x, &Vec_x);
    BVCopyVec((BV)(solver->b_array[0]), 0, Vec_b);

    KSP ksp_solver;
    PC pc;
    KSPCreate(PETSC_COMM_WORLD, &ksp_solver);
    KSPSetType(ksp_solver, KSPCG);
    KSPSetOperators(ksp_solver, (Mat)(solver->A_array[0]), (Mat)(solver->A_array[0]));
    KSPSetTolerances(ksp_solver, solver->tol, solver->tol, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetInitialGuessNonzero(ksp_solver, PETSC_TRUE);
    KSPSetFromOptions(ksp_solver);
    KSPSetUp(ksp_solver);
    KSPSolve(ksp_solver, Vec_b, Vec_x);
    KSPGetResidualNorm(ksp_solver, &(solver->residual));
    KSPGetTotalIterations(ksp_solver, &(solver->iterations));
    KSPDestroy(&ksp_solver);

    BVInsertVec((BV)x, 0, Vec_x);
}

void GCGE_PCG_Solver(AugLinSolver *solver, void **x)
{
    void *A = solver->A_array[0];
    void **b = solver->b_array[0];

    OPS *gcge_ops = solver->gcge_ops;
    int max_iter = MAX_ITS;
    double tol = solver->tol, rate = solver->tol * 1e-8;
    void **cg_space_mv[3];
    cg_space_mv[0] = solver->cg_res[0];
    cg_space_mv[1] = solver->cg_p[0];
    cg_space_mv[2] = solver->cg_w[0];
    double *double_tmp = solver->double_tmp;
    int *int_tmp = solver->int_tmp;
    MultiLinearSolverSetup_BlockPCG(max_iter, rate, tol, "abs",
                                    cg_space_mv, double_tmp, int_tmp, NULL, NULL, gcge_ops);
    int mv_s[2] = {0, 0};
    int mv_e[2] = {1, 1};
    gcge_ops->MultiLinearSolver((void *)A, (void **)b, (void **)x, mv_s, mv_e, gcge_ops);
}

static void Aug_Solver_smooth(AugLinSolver *solver, int its)
{
    void *A = solver->A_array[0];
    void **b = solver->b_array[0];
    void **x = solver->x_array[0];

    OPS *gcge_ops = solver->gcge_ops;
    int max_iter = its;
    double tol = solver->tol, rate = solver->tol * 1e-3;
    void **cg_space_mv[3];
    cg_space_mv[0] = solver->cg_res[0];
    cg_space_mv[1] = solver->cg_p[0];
    cg_space_mv[2] = solver->cg_w[0];
    double *double_tmp = solver->double_tmp;
    int *int_tmp = solver->int_tmp;
    MultiLinearSolverSetup_BlockPCG(max_iter, rate, tol, "abs",
                                    cg_space_mv, double_tmp, int_tmp, NULL, NULL, gcge_ops);
    int mv_s[2] = {0, 0};
    int mv_e[2] = {1, 1};
    gcge_ops->MultiLinearSolver((void *)A, (void **)b, (void **)x, mv_s, mv_e, gcge_ops);
}

static void Aug_Solver_aux_setup(AugLinSolver *solver)
{
    OPS *gcge_ops = solver->gcge_ops;
    PASE_MultiVector aux_b = solver->aux_b;
    PASE_MultiVector aux_x = solver->aux_x;

    int mv_s[2] = {0, 0};
    int mv_e[2] = {1, 1};

    aux_b->b_H = solver->b_array[solver->aux_level];
    gcge_ops->MultiVecInnerProd('N', solver->b_array[0], solver->x_array[0], 0, mv_s, mv_e, aux_b->aux_h, 1, gcge_ops);
    gcge_ops->MultiVecAxpby(0.0, aux_x->b_H, 0.0, aux_x->b_H, mv_s, mv_e, gcge_ops);
    aux_x->aux_h[0] = 1.0;

    /* a in aux_A */
    PASE_Matrix aux_A = solver->aux_A;
    gcge_ops->MatDotMultiVec(solver->A_array[0], solver->x_array[0], solver->cg_p[0], mv_s, mv_e, gcge_ops);
    int idx_level;
    for (idx_level = 1; idx_level <= solver->aux_level; idx_level++)
    {
        BVSetActiveColumns((BV)(solver->cg_p[idx_level - 1]), 0, 1);
        BVSetActiveColumns((BV)(solver->cg_p[idx_level]), 0, 1);
        BVMatMultTranspose((BV)(solver->cg_p[idx_level - 1]), (Mat)(solver->P_array[idx_level - 1]), (BV)(solver->cg_p[idx_level]));
    }
    aux_A->aux_Hh = solver->cg_p[solver->aux_level];
    /* alpha in aux_A */
    double *alpha = aux_A->aux_hh;
    gcge_ops->MultiVecQtAP('N', 'N', solver->x_array[0], solver->A_array[0], solver->x_array[0], 0, mv_s, mv_e,
                           aux_A->aux_hh, 1, solver->cg_res[0], gcge_ops);
}

static void Aug_Solver_aux_solve_gcge_pcg(AugLinSolver *solver, int its)
{
    OPS *pase_ops_to_gcge;
    OPS_Create(&pase_ops_to_gcge);
    GCGE_PASE_SetOps(pase_ops_to_gcge, solver->pase_ops);
    OPS_Setup(pase_ops_to_gcge);

    void *A = (void *)(solver->aux_A);
    void **b = (void **)(solver->aux_b);
    void **x = (void **)(solver->aux_x);

    int max_iter = its * 1e3;
    double tol = solver->tol * 1e-3, rate = solver->tol * 1e-6;
    void **cg_space_mv[3];
    cg_space_mv[0] = (void **)(solver->aux_gcg_mv[0]);
    cg_space_mv[1] = (void **)(solver->aux_gcg_mv[1]);
    cg_space_mv[2] = (void **)(solver->aux_gcg_mv[2]);
    double *double_tmp = solver->double_tmp;
    int *int_tmp = solver->int_tmp;
    MultiLinearSolverSetup_BlockPCG(max_iter, rate, tol, "abs",
                                    cg_space_mv, double_tmp, int_tmp, NULL, NULL, pase_ops_to_gcge);
    int mv_s[2] = {0, 0};
    int mv_e[2] = {1, 1};
    pase_ops_to_gcge->MultiLinearSolver((void *)A, (void **)b, (void **)x, mv_s, mv_e, pase_ops_to_gcge);
}

static void Aug_Solver_aux_prolong(AugLinSolver *solver)
{
    OPS *gcge_ops = solver->gcge_ops;

    PASE_MultiVector aux_x = solver->aux_x;
    int mv_s[2], mv_e[2];
    mv_s[0] = 0;
    mv_e[0] = 1;
    mv_s[1] = 0;
    mv_e[1] = 1;
    int idx_level;
    BVSetActiveColumns((BV)(solver->x_array[0]), 0, 1);
    BVSetActiveColumns((BV)(solver->cg_res[0]), 0, 1);
    BVCopy((BV)(solver->x_array[0]), (BV)(solver->cg_res[0]));
    for (idx_level = solver->aux_level; idx_level > 0; idx_level--)
    {
        BVSetActiveColumns((BV)(solver->x_array[idx_level - 1]), 0, 1);
        BVSetActiveColumns((BV)(solver->x_array[idx_level]), 0, 1);
        BVMatMult((BV)(solver->x_array[idx_level]), (Mat)(solver->P_array[idx_level - 1]), (BV)(solver->x_array[idx_level - 1]));
    }
    gcge_ops->MultiVecAxpby(aux_x->aux_h[0], solver->cg_res[0], 1.0, solver->x_array[0], mv_s, mv_e, gcge_ops);
}

void Aug_Solver(AugLinSolver *solver, void **x, char type[64])
{
    if (strcmp(type, "gcge_pcg") == 0)
    {
        solver->aux_level = GCGEPCGLEVEL(solver->num_levels);
        PetscPrintf(PETSC_COMM_WORLD, "[aux_level] %d\n", solver->aux_level);
    }
    OPS *gcge_ops = solver->gcge_ops;
    PASE_OPS *pase_ops;
    PASE_OPS_Create(&pase_ops, gcge_ops);
    solver->pase_ops = pase_ops;

    //空间准备
    solver->aux_A = (PASE_Matrix)malloc(sizeof(pase_Matrix));
    PASE_Matrix aux_A = solver->aux_A;
    aux_A->num_aux_vec = 1;
    aux_A->A_H = solver->A_array[solver->aux_level];
    aux_A->aux_Hh = solver->cg_res[solver->aux_level];
    aux_A->aux_hh = (double *)malloc(sizeof(double));
    aux_A->factorization = NULL;
    aux_A->AH_inv_aux_Hh = NULL;
    aux_A->aux_hh_inv = NULL;
    aux_A->if_sym = 1;
    aux_A->is_diag = 0;

    solver->aux_b = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    PASE_MultiVector aux_b = solver->aux_b;
    aux_b->num_aux_vec = 1;
    aux_b->num_vec = 1;
    aux_b->b_H = solver->b_array[solver->aux_level];
    aux_b->aux_h = (double *)malloc(sizeof(double));
    aux_b->aux_h_tmp = (double *)malloc(sizeof(double));

    solver->aux_x = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    PASE_MultiVector aux_x = solver->aux_x;
    aux_x->num_aux_vec = 1;
    aux_x->num_vec = 1;
    aux_x->b_H = solver->x_array[solver->aux_level];
    aux_x->aux_h = (double *)malloc(sizeof(double));
    aux_x->aux_h_tmp = (double *)malloc(sizeof(double));

    int aux_level = solver->aux_level;
    solver->aux_gcg_mv[0] = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    solver->aux_gcg_mv[0]->num_aux_vec = 1;
    solver->aux_gcg_mv[0]->num_vec = 1;
    gcge_ops->MultiVecCreateByMat(&(solver->aux_gcg_mv[0]->b_H), 1, solver->A_array[aux_level], gcge_ops);
    solver->aux_gcg_mv[0]->aux_h = (double *)malloc(sizeof(double));
    solver->aux_gcg_mv[0]->aux_h_tmp = (double *)malloc(sizeof(double));

    solver->aux_gcg_mv[1] = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    solver->aux_gcg_mv[1]->num_aux_vec = 1;
    solver->aux_gcg_mv[1]->num_vec = 1;
    gcge_ops->MultiVecCreateByMat(&(solver->aux_gcg_mv[1]->b_H), 1, solver->A_array[aux_level], gcge_ops);
    solver->aux_gcg_mv[1]->aux_h = (double *)malloc(sizeof(double));
    solver->aux_gcg_mv[1]->aux_h_tmp = (double *)malloc(sizeof(double));

    solver->aux_gcg_mv[2] = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    solver->aux_gcg_mv[2]->num_aux_vec = 1;
    solver->aux_gcg_mv[2]->num_vec = 1;
    gcge_ops->MultiVecCreateByMat(&(solver->aux_gcg_mv[2]->b_H), 1, solver->A_array[aux_level], gcge_ops);
    solver->aux_gcg_mv[2]->aux_h = (double *)malloc(sizeof(double));
    solver->aux_gcg_mv[2]->aux_h_tmp = (double *)malloc(sizeof(double));

    int idx;
    for (idx = 0; idx < 2; idx++)
    {
        //前光滑
        // gcge_ops->Printf("\n+++++++++++++++++++++++++++++++++++++++++++[PRE SMOOTH %d]+++++++++++++++++++++++++++++++++++++++++++++++++++\n", idx);
        Aug_Solver_smooth(solver, SMT_ITS);
        //在辅助空间求解
        Aug_Solver_aux_setup(solver);

        // Mat_Output(solver->aux_A->A_H, "A_H");
        // Mat_Output(solver->P_array[solver->aux_level - 1], "I^h_H");
        // MultiVec_Output(solver->x_array[0], 1, "x_hh");
        // MultiVec_Output(solver->aux_A->aux_Hh, 1, "IA_h_x_h");
        // gcge_ops->Printf("x_hTA_hx_h (1, 1) = %e\n", solver->aux_A->aux_hh[0]);
        // MultiVec_Output(solver->aux_x->b_H, 1, "x_H");
        // gcge_ops->Printf("gamma (1, 1) = %e\n", solver->aux_x->aux_h[0]);
        // MultiVec_Output(solver->aux_b->b_H, 1, "b_H");
        // MultiVec_Output(solver->b_array[0], 1, "b_h");
        // gcge_ops->Printf("x_h^Tb_h (1, 1) = %e\n", solver->aux_b->aux_h[0]);

        extern int if_outitr;
        if_outitr=0;
        if (strcmp(type, "gcge_pcg") == 0)
        {
            // gcge_ops->Printf("\n+++++++++++++++++++++++++++++++++++++++++++[AUX %d]+++++++++++++++++++++++++++++++++++++++++++++++++++\n", idx);
            Aug_Solver_aux_solve_gcge_pcg(solver, 100);
        }
        Aug_Solver_aux_prolong(solver);
        if_outitr=1;

        // MultiVec_Output(solver->aux_x->b_H, 1, "x_H_new");
        // gcge_ops->Printf("gamma_new (1, 1) = %e\n", solver->aux_x->aux_h[0]);
        // MultiVec_Output(solver->x_array[0], 1, "x_h_new");

        //后光滑
        // gcge_ops->Printf("\n+++++++++++++++++++++++++++++++++++++++++++[POST SMOOTH %d]+++++++++++++++++++++++++++++++++++++++++++++++++++\n", idx);
        Aug_Solver_smooth(solver, SMT_ITS);
    }
}