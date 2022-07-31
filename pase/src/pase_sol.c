#include "pase_sol.h"
#include "pase_app_gcge.h"

#define PRINT_INFO 1

int info;
int i_one = 1;
double d_zero = 0.0;
double p_one = 1.0;
double m_one = -1.0;
double timer_start, timer_end;

int A_pre = 0;

int PASE_EigenSolver(void *A, void *B, double *eval, void **evec, int nev, PASE_PARAMETER param)
{
    // create gcge_ops
    OPS *gcge_ops;
    OPS_Create(&gcge_ops);
    OPS_SLEPC_Set(gcge_ops);
    OPS_Setup(gcge_ops);

    // create pase_ops
    PASE_OPS *pase_ops;
    PASE_OPS_Create(&pase_ops, gcge_ops);

    // create solver
#if PRINT_INFO
    gcge_ops->Printf("\n============== [ solver creation ] ==============\n");
#endif
    PASE_MG_SOLVER solver = PASE_Mg_solver_create(param, pase_ops);
#if PRINT_INFO
    gcge_ops->Printf("=================================================\n\n");
#endif

    solver->total_time -= gcge_ops->GetWtime();
#if PRINT_INFO
    gcge_ops->Printf("=============== [ solver setup ] ================\n");
#endif
    PASE_Mg_set_up(solver, A, B);
#if PRINT_INFO
    if (solver->print_level > 0)
    {
        PASE_Mg_print_param(solver);
    }
    gcge_ops->Printf("=================================================\n\n");
    gcge_ops->Printf("=================== [ solve ] ===================\n");
#endif
    solver->total_solve_time -= gcge_ops->GetWtime();
    PASE_Mg_solve(solver);
    solver->total_solve_time += gcge_ops->GetWtime();
    solver->total_time += gcge_ops->GetWtime();
#if PRINT_INFO
    gcge_ops->Printf("=================================================\n\n");
#endif
    PASE_Mg_print_timer(solver);
    PASE_Mg_solver_destroy(solver);
}

int PASE_Mg_solve(PASE_MG_SOLVER solver)
{
    OPS *gcge_ops = solver->gcge_ops;
    if (!solver->num_given_eigs)
    {
#if PRINT_INFO
        gcge_ops->Printf("\n-------- initial vecs ----------\n");
#endif
        solver->get_initvec_time -= gcge_ops->GetWtime();
        solver->current_level = solver->initial_level;
        PASE_Direct_solve(solver);
        solver->get_initvec_time += gcge_ops->GetWtime();
#if PRINT_INFO
        gcge_ops->Printf("\n[direct solve] get initial vecs time : %f sec.\n", solver->get_initvec_time);
        gcge_ops->Printf("-------------------------------\n\n");
#endif
    }

    solver->prolong_time -= gcge_ops->GetWtime();
    PASE_Mg_prolong_from_Solution(solver, solver->aux_fine_level);
    solver->current_level = solver->aux_fine_level;
    solver->prolong_time += gcge_ops->GetWtime();

    int idx_cycle;
    for (idx_cycle = 0; idx_cycle < solver->max_cycle_count; idx_cycle++)
    {
        solver->current_cycle = idx_cycle;
#if PRINT_INFO
        gcge_ops->Printf("\n===================== cycle [%d] =====================\n\n", idx_cycle);
#endif
        PASE_Mg_cycle(solver);
#if PRINT_INFO
        gcge_ops->Printf("\n------------------------------\n");
#endif
        solver->error_estimate_time -= gcge_ops->GetWtime();
        PASE_Mg_error_estimate(solver);
        solver->error_estimate_time += gcge_ops->GetWtime();
        if (solver->conv_nev >= solver->user_nev)
        {
            break;
        }
    }
}

int PASE_Mg_error_estimate(PASE_MG_SOLVER solver)
{
    OPS *gcge_ops = solver->gcge_ops;

    int pase_nev = solver->pase_nev;
    int conv_nev = solver->conv_nev;
    int current_level = solver->current_level;
    double atol = solver->atol, rtol = solver->rtol;
    double *eigenvalue = solver->eigenvalues;
    void **solution = solver->solution[current_level];
    void **rhs = solver->cg_rhs[current_level];
    void **res = solver->cg_res[current_level];
    void *A = solver->multigrid->A_array[current_level];
    void *B = solver->multigrid->B_array[current_level];
#if PRINT_INFO
    gcge_ops->Printf("[error estimate] cycle %d\n", solver->current_cycle);
    gcge_ops->Printf("·                conv_nev = %d\n", conv_nev);
#endif

    int flag = 0;
    double *check_multi = (double *)malloc((pase_nev - 1) * sizeof(double));
    double r_norm, u_norm;
    int mv_s[2], mv_e[2];

    int idx_eig;
    for (idx_eig = 0; idx_eig < pase_nev; idx_eig++)
    {
        mv_s[0] = idx_eig;
        mv_e[0] = idx_eig + 1;
        mv_s[1] = 0;
        mv_e[1] = 01;
        gcge_ops->MatDotMultiVec(B, solution, rhs, mv_s, mv_e, gcge_ops);
        gcge_ops->MultiVecInnerProd('D', solution, rhs, 0, mv_s, mv_e, &u_norm, 1, gcge_ops);
        // DefaultMultiVecInnerProd('D', solution, rhs, 0, mv_s, mv_e, &u_norm, 1, gcge_ops);
        u_norm = sqrt(u_norm);
#if PRINT_INFO
        gcge_ops->Printf("·                eigenvalue %d : %e\n", idx_eig, eigenvalue[idx_eig]);
        gcge_ops->Printf("·                | u_norm : %e\n", u_norm);
#endif
        gcge_ops->MatDotMultiVec(A, solution, res, mv_s, mv_e, gcge_ops);
        mv_s[0] = 0;
        mv_e[0] = 1;
        mv_s[1] = 0;
        mv_e[1] = 1;
        gcge_ops->MultiVecAxpby(1.0, res, -eigenvalue[idx_eig], rhs, mv_s, mv_e, gcge_ops);
        gcge_ops->MultiVecInnerProd('N', rhs, rhs, 0, mv_s, mv_e, &r_norm, 1, gcge_ops);
        // DefaultMultiVecInnerProd('N', rhs, rhs, 0, mv_s, mv_e, &r_norm, 1, gcge_ops);
        r_norm = sqrt(r_norm) / u_norm;
#if PRINT_INFO
        gcge_ops->Printf("·                | r_norm : %e\n", r_norm);
#endif
        solver->abs_res_norm[idx_eig] = r_norm;
        if (idx_eig + 1 < pase_nev)
        {
            check_multi[idx_eig] = fabs((eigenvalue[idx_eig] - eigenvalue[idx_eig + 1]) / eigenvalue[idx_eig]);
        }
        if (r_norm < atol || (r_norm / eigenvalue[idx_eig]) < rtol)
        {
            if (flag == 0)
            {
#if PRINT_INFO
                gcge_ops->Printf("·                | 收敛！ : %e\n", r_norm);
#endif
                solver->conv_nev = idx_eig + 1;
            }
        }
        else
        {
            flag = 1;
        }
    }
    if (solver->conv_nev > 1 && solver->conv_nev < pase_nev)
    {
        while (check_multi[solver->conv_nev - 1] < 1e-5 && solver->conv_nev > conv_nev)
        {
#if PRINT_INFO
            gcge_ops->Printf("[error estimate] No.%d 和 No.%d 判定重根！\n", solver->conv_nev - 1, solver->conv_nev);
#endif
            solver->conv_nev--;
            if (solver->conv_nev < 1 && solver->conv_nev >= pase_nev)
            {
                break;
            }
        }
    }
    free(check_multi);
    check_multi = NULL;
    solver->nlock_smooth = solver->conv_nev;
    solver->nlock_direct = 0;

#if PRINT_INFO
    gcge_ops->Printf("[error estimate] final conv_nev = %d\n", solver->conv_nev);
#endif
}

int PASE_Mg_cycle(PASE_MG_SOLVER solver)
{
    OPS *gcge_ops = solver->gcge_ops;
#if PRINT_INFO
    gcge_ops->Printf("(0) pre-smoothing:\n");
    double timer = solver->smooth_time;
#endif
    solver->smooth_time -= gcge_ops->GetWtime();
    PASE_Mg_smoothing(solver);
    solver->smooth_time += gcge_ops->GetWtime();
#if PRINT_INFO
    gcge_ops->Printf("-   pre-smoothing time : %f sec\n\n", solver->smooth_time - timer);
    gcge_ops->Printf("(1) aux space setup:\n");
    timer = solver->build_aux_time;
#endif
    solver->build_aux_time -= gcge_ops->GetWtime();
    PASE_Mg_set_pase_aux_vector(solver);
    PASE_Mg_set_pase_aux_matrix(solver);
    solver->current_level = solver->aux_coarse_level;
    solver->build_aux_time += gcge_ops->GetWtime();
#if PRINT_INFO
    gcge_ops->Printf("\n-   build aux space time : %f sec\n\n", solver->build_aux_time - timer);
    gcge_ops->Printf("(2) aux solver:\n");
    timer = solver->aux_direct_solve_time;
#endif
    solver->aux_direct_solve_time -= gcge_ops->GetWtime();
    PASE_Aux_direct_solve(solver);
    solver->aux_direct_solve_time += gcge_ops->GetWtime();
#if PRINT_INFO
    gcge_ops->Printf("\n-   aux solver time : %f sec\n\n", solver->aux_direct_solve_time - timer);
    gcge_ops->Printf("(3) prolong to fine level:\n");
    timer = solver->prolong_time;
#endif
    solver->prolong_time -= gcge_ops->GetWtime();
    PASE_Mg_prolong_from_pase_aux_Solution(solver, solver->aux_fine_level);
    solver->current_level = solver->aux_fine_level;
    solver->prolong_time += gcge_ops->GetWtime();
#if PRINT_INFO
    gcge_ops->Printf("-   prolong time : %f sec\n\n", solver->prolong_time - timer);
    gcge_ops->Printf("(4) post smoothing:\n");
#endif
    solver->smooth_time -= gcge_ops->GetWtime();
    PASE_Mg_smoothing(solver);
    solver->smooth_time += gcge_ops->GetWtime();
#if PRINT_INFO
    gcge_ops->Printf("-   post-smoothing time : %f sec\n\n", solver->smooth_time - timer);
#endif

    return 0;
}

int PASE_Mg_prolong_from_pase_aux_Solution(PASE_MG_SOLVER solver, int object_level)
{
    OPS *gcge_ops = solver->gcge_ops;

    PASE_MultiVector aux_u = solver->aux_sol;
    void **u = aux_u->b_H;
    double *gamma = aux_u->aux_h;
    int nlock_direct = solver->nlock_direct;
    int conv_nev = solver->conv_nev;
    int pase_nev = solver->pase_nev;
    int current_level = solver->current_level;
    int mv_s[2], mv_e[2];

    if (solver->if_precondition == 1 || solver->if_precondition == 3)
    {
        solver->precondition_time -= gcge_ops->GetWtime();
        PASE_Matrix aux_B = solver->aux_B;
        void **B_inv_b = aux_B->AH_inv_aux_Hh;
        int num = pase_nev * pase_nev;
        dcopy(&num, gamma, &i_one, aux_u->aux_h_tmp, &i_one);
        dscal(&num, &m_one, aux_u->aux_h_tmp, &i_one);
        mv_s[0] = nlock_direct;
        mv_e[0] = pase_nev;
        mv_s[1] = nlock_direct;
        mv_e[1] = pase_nev;
        gcge_ops->MultiVecLinearComb(B_inv_b, u, 0, mv_s, mv_e,
                                     aux_u->aux_h_tmp, pase_nev, &p_one, 0, gcge_ops);
        solver->precondition_time += gcge_ops->GetWtime();
    }

    void **u_H2h = solver->cg_rhs[object_level];
    void **solution = solver->solution[object_level];
    mv_s[0] = nlock_direct;
    mv_e[0] = pase_nev;
    mv_s[1] = nlock_direct;
    mv_e[1] = pase_nev;
    PASE_MULTIGRID_fromItoJ(solver->multigrid,
                            current_level, object_level, mv_s, mv_e, u, u_H2h);
    gcge_ops->MultiVecLinearComb(solution, u_H2h, 0, mv_s, mv_e,
                                 gamma + nlock_direct * pase_nev, pase_nev, &p_one, 0, gcge_ops);
    BVSetActiveColumns((BV)u_H2h, conv_nev, pase_nev);
    BVSetActiveColumns((BV)solution, conv_nev, pase_nev);
    BVCopy((BV)u_H2h, (BV)solution);
}

int PASE_Aux_direct_solve(PASE_MG_SOLVER solver)
{
    void *A = solver->aux_A;
    void *B = solver->aux_B;
    double *eigenvalues = solver->aux_eigenvalues;
    void **eigenvectors = (void **)(solver->aux_sol);
    OPS *gcge_ops = solver->gcge_ops;
    OPS *pase_ops_to_gcge = solver->pase_ops_to_gcge;

    int nevConv = solver->user_nev, multiMax = 1;
    double gapMin = 1e-5;
    int nevGiven = solver->pase_nev;
    int block_size = nevConv / 2, nevMax = 2 * nevConv;
    int nevInit = 2 * nevConv;
    nevInit = nevInit < nevMax ? nevInit : nevMax;
    int max_iter_gcg = solver->max_initial_direct_count;
    double tol_gcg[2] = {solver->aux_rtol, solver->aux_rtol};
    void **gcg_mv_ws[4];
    gcg_mv_ws[0] = (void **)(solver->aux_gcg_mv[0]);
    gcg_mv_ws[1] = (void **)(solver->aux_gcg_mv[1]);
    gcg_mv_ws[2] = (void **)(solver->aux_gcg_mv[2]);
    gcg_mv_ws[3] = (void **)(solver->aux_gcg_mv[3]);
    double *dbl_ws = solver->double_tmp;
    int *int_ws = solver->int_tmp;

    srand(0);
    EigenSolverSetup_GCG(multiMax, gapMin, nevInit, nevMax, block_size,
                         tol_gcg, max_iter_gcg, 0, gcg_mv_ws, dbl_ws, int_ws, pase_ops_to_gcge);

    int check_conv_max_num = 50;

    char initX_orth_method[8] = "mgs";
    int initX_orth_block_size = -1;
    int initX_orth_max_reorth = 2;
    double initX_orth_zero_tol = 2 * DBL_EPSILON; // 1e-12

    char compP_orth_method[8] = "mgs";
    int compP_orth_block_size = -1;
    int compP_orth_max_reorth = 2;
    double compP_orth_zero_tol = 2 * DBL_EPSILON; // 1e-12

    char compW_orth_method[8] = "mgs";
    int compW_orth_block_size = -1;
    int compW_orth_max_reorth = 2;
    double compW_orth_zero_tol = 2 * DBL_EPSILON; // 1e-12
    int compW_bpcg_max_iter = 50;
    double compW_bpcg_rate = 1e-2;
    double compW_bpcg_tol = 1e-14;
    char compW_bpcg_tol_type[8] = "abs";

    int compRR_min_num = -1;
    double compRR_min_gap = gapMin;
    double compRR_tol = 2 * DBL_EPSILON;
#if PRINT_INFO
    gcge_ops->Printf("[aux direct solve] level :\n", solver->aux_coarse_level);
    gcge_ops->Printf("[aux direct solve] gcg设置:\n");
    gcge_ops->Printf("·                  nevConv    = %d \n", nevConv);
    gcge_ops->Printf("·                  block_size = %d \n", block_size);
    gcge_ops->Printf("·                  nevMax     = %d \n", nevMax);
    gcge_ops->Printf("·                  nevGiven   = %d \n", nevGiven);
    gcge_ops->Printf("·                  gapMin     = %e \n", gapMin);
#endif

    EigenSolverSetParameters_GCG(check_conv_max_num,
                                 initX_orth_method, initX_orth_block_size, initX_orth_max_reorth, initX_orth_zero_tol,
                                 compP_orth_method, compP_orth_block_size, compP_orth_max_reorth, compP_orth_zero_tol,
                                 compW_orth_method, compW_orth_block_size, compW_orth_max_reorth, compW_orth_zero_tol,
                                 compW_bpcg_max_iter, compW_bpcg_rate, compW_bpcg_tol, compW_bpcg_tol_type, 0, // without shift
                                 compRR_min_num, compRR_min_gap, compRR_tol,
                                 pase_ops_to_gcge);
    pase_ops_to_gcge->EigenSolver(A, B, eigenvalues, eigenvectors, nevGiven, &nevConv, pase_ops_to_gcge);

#if PRINT_INFO
    gcge_ops->Printf("[aux direct solve] gcg结果:\n");
    gcge_ops->Printf("·                  nevConv    = %d \n", nevConv);
    int idx;
    for (idx = 0; idx < solver->pase_nev; ++idx)
    {
        gcge_ops->Printf("·                  eigs[%2d]   = %e \n", idx, eigenvalues[idx]);
    }
#endif

    memcpy(solver->eigenvalues + solver->conv_nev, solver->aux_eigenvalues + solver->conv_nev, (solver->pase_nev - solver->conv_nev) * sizeof(double));

    return 0;
}

int PASE_Mg_set_pase_aux_vector(PASE_MG_SOLVER solver)
{
    OPS *gcge_ops = solver->gcge_ops;
    int pase_nev = solver->pase_nev;
    int nlock_direct = solver->nlock_direct;
#if PRINT_INFO
    gcge_ops->Printf("-   nlock_direct : %d\n", nlock_direct);
    timer_start = gcge_ops->GetWtime();
#endif

    PASE_MultiVector aux_u = solver->aux_sol;
    void **u = aux_u->b_H;
    double *gamma = aux_u->aux_h;
    int mv_s[2] = {nlock_direct, nlock_direct};
    int mv_e[2] = {pase_nev, pase_nev};

    gcge_ops->MultiVecAxpby(0.0, u, 0.0, u, mv_s, mv_e, gcge_ops);
    int num = pase_nev * (pase_nev - nlock_direct);
    memset(gamma + nlock_direct * pase_nev, 0, pase_nev * (pase_nev - nlock_direct) * sizeof(double));
    int i;
    for (i = 0; i < pase_nev; i++)
    {
        gamma[(i + nlock_direct) * pase_nev + i] = 1.0;
    }
    if (nlock_direct > 0)
    {
        int current_level = solver->current_level;
        int aux_level = solver->aux_coarse_level;
        memset(gamma, 0.0, nlock_direct * pase_nev * sizeof(double));
        mv_s[0] = 0;
        mv_e[0] = nlock_direct;
        mv_s[1] = 0;
        mv_e[1] = nlock_direct;
        PASE_MULTIGRID_fromItoJ(solver->multigrid, current_level, aux_level,
                                mv_s, mv_e, solver->solution[current_level], u);
    }
#if PRINT_INFO
    timer_end = gcge_ops->GetWtime();
    gcge_ops->Printf("-   setup solution time : %f sec\n", timer_end - timer_start);
#endif
}

static int LUfactor_Solve(PASE_MG_SOLVER solver, PASE_Matrix aux_Mat)
{
    KSP ksp = (KSP)(aux_Mat->factorization);
    void **sol = aux_Mat->AH_inv_aux_Hh;
    void **rhs = aux_Mat->aux_Hh;
    Mat solmat, rhsmat;
    int pase_nev = solver->pase_nev;
    BVSetActiveColumns((BV)sol, 0, pase_nev);
    BVSetActiveColumns((BV)rhs, 0, pase_nev);
    BVGetMat((BV)sol, &solmat);
    BVGetMat((BV)rhs, &rhsmat);
    KSPMatSolve(ksp, (Mat)rhsmat, (Mat)solmat);
    BVRestoreMat((BV)sol, &solmat);
    BVRestoreMat((BV)rhs, &rhsmat);
    return 0;
}

int PASE_Mg_set_pase_aux_matrix(PASE_MG_SOLVER solver)
{
    OPS *gcge_ops = solver->gcge_ops;

    int pase_nev = solver->pase_nev;
    int current_level = solver->current_level;
    int idx_level, i, j;
    int nlock_auxmat, mv_s[2], mv_e[2];
    void **solution = solver->solution[current_level];

#if PRINT_INFO
    gcge_ops->Printf("-   nlock_auxmat_A : %d\n", solver->nlock_auxmat_A);
    timer_start = gcge_ops->GetWtime();
#endif
    /* a in aux_A */
    void *Ah = solver->multigrid->A_array[current_level];
    PASE_Matrix aux_A = solver->aux_A;
    void *AH = aux_A->A_H;
    void **a = aux_A->aux_Hh;
    void **rhs = solver->cg_rhs[current_level];
    nlock_auxmat = solver->nlock_auxmat_A;
    mv_s[0] = nlock_auxmat;
    mv_s[1] = nlock_auxmat;
    mv_e[0] = pase_nev;
    mv_e[1] = pase_nev;
    gcge_ops->MatDotMultiVec(Ah, solution, rhs, mv_s, mv_e, gcge_ops);
    for (idx_level = current_level; idx_level < solver->aux_coarse_level; idx_level++)
    {
        PASE_MULTIGRID_fromItoJ(solver->multigrid, idx_level, idx_level + 1, mv_s, mv_e, solver->cg_rhs[idx_level], solver->cg_rhs[idx_level + 1]);
    }
    /* alpha in aux_A */
    mv_s[0] = 0;
    mv_e[0] = pase_nev;
    mv_s[1] = nlock_auxmat;
    mv_e[1] = pase_nev;
    double *alpha = aux_A->aux_hh;
    gcge_ops->MultiVecQtAP('N', 'N', solution, Ah, solution, 0, mv_s, mv_e,
                           alpha + nlock_auxmat * pase_nev, pase_nev,
                           solver->cg_res[current_level], gcge_ops);
    for (i = nlock_auxmat; i < pase_nev; i++)
    {
        for (j = 0; j < nlock_auxmat; j++)
        {
            alpha[j * pase_nev + i] = alpha[i * pase_nev + j];
        }
    }
#if PRINT_INFO
    timer_end = gcge_ops->GetWtime();
    gcge_ops->Printf("-   setup A time : %f sec\n", timer_end - timer_start);
#endif

#if PRINT_INFO
    timer_start = gcge_ops->GetWtime();
    gcge_ops->Printf("-   nlock_auxmat_B : %d\n", solver->nlock_auxmat_B);
#endif
    /* b in aux_B */
    void *Bh = solver->multigrid->B_array[current_level];
    PASE_Matrix aux_B = solver->aux_B;
    void *BH = aux_B->A_H;
    void **b = aux_B->aux_Hh;
    rhs = solver->cg_p[current_level];
    nlock_auxmat = solver->nlock_auxmat_B;
    mv_s[0] = nlock_auxmat;
    mv_s[1] = nlock_auxmat;
    mv_e[0] = pase_nev;
    mv_e[1] = pase_nev;
    gcge_ops->MatDotMultiVec(Bh, solution, rhs, mv_s, mv_e, gcge_ops);
    for (idx_level = current_level; idx_level < solver->aux_coarse_level; idx_level++)
    {
        PASE_MULTIGRID_fromItoJ(solver->multigrid, idx_level, idx_level + 1, mv_s, mv_e, solver->cg_p[idx_level], solver->cg_p[idx_level + 1]);
    }
    /* beta in aux_B */
    mv_s[0] = 0;
    mv_e[0] = pase_nev;
    mv_s[1] = nlock_auxmat;
    mv_e[1] = pase_nev;
    double *beta = aux_B->aux_hh;
    gcge_ops->MultiVecQtAP('N', 'N', solution, Bh, solution, 0, mv_s, mv_e,
                           beta + nlock_auxmat * pase_nev, pase_nev,
                           solver->cg_w[current_level], gcge_ops);
    for (i = nlock_auxmat; i < pase_nev; i++)
    {
        for (j = 0; j < nlock_auxmat; j++)
        {
            beta[j * pase_nev + i] = beta[i * pase_nev + j];
        }
    }
#if PRINT_INFO
    timer_end = gcge_ops->GetWtime();
    gcge_ops->Printf("-   setup B time : %f sec\n", timer_end - timer_start);
#endif

    /* preconditioner in B */
    if (solver->if_precondition == 1 || solver->if_precondition == 3)
    {
        solver->precondition_time -= gcge_ops->GetWtime();
#if PRINT_INFO
        timer_start = gcge_ops->GetWtime();
        gcge_ops->Printf("-   B preconditioner [ON]\n");
#endif
        LUfactor_Solve(solver, aux_B);
        //处理B
        aux_B->is_diag = 1;
        nlock_auxmat = solver->nlock_auxmat_B;
        double *B_inpd = solver->double_tmp;
        void **BH_inv_b = aux_B->AH_inv_aux_Hh;
        mv_s[0] = 0;
        mv_e[0] = pase_nev;
        mv_s[1] = 0;
        mv_e[1] = pase_nev;
        gcge_ops->MultiVecInnerProd('N', b, BH_inv_b, 0, mv_s, mv_e, B_inpd, pase_nev, gcge_ops);
        int num = pase_nev * pase_nev;
        daxpy(&num, &m_one, B_inpd, &i_one, beta, &i_one);
        //处理A
        gcge_ops->MatDotMultiVec(AH, BH_inv_b, b, mv_s, mv_e, gcge_ops);
        mv_s[0] = 0;
        mv_e[0] = pase_nev;
        mv_s[1] = nlock_auxmat;
        mv_e[1] = pase_nev;
        num = pase_nev * (pase_nev - nlock_auxmat);
        gcge_ops->MultiVecInnerProd('N', BH_inv_b, b, 0, mv_s, mv_e, B_inpd, pase_nev, gcge_ops);
        daxpy(&num, &p_one, B_inpd, &i_one, alpha + nlock_auxmat * pase_nev, &i_one);
        gcge_ops->MultiVecInnerProd('N', a, BH_inv_b, 0, mv_s, mv_e, B_inpd, pase_nev, gcge_ops);
        daxpy(&num, &m_one, B_inpd, &i_one, alpha + nlock_auxmat * pase_nev, &i_one);
        gcge_ops->MultiVecInnerProd('N', BH_inv_b, a, 0, mv_s, mv_e, B_inpd, pase_nev, gcge_ops);
        daxpy(&num, &m_one, B_inpd, &i_one, alpha, &i_one);
        for (i = nlock_auxmat; i < pase_nev; i++)
        {
            for (j = 0; j < nlock_auxmat; j++)
            {
                alpha[j * pase_nev + i] = alpha[i * pase_nev + j];
            }
        }
        mv_s[0] = 0;
        mv_e[0] = pase_nev;
        mv_s[1] = 0;
        mv_e[1] = pase_nev;
        gcge_ops->MultiVecAxpby(-1.0, b, 1.0, a, mv_s, mv_e, gcge_ops);
        //处理u
        gcge_ops->MultiVecAxpby(1.0, BH_inv_b, 0.0, solver->aux_sol->b_H, mv_s, mv_e, gcge_ops);
        solver->precondition_time += gcge_ops->GetWtime();
#if PRINT_INFO
        timer_end = gcge_ops->GetWtime();
        gcge_ops->Printf("-   B precondition time : %f sec\n", timer_end - timer_start);
#endif
    }
#if PRINT_INFO
    else
    {
        gcge_ops->Printf("-   B preconditioner [OFF]\n");
    }
#endif

    /* preconditioner in A */
    if (solver->if_precondition == 2 || solver->if_precondition == 3)
    {
        solver->precondition_time -= gcge_ops->GetWtime();
#if PRINT_INFO
        timer_start = gcge_ops->GetWtime();
        gcge_ops->Printf("-   A preconditioner [ON]\n");
#endif
        LUfactor_Solve(solver, aux_A);
        //处理A
        nlock_auxmat = solver->nlock_auxmat_A;
        double *alpha_inv = aux_A->aux_hh_inv;
        void **AH_inv_a = aux_A->AH_inv_aux_Hh;
        mv_s[0] = 0;
        mv_e[0] = pase_nev;
        mv_s[1] = 0;
        mv_e[1] = pase_nev;
        gcge_ops->MultiVecInnerProd('N', a, AH_inv_a, 0, mv_s, mv_e, alpha_inv, pase_nev, gcge_ops);
        // DefaultMultiVecInnerProd('N', a, AH_inv_a, 0, mv_s, mv_e, alpha_inv, pase_nev, gcge_ops);
        int num = pase_nev * pase_nev;
        dscal(&num, &m_one, alpha_inv, &i_one);
        daxpy(&num, &p_one, alpha, &i_one, alpha_inv, &i_one);
        int *ipiv = solver->int_tmp;
        double *WORK = solver->double_tmp;
        dgetrf(&pase_nev, &pase_nev, alpha_inv, &pase_nev, ipiv, &info);
        assert(info == 0);
        dgetri(&pase_nev, alpha_inv, &pase_nev, ipiv, WORK, &num, &info);
        assert(info == 0);
        MPI_Bcast(alpha_inv, pase_nev * pase_nev, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        solver->precondition_time += gcge_ops->GetWtime();
        A_pre = 1;
#if PRINT_INFO
        timer_end = gcge_ops->GetWtime();
        gcge_ops->Printf("-   A precondition time : %f sec\n", timer_end - timer_start);
#endif
    }
#if PRINT_INFO
    else
    {
        gcge_ops->Printf("-   A preconditioner [OFF]\n");
    }
#endif
}

static void smoothing_bcg(PASE_MG_SOLVER solver)
{
    OPS *gcge_ops = solver->gcge_ops;
    int current_level = solver->current_level;
    void *A = solver->multigrid->A_array[current_level];
    void *B = solver->multigrid->B_array[current_level];
    void **solution = solver->solution[current_level];
    void **rhs = solver->cg_rhs[current_level];
    int pase_nev = solver->pase_nev;
    int nlock_smooth = solver->nlock_smooth;
    double *eigenvalues = solver->eigenvalues;

    int max_iter = solver->max_pre_count;
    int max_bcg_iter = 10 * max_iter;
    double rate = 1e-2, tol = 1e-12;
    double *dbl_ws = solver->double_tmp;
    int *int_ws = solver->int_tmp;
    void **cg_space_mv[3];
    cg_space_mv[0] = solver->cg_res[current_level];
    cg_space_mv[1] = solver->cg_p[current_level];
    cg_space_mv[2] = solver->cg_w[current_level];
    MultiLinearSolverSetup_BlockPCG(max_bcg_iter, rate, tol, "abs",
                                    cg_space_mv, dbl_ws, int_ws, NULL, NULL, gcge_ops);

    int mv_s[2] = {nlock_smooth, nlock_smooth};
    int mv_e[2] = {pase_nev, pase_nev};
    int niter;
    for (niter = 0; niter < max_iter; niter++)
    {
        assert(mv_e[0] - mv_s[0] == mv_e[1] - mv_s[1]);
        gcge_ops->MatDotMultiVec(B, solution, rhs, mv_s, mv_e, gcge_ops);
        gcge_ops->MultiVecLinearComb(NULL, rhs, 0, mv_s, mv_e, NULL, 0, eigenvalues + mv_s[0], 1, gcge_ops);
        gcge_ops->MultiLinearSolver(A, rhs, solution, mv_s, mv_e, gcge_ops);
    }
}

int PASE_Mg_smoothing(PASE_MG_SOLVER solver)
{
    OPS *gcge_ops = solver->gcge_ops;
#if PRINT_INFO
    gcge_ops->Printf("-   smooth_level : %d\n", solver->current_level);
    gcge_ops->Printf("-   max_pre_count : %d\n", solver->max_pre_count);
    gcge_ops->Printf("-   nlock_smooth : %d\n", solver->nlock_smooth);
#endif

    if (0 == strcmp("bcg", solver->smooth_type))
    {
#if PRINT_INFO
        gcge_ops->Printf("-   smooth_type : BCG\n");
#endif
        smoothing_bcg(solver);
    }

    return 0;
}

int PASE_Mg_prolong_from_Solution(PASE_MG_SOLVER solver, int object_level)
{
    int current_level = solver->current_level;
    if (current_level == object_level)
    {
        return 0;
    }

    int conv_nev = solver->conv_nev;
    int pase_nev = solver->pase_nev;
    int mv_s[2] = {conv_nev, conv_nev};
    int mv_e[2] = {pase_nev, pase_nev};

    PASE_MULTIGRID_fromItoJ(solver->multigrid, current_level, object_level, mv_s, mv_e, solver->solution[current_level], solver->solution[object_level]);
}

int PASE_Direct_solve(PASE_MG_SOLVER solver)
{
    OPS *gcge_ops = solver->gcge_ops;
    int initial_level = solver->initial_level;
    void *A = solver->multigrid->A_array[initial_level];
    void *B = solver->multigrid->B_array[initial_level];
    double *eigenvalues = solver->eigenvalues;
    void **eigenvectors = solver->solution[initial_level];

    /* gcg 设置 */
    int nevConv = solver->user_nev;
    int multiMax = 1;
    double gapMin = 1e-5;
    int nevGiven = solver->num_given_eigs;
    int block_size = nevConv / 2, nevMax = 2 * nevConv;
    int nevInit = 2 * nevConv;
    nevInit = nevInit < nevMax ? nevInit : nevMax;
    int max_iter_gcg = solver->max_initial_direct_count;
    double tol_gcg[2] = {solver->initial_rtol, solver->initial_rtol};
#if PRINT_INFO
    gcge_ops->Printf("[direct solve] level :\n", initial_level);
    gcge_ops->Printf("[direct solve] gcg设置:\n");
    gcge_ops->Printf("·              nevConv    = %d \n", nevConv);
    gcge_ops->Printf("·              block_size = %d \n", block_size);
    gcge_ops->Printf("·              nevMax     = %d \n", nevMax);
    gcge_ops->Printf("·              nevGiven   = %d \n", nevGiven);
    gcge_ops->Printf("·              gapMin     = %e \n", gapMin);
#endif
    if (eigenvalues == NULL)
    {
        eigenvalues = (double *)calloc(nevMax, sizeof(double));
    }
    if (eigenvectors == NULL)
    {
        gcge_ops->MultiVecCreateByMat(&eigenvectors, nevMax, A, gcge_ops);
        gcge_ops->MultiVecSetRandomValue(eigenvectors, 0, nevMax, gcge_ops);
    }
    void **gcg_mv_ws[4];
    double *dbl_ws;
    int *int_ws;
    gcge_ops->MultiVecCreateByMat(&gcg_mv_ws[0], nevMax + 2 * block_size, A, gcge_ops);
    gcge_ops->MultiVecSetRandomValue(gcg_mv_ws[0], 0, nevMax + 2 * block_size, gcge_ops);
    gcge_ops->MultiVecCreateByMat(&gcg_mv_ws[1], block_size, A, gcge_ops);
    gcge_ops->MultiVecSetRandomValue(gcg_mv_ws[1], 0, block_size, gcge_ops);
    gcge_ops->MultiVecCreateByMat(&gcg_mv_ws[2], block_size, A, gcge_ops);
    gcge_ops->MultiVecSetRandomValue(gcg_mv_ws[2], 0, block_size, gcge_ops);
    gcge_ops->MultiVecCreateByMat(&gcg_mv_ws[3], block_size, A, gcge_ops);
    gcge_ops->MultiVecSetRandomValue(gcg_mv_ws[3], 0, block_size, gcge_ops);

    dbl_ws = solver->double_tmp;
    int_ws = solver->int_tmp;
    srand(0);
    EigenSolverSetup_GCG(multiMax, gapMin, nevInit, nevMax, block_size,
                         tol_gcg, max_iter_gcg, 0, gcg_mv_ws, dbl_ws, int_ws, gcge_ops);
    int check_conv_max_num = 50;

    char initX_orth_method[8] = "mgs";
    int initX_orth_block_size = -1;
    int initX_orth_max_reorth = 2;
    double initX_orth_zero_tol = 2 * DBL_EPSILON; // 1e-12

    char compP_orth_method[8] = "mgs";
    int compP_orth_block_size = -1;
    int compP_orth_max_reorth = 2;
    double compP_orth_zero_tol = 2 * DBL_EPSILON; // 1e-12

    char compW_orth_method[8] = "mgs";
    int compW_orth_block_size = -1;
    int compW_orth_max_reorth = 2;
    double compW_orth_zero_tol = 2 * DBL_EPSILON; // 1e-12
    int compW_bpcg_max_iter = 30;
    double compW_bpcg_rate = 1e-2;
    double compW_bpcg_tol = 1e-14;
    char compW_bpcg_tol_type[8] = "abs";

    int compRR_min_num = -1;
    double compRR_min_gap = gapMin;
    double compRR_tol = 2 * DBL_EPSILON;
    EigenSolverSetParameters_GCG(check_conv_max_num,
                                 initX_orth_method, initX_orth_block_size, initX_orth_max_reorth, initX_orth_zero_tol,
                                 compP_orth_method, compP_orth_block_size, compP_orth_max_reorth, compP_orth_zero_tol,
                                 compW_orth_method, compW_orth_block_size, compW_orth_max_reorth, compW_orth_zero_tol,
                                 compW_bpcg_max_iter, compW_bpcg_rate, compW_bpcg_tol, compW_bpcg_tol_type, 1, // without shift
                                 compRR_min_num, compRR_min_gap, compRR_tol,
                                 gcge_ops);
    gcge_ops->EigenSolver(A, B, eigenvalues, eigenvectors, nevGiven, &nevConv, gcge_ops);

#if PRINT_INFO
    gcge_ops->Printf("[direct solve] gcg结果:\n");
    gcge_ops->Printf("·              nevConv    = %d \n", nevConv);
    int idx;
    for (idx = 0; idx < solver->pase_nev; ++idx)
    {
        gcge_ops->Printf("·              eigs[%2d]   = %e \n", idx, eigenvalues[idx]);
    }
#endif

    gcge_ops->MultiVecDestroy(&gcg_mv_ws[0], nevMax + 2 * block_size, gcge_ops);
    gcge_ops->MultiVecDestroy(&gcg_mv_ws[1], block_size, gcge_ops);
    gcge_ops->MultiVecDestroy(&gcg_mv_ws[2], block_size, gcge_ops);
    gcge_ops->MultiVecDestroy(&gcg_mv_ws[3], block_size, gcge_ops);

    memcpy(solver->aux_eigenvalues, solver->eigenvalues, solver->pase_nev * sizeof(double));

    return 0;
}

int PASE_Mg_print_param(PASE_MG_SOLVER solver)
{
    int i = 0;
    int num_levels = solver->num_levels;
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "|=================== PARAMETERS ====================|\n");
    PetscPrintf(PETSC_COMM_WORLD, "|------------------ [eigenvalue] -------------------|\n");
    PetscPrintf(PETSC_COMM_WORLD, "| user_nev                 = %12d           |\n", solver->user_nev);
    PetscPrintf(PETSC_COMM_WORLD, "| pase_nev                 = %12d           |\n", solver->pase_nev);
    PetscPrintf(PETSC_COMM_WORLD, "| num_given_eigs           = %12d           |\n", solver->num_given_eigs);
    PetscPrintf(PETSC_COMM_WORLD, "|------------------ [converging] -------------------|\n");
    PetscPrintf(PETSC_COMM_WORLD, "| rtol                     = %12e           |\n", solver->rtol);
    PetscPrintf(PETSC_COMM_WORLD, "| atol                     = %12e           |\n", solver->atol);
    PetscPrintf(PETSC_COMM_WORLD, "| aux_rtol                 = %12e           |\n", solver->aux_rtol);
    PetscPrintf(PETSC_COMM_WORLD, "| initial_rtol             = %12e           |\n", solver->initial_rtol);
    PetscPrintf(PETSC_COMM_WORLD, "|------------------ [multigrid] --------------------|\n");
    PetscPrintf(PETSC_COMM_WORLD, "| num_levels               = %12d           |\n", solver->num_levels);
    PetscPrintf(PETSC_COMM_WORLD, "| initial_level            = %12d           |\n", solver->initial_level);
    PetscPrintf(PETSC_COMM_WORLD, "| aux_coarse_level         = %12d           |\n", solver->aux_coarse_level);
    PetscPrintf(PETSC_COMM_WORLD, "| aux_fine_level           = %12d           |\n", solver->aux_fine_level);
    PetscPrintf(PETSC_COMM_WORLD, "|------------------ [iterations] -------------------|\n");
    PetscPrintf(PETSC_COMM_WORLD, "| max_cycle_count          = %12d           |\n", solver->max_cycle_count);
    PetscPrintf(PETSC_COMM_WORLD, "| max_pre_count            = %12d           |\n", solver->max_pre_count);
    PetscPrintf(PETSC_COMM_WORLD, "| max_post_count           = %12d           |\n", solver->max_post_count);
    PetscPrintf(PETSC_COMM_WORLD, "| max_direct_count         = %12d           |\n", solver->max_direct_count);
    PetscPrintf(PETSC_COMM_WORLD, "| max_initial_direct_count = %12d           |\n", solver->max_initial_direct_count);
    PetscPrintf(PETSC_COMM_WORLD, "|===================================================|\n\n");
    return 0;
}

int PASE_Mg_set_up(PASE_MG_SOLVER solver, void *A, void *B)
{
    OPS *gcge_ops = solver->gcge_ops;
    solver->solver_setup_time -= gcge_ops->GetWtime();
#if PRINT_INFO
    double multigrid_setup_time = 0.0;
    double aux_space_setup_time = 0.0;
#endif
    /* aux space ops */
    OPS *pase_ops_to_gcge;
    OPS_Create(&pase_ops_to_gcge);
    GCGE_PASE_SetOps(pase_ops_to_gcge, solver->pase_ops);
    OPS_Setup(pase_ops_to_gcge);
    solver->pase_ops_to_gcge = pase_ops_to_gcge;

    /* 特征值所需空间 */
    int max_nev = solver->max_nev;
    solver->eigenvalues = (double *)calloc(max_nev, sizeof(double));
    solver->aux_eigenvalues = (double *)calloc(max_nev, sizeof(double));
    solver->abs_res_norm = (double *)calloc(max_nev, sizeof(double));

#if PRINT_INFO
    if (solver->print_level > 1)
    {
        multigrid_setup_time -= gcge_ops->GetWtime();
    }
#endif
    /* multigrid setup */
    PASE_MULTIGRID_Create(&(solver->multigrid), solver->num_levels, A, B, solver->pase_nev, gcge_ops);
    solver->num_levels = solver->multigrid->num_levels;
    solver->mg_coarsest_level = solver->num_levels - 1;
    solver->mg_finest_level = 0;
    solver->aux_coarse_level = solver->multigrid->pase_coarse_level;
    solver->aux_fine_level = solver->multigrid->pase_fine_level;
    solver->initial_level = (solver->num_given_eigs == 0) ? solver->aux_coarse_level : solver->initial_level;
    solver->current_level = solver->initial_level;
    solver->current_cycle = 0;
#if PRINT_INFO
    if (solver->print_level > 1)
    {
        multigrid_setup_time += gcge_ops->GetWtime();
    }
#endif

    /* 向量矩阵所需空间 */
    solver->solution = solver->multigrid->solution;
    solver->cg_rhs = solver->multigrid->cg_rhs;
    solver->cg_res = solver->multigrid->cg_res;
    solver->cg_p = solver->multigrid->cg_p;
    solver->cg_w = solver->multigrid->cg_w;
    solver->double_tmp = solver->multigrid->double_tmp;
    solver->int_tmp = solver->multigrid->int_tmp;

#if PRINT_INFO
    if (solver->print_level > 1)
    {
        aux_space_setup_time -= gcge_ops->GetWtime();
    }
#endif
    /* aux space */
    PASE_Mg_pase_aux_matrix_create(solver);
    PASE_Mg_pase_aux_vector_create(solver);
    if (solver->if_precondition)
    {
        solver->precondition_time -= gcge_ops->GetWtime();
        PASE_Mg_pase_aux_matrix_factorization(solver);
        solver->precondition_time += gcge_ops->GetWtime();
    }
    /* 创建aux_gcg需要的空间 */
    int pase_nev = solver->pase_nev;
    int user_nev = solver->user_nev;
    int aux_coarse_level = solver->aux_coarse_level;
    solver->aux_gcg_mv[0] = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    solver->aux_gcg_mv[0]->num_aux_vec = pase_nev;
    solver->aux_gcg_mv[0]->num_vec = pase_nev + user_nev;
    gcge_ops->MultiVecCreateByMat(&(solver->aux_gcg_mv[0]->b_H), pase_nev + user_nev, solver->multigrid->A_array[aux_coarse_level], gcge_ops);
    solver->aux_gcg_mv[0]->aux_h = (double *)malloc(pase_nev * (pase_nev + user_nev) * sizeof(double));
    solver->aux_gcg_mv[0]->aux_h_tmp = (double *)malloc(pase_nev * (pase_nev + user_nev) * sizeof(double));

    solver->aux_gcg_mv[1] = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    solver->aux_gcg_mv[1]->num_aux_vec = pase_nev;
    solver->aux_gcg_mv[1]->num_vec = user_nev;
    gcge_ops->MultiVecCreateByMat(&(solver->aux_gcg_mv[1]->b_H), user_nev, solver->multigrid->A_array[aux_coarse_level], gcge_ops);
    solver->aux_gcg_mv[1]->aux_h = (double *)malloc(pase_nev * user_nev * sizeof(double));
    solver->aux_gcg_mv[1]->aux_h_tmp = (double *)malloc(pase_nev * user_nev * sizeof(double));

    solver->aux_gcg_mv[2] = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    solver->aux_gcg_mv[2]->num_aux_vec = pase_nev;
    solver->aux_gcg_mv[2]->num_vec = user_nev;
    gcge_ops->MultiVecCreateByMat(&(solver->aux_gcg_mv[2]->b_H), user_nev, solver->multigrid->A_array[aux_coarse_level], gcge_ops);
    solver->aux_gcg_mv[2]->aux_h = (double *)malloc(pase_nev * user_nev * sizeof(double));
    solver->aux_gcg_mv[2]->aux_h_tmp = (double *)malloc(pase_nev * user_nev * sizeof(double));

    solver->aux_gcg_mv[3] = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    solver->aux_gcg_mv[3]->num_aux_vec = pase_nev;
    solver->aux_gcg_mv[3]->num_vec = user_nev;
    gcge_ops->MultiVecCreateByMat(&(solver->aux_gcg_mv[3]->b_H), user_nev, solver->multigrid->A_array[aux_coarse_level], gcge_ops);
    solver->aux_gcg_mv[3]->aux_h = (double *)malloc(pase_nev * user_nev * sizeof(double));
    solver->aux_gcg_mv[3]->aux_h_tmp = (double *)malloc(pase_nev * user_nev * sizeof(double));

    solver->solver_setup_time += gcge_ops->GetWtime();
#if PRINT_INFO
    if (solver->print_level > 1)
    {
        aux_space_setup_time += gcge_ops->GetWtime();
        gcge_ops->Printf("------ timer -------\n");
        gcge_ops->Printf("[solver setup] 总时间: %f\n", solver->solver_setup_time);
        gcge_ops->Printf("(1) multigrid setup : %f (%g%%)\n", multigrid_setup_time,
                         100 * multigrid_setup_time / solver->solver_setup_time);
        gcge_ops->Printf("(2) aux space setup : %f (%g%%)\n", aux_space_setup_time,
                         100 * aux_space_setup_time / solver->solver_setup_time);
        gcge_ops->Printf("------ timer -------\n\n");
    }
#endif

    return 0;
}

int PASE_Mg_pase_aux_vector_create(PASE_MG_SOLVER solver)
{
    solver->aux_sol = (PASE_MultiVector)malloc(sizeof(pase_MultiVector));
    PASE_MultiVector aux_sol = solver->aux_sol;
    aux_sol->num_aux_vec = solver->pase_nev;
    aux_sol->num_vec = solver->pase_nev;
    aux_sol->b_H = solver->solution[solver->aux_coarse_level];
    aux_sol->aux_h = (double *)malloc(solver->pase_nev * solver->pase_nev * sizeof(double));
    aux_sol->aux_h_tmp = (double *)malloc(solver->pase_nev * solver->pase_nev * sizeof(double));
}

int PASE_Mg_pase_aux_matrix_factorization(PASE_MG_SOLVER solver)
{
    if (solver->if_precondition == 2 || solver->if_precondition == 3)
    {
        PASE_Matrix aux_A = solver->aux_A;
        KSP A_ksp;
        PC A_factor;
        KSPCreate(PETSC_COMM_WORLD, &A_ksp);
        KSPSetOperators(A_ksp, aux_A->A_H, aux_A->A_H);
        KSPSetType(A_ksp, KSPPREONLY);
        KSPGetPC(A_ksp, &A_factor);
        PCSetType(A_factor, PCLU);
        PCFactorSetMatSolverType(A_factor, MATSOLVERSUPERLU_DIST);
        PCFactorSetUpMatSolverType(A_factor);
        KSPSetFromOptions(A_ksp);
        KSPSetUp(A_ksp);
        aux_A->factorization = (void *)A_ksp;
    }

    if (solver->if_precondition == 1 || solver->if_precondition == 3)
    {
        PASE_Matrix aux_B = solver->aux_B;
        KSP B_ksp;
        PC B_factor;
        KSPCreate(PETSC_COMM_WORLD, &B_ksp);
        KSPSetOperators(B_ksp, aux_B->A_H, aux_B->A_H);
        KSPSetType(B_ksp, KSPPREONLY);
        KSPGetPC(B_ksp, &B_factor);
        PCSetType(B_factor, PCLU);
        PCFactorSetMatSolverType(B_factor, MATSOLVERSUPERLU_DIST);
        PCFactorSetUpMatSolverType(B_factor);
        KSPSetFromOptions(B_ksp);
        KSPSetUp(B_ksp);
        aux_B->factorization = (void *)B_ksp;
    }
}

int PASE_Mg_pase_aux_matrix_create(PASE_MG_SOLVER solver)
{
    int aux_level = solver->aux_coarse_level;
    int pase_nev = solver->pase_nev;
    int if_prc = solver->if_precondition;

    solver->aux_A = (PASE_Matrix)malloc(sizeof(pase_Matrix));
    PASE_Matrix aux_A = solver->aux_A;
    aux_A->num_aux_vec = pase_nev;
    aux_A->A_H = solver->multigrid->A_array[aux_level];
    aux_A->aux_Hh = solver->cg_rhs[aux_level];
    aux_A->aux_hh = (double *)malloc(pase_nev * pase_nev * sizeof(double));
    if (if_prc == 2 || if_prc == 3)
    {
        aux_A->factorization = NULL;
        aux_A->AH_inv_aux_Hh = solver->cg_res[aux_level];
        aux_A->aux_hh_inv = (double *)malloc(pase_nev * pase_nev * sizeof(double));
    }
    else
    {
        aux_A->factorization = NULL;
        aux_A->AH_inv_aux_Hh = NULL;
        aux_A->aux_hh_inv = NULL;
    }
    aux_A->if_sym = 1;
    aux_A->is_diag = 0;

    solver->aux_B = (PASE_Matrix)malloc(sizeof(pase_Matrix));
    PASE_Matrix aux_B = solver->aux_B;
    aux_B->num_aux_vec = pase_nev;
    aux_B->A_H = solver->multigrid->B_array[aux_level];
    aux_B->aux_Hh = solver->cg_p[aux_level];
    aux_B->aux_hh = (double *)malloc(pase_nev * pase_nev * sizeof(double));
    if (if_prc == 1 || if_prc == 3)
    {
        aux_B->factorization = NULL;
        aux_B->AH_inv_aux_Hh = solver->cg_w[aux_level];
        aux_B->aux_hh_inv = (double *)malloc(pase_nev * pase_nev * sizeof(double));
    }
    else
    {
        aux_B->factorization = NULL;
        aux_B->AH_inv_aux_Hh = NULL;
        aux_B->aux_hh_inv = NULL;
    }
    aux_B->if_sym = 1;
    aux_B->is_diag = 0;
}

PASE_MG_SOLVER PASE_Mg_solver_create(PASE_PARAMETER param, PASE_OPS *pase_ops)
{
    PASE_MG_SOLVER solver = (PASE_MG_SOLVER)malloc(sizeof(PASE_MG_SOLVER_PRIVATE));

    solver->multigrid_type = param->multigrid_type;
#if PRINT_INFO
    if (param->print_level)
    {
        if (solver->multigrid_type == 0)
        {
            pase_ops->gcge_ops->Printf("\n[solver_create] Multigrid Type is == AMG ==\n");
        }
        if (solver->multigrid_type == 1)
        {
            pase_ops->gcge_ops->Printf("\n[solver_create] Multigrid Type is == GMG ==\n");
        }
    }
#endif
    solver->if_precondition = param->if_precondition;
#if PRINT_INFO
    if (param->print_level)
    {
        if (solver->if_precondition)
        {
            pase_ops->gcge_ops->Printf("[solver_create] Preconditioner(%d) On!\n", solver->if_precondition);
        }
        else
        {
            pase_ops->gcge_ops->Printf("[solver_create] Preconditioner Off!\n");
        }
    }
#endif
    solver->print_level = param->print_level;
    strcpy(solver->smooth_type, "bcg");

    solver->multigrid = NULL;
    solver->num_levels = param->num_levels;
    solver->initial_level = param->initial_level;
    solver->mg_coarsest_level = param->mg_coarsest_level;
    solver->mg_finest_level = param->mg_finest_level;
    solver->aux_coarse_level = param->aux_coarse_level;
    solver->aux_fine_level = param->aux_fine_level;
    solver->current_level = solver->initial_level;
    solver->current_cycle = 0;

    solver->max_cycle_count = param->max_cycle_count;
    solver->max_pre_count = param->max_pre_count;
    solver->max_post_count = param->max_post_count;
    solver->max_direct_count = param->max_direct_count;
    solver->max_initial_direct_count = param->max_initial_direct_count;

    solver->user_nev = param->nev;
    solver->more_nev = param->more_nev;
    solver->max_nev =
        ((2 * param->nev) < (param->nev + param->more_nev)) ? (2 * param->nev) : (param->nev + param->more_nev);
#if PRINT_INFO
    if (solver->print_level)
    {
        pase_ops->gcge_ops->Printf("[CREATING] user_nev = %d | pase_nev = %d\n", solver->user_nev, solver->max_nev);
    }
#endif
#if PRINT_INFO
    if (solver->print_level)
    {
        pase_ops->gcge_ops->Printf("[CREATING] init given eigs num: %d\n\n", param->num_given_eigs);
    }
#endif
    solver->pase_nev = solver->max_nev;
    solver->num_given_eigs = param->num_given_eigs;
    if (solver->num_given_eigs > solver->user_nev)
    {
        solver->pase_nev = solver->num_given_eigs;
#if PRINT_INFO
        if (solver->print_level)
        {
            pase_ops->gcge_ops->Printf("[CREATING] 修正 pase_nev = %d\n\n", solver->pase_nev);
        }
#endif
    }
    solver->conv_nev = 0;
    solver->nlock_smooth = 0;
    solver->nlock_direct = 0;
    solver->nlock_auxmat_A = 0;
    solver->nlock_auxmat_B = 0;

    solver->rtol = param->rtol;
    solver->atol = param->atol;
    solver->initial_rtol = param->initial_rtol;
    solver->aux_rtol = param->aux_rtol;
    solver->abs_res_norm = NULL;

    solver->eigenvalues = NULL;
    solver->solution = NULL;
    solver->cg_rhs = NULL;
    solver->cg_res = NULL;
    solver->cg_p = NULL;
    solver->cg_w = NULL;
    solver->double_tmp = NULL;
    solver->int_tmp = NULL;
    solver->aux_gcg_mv[0] = NULL;
    solver->aux_gcg_mv[1] = NULL;
    solver->aux_gcg_mv[2] = NULL;
    solver->aux_gcg_mv[3] = NULL;

    solver->gcge_ops = pase_ops->gcge_ops;
    solver->pase_ops = pase_ops;
    solver->pase_ops_to_gcge = NULL;

    solver->aux_sol = NULL;
    solver->aux_A = NULL;
    solver->aux_B = NULL;
    solver->aux_eigenvalues = NULL;

    solver->solver_setup_time = 0.0;
    solver->get_initvec_time = 0.0;
    solver->smooth_time = 0.0;
    solver->build_aux_time = 0.0;
    solver->prolong_time = 0.0;
    solver->error_estimate_time = 0.0;
    solver->aux_direct_solve_time = 0.0;
    solver->precondition_time = 0.0;
    solver->total_solve_time = 0.0;
    solver->total_time = 0.0;

    return solver;
}

int PASE_Matrix_destroy_sub(PASE_Matrix *aux_A)
{
    free((*aux_A)->aux_hh);
    (*aux_A)->aux_hh = NULL;
    if ((*aux_A)->aux_hh_inv != NULL)
    {
        free((*aux_A)->aux_hh_inv);
        (*aux_A)->aux_hh_inv = NULL;
    }
    free(*aux_A);
    *aux_A = NULL;
    return 0;
}

int PASE_MultiVector_destroy_sub(PASE_MultiVector *aux_sol)
{
    free((*aux_sol)->aux_h);
    (*aux_sol)->aux_h = NULL;
    free((*aux_sol)->aux_h_tmp);
    (*aux_sol)->aux_h_tmp = NULL;
    free(*aux_sol);
    *aux_sol = NULL;
    return 0;
}

int PASE_Mg_solver_destroy(PASE_MG_SOLVER solver)
{
    if (solver->eigenvalues != NULL)
    {
        free(solver->eigenvalues);
        solver->eigenvalues = NULL;
    }
    if (solver->abs_res_norm != NULL)
    {
        free(solver->abs_res_norm);
        solver->abs_res_norm = NULL;
    }
    if (solver->multigrid != NULL)
    {
        PASE_MULTIGRID_Destroy(&(solver->multigrid));
    }

    //-----------------------------------------------------
    //释放 aux_sol空间
    if (solver->aux_eigenvalues != NULL)
    {
        free(solver->aux_eigenvalues);
        solver->aux_eigenvalues = NULL;
    }

    PASE_Matrix_destroy_sub(&(solver->aux_A));
    PASE_Matrix_destroy_sub(&(solver->aux_B));
    PASE_MultiVector_destroy_sub(&(solver->aux_sol));
    BVDestroy((BV *)(&(solver->aux_gcg_mv[0]->b_H)));
    PASE_MultiVector_destroy_sub(&(solver->aux_gcg_mv[0]));
    PASE_MultiVector_destroy_sub(&(solver->aux_gcg_mv[1]));
    PASE_MultiVector_destroy_sub(&(solver->aux_gcg_mv[2]));
    PASE_MultiVector_destroy_sub(&(solver->aux_gcg_mv[3]));

    if (solver->pase_ops != NULL)
        PASE_OPS_Destroy(&(solver->pase_ops));
    if (solver->gcge_ops != NULL)
        OPS_Destroy(&(solver->gcge_ops));
    if (solver->pase_ops_to_gcge != NULL)
        OPS_Destroy(&(solver->pase_ops_to_gcge));

    free(solver);
    return 0;
}

int PASE_Mg_print_timer(PASE_MG_SOLVER solver)
{
    OPS *gcge_ops = solver->gcge_ops;
    int idx_eigen = 0;
    if (solver->print_level > 1)
    {
        gcge_ops->Printf("\n");
        gcge_ops->Printf("=============================================================\n");
        for (idx_eigen = 0; idx_eigen < solver->pase_nev; idx_eigen++)
        {
            gcge_ops->Printf("%3d-th eig=%.8e, abs_res = %.8e\n", idx_eigen, solver->eigenvalues[idx_eigen], solver->abs_res_norm[idx_eigen]);
        }
    }
    if (solver->print_level > 0)
    {
        gcge_ops->Printf("=============================================================\n");
        gcge_ops->Printf("· solver setup time         = %f seconds\n", solver->solver_setup_time);
        gcge_ops->Printf("· get initvec time          = %f seconds\n", solver->get_initvec_time);
        gcge_ops->Printf("· smooth time               = %f seconds\n", solver->smooth_time);
        gcge_ops->Printf("· build aux time            = %f seconds\n", solver->build_aux_time);
        gcge_ops->Printf("· prolong time              = %f seconds\n", solver->prolong_time);
        gcge_ops->Printf("· error estimate time       = %f seconds\n", solver->error_estimate_time);
        gcge_ops->Printf("· aux direct solve time     = %f seconds\n", solver->aux_direct_solve_time);
        gcge_ops->Printf("· precondition time         = %f seconds\n", solver->precondition_time);
        gcge_ops->Printf("· total solve time          = %f seconds\n", solver->total_solve_time);
        gcge_ops->Printf("· total time                = %f seconds\n", solver->total_time);
        gcge_ops->Printf("=============================================================\n");
    }
}