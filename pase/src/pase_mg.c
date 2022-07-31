#include "pase_mg.h"
#include "pase_param.h"

#define PRINT_INFO 1
#define PRINT_TIME_INFO 1
#define PETSC_GAMG_USE 1

static void get_amg_array_by_petsc_gamg(PASE_MULTIGRID *multi_grid, void *A, void *B, OPS *gcge_ops)
{
#if PRINT_TIME_INFO
    double get_amg_array_time = 0.0;
    double get_A_P_array_time = 0.0;
    double get_B_array_time = 0.0;
    get_amg_array_time -= gcge_ops->GetWtime();
#endif
#if PRINT_INFO
    gcge_ops->Printf("\n[petsc_gamg] multigrid 初始层数: %d levels\n", (*multi_grid)->num_levels);
#endif
#if PRINT_TIME_INFO
    get_A_P_array_time -= gcge_ops->GetWtime();
#endif
    PetscErrorCode ierr;
    Mat petsc_A = (Mat)A, petsc_B = (Mat)B;
    Mat *A_array = NULL, *B_array = NULL, *P_array = NULL;
    Mat *Aarr = NULL, *Parr = NULL;

    PC pc;
    PCCreate(PETSC_COMM_WORLD, &pc);
    PCSetOperators(pc, petsc_A, petsc_A);
    PCSetType(pc, PCGAMG);
    PCGAMGSetNlevels(pc, (*multi_grid)->num_levels);
    PCGAMGSetType(pc, PCGAMGAGG);
    PCGAMGSetRepartition(pc, PETSC_TRUE);
    PCSetUp(pc);
    PCGetCoarseOperators(pc, &((*multi_grid)->num_levels), &Aarr);
    PCGetInterpolations(pc, &((*multi_grid)->num_levels), &Parr);
    // size of Parr = num_levels-1
    if ((*multi_grid)->coarsest_level > (*multi_grid)->num_levels - 1)
    {
        (*multi_grid)->coarsest_level = (*multi_grid)->num_levels - 1;
    }
#if PRINT_INFO
    gcge_ops->Printf("[petsc_gamg] multigrid 层数更新: %d levels\n", (*multi_grid)->num_levels);
#endif

    A_array = (Mat *)malloc(((*multi_grid)->num_levels) * sizeof(Mat));
    P_array = (Mat *)malloc(((*multi_grid)->num_levels - 1) * sizeof(Mat));
    (*multi_grid)->array_size_g = (int *)malloc(((*multi_grid)->num_levels) * sizeof(int));
    (*multi_grid)->array_size_l = (int *)malloc(((*multi_grid)->num_levels) * sizeof(int));

    A_array[0] = petsc_A;
    int idx_level, m, n;
    MatGetSize(A_array[0], &m, &n);
    (*multi_grid)->array_size_g[0] = m;
#if PRINT_INFO
    gcge_ops->Printf("·            size of level %d: %d x %d\n", 0, m, n);
#endif
    MatGetLocalSize(A_array[0], &m, &n);
    (*multi_grid)->array_size_l[0] = m;
    for (idx_level = 1; idx_level < ((*multi_grid)->num_levels); idx_level++)
    {
        P_array[idx_level - 1] = Parr[((*multi_grid)->num_levels) - idx_level - 1];
#if PRINT_INFO
        MatGetSize(P_array[idx_level - 1], &m, &n);
        gcge_ops->Printf("·             |  interpolation %d: %d x %d\n", idx_level - 1, m, n);
#endif
        A_array[idx_level] = Aarr[((*multi_grid)->num_levels) - idx_level - 1];
        MatGetSize(A_array[idx_level], &m, &n);
        (*multi_grid)->array_size_g[idx_level] = m;
#if PRINT_INFO
        gcge_ops->Printf("·            size of level %d: %d x %d\n", idx_level, m, n);
#endif
        MatGetLocalSize(A_array[idx_level], &m, &n);
        (*multi_grid)->array_size_l[idx_level] = m;
    }
    PetscFree(Aarr);
    PetscFree(Parr);
    PCDestroy(&pc);
#if PRINT_TIME_INFO
    get_A_P_array_time += gcge_ops->GetWtime();
    get_B_array_time -= gcge_ops->GetWtime();
#endif
    B_array = malloc((*multi_grid)->num_levels * sizeof(Mat));
    B_array[0] = petsc_B;
    for (idx_level = 1; idx_level < ((*multi_grid)->num_levels); idx_level++)
    {
        MatPtAP(B_array[idx_level - 1], P_array[idx_level - 1],
                MAT_INITIAL_MATRIX, PETSC_DEFAULT, &(B_array[idx_level]));
    }

    (*multi_grid)->A_array = (void **)A_array;
    (*multi_grid)->B_array = (void **)B_array;
    (*multi_grid)->P_array = (void **)P_array;

    int num_levels = (*multi_grid)->num_levels;
    BV *solution = (BV *)malloc(num_levels * sizeof(BV));
    BV *cg_rhs = (BV *)malloc(num_levels * sizeof(BV));
    BV *cg_res = (BV *)malloc(num_levels * sizeof(BV));
    BV *cg_p = (BV *)malloc(num_levels * sizeof(BV));
    BV *cg_w = (BV *)malloc(num_levels * sizeof(BV));
    for (idx_level = 0; idx_level < num_levels; idx_level++)
    {
        gcge_ops->MultiVecCreateByMat((void ***)(&(solution[idx_level])), (*multi_grid)->pase_nev,
                                      (void **)A_array[idx_level], gcge_ops);
        gcge_ops->MultiVecCreateByMat((void ***)(&(cg_rhs[idx_level])), (*multi_grid)->step_size,
                                      (void **)A_array[idx_level], gcge_ops);
        gcge_ops->MultiVecCreateByMat((void ***)(&(cg_res[idx_level])), (*multi_grid)->step_size,
                                      (void **)A_array[idx_level], gcge_ops);
        gcge_ops->MultiVecCreateByMat((void ***)(&(cg_p[idx_level])), (*multi_grid)->step_size,
                                      (void **)A_array[idx_level], gcge_ops);
        gcge_ops->MultiVecCreateByMat((void ***)(&(cg_w[idx_level])), (*multi_grid)->step_size,
                                      (void **)A_array[idx_level], gcge_ops);
    }
    (*multi_grid)->solution = (void ***)solution;
    (*multi_grid)->cg_rhs = (void ***)cg_rhs;
    (*multi_grid)->cg_res = (void ***)cg_res;
    (*multi_grid)->cg_p = (void ***)cg_p;
    (*multi_grid)->cg_w = (void ***)cg_w;

    int size = 2 * (*multi_grid)->pase_nev;
    size = 3 * size * size + 11 * size + size;
    double *double_tmp = (double *)malloc(size * sizeof(double));
    (*multi_grid)->double_tmp = double_tmp;
    int *int_tmp = (int *)malloc(size * sizeof(int));
    (*multi_grid)->int_tmp = int_tmp;

#if PRINT_TIME_INFO
    get_B_array_time += gcge_ops->GetWtime();
    get_amg_array_time += gcge_ops->GetWtime();
    gcge_ops->Printf("\n------ timer -------\n");
    gcge_ops->Printf("[get_amg_array] 总时间: %f\n", get_amg_array_time);
    gcge_ops->Printf("(1): 调用GAMG获取A_array和P_array时间:%f\n", get_A_P_array_time);
    gcge_ops->Printf("(2): 计算B_array时间:                 %f\n", get_B_array_time);
    gcge_ops->Printf("------ timer -------\n\n");
#endif
    return;
}

int PASE_MULTIGRID_Create(PASE_MULTIGRID *multi_grid, int num_levels, void *A, void *B, int pase_nev, OPS *gcge_ops)
{

    (*multi_grid) = (PASE_MULTIGRID)malloc(sizeof(pase_MultiGrid));
    (*multi_grid)->num_levels = num_levels;
    (*multi_grid)->coarsest_level = num_levels - 1;
    (*multi_grid)->pase_coarse_level = num_levels - 1;
    (*multi_grid)->pase_fine_level = 0;
    (*multi_grid)->gcge_ops = gcge_ops;
    (*multi_grid)->array_size_g = NULL;
    (*multi_grid)->array_size_l = NULL;
    (*multi_grid)->A_array = NULL;
    (*multi_grid)->B_array = NULL;
    (*multi_grid)->P_array = NULL;
    (*multi_grid)->pase_nev = pase_nev;
    (*multi_grid)->step_size = pase_nev;
    (*multi_grid)->solution = NULL;
    (*multi_grid)->cg_rhs = NULL;
    (*multi_grid)->cg_res = NULL;
    (*multi_grid)->cg_p = NULL;
    (*multi_grid)->cg_w = NULL;
    (*multi_grid)->double_tmp = NULL;
    (*multi_grid)->int_tmp = NULL;

#if PETSC_GAMG_USE
    get_amg_array_by_petsc_gamg(multi_grid, A, B, gcge_ops);
#endif

    PASE_MULTIGRID_TwoGirdLevel(*multi_grid);

    return 0;
}

int PASE_MULTIGRID_TwoGirdLevel(PASE_MULTIGRID multi_grid)
{
    int idx_level, coarse_level = 0;
    int *array_size_l = multi_grid->array_size_l;
    for (idx_level = multi_grid->coarsest_level; idx_level > 0; idx_level--)
    {
        if (array_size_l[idx_level] > 32)
        {
            coarse_level = idx_level;
            break;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &coarse_level, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    multi_grid->pase_coarse_level = coarse_level;
#if PRINT_INFO
    multi_grid->gcge_ops->Printf("[TwoGird] 选择粗层为: Level %d\n", coarse_level);
    multi_grid->gcge_ops->Printf("[TwoGird] 选择细层为: Level %d\n\n", multi_grid->pase_fine_level);
#endif
}

int PASE_MULTIGRID_Destroy(PASE_MULTIGRID *multi_grid)
{
#if PETSC_GAMG_USE
    PetscErrorCode ierr;
    int level;
    for (level = 0; level < (*multi_grid)->num_levels; ++level)
    {
        ierr = MatDestroy((Mat *)(&((*multi_grid)->A_array[level])));
        ierr = MatDestroy((Mat *)(&((*multi_grid)->B_array[level])));
        ierr = BVDestroy((BV *)(&((*multi_grid)->solution[level])));
        ierr = BVDestroy((BV *)(&((*multi_grid)->cg_rhs[level])));
        ierr = BVDestroy((BV *)(&((*multi_grid)->cg_res[level])));
        ierr = BVDestroy((BV *)(&((*multi_grid)->cg_p[level])));
        ierr = BVDestroy((BV *)(&((*multi_grid)->cg_w[level])));
    }
    for (level = 0; level < (*multi_grid)->num_levels - 1; ++level)
    {
        ierr = MatDestroy((Mat *)(&((*multi_grid)->P_array[level])));
    }
#endif
    free((*multi_grid)->A_array);
    free((*multi_grid)->B_array);
    free((*multi_grid)->P_array);
    free((*multi_grid)->solution);
    free((*multi_grid)->cg_rhs);
    free((*multi_grid)->cg_res);
    free((*multi_grid)->cg_p);
    free((*multi_grid)->cg_w);
    free((*multi_grid)->double_tmp);
    free((*multi_grid)->int_tmp);
    free((*multi_grid)->array_size_g);
    free((*multi_grid)->array_size_l);

    free(*multi_grid);
    *multi_grid = NULL;

    return 0;
}

static int Prolong_from_I_to_J(PASE_MULTIGRID multigrid, int level_i, int level_j,
                               int *mv_s, int *mv_e, void **pvx_i, void **pvx_j)
{
    int idx_level;
    OPS *gcge_ops = multigrid->gcge_ops;
    void **from_mv, **to_mv;
    int start[2], end[2];
    void ***tmp = multigrid->cg_res;
    void **P_array = multigrid->P_array;
    for (idx_level = level_i; idx_level > level_j; idx_level--)
    {
        if (idx_level == level_i)
        {
            from_mv = pvx_i;
            start[0] = mv_s[0];
            end[0] = mv_e[0];
        }
        else
        {
            from_mv = tmp[idx_level];
            start[0] = 0;
            end[0] = mv_e[0] - mv_s[0];
        }
        if (idx_level == level_j + 1)
        {
            to_mv = pvx_j;
            start[1] = mv_s[1];
            end[1] = mv_e[1];
        }
        else
        {
            to_mv = tmp[idx_level - 1];
            start[1] = 0;
            end[1] = mv_e[0] - mv_s[0];
        }
        gcge_ops->MatDotMultiVec(P_array[idx_level - 1], from_mv, to_mv, start, end, gcge_ops);
    }
    return 0;
}

static int Restrict_from_I_to_J(PASE_MULTIGRID multigrid, int level_i, int level_j,
                                int *mv_s, int *mv_e, void **pvx_i, void **pvx_j)
{
    int idx_level;
    OPS *gcge_ops = multigrid->gcge_ops;
    void **from_mv, **to_mv;
    int start[2], end[2];
    void ***tmp = multigrid->cg_res;
    void **P_array = multigrid->P_array;
    for (idx_level = level_i; idx_level < level_j; idx_level++)
    {
        if (idx_level == level_i)
        {
            from_mv = pvx_i;
            start[0] = mv_s[0];
            end[0] = mv_e[0];
        }
        else
        {
            from_mv = tmp[idx_level];
            start[0] = 0;
            end[0] = mv_e[0] - mv_s[0];
        }

        if (idx_level == level_j - 1)
        {
            to_mv = pvx_j;
            start[1] = mv_s[1];
            end[1] = mv_e[1];
        }
        else
        {
            to_mv = tmp[idx_level + 1];
            start[1] = 0;
            end[1] = mv_e[0] - mv_s[0];
        }
        gcge_ops->MatTransDotMultiVec(P_array[idx_level], from_mv, to_mv, start, end, gcge_ops);
    }
    return 0;
}

int PASE_MULTIGRID_fromItoJ(PASE_MULTIGRID multi_grid, int level_i, int level_j,
                            int *mv_s, int *mv_e, void **pvx_i, void **pvx_j)
{
#if PRINT_INFO
    multi_grid->gcge_ops->Printf("···· [ Level %d ] --> [ Level %d ] ····\n", level_i, level_j);
#endif
    if (level_i < level_j)
    {
        Restrict_from_I_to_J(multi_grid, level_i, level_j, mv_s, mv_e, pvx_i, pvx_j);
    }
    else if (level_i > level_j)
    {
        Prolong_from_I_to_J(multi_grid, level_i, level_j, mv_s, mv_e, pvx_i, pvx_j);
    }
    else if (level_i == level_j)
    {
        multi_grid->gcge_ops->MultiVecAxpby(1.0, pvx_i, 0.0, pvx_j, mv_s, mv_e, multi_grid->gcge_ops);
    }
    return 0;
}