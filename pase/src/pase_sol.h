#ifndef __PASE_SOL_H__
#define __PASE_SOL_H__

#include "pase_mg.h"
#include "pase_param.h"
#include "ops_eig_sol_gcg.h"
#include "ops.h"
#include "petscmat.h"
#include "pase_ops.h"

typedef struct PASE_MG_SOLVER_PRIVATE_
{
    //============================================
    /* multigrid类型 0-AMG 1-GMG */
    int multigrid_type;
    /* 预处理与否 0不做 1只做B 2只做A 3都做 */
    int if_precondition;
    /* 信息打印 0将打印关键步骤起始信息 1将打印重要结果 2将打印时间等更详细内容*/
    int print_level;
    /* 光滑解法 */
    char smooth_type[8];
    //======================= 网格体系 =====================
    /* 多重网格结构 */
    PASE_MULTIGRID multigrid;
    /* 层数 */
    int num_levels;
    /* 初始层 */
    int initial_level;
    /* 多重网格最粗层 */
    int mg_coarsest_level;
    /* 多重网格最细层 */
    int mg_finest_level;
    /* 辅助空间粗层 */
    int aux_coarse_level;
    /* 辅助空间细层 */
    int aux_fine_level;
    /* 当前层数 */
    int current_level;
    /* 当前cycle数 */
    int current_cycle;
    //======================= 迭代次数 =====================
    /* 最多的cycle数 */
    int max_cycle_count;
    /* pre smooth的迭代参数 */
    int max_pre_count;
    /* post smooth的迭代参数 */
    int max_post_count;
    /* 辅助空间gcg的迭代次数参数 */
    int max_direct_count;
    /* 初始求解gcg的迭代次数参数 */
    int max_initial_direct_count;
    //======================= 个数控制 =====================
    /* 用户要求的特征对个数 */
    int user_nev;
    /* 多求得特征对个数 */
    int more_nev;
    /* 给定的特征对个数 */
    int num_given_eigs;
    /* 最多算的特征对个数 */
    int max_nev;
    /* 在实际计算的时候求解的特征对个数（在本算法中一般设置为 nev+5 和2*nev的较小值）*/
    int pase_nev;
    /* 已经收敛的特征对个数 */
    int conv_nev;
    /* smooth时的lock的个数 */
    int nlock_smooth;
    /* aux_direct时的lock的个数 */
    int nlock_direct;
    /* aux_A->b_H中lock住的向量列数 */
    int nlock_auxmat_A;
    /* aux_B->b_H中lock住的向量列数 */
    int nlock_auxmat_B;
    //======================= 收敛准则 =====================
    /* 相对残差收敛准则 */
    /* ||A*x-\lambda*B*x||_2/(\lambda*||x||_2) < rtol */
    double rtol;
    /* 绝对残差收敛准则 */
    /* ||A*x-\lambda*B*x||_2/||x||_2 < atol */
    double atol;
    /* 初始求解gcg的残差收敛准则 */
    double initial_rtol;
    /* 辅助空间gcg的残差收敛准则 */
    double aux_rtol;
    /* 存储特征对的绝对残差 */
    double *abs_res_norm;
    //======================= 特征对 =======================
    /* 特征值 */
    double *eigenvalues;
    /* 每一层的向量空间 */
    void ***solution;
    //======================= 临时空间 =====================
    /* BCG/BMG 所需向量空间 */
    void ***cg_rhs;
    void ***cg_res;
    void ***cg_p;
    void ***cg_w;
    /* double型临时空间 */
    double *double_tmp;
    /* int型临时空间 */
    int *int_tmp;
    /* direct aux solver要用到的空间 */
    PASE_MultiVector aux_gcg_mv[4];
    //======================= ops ===========================
    /* 普通矩阵向量操作 */
    OPS *gcge_ops;
    /* pase aux 复合矩阵的矩阵向量操作 */
    PASE_OPS *pase_ops;
    /* pase_ops接到gcge的操作 */
    OPS *pase_ops_to_gcge;
    //======================= aux space =====================
    /* pase aux 复合向量 */
    PASE_MultiVector aux_sol;
    /* pase aux 复合矩阵 */
    PASE_Matrix aux_A;
    /* pase aux 复合矩阵 */
    PASE_Matrix aux_B;
    /* pase aux 复合矩阵特征值 */
    double *aux_eigenvalues;
    //======================= 统计时间 =====================
    double solver_setup_time;
    double get_initvec_time;
    double smooth_time;
    double build_aux_time;
    double prolong_time;
    double error_estimate_time;
    double aux_direct_solve_time;
    double precondition_time;
    double total_solve_time;
    double total_time;
    //======================================================
} PASE_MG_SOLVER_PRIVATE;
typedef PASE_MG_SOLVER_PRIVATE *PASE_MG_SOLVER;

int PASE_EigenSolver(void *A, void *B, double *eval, void **evec, int nev, PASE_PARAMETER param);
int PASE_Mg_set_up(PASE_MG_SOLVER solver, void *A, void *B);
int PASE_Mg_solve(PASE_MG_SOLVER solver);
int PASE_Direct_solve(PASE_MG_SOLVER solver);
int PASE_Mg_cycle(PASE_MG_SOLVER solver);
int PASE_Mg_smoothing(PASE_MG_SOLVER solver);
int PASE_Aux_direct_solve(PASE_MG_SOLVER solver);
int PASE_Mg_error_estimate(PASE_MG_SOLVER solver);

int PASE_Mg_prolong_from_Solution(PASE_MG_SOLVER solver, int object_level);
int PASE_Mg_prolong_from_pase_aux_Solution(PASE_MG_SOLVER solver, int object_level);
int PASE_Mg_set_pase_aux_vector(PASE_MG_SOLVER solver);
int PASE_Mg_set_pase_aux_matrix(PASE_MG_SOLVER solver);

PASE_MG_SOLVER PASE_Mg_solver_create(PASE_PARAMETER param, PASE_OPS *pase_ops);
int PASE_Mg_pase_aux_vector_create(PASE_MG_SOLVER solver);
int PASE_Mg_pase_aux_matrix_create(PASE_MG_SOLVER solver);
int PASE_Mg_pase_aux_matrix_factorization(PASE_MG_SOLVER solver);

int PASE_Matrix_destroy_sub(PASE_Matrix *aux_A);
int PASE_MultiVector_destroy_sub(PASE_MultiVector *aux_sol);
int PASE_Mg_solver_destroy(PASE_MG_SOLVER solver);

int PASE_Mg_print_param(PASE_MG_SOLVER solver);
int PASE_Mg_print_timer(PASE_MG_SOLVER solver);

#endif