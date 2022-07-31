#ifndef _PASE_PARAM_H_
#define _PASE_PARAM_H_

#include <stdbool.h>

typedef struct PASE_PARAMETER_PRIVATE_
{
    //============================================
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
    //============================================
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
    //============================================
    /* 要求的特征对个数 */
    int nev;
    /* 多求得特征对个数 */
    int more_nev;
    /* 给定的特征对个数 */
    int num_given_eigs;
    //============================================
    /* 相对残差收敛准则 */
    double rtol;
    /* 绝对残差收敛准则 */
    double atol;
    /* 初始求解gcg的残差收敛准则 */
    double initial_rtol;
    /* 辅助空间gcg的残差收敛准则 */
    double aux_rtol;
    //============================================
    /* multigrid类型 0-AMG 1-GMG */
    int multigrid_type;
    /* 预处理与否 */
    int if_precondition;
    /* 信息打印 */
    int print_level;
    //============================================
} PASE_PARAMETER_PRIVATE;

typedef PASE_PARAMETER_PRIVATE *PASE_PARAMETER;

void PASE_PARAMETER_Create(PASE_PARAMETER *param, int num_levels, int nev);

void PASE_PARAMETER_Destroy(PASE_PARAMETER *param);

void PASE_PARAMETER_Get_from_command_line(PASE_PARAMETER param, int argc, char *argv[]);

#endif
