#ifndef _pase_mg_h_
#define _pase_mg_h_

#include "pase_ops.h"

typedef struct pase_MultiGrid_struct
{
    /* level */
    int num_levels;
    int coarsest_level;
    int pase_coarse_level;
    int pase_fine_level;
    /* size */
    int *array_size_g;
    int *array_size_l;
    /* matrix */
    void **A_array;
    void **B_array;
    void **P_array;
    /* vector */
    int pase_nev;
    int step_size;
    void ***solution;
    void ***cg_rhs;
    void ***cg_res;
    void ***cg_p;
    void ***cg_w;
    /* tmp */
    double *double_tmp;
    int *int_tmp;
    /* ops */
    OPS *gcge_ops;
} pase_MultiGrid;
typedef struct pase_MultiGrid_struct *PASE_MULTIGRID;

int PASE_MULTIGRID_Create(PASE_MULTIGRID *multi_grid, int num_levels, void *A, void *B, int pase_nev, OPS *gcge_ops);
int PASE_MULTIGRID_Destroy(PASE_MULTIGRID *multi_grid);

int PASE_MULTIGRID_TwoGirdLevel(PASE_MULTIGRID multi_grid);
int PASE_MULTIGRID_fromItoJ(PASE_MULTIGRID multi_grid, int level_i, int level_j,
                            int *mv_s, int *mv_e, void **pvx_i, void **pvx_j);

#endif