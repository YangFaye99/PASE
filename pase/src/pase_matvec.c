#include "pase_matvec.h"
#include "ops.h"
#include "app_slepc.h"
#include "ops_eig_sol_gcg.h"

void Mat_Output(void *A, char *name)
{
    Mat mat = (Mat)A;
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, name, FILE_MODE_WRITE, &viewer);
    MatView(mat, viewer);
    PetscViewerDestroy(&viewer);
}

void MultiVec_Output(void **MV, int vec_num, char *name)
{
    BV mv = (BV)MV;
    Mat tmp;
    BVSetActiveColumns(mv, 0, vec_num);
    BVGetMat(mv, &tmp);
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, name, FILE_MODE_WRITE, &viewer);
    MatView(tmp, viewer);
    PetscViewerDestroy(&viewer);
    BVRestoreMat(mv, &tmp);
}

void PASE_Mat_Output(PASE_Matrix aux_A, char *name)
{
    int num = aux_A->num_aux_vec;

    char mat_name[20], end_name[10];
    strcpy(mat_name, name);
    strcpy(end_name, "_H");
    strcat(mat_name, end_name);
    PetscViewer viewer_H;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, mat_name, FILE_MODE_WRITE, &viewer_H);
    MatView(aux_A->A_H, viewer_H);
    PetscViewerDestroy(&viewer_H);
    /*辅助空间多向量*/
    char vec_name[20];
    strcpy(vec_name, name);
    strcpy(end_name, "_BV");
    strcat(vec_name, end_name);
    Mat tmp;
    BVSetActiveColumns((BV)(aux_A->aux_Hh), 0, num);
    BVGetMat((BV)(aux_A->aux_Hh), &tmp);
    PetscViewer viewer_BV;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, vec_name, FILE_MODE_WRITE, &viewer_BV);
    MatView(tmp, viewer_BV);
    PetscViewerDestroy(&viewer_BV);
    BVRestoreMat((BV)(aux_A->aux_Hh), &tmp);
    /*辅助空间并行部分*/
    int idx, idy;
    for (idx = 0; idx < num; idx++)
    {
        for (idy = 0; idy < num; idy++)
        {
            DefaultPrintf("%s_dense(%d,%d)=%2.18f;\n", name, idy + 1, idx + 1,
                          aux_A->aux_hh[idx * num + idy]);
        }
    }
    return;
}

void PASE_MultiVec_Output(PASE_MultiVector MV, int vec_num, char *name)
{
    char mat_name[20], end_name[10];
    strcpy(mat_name, name);
    strcpy(end_name, "_H");
    strcat(mat_name, end_name);
    MultiVec_Output(MV->b_H, vec_num, mat_name);
    int idx, idy;
    int num_vec = MV->num_vec;
    int num_aux_vec = MV->num_aux_vec;
    for (idx = 0; idx < vec_num; idx++)
    {
        for (idy = 0; idy < num_aux_vec; idy++)
        {
            DefaultPrintf("%s_dense(%d,%d)=%2.18f;\n", name, idy + 1, idx + 1,
                          MV->aux_h[idx * num_aux_vec + idy]);
        }
    }
}