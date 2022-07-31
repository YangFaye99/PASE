#include "aug_lin.h"
#include "app_slepc.h"

static char help[] = "Test with Andrews Mat.\n";

void MatrixRead(Mat *A);

int main(int argc, char *argv[])
{
    PetscErrorCode ierr;
    SlepcInitialize(&argc, &argv, (char *)0, help);
    srand(1);

    Mat A;
    MatrixRead(&A);

    OPS *gcge_ops;
    OPS_Create(&gcge_ops);
    OPS_SLEPC_Set(gcge_ops);
    OPS_Setup(gcge_ops);

    BV b, x;
    gcge_ops->MultiVecCreateByMat((void ***)(&b), 1, (void *)A, gcge_ops);
    gcge_ops->MultiVecCreateByMat((void ***)(&x), 1, (void *)A, gcge_ops);
    gcge_ops->MultiVecSetRandomValue((void **)b, 0, 1, gcge_ops);
    BVNormalize(b, NULL);

    AugLinSolver *solver;
    AugLinSolver_Create(&solver, (void *)A, (void **)b, 10, 1e-8);
    // AugLinSolve(solver, (void **)x, "gcge_pcg");
    AugLinSolve(solver, (void **)x, "aug_gcge_pcg");
    AugLinSolve(solver, solver->x_array[0], "gcge_pcg");
    AugLinSolver_Free(&solver);

    ierr = SlepcFinalize();
    return 0;
}

void MatrixRead(Mat *A)
{
#if 1
    Mat matA;
    char filename[128] = "${HOME}/matdata/gridgena.petsc.bin";
    PetscInt nrows, ncols, row_s, row_e;
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer);
    MatCreate(PETSC_COMM_WORLD, &matA);
    MatLoad(matA, viewer);
    PetscViewerDestroy(&viewer);
    MatGetSize(matA, &nrows, &ncols);
    MatGetOwnershipRange(matA, &row_s, &row_e);
    *A = matA;

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0)
    {
        printf("[MatrixRead] matrix size: %d x %d\n", nrows, ncols);
    }
#else
    int n = 300;
    int N = n * n;

    Mat matA;
    MatCreate(PETSC_COMM_WORLD, &matA);
    MatSetSizes(matA, PETSC_DECIDE, PETSC_DECIDE, N, N);
    MatSetFromOptions(matA);
    MatSetUp(matA);
    int II, Istart, Iend, i, j;
    MatGetOwnershipRange(matA, &Istart, &Iend);
    for (II = Istart; II < Iend; II++)
    {
        i = II / n;
        j = II - i * n;
        if (i > 0)
        {
            MatSetValue(matA, II, II - n, -1.0, INSERT_VALUES);
        }
        if (i < n - 1)
        {
            MatSetValue(matA, II, II + n, -1.0, INSERT_VALUES);
        }
        if (j > 0)
        {
            MatSetValue(matA, II, II - 1, -1.0, INSERT_VALUES);
        }
        if (j < n - 1)
        {
            MatSetValue(matA, II, II + 1, -1.0, INSERT_VALUES);
        }
        MatSetValue(matA, II, II, 4.0, INSERT_VALUES);
    }
    MatAssemblyBegin(matA, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matA, MAT_FINAL_ASSEMBLY);
    for (II = Istart; II < Iend; II++)
    {
        i = II / n;
        j = II - i * n;
        if (i > 0 && i < n - 1 && j > 0 && j < n - 1)
        {
        }
        else
        {
            if (i > 0)
            {
                MatSetValue(matA, II, II - n, 0.0, INSERT_VALUES);
                MatSetValue(matA, II - n, II, 0.0, INSERT_VALUES);
            }
            if (i < n - 1)
            {
                MatSetValue(matA, II, II + n, 0.0, INSERT_VALUES);
                MatSetValue(matA, II + n, II, 0.0, INSERT_VALUES);
            }
            if (j > 0)
            {
                MatSetValue(matA, II, II - 1, 0.0, INSERT_VALUES);
                MatSetValue(matA, II - 1, II, 0.0, INSERT_VALUES);
            }
            if (j < n - 1)
            {
                MatSetValue(matA, II, II + 1, 0.0, INSERT_VALUES);
                MatSetValue(matA, II + 1, II, 0.0, INSERT_VALUES);
            }
            MatSetValue(matA, II, II, 1.0, INSERT_VALUES);
        }
    }
    MatAssemblyBegin(matA, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matA, MAT_FINAL_ASSEMBLY);
    *A = matA;
#endif
}