#include "pase.h"
#include "app_slepc.h"

static char help[] = "Test with Ga3As3H12 Mat.\n";

void MatrixRead(Mat *A, Mat *B);

int main(int argc, char *argv[])
{
    PetscErrorCode ierr;
    SlepcInitialize(&argc, &argv, (char *)0, help);
    srand(1);

    Mat A, B;
    MatrixRead(&A, &B);

    int nev = 10, num_levels = 5;
    double *eval = NULL;
    void **evec = NULL;
    PASE_PARAMETER param;
    PASE_PARAMETER_Create(&param, num_levels, nev);
    
    PASE_EigenSolver((void *)A, (void *)B, eval, evec, nev, param);

    PASE_PARAMETER_Destroy(&param);

    ierr = SlepcFinalize();
    return 0;
}

void MatrixRead(Mat *A, Mat *B)
{
    Mat matA;
    char filename[128] = "../data/Ga3As3H12.petsc.bin";
    PetscInt nrows, ncols, row_s, row_e;
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer);
    MatCreate(PETSC_COMM_WORLD, &matA);
    MatLoad(matA, viewer);
    PetscViewerDestroy(&viewer);
    MatGetSize(matA, &nrows, &ncols);
    MatGetOwnershipRange(matA, &row_s, &row_e);
    *A = matA;

    Mat matB;
    MatCreate(PETSC_COMM_WORLD, &matB);
    MatSetSizes(matB, row_e - row_s, row_e - row_s, PETSC_DECIDE, PETSC_DECIDE);
    MatSetFromOptions(matB);
    MatSetUp(matB);
    int II;
    for (II = row_s; II < row_e; II++)
    {
        MatSetValue(matB, II, II, 1.0, INSERT_VALUES);
    }
    MatAssemblyBegin(matB, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matB, MAT_FINAL_ASSEMBLY);
    *B = matB;

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0)
    {
        printf("*************************************************\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("\t[Ga3As3H12] matrix size: %d x %d\n", nrows, ncols);
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("*************************************************\n\n");
    }
}