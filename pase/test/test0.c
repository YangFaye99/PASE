#include "pase.h"
#include "app_slepc.h"

static char help[] = "Test with Laplace Mat.\n";

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

    Mat matB;
    MatCreate(PETSC_COMM_WORLD, &matB);
    MatSetSizes(matB, PETSC_DECIDE, PETSC_DECIDE, N, N);
    MatSetFromOptions(matB);
    MatSetUp(matB);
    for (II = Istart; II < Iend; II++)
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
        printf("[MatrixRead] matrix size: %d x %d\n", N, N);
    }
}