# Project: TestOPS
#
# Makefile created by Dev-C++ 5.11
PETSCFLAGS = -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden 
PETSCINC   = -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
LIBPETSC   = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/local/lib -L/usr/local/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lpetsc -lflapack -lfblas -lpthread -lm -lX11 -lstdc++ -ldl -lmpifort -lmpi -lgfortran -lgcc_s -lquadmath 

SLEPCFLAGS = -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden 
SLEPCINC   = -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
LIBSLEPC   = -Wl,-rpath,${SLEPC_DIR}/${PETSC_ARCH}/lib -L${SLEPC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L/usr/local/lib -L/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu  -Wl,-rpath,/usr/local/lib -L/usr/local/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lslepc -lpetsc -lflapack -lfblas -lpthread -lm -lX11 -lstdc++ -ldl -lmpifort -lmpi -lgfortran -lgcc_s -lquadmath

OBJ      = ../src/ops.o ../src/ops_eig_sol_gcg.o ../src/ops_lin_sol.o ../src/ops_multi_vec.o ../src/ops_orth.o ../app/app_lapack.o ../app/app_ccs.o ../app/app_slepc.o
LINKOBJ  = ../src/ops.o ../src/ops_eig_sol_gcg.o ../src/ops_lin_sol.o ../src/ops_multi_vec.o ../src/ops_orth.o ../app/app_lapack.o ../app/app_ccs.o ../app/app_slepc.o

LIBBLASLAPACK = -lblas -llapack
LIB = -lpthread -lmpi -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl
LIBS = $(LIBBLASLAPACK) $(LIB) -lslepc -lpetsc

INCS     = -I../src -I../app
CXXINCS  = 
CXXFLAGS = $(CXXINCS) -O2 -Wall -fPIC -g  
CFLAGS   = $(INCS) -O2 -Wall -fPIC -g  -fopenmp #-qopenmp -DMEMWATCH -DMW_STDIO ##
RM       = rm -f

CPP      = mpic++ -D__DEBUG__
CC       = mpicc  -D__DEBUG__

clean: 
	rm -f *.om

all:$(OBJ)
	$(CC) -shared -o libgcge.so  $(OBJ) $(LIBS) $(LIBPETSC) $(LIBSLEPC) $(SLEPCINC) $(SLEPCFLAGS)


../src/ops.o: ../src/ops.c
	$(CC) -c ../src/ops.c -o ../src/ops.o $(CFLAGS)

../src/ops_eig_sol_gcg.o: ../src/ops_eig_sol_gcg.c
	$(CC) -c ../src/ops_eig_sol_gcg.c -o ../src/ops_eig_sol_gcg.o $(CFLAGS)

../src/ops_lin_sol.o: ../src/ops_lin_sol.c
	$(CC) -c ../src/ops_lin_sol.c -o ../src/ops_lin_sol.o $(CFLAGS)

../src/ops_multi_vec.o: ../src/ops_multi_vec.c
	$(CC) -c ../src/ops_multi_vec.c -o ../src/ops_multi_vec.o $(CFLAGS)

../src/ops_orth.o: ../src/ops_orth.c
	$(CC) -c ../src/ops_orth.c -o ../src/ops_orth.o $(CFLAGS)

../app/app_ccs.o: ../app/app_ccs.c
	$(CC) -c ../app/app_ccs.c -o ../app/app_ccs.o $(CFLAGS) $(LIBBLASLAPACK)

../app/app_lapack.o: ../app/app_lapack.c
	$(CC) -c ../app/app_lapack.c -o ../app/app_lapack.o $(CFLAGS)  $(LIBBLASLAPACK)

../app/app_slepc.o: ../app/app_slepc.c
	$(CC) -c ../app/app_slepc.c -o ../app/app_slepc.o $(CFLAGS)  $(LIBSLEPC) $(SLEPCINC) $(SLEPCFLAGS)
