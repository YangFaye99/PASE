HYPREFLAGS = -g -Wall -DHAVE_CONFIG_H -DHYPRE_TIMING
HYPREINC   = -I${HYPRE_DIR}/include
LIBHYPRE   = -L${HYPRE_DIR}/lib -lHYPRE -lm -lstdc++

PETSCFLAGS = -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden 
PETSCINC   = -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
LIBPETSC   = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/local/lib -L/usr/local/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lpetsc -lflapack -llapack -lfblas -lpthread -lm -lstdc++ -ldl -lmpifort -lmpi -lgfortran -lgcc_s -lquadmath 

SLEPCFLAGS = -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden 
SLEPCINC   = -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
LIBSLEPC   = -Wl,-rpath,${SLEPC_DIR}/${PETSC_ARCH}/lib -L${SLEPC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L/usr/local/lib -L/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu  -Wl,-rpath,/usr/local/lib -L/usr/local/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lslepc -lpetsc -lflapack -llapack -lfblas -lpthread -lm -lstdc++ -ldl -lmpifort -lmpi -lgfortran -lgcc_s -lquadmath

PHGFLAGS   = 
PHGINC     = 
LIBPHG     = 

#MKLINC     = -I${MKL_ROOT}/include
#LIBMKL     = -L${MKL_ROOT}/lib/intel64 -Wl,-rpath,${MKL_ROOT}/lib/intel64 -lmkl_intel_lp64 -L/home/lyf/intel/oneapi/compiler/2022.0.2/linux/lib -Wl,-rpath,/home/lyf/intel/oneapi/compiler/2022.0.2/linux/lib -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

CPP      = mpic++ -D__DEBUG__
CC       = mpicc 
LDFLAGS  = 

vpath %.h $(GCGEHOME)/src
vpath %.h $(GCGEHOME)/app
vpath %.c $(GCGEHOME)/src
vpath %.c $(GCGEHOME)/app