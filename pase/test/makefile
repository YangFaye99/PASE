############################################################
GCGEHOME = ../..
PASEHOME = ..
include $(GCGEHOME)/config/make.MPI.inc
SOURCE = $(wildcard $(GCGEHOME)/app/*.c) $(wildcard $(GCGEHOME)/src/*.c) $(wildcard $(PASEHOME)/src/*.c)
OBJECTS = $(SOURCE:.c=.o)

MPIRUN = mpirun -np 10
DEBUGEXEC = 
DEBUGFLAGS = 

############################################################
OBJS0  = test0.o $(OBJECTS)
OBJS1  = test1.o $(OBJECTS)
OBJS2  = test2.o $(OBJECTS)

EXE0   = ./test0.exe 
EXE1   = ./test1.exe 
EXE2   = ./test2.exe 

#INC    = -I$(GCGEINC)
#LIBS   = $(LIBPASE) $(LIBGCGECSR) $(LIBGCGE) $(LIBLAPACK) $(LIBBLAS) $(LIBF)

INC    = -I$(GCGEHOME)/app -I$(GCGEHOME)/src -I$(PASEHOME)/src $(PETSCINC) $(SLEPCINC) $(HYPREINC) 
LIBS   = $(LIBPHG) $(LIBSLEPC) $(LIBPETSC) $(LIBHYPRE) $(LIBBLASLAPACK) 

############################################################

all:	help

ex0: $(EXE0)
ex1: $(EXE1)
ex2: $(EXE2)

clean:
	@$(RM) $(OBJECTS) *.o *.exe

help:
	@echo " "
	@echo "make {exe|run|clean|help}"
	@echo " "
	@echo "   xx   - 	create exe"
	@echo " "
	@echo "   runxx    - test_xx"
	@echo " "
	@echo "   clean - remove temporary files"
	@echo "   help  - print this message"
	@echo " "

run0:
	$(MPIRUN) $(DEBUGEXEC) $(DEBUGFLAGS) $(EXE0)
run1:
	$(MPIRUN) $(DEBUGEXEC) $(DEBUGFLAGS) $(EXE1)
run2:
	$(MPIRUN) $(DEBUGEXEC) $(DEBUGFLAGS) $(EXE2)

############################################################

.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	$(CC)  $(CFLAGS) -c $*.c -o $*.o $(INC)

############################################################

$(EXE0): $(OBJS0)
	@echo ''
	@echo 'example: Laplace'
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo ''

$(EXE1): $(OBJS1)
	@echo ''
	@echo 'example: Anderews'
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo ''

$(EXE2): $(OBJS2)
	@echo ''
	@echo 'example: Ga3As3H12'
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo ''

############################################################
