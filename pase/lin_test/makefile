############################################################
GCGEHOME = ../..
PASEHOME = ..
include $(GCGEHOME)/config/make.MPI.inc
SOURCE = $(wildcard $(GCGEHOME)/app/*.c) $(wildcard $(GCGEHOME)/src/*.c) $(wildcard $(PASEHOME)/src/*.c) $(PASEHOME)/lin_test/aug_lin.c
OBJECTS = $(SOURCE:.c=.o)

MPIRUN = mpirun -np 8
DEBUGEXEC = 
DEBUGFLAGS = 

############################################################
OBJS0  = test.o $(OBJECTS)

EXE0   = ./test.exe 

#INC    = -I$(GCGEINC)
#LIBS   = $(LIBPASE) $(LIBGCGECSR) $(LIBGCGE) $(LIBLAPACK) $(LIBBLAS) $(LIBF)

INC    = -I$(GCGEHOME)/app -I$(GCGEHOME)/src -I$(PASEHOME)/src -I$(PASEHOME)/lin_test $(PETSCINC) $(SLEPCINC) $(HYPREINC) $(MKLINC)
LIBS   = $(LIBPHG) $(LIBSLEPC) $(LIBPETSC) $(LIBHYPRE) $(LIBMKL) 

############################################################

all:	help

ex0: $(EXE0)

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

############################################################

.c.o:
	@echo '$(CC)  -c $(CFLAGS) ' $*.c
	$(CC)  $(CFLAGS) -c $*.c -o $*.o $(INC)

############################################################

$(EXE0): $(OBJS0)
	@echo ''
	@echo 'example: Andrews'
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo ''

############################################################
