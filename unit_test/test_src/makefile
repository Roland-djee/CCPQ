#------------- Generic makefile for fortran codes -------------#

# compiler flags

COMP = ifort
#COMP = gfortran

# ifort FLAGS

#WFLAGS   = -warn argument checking -check underflow
#WFLAGS   = -C -check all -warn all,nodec,interfaces -gen_interfaces -traceback # available with ifort 
#WFLAGS   = -check all -warn all -traceback
WFLAGS     = -nowarn
#CFLAGS   = -c -g -heap-arrays # available with ifort
#LAPACK   = -L/usr/atlas/lib -llapack -lblas -lgfortran
#LAPACK   = -L/usr/lib -llapack -lblas
#OPTFLAGS = -O3 -mtune=core2
#PARFLAGS = -openmp
OPTFLAGS = -O3 
#DEBUGFLAGS = -g # available with ifort

# gfortran FLAGS

#WFLAGS     = -Wall -Wextra
#DEBUGFLAGS = -g -fbacktrace
#OPTFLAGS   = -O3

#LAPACK = -L/usr/lib  -llapack -lblas

# Common FLAGS

#CFLAGS  = -c -g
CFLAGS = -c
LDFLAGS = -o
#PROFLAGS = -pg

# Suffixes
SUFFIXES = .exe .o .f90 
.SUFFIXES: $(SUFFIXES) .

# executable

EXEC1 = test_read
EXEC2 = test_crystal
EXEC3 = test_interactions
EXEC4 = test_decoherence
#EXEC = compare_decoherence

# modules lists

MOD = ../../lib/Types.f90\
	../../lib/Constants.f90\
	../../lib/Read.f90\
	../../lib/Crystal.f90\
	../../lib/Write.f90\
	../../lib/Interactions.f90\
	../../lib/Deco.f90

MOD2 = ./Types.o\
	./Constants.o\
	./Read.o\
	./Crystal.o\
	./Write.o\
	./Interactions.o\
	./Deco.o

# source files

SRC1  = test_read.f90
SRC2  = test_crystal.f90 
SRC3  = test_interactions.f90 
SRC4  = test_decoherence.f90 

OBJS1  = $(MOD2)\
	$(SRC1:.f90=.o)
OBJS2  = $(MOD2)\
	$(SRC2:.f90=.o)
OBJS3  = $(MOD2)\
	$(SRC3:.f90=.o)
OBJS4  = $(MOD2)\
	$(SRC4:.f90=.o)

all: $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4) Class Mrproper

# compiling modules

Modules:
	@$(COMP) $(CFLAGS) $(PARFLAGS) $(WFLAGS) $(MOD) $(PROFLAGS)

# compiling sources

Sources1:
	@$(COMP) $(PARFLAGS) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $(SRC1) $(PROFLAGS)
Sources2:
	@$(COMP) $(PARFLAGS) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $(SRC2) $(PROFLAGS)
Sources3:
	@$(COMP) $(PARFLAGS) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $(SRC3) $(PROFLAGS)
Sources4:
	@$(COMP) $(PARFLAGS) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $(SRC4) $(PROFLAGS)

# linking

test_read: Modules Sources1 
	@$(COMP) $(DEBUGFLAGS) $(OBJS1) $(LAPACK) $(PARFLAGS) $(LDFLAGS) $(EXEC1).exe $(PROFLAGS)
	@echo "--------Compilation test_read done--------"
test_crystal: Sources2 
	@$(COMP) $(DEBUGFLAGS) $(OBJS2) $(LAPACK) $(PARFLAGS) $(LDFLAGS) $(EXEC2).exe $(PROFLAGS)
	@echo "--------Compilation test_crystal done--------"
test_interactions: Sources3 
	@$(COMP) $(DEBUGFLAGS) $(OBJS3) $(LAPACK) $(PARFLAGS) $(LDFLAGS) $(EXEC3).exe $(PROFLAGS)
	@echo "--------Compilation test_interactions done--------"
test_decoherence: Sources4 
	@$(COMP) $(DEBUGFLAGS) $(OBJS4) $(LAPACK) $(PARFLAGS) $(LDFLAGS) $(EXEC4).exe $(PROFLAGS)
	@echo "--------Compilation test_decoherence done--------"

# sorting files

Class:
	@mv *.o ../items/.
	@mv *.exe ../
#	@mv *.f90~ ../Save/.
	@echo "--------Data classification done----"

# cleaning

Clean:
	@rm -rf *.mod
	@rm -rf *genmod.f90
	@rm -rf *.f90~

Mrproper: Clean
	@rm -rf $(EXEC1)
	@echo "--------Data cleaning done-------"

Data:
	@rm -rf ../*.dat
	@echo "--------Data cleaning done-------"

