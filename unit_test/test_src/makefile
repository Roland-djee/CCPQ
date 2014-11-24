#############################################################################
#                                                                           #
#	MAKEFILE GENERIQUE POUR LA COMPILATION DES CODES FORTRAN            #
#                                                                           #
#############################################################################

# Indiquer quel compilateur est à utiliser

COMP = ifort
#COMP = gfortran

# Spécifier les options du compilateur CFLAG = Compilation FLAGS,
# WFLAGS = Warning FLAGS, LDFLAGS = Link eDition FLAGS, OPTFLAGS = OPTimizing
# FLAGS, DEBUGFLAGS = DEBUGging FLAGS, LAPACK = édition de liens avec la 
# librairie LAPACKS, PROFLAGS = PROfiling FLAGS

# ifort FLAGS

#WFLAGS   = -warn argument checking -check underflow
WFLAGS   = -C -check all -warn all,nodec,interfaces -gen_interfaces -traceback # available with ifort 
#WFLAGS   = -check all -warn all -traceback
#CFLAGS   = -c -g -heap-arrays # available with ifort
#LAPACK   = -L/usr/atlas/lib -llapack -lblas -lgfortran
#LAPACK   = -L/usr/lib -llapack -lblas
#OPTFLAGS = -O3 -mtune=core2
PARFLAGS = -openmp
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

# Reconnaître les Extension de nom de fichier *.out, *.o, .f90 et .h
# comme suffixe

SUFFIXES = .exe .o .f90 
.SUFFIXES: $(SUFFIXES) .

# Nom de l'exécutable

EXEC1 = test_read
EXEC2 = test_crystal
EXEC3 = test_interactions
EXEC4 = test_decoherence
#EXEC = compare_decoherence

# Liste des fichiers modules nécessaires pour le programme final

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

# Listes de fichiers sources nécessaires pour le programme final

SRC1  = test_read.f90
SRC2  = test_crystal.f90 
SRC3  = test_interactions.f90 
SRC4  = test_decoherence.f90 
#SRC  = compare_decoherence.f90 

# Liste de fichiers objets nécessaires pour le programme final.
# Il s'agit des mêmes objets que ceux source avec une extension
# .o.

OBJS1  = $(MOD2)\
	$(SRC1:.f90=.o)
OBJS2  = $(MOD2)\
	$(SRC2:.f90=.o)
OBJS3  = $(MOD2)\
	$(SRC3:.f90=.o)
OBJS4  = $(MOD2)\
	$(SRC4:.f90=.o)

all: $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4) Class Mrproper
 
# Etape de compilation de tous les fichiers modules .f90 en .o

Modules:
	@$(COMP) $(CFLAGS) $(PARFLAGS) $(WFLAGS) $(MOD) $(PROFLAGS)

# Etape de compilation de tous les fichiers sources .f90 en .o

Sources1:
	@$(COMP) $(PARFLAGS) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $(SRC1) $(PROFLAGS)
Sources2:
	@$(COMP) $(PARFLAGS) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $(SRC2) $(PROFLAGS)
Sources3:
	@$(COMP) $(PARFLAGS) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $(SRC3) $(PROFLAGS)
Sources4:
	@$(COMP) $(PARFLAGS) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $(SRC4) $(PROFLAGS)

# Etape d'édition de liens

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

# Classement des fichiers

Class:
	@mv *.o ../items/.
	@mv *.exe ../
#	@mv *.f90~ ../Save/.
	@echo "----Data classification done----"

# Nettoyage

Clean:
	@rm -rf *.mod
	@rm -rf *genmod.f90
	@rm -rf *.f90~

Mrproper: Clean
	@rm -rf $(EXEC1)
	@echo "-------Data cleaning done-------"

Data:
	@rm -rf ../*.dat
	@echo "-------Data cleaning done-------"

