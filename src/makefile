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

EXEC = Decoherence

# Liste des fichiers modules nécessaires pour le programme final

MOD = ../lib/Types.f90\
	../lib/Constants.f90\
	../lib/Read.f90\
	../lib/Crystal.f90\
	../lib/Write.f90\
	../lib/Interactions.f90\
	../lib/Deco.f90

MOD2 = ./Types.o\
	./Constants.o\
	./Read.o\
	./Crystal.o\
	./Write.o\
	./Interactions.o\
	./Deco.o

# Listes de fichiers sources nécessaires pour le programme final

SRC  = Decoherence.f90 

# Liste de fichiers objets nécessaires pour le programme final.
# Il s'agit des mêmes objets que ceux source avec une extension
# .o.

OBJS  = $(MOD2)\
	$(SRC:.f90=.o)

all: $(EXEC) Class Mrproper
 
# Etape de compilation de tous les fichiers modules .f90 en .o

Modules:
	@$(COMP) $(CFLAGS) $(PARFLAGS) $(WFLAGS) $(MOD) $(PROFLAGS)

# Etape de compilation de tous les fichiers sources .f90 en .o

Sources:
	@$(COMP) $(PARFLAGS) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $(SRC) $(PROFLAGS)

# Etape d'édition de liens

Decoherence: Modules Sources 
	@$(COMP) $(DEBUGFLAGS) $(OBJS) $(LAPACK) $(PARFLAGS) $(LDFLAGS) $(EXEC).exe $(PROFLAGS)
	@echo "--------Compilation done--------"

# Classement des fichiers

Class:
	@mv *.o ../items/.
	@mv *.exe ../bin/.
#	@mv *.f90~ ../Save/.
	@echo "----Data classification done----"

# Nettoyage

Clean:
	@rm -rf *.mod

Mrproper: Clean
	@rm -rf $(EXEC)
	@echo "-------Data cleaning done-------"

