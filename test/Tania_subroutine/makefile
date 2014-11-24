#!/bin/make

###### Generic Makefile for fortran codes #####


# Compiler flag

COMP = ifort
#COMP = gfortran

# ifort flags

#WFLAGS   = -warn argument checking -check underflow
#WFLAGS   = -C -check all -warn all,nodec,interfaces -gen_interfaces -traceback # available with ifort 
#WFLAGS   = -check all -warn all -traceback
WFLAGS   = -warn none
#CFLAGS   = -c -g -heap-arrays # available with ifort
#OPTFLAGS = -O3 -mtune=core2
#PARFLAGS = -openmp
#DEBUGFLAGS = -g # available with ifort

# gfortran flags

#WFLAGS     = -Wall -Wextra
#DEBUGFLAGS = -g -fbacktrace
#OPTFLAGS   = -O3

# common flags

#CFLAGS  = -c -g
CFLAGS = -c
LDFLAGS = -o
#PROFLAGS = -pg

# Rules for compilation
all: ROLAND.DD mrproper

# compiling the sources

ROLAND.DD:
	@$(COMP) ROLAND.DD.f90 $(LDFLAGS) ROLAND.DD.out
	@echo "--------Compilation ROLAND.DD done-----------"

mrproper: clean data

clean:
	@rm -rf *.mod
	@rm -rf *.o
	@rm -rf *genmod.f90
	@rm -rf *.f90~
data:
	@rm ./*.dat
	@echo "--------Data cleaning done-------------------"

