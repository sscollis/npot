#==============================================================================
#
#  Makefile for npot (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 2-13-2020 
#
#==============================================================================
NAME   = npot 
DEBUG  = 
FFLAGS = -cpp -fdefault-real-8 -c $(DEBUG)
OPT    = -O2 -fopenmp
OFLAGS = $(OPT) $(DEBUG) -o $(NAME)
LIB    = -L$(HOME)/local/OpenBLAS/lib -lopenblas
COMP   = gfortran 
F77    = gfortran
CC     = gcc-9
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
# Modules
#
MODS = global.o consts.o stencil.o
#
#  These objects depend on global.f90
#
OBJS = npot.o error.o grad.o grad2.o dtcfl.o smoother.o itrbc.o \
       rhsbc.o resstat.o initial.o input.o lhs1.o lhs2.o
OBJ2 = penta1bc.o penta2bc.o penta1p_blk.o penta2p_blk.o \
       penta1bc_blk.o penta2bc_blk.o solve.o
#
#  End of objects
#
$(NAME): $(MODS) $(OBJS) $(OBJ2) 
	 $(COMP) $(OPT) $(MODS) $(OFLAGS) $(OBJS) $(OBJ2) $(LIB)

clean:
	$(RM) *.o *.mod

$(OBJ2):

$(OBJS): $(MODS) 

$(MODS):

.f90.o:
	$(COMP) $(OPT) $(FFLAGS) $*.f90 

.f.o:
	$(F77) $(OPT) $(FFLAGS) -ffixed-line-length-120 -std=legacy $*.f
