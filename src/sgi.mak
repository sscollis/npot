#==============================================================================
#
#  Makefile for npot (for SGI)
#
#  Author:  Scott Collis
#
#  Revised: 3-17-96
#
#==============================================================================
NAME   = npot 
DEBUG  =
FFLAGS = -ftpp -r8 -c $(DEBUG)
OPT    = -O2 -OPT:Olimit=0 -mp -MP:old_mp=OFF
#OPT   = -g
OFLAGS = -r8 $(DEBUG) -o $(NAME)
LIB    =
COMP   = f90
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
	/bin/rm *.o *.mod

$(OBJ2):

$(OBJS): $(MODS) 

$(MODS):

.f90.o:
	$(COMP) $(OPT) $(FFLAGS) $*.f90 

.f.o:
	$(COMP) $(OPT) -extend_source $(FFLAGS) $*.f
