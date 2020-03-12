This is a lifting airfoil so there must be circulation:

On Intel Mac with gfortran:

env GFORTRAN_CONVERT_UNIT='swap' ./npot -wc < input.dat

env GFORTRAN_CONVERT_UNIT='swap' ./npot -wc < input.res

The mesh was clearly made with chimera/hypgen by I haven't yet been 
able to duplicate it.

S. Scott Collis
Sat Feb 15 13:49:45 MST 2020
