# NPOT 

Nonconservative compressible potential flow solver.

    confpc -y1
      dymin = 0.02
      nx, ny = 256 63
      xmax, ymax = 1500 1500

To build do the following:

    ln -s <platform>.mak Makefile
    make

Where, for example, `<platform>` is `gcc` when using `gfortran` on Mac OS X

S. Scott Collis\
