## NPOT

Solves compressible potential flow equations in nonconservative form 
using a high-order finite-difference method.

Directory   |   Description 
------------|-------------------------------
`src`       | Main source files
`test`      | Examples and test cases

To build:

    cd src
    ln -s <platform>.mak Makfile
    make

where, for example, <platform> is `gcc` when using `gfortran` on Mac OS X.

Notes:

  1. You may need to copy the `<platform>.mak` that is closest to your
     setup to another `name.mak` and edit it before proceeding.
  2. When running the `test` cases, it is recommended that you experiment
     with different values of `sigma`, the pseudo-timestep on multiple
     restarts as a means of accelerating convergence to steady-state.
  3. While it would be nice and relatively easy to automate this, that
     remains a todo item.

S. Scott Collis \
