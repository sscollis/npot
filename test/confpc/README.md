## Conformal mesh for a parabolic cylinder

1. Check the paths in `run.sh`
2. Note that you need the `confpc` mesh generator
3. Also note that `confpc` uses the old `ji` ordering so
   that you need to use the `-ms` option when running `npot`
4. `./run.sh`

Note that to converge this it is helpful to restart with
several different values of `sigma` 

To use paraview you may want to do:

    ln -s grid.dat grid.xyz

S. Scott Collis \
Wed Feb 19 06:22:30 MST 2020