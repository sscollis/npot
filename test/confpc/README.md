## Conformal mesh for a parabolic cylinder

Mach 0.4 potential flow past a parabolic cylinder

### Steps to run

1. Check the paths in `run.sh`
2. Note that you need the `confpc` mesh generator
3. Also note that `confpc` uses the old `ji` ordering so
   that you need to use the `-ms` option when running `npot`
4. `./run.sh`

Note that to converge this it is helpful to restart with
several different values of `sigma` 

### Visualization

To use paraview you may want to do:

    ln -s grid.dat grid.xyz

### Sample Results

#### Contours of streamwise velocity

![Streamwise velocity](https://github.com/sscollis/npot/blob/master/test/confpc/u.png)

#### Computational Mesh

![Mesh](https://github.com/sscollis/npot/blob/master/test/confpc/mesh.png)

#### Close-up of Computational Mesh

![Mesh](https://github.com/sscollis/npot/blob/master/test/confpc/mesh-cu.png)

S. Scott Collis \
Wed Feb 19 06:22:30 MST 2020
