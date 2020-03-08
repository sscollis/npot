### O-grid for a circular cylinder

Mach 0.4 compressible, potential flow past a circular cylinder using an O-grid.

![Streamwise velocity](https://github.com/sscollis/npot/blob/master/test/ogrid/ogrid.png)

![Mesh](https://github.com/sscollis/npot/blob/master/test/ogrid/mesh.png)

### Steps to make the run

1. Check the paths in `run.sh`
2. Note that you need the cyl mesh generator and
   pre pre-processors from `lns3d`
3. `./run.sh`

To use paraview you may want to do:

    ln -s grid.dat grid.xyz

S. Scott Collis\
Tue Feb 18 06:38:00 MST 2020
