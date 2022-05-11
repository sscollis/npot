### O-grid for a circular cylinder

Mach 0.4 compressible, potential flow past a circular cylinder using an O-grid.

### Steps to make the run

1. Check the paths in `run.sh`
2. Note that you need the cyl mesh generator and
   pre pre-processors from `lns3d`
3. `./run.sh`

### Visualization

To use paraview you may want to do:

    ln -s grid.dat grid.xyz

### Sample Results

#### Contours of streamwise velocity

![Streamwise velocity](https://github.com/sscollis/npot/blob/master/test/ogrid/u.png)

#### Closeup of mesh

![Mesh](https://github.com/sscollis/npot/blob/master/test/ogrid/mesh.png)

S. Scott Collis\
