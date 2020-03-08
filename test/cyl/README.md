### Level-set grid for 1/2 circular cylinder

Compressible potential flow at Mach 0.4 over a circular cylinder using a 
1/2 mesh with symmetry conditions.

![Streamwise velocity](https://github.com/sscollis/npot/blob/master/test/cyl/u.png)

![Mesh](https://github.com/sscollis/npot/blob/master/test/cyl/mesh.png)

### Steps to make the run
1. Check the paths in `run.sh`
2. Note that you need the ogrid mesh generator from lns3d
3. `./run.sh`

To use `paraview` you may want to do:

    ln -s grid.dat grid.xyz

S. Scott Collis\
Sun Mar  8 08:10:46 MDT 2020
