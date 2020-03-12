### NACA airfoil

This is a lifting airfoil so there must be circulation:

On Intel Mac with gfortran:

    env GFORTRAN_CONVERT_UNIT='swap' ./npot -wc < input.dat

    env GFORTRAN_CONVERT_UNIT='swap' ./npot -wc < input.res

The mesh was made with chimera/hypgen and I have included the input file
`hypgen.inp` for those that are able to access and use that software.
Otherewise, the `grid.dat` and `metric.dat` are provided, although you
can make the `metric.dat` using `lns3d/pre`.

Note that this is a lifting airfoil so that there is a wake cut and the
specification of circulation.  This should actually be computed as part of
the solution process.

### Sample Results

#### Contours of streamwise velocity

![Streamwise velocity](https://github.com/sscollis/npot/blob/master/test/airfoil/u.png)

#### Computational Mesh

![Mesh](https://github.com/sscollis/npot/blob/master/test/airfoil/mesh.png)

#### Computational Mesh close-up

![Mesh](https://github.com/sscollis/npot/blob/master/test/airfoil/mesh-cu.png)

S. Scott Collis\
Thu Mar 12 07:35:17 MDT 2020
