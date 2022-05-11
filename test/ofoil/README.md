## Ofoil

### O-grid for an airfoil

This is Mach 0.3 potential flow over a simple analytical airfoil of the form

    y = t/0.2 * ( 0.2969 * sqrt(x) - 0.1281 * x - 0.3516 * x**2 + &
        0.2843 * x**3 - 0.1015 * x**4 )

Where the thickness, `t=0.1`.  Note that the trailing edge has been rounded
using a radius `0.00255` circle.

### Sample Results

#### Contours of density

![Streamwise velocity](https://github.com/sscollis/npot/blob/master/test/ofoil/rho.png)

#### Computational Mesh

![Mesh](https://github.com/sscollis/npot/blob/master/test/ofoil/mesh.png)

#### Computational Mesh close-up with density contours

![Mesh](https://github.com/sscollis/npot/blob/master/test/ofoil/rho-mesh.png)

S. Scott Collis \
