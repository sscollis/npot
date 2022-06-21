## Conformal mesh for a parabolic cylinder

Mach 0.4 potential flow past a parabolic cylinder

### Steps to run

1. Check the paths in `run.sh`
2. Note that you need the `confpc` mesh generator
3. ~~Also note that `confpc` uses the old `ji` ordering so
   that you need to use the `-ms` option when running `npot`~~
4. NB:  I have updated `confpc` so that it outputs IJ ordering so `-ms` is
   no longer needed.
4. `./run.sh 0`  where `0` indicates the case number to be run which amount
   to different versions of the input file to `confpc`.

Note that to enhance convergence, this script restarts with
several different values of `sigma` 

### Visualization

To use paraview you may want to do:

    ln -s grid.dat grid.xyz

Wall values (like $C_p$ and $\partial p/\partial s$ are available in 
`wall.dat.0` (replace `0` with other case numbers as needed).

You can plot the wall pressure gradient using `gnuplot` with
```bash
gnuplot
load "dpds.com" 
```

### Sample Results

#### Contours of streamwise velocity

![Streamwise velocity](https://github.com/sscollis/npot/blob/master/test/confpc/u.png)

#### Computational Mesh

![Mesh](https://github.com/sscollis/npot/blob/master/test/confpc/mesh.png)

#### Close-up of Computational Mesh

![Mesh](https://github.com/sscollis/npot/blob/master/test/confpc/mesh-cu.png)

S. Scott Collis \
