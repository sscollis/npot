#!/bin/bash
LNS3D_DIR=../../../lns3d
MESH_DIR=$LNS3D_DIR/mesh
NPOT=../../src/npot
$MESH_DIR/ogrid -a < ogrid.inp > ogrid.log
ln -fs foil.xyz grid.dat
$LNS3D_DIR/pre/src_ij/pre < pre.inp > pre.log
$NPOT < npot.inp | tee npot.log
$NPOT < res1.inp | tee npot.log
$NPOT < res2.inp | tee npot.log
$NPOT < res3.inp | tee npot.log
$NPOT < res4.inp | tee npot.log
#$NPOT < res3.inp | tee npot.log
#$NPOT < res2.inp | tee npot.log
$NPOT < res1.inp | tee npot.log
exit 0
