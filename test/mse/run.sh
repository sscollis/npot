#!/bin/bash
LNS3D_DIR=../../../lns3d
NPOT=../../src/npot
$LNS3D_DIR/mesh/mse < mse.inp > mse.log
$NPOT < npot.inp | tee npot.log
$NPOT < res-1.inp | tee npot.log
$NPOT < res-2.inp | tee npot.log
$NPOT < res-3.inp | tee npot.log
$NPOT < res-1.inp | tee npot.log
exit 0
