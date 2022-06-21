#!/bin/bash
LNS3D_DIR=../../../lns3d
NPOT=../../src/npot
$LNS3D_DIR/mesh/confpc $2 < confpc.inp.$1 > confpc.log
$NPOT < npot.inp | tee npot.log
$NPOT < restart-1.inp | tee npot.log
$NPOT < restart-2.inp | tee npot.log
$NPOT < restart-3.inp | tee npot.log
$NPOT < restart-4.inp | tee npot.log
$NPOT < restart-5.inp | tee npot.log
$NPOT < restart-4.inp | tee npot.log
$NPOT < restart-3.inp | tee npot.log
$LNS3D_DIR/util/npost lns.dat
\cp wall.dat wall.dat.$1
exit 0
