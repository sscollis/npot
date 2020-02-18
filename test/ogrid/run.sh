#!/bin/bash
LNS3DDIR=../../../lns3d
NPOT=../../src/npot
$LNS3DDIR/mesh/ogrid < ogrid.inp > ogrid.log
$LNS3DDIR/pre/src_ij/pre < pre.inp > pre.log
mv circ.xyz grid.dat
$NPOT < npot.inp | tee npot.log
exit 0
