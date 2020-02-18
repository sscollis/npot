#!/bin/bash
LNS3DDIR=../../../lns3d
NPOTDIR=../..
NPOT=$NPOTDIR/src/npot
$LNS3DDIR/mesh/cyl < cyl.inp > cyl.log
$NPOT < npot.inp | tee npot.log
exit 0
