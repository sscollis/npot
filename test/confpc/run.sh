#!/bin/bash
LNS3D_DIR=../../../lns3d
NPOT=../../src/npot
$LNS3D_DIR/mesh/confpc < confpc.inp > confpc.log
$NPOT -ms < npot.inp | tee npot.log
exit 0
