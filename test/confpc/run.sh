#!/bin/bash
LNS3D_DIR=../../../lns3d
NPOT=../../src/npot
$LNS3D_DIR/mesh/confpc < confpc.inp > confpc.log
$NPOT -ms < npot.inp | tee npot.log
$NPOT -ms < restart.inp | tee npot.log
$NPOT -ms < restart-2.inp | tee npot.log
$NPOT -ms < restart-3.inp | tee npot.log
$NPOT -ms < restart-4.inp | tee npot.log
$NPOT -ms < restart-5.inp | tee npot.log
$NPOT -ms < restart-4.inp | tee npot.log
$NPOT -ms < restart-3.inp | tee npot.log
exit 0
