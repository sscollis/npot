#!/bin/bash
LNS3D_DIR=../../../lns3d
NPOT=../../src/npot
if [[ $# -lt 1 ]]; then
  echo "Usage:  run.sh case args" 
  exit 2
else
  case=$1
  shift
fi
echo "Running case $case with args: $@"
#set -x
$LNS3D_DIR/mesh/confpc $@ < confpc.inp.$case > confpc.log
$NPOT < npot.inp.0 | tee npot.log
$NPOT < npot.inp.1 | tee npot.log
$NPOT < npot.inp.2 | tee npot.log
$NPOT < npot.inp.3 | tee npot.log
$NPOT < npot.inp.4 | tee npot.log
$NPOT < npot.inp.5 | tee npot.log
$NPOT < npot.inp.4 | tee npot.log
$NPOT < npot.inp.3 | tee npot.log
$LNS3D_DIR/util/npost lns.dat
\cp wall.dat wall.dat.$case
\cp lns.dat.q lns.q.$case
exit 0
