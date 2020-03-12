#!/bin/bash
LNS3D_DIR=../../../lns3d
NPOT=../../src/npot
env GFORTRAN_CONVERT_UNIT='swap' $NPOT -wc < npot.inp | tee npot.log
env GFORTRAN_CONVERT_UNIT='swap' $NPOT -wc < res1.inp | tee npot.log
env GFORTRAN_CONVERT_UNIT='swap' $NPOT -wc < res2.inp | tee npot.log
#env GFORTRAN_CONVERT_UNIT='swap' $NPOT -wc < res1.inp | tee npot.log
#env GFORTRAN_CONVERT_UNIT='swap' $NPOT -wc < res2.inp | tee npot.log
exit 0
