#!/bin/bash
#
# Build a fresh version of npot
#
# Revised:  8/15/22
# Author:   S.Scott Collis 
#
set -e
cd src
\ln -fs gcc.mak Makefile
make clean && make $@
exit $? 
