#!/bin/bash
#
# Build a fresh version
#
# Revised:  8/15/22
# Author:   S.Scott Collis 
#
set -e
cd src
\ln -fs gcc.mak Makefile
make clean && make $@
exit $? 
