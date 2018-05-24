#!/bin/bash

rm *.so
rm *.o
export MY_LOCAL_R_PATH=$(R RHOME)
export OMP_NUM_THREADS=4
R CMD SHLIB *.cpp -o libpacs.so
Rscript ../R/gogo.r
