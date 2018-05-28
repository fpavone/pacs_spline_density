#!/bin/bash

rm *.so
rm *.o
export OMP_NUM_THREADS=4
R CMD SHLIB *.cpp -o libpacs.so
Rscript ../R/gogo.r
