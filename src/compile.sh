#!/bin/bash

rm *.so
rm *.o
export MY_LOCAL_R_PATH=$(R RHOME)
R CMD SHLIB *.cpp -o libpacs.so
Rscript ../R/gogo.r
