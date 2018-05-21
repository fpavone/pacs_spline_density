#!/bin/bash

rm *.so
rm *.o
R CMD SHLIB *.cpp -o libpacs.so
Rscript ../R/gogo2.r
