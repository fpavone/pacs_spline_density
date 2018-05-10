#!/bin/bash

#rm *.so
#rm *.o
R CMD SHLIB *.cpp -o libpacs.so
R
