#!/bin/bash

rm *.so
rm *.o
R CMD SHLIB *.cpp -o libpacs.so
g++ -std=c++14 -L. main.o -o main -lpacs
