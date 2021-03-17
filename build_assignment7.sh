#!/bin/bash

# build object files
gcc -I ./random.h -c ./random.c

# build assignment files
gcc -I include/ ./assignment7_5.c -o ./assignment7_5.exe random.o
gcc -I include/ ./final.c -o ./final.exe random.o

# build data files
./assignment7_5.exe > ./data7/7_5.dat 
./final.exe > ./data7/gauss.dat
