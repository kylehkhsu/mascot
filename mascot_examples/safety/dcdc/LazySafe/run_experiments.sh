#!/bin/bash

make clean
make

for((numAbs=1;numAbs<=10;numAbs++));
do
    ./dcdc numAbs
done
