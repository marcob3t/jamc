#!/bin/bash

make
cd examples
../ljmd-serial.x < argon_2916.inp
gprof ../ljmd-serial.x > serial.txt
cd ..
make clean
