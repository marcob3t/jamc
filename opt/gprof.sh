#!/bin/bash

make serial
cd examples
rm gmon.out
../ljmd-serial.x < argon_108.inp
gprof ../ljmd-serial.x > serial.txt
cd ..
make clean
