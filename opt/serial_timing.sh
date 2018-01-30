#! /bin/bash

make serial
cd examples
for j in {1..5}
do
../ljmd-serial.x < argon_108.inp
done
cd ..
make clean
