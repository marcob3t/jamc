#! /bin/bash

make
cd examples
for i in {1..5}
do
export OMP_NUM_THREADS=$i
echo $OMP_NUM_THREADS
../ljmd-serial.x < argon_108.inp
../ljmd-serial.x < argon_2916.inp
done
cd ..
make clean
