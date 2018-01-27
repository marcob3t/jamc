#! /bin/bash

cd ..
make openmp
cd examples
for i in {1..10}
do
export OMP_NUM_THREADS=$i
echo $OMP_NUM_THREADS
for j in {1..5}
do
#../ljmd-openmp.x < argon_108.inp
../ljmd-openmp.x < argon_2916.inp
done
done
cd ..
make clean
