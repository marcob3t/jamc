# Script to plot scaling

## omp agg --------------------------------
set title 'LJMD with OpenMP - argon 2916 (aggressive truncation) - Ulysses full node'
set xlabel 'omp threads'
set ylabel 'speedup'

plot 'omp_agg.txt' u 3:2 w lp pt 7 ps 3 title ''

set terminal png
set output 'omp_agg.png'
rep
unset terminal

## omp aggNewtonAtomic --------------------------
set title 'LJMD with OpenMP - argon 2916 (Newton + Atomic) - Ulysses full node'
set xlabel 'omp threads'
set ylabel 'speedup'

plot 'omp_aggAtomic.txt' u 3:2 w lp pt 7 title ''

set terminal png
set output 'omp_aggAtomic.png'
rep
unset terminal

## omp aggNewtonReplicatedMemLINLOOPSIndexArray --------------------------
set title 'LJMD with OpenMP - argon 2916 (Newton + Linearized Loops + IndexArray) - Ulysses full node'
set xlabel 'omp threads'
set ylabel 'speedup'

plot 'omp_IndexArray.txt' u 3:2 w lp pt 7 title ''

set terminal png
set output 'omp_IndexArray.png'
rep
unset terminal

## omp aggNewtonReplicatedMemLinearLoopInPlaceIndex --------------------------
set title 'LJMD with OpenMP - argon 2916 (aggressive truncation + Atomic) - Ulysses full node'
set xlabel 'omp threads'
set ylabel 'speedup'

plot 'omp_InPlace.txt' u 3:2 w lp pt 7 title ''

set terminal png
set output 'omp_InPlace.png'
rep
unset terminal

## mpi ------------------------------------
set title 'LJMD with MPI - argon 2916'
set xlabel 'MPI processes'
set ylabel 'speedup'
set grid ytics
set grid xtics

plot 'mpi.txt' u 3:2 w lp pt 7 ps 2 title ''

set terminal png
set output 'mpi.png'
rep
unset terminal