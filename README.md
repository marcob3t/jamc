# jamc report
P1.6 group assignment (GROUP 1): Lennard-Jones Molecular Dynamics

## Collaborators:

             - Jiaxin Wang         ---> GitHub: ricphy 
             
             - Marco Bettiol       ---> GitHub: marcob3t
             
             - Carolina Bonivento  ---> GitHub: carolinabonivento
             
             - Alejandra Foggia    ---> GitHub: amfoggia

## Contributions:

#### Jiaxin Wang:
* unit test for kinetic energy: [test_ekin](./test/test_ekin.c)
* calculation optimization: report below in "OPT log" section
* cell list module in cpp: [cell](./src/cell.cc)
* unit test for cell list: [test_cell](./test/test_cell.c)

#### Marco Bettiol:
* unit test for integration: [test_velverlet_1](./test/test_velverlet_1.c) [test_velverlet_2](./test/test_velverlet_2.c)
* multi-threading: [force](./src/force.c)

#### Carolina Bonivento:
* unit test for input/output: [test_in](./test/test_in.c) [test_out](./test/test_out.c)
* python interface:

#### Alejandra Foggia:
* unit test for force calculation: [test_force](./test/test_force.c)
* MPI


## OPT log:
(serial code run with i7 3.1GHz, parallel code run with Ulysses)

* profiling original code


* aggressive truncation, testing with argon_2916

velverlet function accumulated timing:

|time (s)|feature|
|-----------------|-------|
|~95 |original|
|~75 |+aggressive truncate|

aggressive truncation means:

```
/* get distance between particle i and j */
rx=pbc(sys->rx[i] - sys->rx[j], boxby2);
rsq = rx*rx;
if(rsq>rcutsq) continue; // 1D pre truncate
ry=pbc(sys->ry[i] - sys->ry[j], boxby2);
rsq += ry*ry;
if(rsq>rcutsq) continue; // 1D pre truncate
rz=pbc(sys->rz[i] - sys->rz[j], boxby2);
rsq += rz*rz;
if(rsq>rcutsq) continue; // 1D pre truncate
```
which is a temporary solution before implementing cell-list

* in further development we separate velverlet into velverlet_1,force,velverlet_2, in which we found force function dominates computing time, thus in following optimization we focus on force function

force function accumulated timing:

|argon_108 (ms)|argon_2916 (ms)|feature|
|--------------|---------------|-------|
|2999.177|99288.438|original|
|2477.652|76481.576|+agg. trunc|
|1384.265|42099.080|+Newton|
|1198.499|37016.217|+Newton +agg. trunc|

comments: aggressive truncation is always benefitial, so we keep this algorithm as standard

* profiling with above modifications



* OMP and MPI performances

force function accumulated timing:

|argon_2916 (ms)|speedup|MPI_procs/OMP_threads|feature|
|--------------|--------|---------------------|-------|
|73217.0729|1.00|1/1|omp +agg. trunc|
|37895.96006|1.93|1/2|omp +agg. trunc|
|26387.43521|2.77|1/3|omp +agg. trunc|
|20478.64648|3.58|1/4|omp +agg. trunc|
|16802.16675|4.36|1/5|omp +agg. trunc|
|13997.1981|5.23|1/6|omp +agg. trunc|
|12071.90107|6.07|1/7|omp +agg. trunc|
|10585.75518|6.92|1/8|omp +agg. trunc|
|9361.658057|7.82|1/9|omp +agg. trunc|
|8516.886572|8.60|1/10|omp +agg. trunc|

|argon_2916 (ms)|speedup|MPI_procs/OMP_threads|feature|
|--------------|--------|---------------------|-------|
|||1/2|omp +agg. +Newton|
|||1/3|omp +agg. +Newton|
|||1/4|omp +agg. +Newton|
|||1/5|omp +agg. +Newton|
|||1/6|omp +agg. +Newton|
|||1/7|omp +agg. +Newton|
|||1/8|omp +agg. +Newton|
|||1/9|omp +agg. +Newton|
|||1/10|omp +agg. +Newton|



* applying cell-list

force function accumulated timing:

|argon_108 (ms)|argon_2916 (ms)|feature|
|--------------|---------------|-------|
|              |               |serial cell-list|

* profiling with cell list



