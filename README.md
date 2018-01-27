# jamc report
P1.6 group assignment (GROUP 1): Lennard-Jones Molecular Dynamics

## Collaborators:

             - Jiaxin Wang         ---> GitHub: ricphy 
             
             - Marco Bettiol       ---> GitHub: marcob3t
             
             - Carolina Bonivento  ---> GitHub: carolinabonivento
             
             - Alejandra Foggia    ---> GitHub: amfoggia

## Contributions:

#### Jiaxin Wang:
* unit test for kinetic energy [this file](./test/test_ekin.c)
* calculation optimization (report below in "OPT log" section)
* cell list module in cpp [this file](./src/cell.cc)
* unit test for cell list [this file](./test/test_cell.c)

#### Marco Bettiol:
* unit test for integration [this file](./test/test_velverlet_1.c) [this file](./test/test_velverlet_2.c)
* multi-threading

#### Carolina Bonivento:
* unit test for input/output
* python interface

#### Alejandra Foggia:
* unit test for force calculation
* MPI


## OPT log:
(serial code run with i7 3.1GHz, parallel code run with Ulysses)

* profiling original code


* tiny modifications, testing with argon_2916

force function single call timing:
|timing (ms)|feature|
|-------------------|-------|
|120-130|original|
|115-120|+ 1D pre truncate|
|110-120|+ use r-square|
|95-105 |+ inline pbc|
|93-100 |+ replace i==j|
|93-100 |+ rsq_inv|
(in old gcc std, we can only use "static inline")

* aggressive truncation, testing with argon_2916

velverlet function accumulated timing:
|time (s)|feature|
|-----------------|-------|
|~95 |original|
|~75 |+ tearing down spherical truncate into 1D aggressive truncate|

1D aggressive truncate means:

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
|2477.652|76481.576|+agg. truncate|
|1384.265|42099.080|+Newton|
|1198.499|37016.217|+Newton +agg. truncate|

* profiling with above modifications



* OMP and MPI performances

force function accumulated timing:
|argon_108 (ms)|argon_2916 (ms)|feature|
|||omp +agg. truncate|
|||omp +Newton +agg. truncate|
|||mpi +Newton +agg. truncate|



* applying cell-list

force function accumulated timing:
|argon_108 (ms)|argon_2916 (ms)|feature|
|--------------|---------------|-------|
|              |               |serial cell-list|

* profiling with cell list



