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
* cell list module: [cell](./src/cell.cc)
* unit test for cell list: [test_cell](./test/test_cell.c)

#### Marco Bettiol:
* unit test for integration: [test_velverlet_1](./test/test_velverlet_1.c) [test_velverlet_2](./test/test_velverlet_2.c)
* multi-threading: [force](./src/force.c)

#### Carolina Bonivento:
* unit test for input/output: [test_in](./test/test_in.c) [test_out](./test/test_out.c)
* python interface:

#### Alejandra Foggia:
* unit test for force calculation: [test_force](./test/test_force.c)
* MPI: 
* apply cell list: 


## OPT log:
(related source files in [here](./opt/) )

* profiling original code

```
Each sample counts as 0.01 seconds.
%   cumulative   self              self     total
time   seconds   seconds    calls  us/call  us/call  name
68.75      4.85     4.85    10001   485.34   684.99  force(_mdsys*)
28.00      6.83     1.98 346714668     0.01     0.01  pbc(double, double)
0.57      6.87     0.04    10001     4.00     4.00  ekin(_mdsys*)
0.57      6.91     0.04    10000     4.00     4.00  velverlet_1(_mdsys*)
0.28      6.93     0.02    30006     0.67     0.67  azzero(double*, int)
0.07      6.94     0.01    10000     0.50     0.50  velverlet_2(_mdsys*)
0.00      6.94     0.00      101     0.00     0.00  output(_mdsys*, _IO_FILE*, _IO_FILE*)
0.00      6.94     0.00       12     0.00     0.00  get_a_line(_IO_FILE*, char*)
```

* serial code optimization without cell list

aggressive truncation is better than spherical truncation:

```
/* get distance between particle i and j */
rx=pbc(sys->rx[i] - sys->rx[j], boxby2);
rsq = rx*rx;
if(rsq>rcutsq) continue; 
ry=pbc(sys->ry[i] - sys->ry[j], boxby2);
rsq += ry*ry;
if(rsq>rcutsq) continue; 
rz=pbc(sys->rz[i] - sys->rz[j], boxby2);
rsq += rz*rz;
if(rsq>rcutsq) continue; 
```
which is a temporary solution before implementing cell-list

profiling with using aggressive truncatioin

```
Each sample counts as 0.01 seconds.
%   cumulative   self              self     total
time   seconds   seconds    calls  us/call  us/call  name
73.49      5.50     5.50    10001   550.36   737.48  force(_mdsys*)
24.99      7.38     1.87 319175220     0.01     0.01  pbc(double, double)
1.00      7.45     0.08    10001     7.50     7.50  ekin(_mdsys*)
0.33      7.48     0.03                             stamp()
0.13      7.49     0.01    10000     1.00     1.00  velverlet_1(_mdsys*)
0.13      7.50     0.01    10000     1.00     1.00  velverlet_2(_mdsys*)
0.00      7.50     0.00    30006     0.00     0.00  azzero(double*, int)
0.00      7.50     0.00      101     0.00     0.00  output(_mdsys*, _IO_FILE*, _IO_FILE*)
0.00      7.50     0.00       12     0.00     0.00  get_a_line(_IO_FILE*, char*)
```

Newton 3rd law:

```
for(i=0; i < (sys->natoms); ++i) {
    for(j=i+1; j < (sys->natoms); ++j) {
        // get distance between particle i and j
        rx=pbc(sys->rx[i] - sys->rx[j], boxby2);
        rsq = rx*rx;
        if(rsq>rcutsq) continue;
        ry=pbc(sys->ry[i] - sys->ry[j], boxby2);
        rsq += ry*ry;
        if(rsq>rcutsq) continue;
        rz=pbc(sys->rz[i] - sys->rz[j], boxby2);
        rsq += rz*rz;
        if(rsq>rcutsq) continue;

        rsq_inv = 1.0/rsq;
        r6 = rsq_inv*rsq_inv*rsq_inv;
        ffac = (48*c12*r6-24*c6)*r6*rsq_inv;
        sys->epot += 4*r6*(c12*r6-c6);
        sys->fx[i] += rx*ffac; sys->fx[j] -= rx*ffac;
        sys->fy[i] += ry*ffac; sys->fy[j] -= ry*ffac;
        sys->fz[i] += rz*ffac; sys->fz[j] -= rz*ffac;
    }
}
```

profiling with Newton 3rd law and aggressive truncation

```
Each sample counts as 0.01 seconds.
%   cumulative   self              self     total
time   seconds   seconds    calls  us/call  us/call  name
73.59      3.00     3.00    10001   300.21   401.78  force(_mdsys*)
24.65      4.01     1.01 159587610     0.01     0.01  pbc(double, double)
1.23      4.06     0.05    10001     5.00     5.00  ekin(_mdsys*)
0.25      4.07     0.01    30006     0.33     0.33  azzero(double*, int)
0.25      4.08     0.01    10000     1.00     1.00  velverlet_1(_mdsys*)
0.12      4.08     0.01                             stamp()
0.00      4.08     0.00    10000     0.00     0.00  velverlet_2(_mdsys*)
0.00      4.08     0.00      101     0.00     0.00  output(_mdsys*, _IO_FILE*, _IO_FILE*)
0.00      4.08     0.00       12     0.00     0.00  get_a_line(_IO_FILE*, char*)
```

force function accumulated timing:

|argon_108 (ms)|argon_2916 (ms)|feature|
|--------------|---------------|-------|
|2999.177|99288.438|original|
|2477.652|76481.576|+agg. trunc|
|1384.265|42099.080|+Newton|
|1198.499|37016.217|+Newton +agg. trunc|

comments: aggressive truncation is benefitial, so we keep it as standard *until* cell list


* OMP and MPI performances

force function accumulated timing:

|argon_2916 (ms)|speedup|MPI_procs/OMP_threads|feature|
|--------------|--------|---------------------|-------|
|73217.0729|1.00|1/1|+agg. trunc|
|37895.96006|1.93|1/2|+agg. trunc|
|26387.43521|2.77|1/3|+agg. trunc|
|20478.64648|3.58|1/4|+agg. trunc|
|16802.16675|4.36|1/5|+agg. trunc|
|13997.1981|5.23|1/6|+agg. trunc|
|12071.90107|6.07|1/7|+agg. trunc|
|10585.75518|6.92|1/8|+agg. trunc|
|9361.658057|7.82|1/9|+agg. trunc|
|8516.886572|8.60|1/10|+agg. trunc|

|argon_2916 (ms)|speedup|MPI_procs/OMP_threads|feature|
|--------------|--------|---------------------|-------|
|||1/10|+Newton +agg.|
|||2/10|+Newton +agg.|
|||3/10|+Newton +agg.|
|||4/10|+Newton +agg.|


* applying cell-list

profiling with cell list

force function accumulated timing:

|argon_108 (ms)|argon_2916 (ms)|feature|
|--------------|---------------|-------|
|              |               |serial cell-list|





