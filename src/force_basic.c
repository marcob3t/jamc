#include "ljmd.h"

void force_basic(mdsys_t *sys)
{
    double rsq,rsq_inv,r6,ffac;
    double rx,ry,rz;
    int i,j;
    double epot=0.0; // needed for reduction with openmp
    
    // zero forces
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    double boxby2 = 0.5*sys->box;// pre-calculate
    double rcutsq = sys->rcut*sys->rcut;// pre-calculate, take square
    double c6 = sys->epsilon*pow(sys->sigma,6);
    double c12 = sys->epsilon*pow(sys->sigma,12);

    int nprocs, rank, local_niter, lower_bound, upper_bound;
#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rest = sys->natoms%nprocs;
    if (rank < rest) {
        local_niter = sys->natoms/nprocs + 1;
        lower_bound = local_niter * rank;
        upper_bound = (rank + 1) * local_niter;
    }
    else {
        local_niter = sys->natoms/nprocs;
        lower_bound = local_niter * rank + rest;
        upper_bound = (rank + 1) * local_niter + rest;
    }
    // Broadcast the positions to all processes
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    nprocs = 1;
    rank = 0;
    local_niter= sys->natoms;
    lower_bound = 0;
    upper_bound = sys->natoms;
#endif /* USE_MPI */    
    
#ifdef _OPENMP
#pragma omp parallel private(i,j,rx,ry,rz,rsq,rsq_inv,r6,ffac) reduction(+:epot)
#endif
    {
#ifdef _OPENMP
#ifndef CHUNKSIZE
#define CHUNKSIZE=9
#endif
#pragma omp for schedule(dynamic,CHUNKSIZE)
#endif
        for(i=lower_bound; i < upper_bound; ++i) {
            for(j=0; j < i; ++j) {
                // particles have no interactions with themselves
                
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
                epot += 2*r6*(c12*r6-c6);
                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;
                
            }
            // particles have no interactions with themselves
            // rm (i==j) decision with separated loop
            for(j=i+1; j < sys->natoms; ++j) {
                // particles have no interactions with themselves
                
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
                epot += 2*r6*(c12*r6-c6);
                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;
            }
        }
    }
    sys->epot = epot;
    
#ifdef USE_MPI
    // Communicate forces
    if (rank == 0){
        MPI_Reduce(MPI_IN_PLACE, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Reduce(sys->fx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(sys->fy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(sys->fz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    
    // Communicate epot
    if (rank == 0)
        MPI_Reduce(MPI_IN_PLACE, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    else
        MPI_Reduce(&sys->epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif /* USE_MPI */
}
