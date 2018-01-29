#include "ljmd.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/* WE INCLUDE ALL HISTORICAL VERSIONS OF THIS FUNCTION */

/* original */
/*
 void force(mdsys_t *sys)
 {
 double r,ffac;
 double rx,ry,rz;
 int i,j;
 
 // zero energy and forces
 sys->epot=0.0;
 azzero(sys->fx,sys->natoms);
 azzero(sys->fy,sys->natoms);
 azzero(sys->fz,sys->natoms);
 
 for(i=0; i < (sys->natoms); ++i) {
 for(j=0; j < (sys->natoms); ++j) {
 
 // particles have no interactions with themselves
 if (i==j) continue;
 
 // get distance between particle i and j
 rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
 ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
 rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
 r = sqrt(rx*rx + ry*ry + rz*rz);
 
 // compute force and energy if within cutoff
 if (r < sys->rcut) {
 ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
 +6*pow(sys->sigma/r,6.0)/r);
 
 sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
 -pow(sys->sigma/r,6.0));
 
 sys->fx[i] += rx/r*ffac;
 sys->fy[i] += ry/r*ffac;
 sys->fz[i] += rz/r*ffac;
 }
 }
 }
 }
 */

/* aggressive truncate */
/*
void force(mdsys_t *sys)
{
    double rsq,rsq_inv,r6,ffac;
    double rx,ry,rz;
    int i,j;
    double epot=0.0; // needed for reduction with openmp
    
    // zero energy and forces
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    double boxby2 = 0.5*sys->box;// pre-calculate
    double rcutsq = sys->rcut*sys->rcut;// pre-calculate, take square
    double c6 = sys->epsilon*pow(sys->sigma,6);
    double c12 = sys->epsilon*pow(sys->sigma,12);
    
    for(i=0; i < (sys->natoms); ++i) {
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
            sys->epot += 2*r6*(c12*r6-c6);
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
            sys->epot += 2*r6*(c12*r6-c6);
            sys->fx[i] += rx*ffac;
            sys->fy[i] += ry*ffac;
            sys->fz[i] += rz*ffac;
            
        }
    }
}
*/


/* Newton */
/*
void force(mdsys_t *sys)
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;
    
    // zero energy and forces
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    
    for(i=0; i < (sys->natoms); ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {
            
            // get distance between particle i and j
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r = sqrt(rx*rx + ry*ry + rz*rz);
            
            // compute force and energy if within cutoff
            if (r < sys->rcut) {
                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                          +6*pow(sys->sigma/r,6.0)/r);
                // a factor of two is necessary
                sys->epot += 2*0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                                   -pow(sys->sigma/r,6.0));
                
                sys->fx[i] += rx/r*ffac; sys->fx[j] -= rx/r*ffac;
                sys->fy[i] += ry/r*ffac; sys->fy[j] -= ry/r*ffac;
                sys->fz[i] += rz/r*ffac; sys->fz[j] -= rz/r*ffac;
            }
        }
    }
}
*/

/* aggressive Newton */
/*
void force(mdsys_t *sys)
{
    double rsq,rsq_inv,r6,ffac;
    double rx,ry,rz;
    int i,j;
    
    // zero energy and forces
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    double boxby2 = 0.5*sys->box;// pre-calculate
    double rcutsq = sys->rcut*sys->rcut;// pre-calculate, take square
    double c6 = sys->epsilon*pow(sys->sigma,6);
    double c12 = sys->epsilon*pow(sys->sigma,12);
    
    for(i=0; i < (sys->natoms); ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {
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
            sys->epot += 4*r6*(c12*r6-c6);
            sys->fx[i] += rx*ffac; sys->fx[j] -= rx*ffac;
            sys->fy[i] += ry*ffac; sys->fy[j] -= ry*ffac;
            sys->fz[i] += rz*ffac; sys->fz[j] -= rz*ffac;
        }
    }
}
*/


/* omp aggressive */
/*
void force(mdsys_t *sys)
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
        for(i=0; i < (sys->natoms); ++i) {
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
}
*/


/* newton omp aggressive with atomic */
/*
void force(mdsys_t *sys)
{
    double rsq,rsq_inv,r6,ffac;
    double rx,ry,rz;
    int i,j;
    double epot=0.0; // needed for reduction with openmp

    // zero energy and forces
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    double boxby2 = 0.5*sys->box;// pre-calculate
    double rcutsq = sys->rcut*sys->rcut;// pre-calculate, take square
    double c6 = sys->epsilon*pow(sys->sigma,6);
    double c12 = sys->epsilon*pow(sys->sigma,12);

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
        for(i=0; i < (sys->natoms); ++i) {
            for(j=i+1; j < (sys->natoms); ++j) {
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
                epot += 4*r6*(c12*r6-c6);
#ifdef _OPENMP
#pragma omp atomic update
#endif
                sys->fx[i] += rx*ffac;
#ifdef _OPENMP
#pragma omp atomic update
#endif
                sys->fx[j] -= rx*ffac;
#ifdef _OPENMP
#pragma omp atomic update
#endif
                sys->fy[i] += ry*ffac;
#ifdef _OPENMP
#pragma omp atomic update
#endif
                sys->fy[j] -= ry*ffac;
#ifdef _OPENMP
#pragma omp atomic update
#endif
                sys->fz[i] += rz*ffac;
#ifdef _OPENMP
#pragma omp atomic update
#endif
                sys->fz[j] -= rz*ffac;
            }
        }
    }
    sys->epot = epot;
}
*/


/* newton omp aggressive with replicated memory */
void force(mdsys_t *sys)
{
    double rsq,rsq_inv,r6,ffac;
    double rx,ry,rz;
    int n,i,j,offset;
    double * fx, * fy, * fz; // local pointers for openmp
    double epot=0.0; // needed for reduction with openmp
    int niters = (sys->natoms) * (sys->natoms - 1) / 2; // to linearize (and balance) the loop with openmp / mpi
    // TO DO: niters should be divided by the number of processes
    // accordingly to MPI policies
    int * indexes_i, * indexes_j; // to store indexes

    // allocate and compute indexes
    indexes_i=(int *)malloc(niters*sizeof(int));
    indexes_j=(int *)malloc(niters*sizeof(int));
    for (i=0, n=0; i<sys->natoms; ++i) {
      for (j=i+1; j<sys->natoms; ++j, ++n) {
        indexes_i[n] = i;
        indexes_j[n] = j;
      }
    }

    double boxby2 = 0.5*sys->box;// pre-calculate
    double rcutsq = sys->rcut*sys->rcut;// pre-calculate, take square
    double c6 = sys->epsilon*pow(sys->sigma,6);
    double c12 = sys->epsilon*pow(sys->sigma,12);

#ifdef _OPENMP
#pragma omp parallel private(fx,fy,fz,i,j,rx,ry,rz,rsq,rsq_inv,r6,ffac,offset) reduction(+:epot)
#endif
    {

        // who is who
        int thid; // thread id
        int nthds; // number of threads

#ifdef _OPENMP
        thid = omp_get_thread_num();
        nthds = omp_get_num_threads();
#else
        thid = 0;
        nthds = 1;
#endif

        // assign slices of force array to threads
        fx = sys->fx + (thid * sys->natoms);
        fy = sys->fy + (thid * sys->natoms);
        fz = sys->fz + (thid * sys->natoms);

        // zero energy and forces
        azzero(fx,sys->natoms);
        azzero(fy,sys->natoms);
        azzero(fz,sys->natoms);

        // main loop
#ifdef _OPENMP
#ifndef CHUNKSIZE
#define CHUNKSIZE=1
#endif
#pragma omp for schedule(dynamic,CHUNKSIZE)
#endif
        for(n=0; n<niters; ++n) {

            i=indexes_i[n]; // obtain original i index
            j=indexes_j[n]; // obtain original j index

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
            epot += 4*r6*(c12*r6-c6);
            fx[i] += rx*ffac;
            fx[j] -= rx*ffac;
            fy[i] += ry*ffac;
            fy[j] -= ry*ffac;
            fz[i] += rz*ffac;
            fz[j] -= rz*ffac;
        }

        // after a work sharing construct an omp barrier is implied
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (i=0; i<sys->natoms; ++i) {
            for (offset=sys->natoms; offset<nthds*sys->natoms; offset+=sys->natoms) {
              sys->fx[i] += sys->fx[offset + i];
              sys->fy[i] += sys->fy[offset + i];
              sys->fz[i] += sys->fz[offset + i];
            }
        }

    }

    sys->epot = epot;

    // free memory
    free(indexes_i);
    free(indexes_j);

}
