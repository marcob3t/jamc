#include "ljmd.h"
#include "cell.h"
#ifdef _OPENMP
#include "omp.h"
#endif

/* main */
int main(int argc, char **argv) 
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    mdsys_t sys; FILE *fp;
#ifdef _OPENMP
    int nthds = omp_get_max_threads(); // openmp
#else
    int nthds = 1;
#endif
#ifdef TIMING
    double timer = 0; // for recording time
#endif

    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    sys.natoms=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.redmass = atof(line)*mvsq2e; // reduced mass
    if(get_a_line(stdin,line)) return 1;
    sys.epsilon=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.sigma=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.rcut=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.box=atof(line);
    if(get_a_line(stdin,restfile)) return 1;
    if(get_a_line(stdin,trajfile)) return 1;
    if(get_a_line(stdin,ergfile)) return 1;
    if(get_a_line(stdin,line)) return 1;
    sys.nsteps=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.dt=atof(line);
    if(get_a_line(stdin,line)) return 1;
    nprint=atoi(line);
    
    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*nthds*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*nthds*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*nthds*sizeof(double));

    /* read restart */
    fp=fopen(restfile,"r");
    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
        }
        fclose(fp);
        azzero(sys.fx, sys.natoms);
        azzero(sys.fy, sys.natoms);
        azzero(sys.fz, sys.natoms);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys);
    ekin(&sys);
    
#ifndef TIMING
    FILE *traj,*erg;
    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");
    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);
#endif
    

    /**************************************************/
    /* main MD loop */
    
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
        
        /* propagate system and recompute energies */
#ifdef TIMING
        velverlet_1(&sys);
        timer -= stamp();
        force(&sys);
        timer += stamp();
        velverlet_2(&sys);
#else
        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);
        velverlet_1(&sys);
        force(&sys);
        velverlet_2(&sys);
#endif
        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
#ifdef TIMING
    printf("\n timing \n\t %f ms \n\n",timer);
#else
    printf("Simulation Done.\n");
    fclose(erg);
    fclose(traj);
#endif

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);

    return 0;
}
