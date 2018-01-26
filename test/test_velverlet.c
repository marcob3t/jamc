// reference output of this test should be test_velverlet.dat

#include "ljmd.h"



/* append data to output - simple version. */
void simple_output(mdsys_t *sys, FILE *ofp) {

    int i;

    for (i=0; i<sys->natoms; ++i) {
        fprintf(ofp, "%20.8f %20.8f %20.8f\n%20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i], sys->rz[i], sys->vx[i], sys->vy[i], sys->vz[i]);
    }

}



/* main */
int main(int argc, char **argv) {
  
    int i;
    char line[BLEN], restfile[BLEN];
    /* set output file (hard-coded) */
    char* outfile = "test_velverlet.dat";
    /* set force file from command line argument */
    char* forcefile = argv[1];

    FILE *ifp, *ffp, *ofp;
    mdsys_t sys;

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
    if(get_a_line(stdin,line)) return 1; // disregarding trajfile
    if(get_a_line(stdin,line)) return 1; // disregarding ergfile
    if(get_a_line(stdin,line)) return 1;
    sys.nsteps=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.dt=atof(line);
    if(get_a_line(stdin,line)) return 1;
    // nprint=atoi(line);

    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*sizeof(double));

    /* read restart */
    ifp=fopen(restfile,"r");
    if(ifp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(ifp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        for (i=0; i<sys.natoms; ++i) {
            fscanf(ifp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
        }
        fclose(ifp);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    /* initialize forces.*/
    // force(&sys);
    ffp=fopen(forcefile,"r");
    if(ffp) {
        int buf;
        for (i=0; i<sys.natoms; ++i) {
            fscanf(ffp,"%d%lf%lf%lf\n", &buf, sys.fx+i, sys.fy+i, sys.fz+i);
        }
        fclose(ffp);
    } else {
        perror("cannot read force file");
        return 3;
    }    

    /* initialize other members (not needed actually) */
    sys.nfi=0;
    sys.ekin=1.;
    sys.epot=1.;
    sys.temp=1.;

    /* open output file */
    ofp=fopen(outfile,"w");

    /* output initial state */
    simple_output(&sys, ofp);

    /* propagate system with velverlet */
    velverlet(&sys);

    /* output updated state */
    simple_output(&sys, ofp);

    /* clean up: close files, free memory */
    fclose(ofp);

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
