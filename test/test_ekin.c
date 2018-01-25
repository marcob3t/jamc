// reference output of this test should be test_ekin.dat
#include "ljmd.h"

/* test main */
int main(int argc, char **argv){
    int i;
    char restfile[BLEN], line[BLEN];
    FILE *fp,*fn;
    mdsys_t sys;
    
    /* read input file */
    // number of atoms
    if(get_a_line(stdin,line)) return 1;
    sys.natoms=atoi(line);
    // mass of atoms
    if(get_a_line(stdin,line)) return 1;
    sys.mass=atof(line);
    // keep these trash just to read restfile name
    if(get_a_line(stdin,line)) return 1;
    if(get_a_line(stdin,line)) return 1;
    if(get_a_line(stdin,line)) return 1;
    if(get_a_line(stdin,line)) return 1;
    // initial setting file name
    if(get_a_line(stdin,restfile)) return 1;
    // alloc velocity field
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    // catching un-wanted information
    double trash[3];
    
    /* read restart */
    fp=fopen(restfile,"r");
    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",trash, trash+1, trash+2);
        }
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
        }
        fclose(fp);
    } else {
        perror("cannot read restart file");
        return 3;
    }
    
    ekin(&sys);
    fn = fopen("../reference/test_ekin.dat","w");
    // print out to reference file
    fprintf(fn,"% 20.8f\n",sys.ekin);
    
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    return 0;
}
