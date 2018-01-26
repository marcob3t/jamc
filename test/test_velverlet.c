// reference output of this test should be test_velverlet.dat

#include "ljmd.h"



/* main */
int main(int argc, char **argv) {
  
    mdsys_t sys;
    int max_test_sample = 4; // max number of atoms used for testing
    int i;

    /* allocate memory */
    sys.rx=(double *)malloc(max_test_sample * sizeof(double));
    sys.ry=(double *)malloc(max_test_sample * sizeof(double));
    sys.rz=(double *)malloc(max_test_sample * sizeof(double));
    sys.vx=(double *)malloc(max_test_sample * sizeof(double));
    sys.vy=(double *)malloc(max_test_sample * sizeof(double));
    sys.vz=(double *)malloc(max_test_sample * sizeof(double));
    sys.fx=(double *)malloc(max_test_sample * sizeof(double));
    sys.fy=(double *)malloc(max_test_sample * sizeof(double));
    sys.fz=(double *)malloc(max_test_sample * sizeof(double));



    /* TEST # 1: 3 particles far away form each other */

    /* initialize system */
    sys.natoms = 3;
    sys.mass = 39.948;
    sys.epsilon = 0.2379;
    sys.sigma = 3.405;
    sys.rcut = 8.5;
    sys.box = 17.1580;
    sys.nsteps = 1;
    sys.dt = 5.0;

    /* initialize other members (not needed actually) */
    sys.nfi=0;
    sys.ekin=1.;
    sys.epot=1.;
    sys.temp=1.;

    /* initialize system positions */
    sys.rx[0] = 0.;
    sys.rx[0] = 0.;
    sys.rx[0] = 0.;

    sys.rx[1] = 8.54;
    sys.rx[1] = 0.;
    sys.rx[1] = 0.;

    sys.rx[2] = 0.;
    sys.rx[2] = 8.54;
    sys.rx[2] = 0.;

    /* initialize system velocities */

    /* initialize system forces */

    /* propagate system with velverlet and check */
    // velverlet(&sys);

    /* compute and print forces (TO BE REMOVED) */
    force(&sys);
    for (i=0; i < sys.natoms; ++i) {
      printf("%20.8f%20.8f%20.8f\n", sys.fx[i], sys.fy[i], sys.fz[i]);
    }

    /* ************************************************* */



    /* TEST # 2: 3 particles, 2 close to each other */

    /* initialize system */

    /* initialize system positions */
    sys.rx[0] = 0.;
    sys.rx[0] = 0.;
    sys.rx[0] = 0.;

    sys.rx[1] = 1.;
    sys.rx[1] = 0.;
    sys.rx[1] = 0.;

    sys.rx[2] = 0.;
    sys.rx[2] = 8.54;
    sys.rx[2] = 0.;

    /* initialize system velocities */

    /* initialize system forces */

    /* propagate system with velverlet and check */
    // velverlet(&sys);

    /* compute and print forces (TO BE REMOVED) */
    force(&sys);
    for (i=0; i < sys.natoms; ++i) {
      printf("%20.8f%20.8f%20.8f\n", sys.fx[i], sys.fy[i], sys.fz[i]);
    }

    /* ************************************************* */



    /* TEST # 3: 3 particles, 3 close to each other */

    /* initialize system */

    /* initialize system positions */
    sys.rx[0] = 0.;
    sys.rx[0] = 0.;
    sys.rx[0] = 0.;

    sys.rx[1] = 1.;
    sys.rx[1] = 0.;
    sys.rx[1] = 0.;

    sys.rx[2] = -1.;
    sys.rx[2] = 0.;
    sys.rx[2] = 0.;

    /* initialize system velocities */

    /* initialize system forces */

    /* propagate system with velverlet and check */
    // velverlet(&sys);

    /* compute and print forces (TO BE REMOVED) */
    force(&sys);
    for (i=0; i < sys.natoms; ++i) {
      printf("%20.8f%20.8f%20.8f\n", sys.fx[i], sys.fy[i], sys.fz[i]);
    }

    /* ************************************************* */



    /* TEST # 4: 4 particles, 3 close to each other */

    /* initialize system */

    sys.natoms = 4;

    /* initialize system positions */
    sys.rx[0] = 0.;
    sys.rx[0] = 0.;
    sys.rx[0] = 0.;

    sys.rx[1] = 1.;
    sys.rx[1] = 0.;
    sys.rx[1] = 0.;

    sys.rx[2] = -1.;
    sys.rx[2] = 0.;
    sys.rx[2] = 0.;

    sys.rx[3] = 0.;
    sys.rx[3] = 8.54;
    sys.rx[3] = 0.;

    /* initialize system velocities */

    /* initialize system forces */

    /* propagate system with velverlet and check */
    // velverlet(&sys);

    /* compute and print forces (TO BE REMOVED) */
    force(&sys);
    for (i=0; i < sys.natoms; ++i) {
      printf("%20.8f%20.8f%20.8f\n", sys.fx[i], sys.fy[i], sys.fz[i]);
    }

    /* ************************************************* */



    /* clean up: free memory */
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
