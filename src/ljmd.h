/* 
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#ifndef JLMD_H
#define JLMD_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information 
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
};

typedef struct _mdsys mdsys_t;

/* helper function: read a line and then return
   the first string with whitespace stripped off */
int get_a_line(FILE *fp, char *buf);

/* helper function: zero out an array */
void azzero(double *d, const int n);

/* helper function: apply minimum image convention */
double pbc(double x, const double boxby2);

/* compute kinetic energy */
void ekin(mdsys_t *sys);

/* compute forces */
void force(mdsys_t *sys);

/* velocity verlet */
void velverlet(mdsys_t *sys);

/* append data to output. */
void output(mdsys_t *sys, FILE *erg, FILE *traj);

#endif