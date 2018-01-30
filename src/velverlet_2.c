#include "ljmd.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* velocity verlet */
void velverlet_2(mdsys_t *sys)
{
    int i;
    double coef = 0.5*sys->dt / sys->redmass;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(i)
#endif
    /* second part: propagate velocities by another half step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += coef * sys->fx[i];
        sys->vy[i] += coef * sys->fy[i];
        sys->vz[i] += coef * sys->fz[i];
    }
}
