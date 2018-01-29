#include <sys/time.h>
#include <vector>
#include "ljmd.h"

/* helper function: zero out an array */
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

/* compute kinetic energy */
void ekin(mdsys_t *sys)
{   
    int i;
    
    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*sys->redmass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}

double pbc(double x, const double boxby2) {
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

/* timing */
double stamp(){
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*1e+3 + tv.tv_usec*1e-3;
}

// pair up cells in sys
void pair(mdsys_t *sys) {
    int i,j,k;
    for(i=0;i<sys->cn;++i){
        for(j=0;j<sys->cn;++j){
            for(k=0;k<sys->cn;++k){
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i-1,j-1,k-1));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i-1,j-1,k));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i-1,j-1,k+1));//1 line
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i-1,j,k-1));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i-1,j,k));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i-1,j,k+1));//2 line
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i-1,j+1,k-1));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i-1,j+1,k));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i-1,j+1,k+1));//1 facard
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i,j-1,k-1));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i,j-1,k));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i,j-1,k+1));//1 line
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i,j,k-1));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i,j,k+1));//2 line
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i,j+1,k-1));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i,j+1,k));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i,j+1,k+1));//1 facard
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i+1,j-1,k-1));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i+1,j-1,k));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i+1,j-1,k+1));//1 line
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i+1,j,k-1));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i+1,j,k));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i+1,j,k+1));//2 line
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i+1,j+1,k-1));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i+1,j+1,k));
                
                sys->pair.push_back(index3d(sys,i,j,k));
                sys->pair.push_back(index3d(sys,i+1,j+1,k+1));//1 facard
            }
        }
    }
}


// translate cell (i,j,k) to cel position but with boundary fixing
int index3d(mdsys_t *sys, int i, int j, int k) {
    while (i>sys->cn-1) i -= sys->cn;
    while (i<0) i += sys->cn;
    while (j>sys->cn-1) j -= sys->cn;
    while (j<0) j += sys->cn;
    while (k>sys->cn-1) k -= sys->cn;
    while (k<0) k += sys->cn;
    return k+sys->cn*(j+i*sys->cn);
}
