#include <stdlib.h>
#include <math.h>
#include <vector>
#include "ljmd.h"

// sort particles into cells
void sort(mdsys_t *sys,cell_t *cel) {
    double time = stamp();
    for(int it=0;it<sys->cn*sys->cn*sys->cn;++it) {
        cel[it].idx.clear();
    }
    int pos[3];
    for(int i=0;i<sys->natoms;++i) {// loop through atoms
        pos[0] = floor((sys->rx[i]+0.5*sys->box)/sys->cl);
        pos[1] = floor((sys->ry[i]+0.5*sys->box)/sys->cl);
        pos[2] = floor((sys->rz[i]+0.5*sys->box)/sys->cl);
        // apply pierodic boundary condition
        int cidx = index3d(sys,pos[0],pos[1],pos[2]);// which cell
        cel[cidx].idx.push_back(i);
    }
    printf("sort timing %f\n",stamp()-time);
}

/* aggressive cell list */
void cell_force(mdsys_t *sys,cell_t *cell){
    double rsq,rsq_inv,r6,ffac;
    double rx,ry,rz;
    
    // zero energy and forces
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    double boxby2 = 0.5*sys->box;// pre-calculate
    double rcutsq = sys->rcut*sys->rcut;// pre-calculate, take square
    double c6 = sys->epsilon*pow(sys->sigma,6);
    double c12 = sys->epsilon*pow(sys->sigma,12);
    printf("pair size %lu \t",sys->pair.size());
    double time = stamp();
    for(int it=0;it<sys->pair.size()/2;++it) {// loop through cell pairs
        int c1 = sys->pair[2*it];
        int c2 = sys->pair[2*it+1];
        
        for(int ii=0;ii<cell[c1].idx.size();++ii){
            int i = cell[c1].idx[ii];
            for(int jj=0;jj<cell[c2].idx.size();++jj){
                int j = cell[c2].idx[jj];
                
                
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
            }// c1
        }// c2
    }// cell pairs
    printf("force timing %f\n",stamp()-time);
}
