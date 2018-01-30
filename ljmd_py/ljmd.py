#!/usr/bin/env python

from ctypes import *
from math import *
import sys

class MDSYS(Structure):
   _parms_ = [ ("natoms", c_int), ("nfi", c_int), ("nsteps", c_int), ("dt", c_double), ("redmass", c_double), ("epsilon", c_double), ("sigma", c_double), ("box", c_double), ("rcut", c_double), ("ekin", c_double), ("epot", c_double), ("temp", c_double), ("rx", POINTER(c_double)), ("ry", POINTER(c_double)), ("rz", POINTER(c_double)),("vx", POINTER(c_double)),("vy", POINTER(c_double)), ("vz", POINTER(c_double)), ("fx", POINTER(c_double)),("fy", POINTER(c_double)), ("fz", POINTER(c_double))]

def read():

    my_list = sys.stdin.readlines()
    x = []
    for i in my_list:
        tmp=i.split()
        dat=tmp[0]
        x.append(dat)
    
    return (x)


if __name__ == "__main__":
    argv = sys.argv[:]

    dso = CDLL("../ljmd.so")
    
    list=read()
    
    mdsys=MDSYS()
    mvsq2e=2390.05736153349
    
    mdsys.natoms=int(list[0])
    mdsys.redmass =float(list[1]) *mvsq2e
    mdsys.epsilon=float(list[2])
    mdsys.sigma=float(list[3])
    mdsys.rcut=float(list[4])
    mdsys.box=float(list[5])
    restfile=(list[6])
    trajfile=(list[7])
    ergfile=(list[8])
    mdsys.nsteps=int(list[9])
    mdsys.dt=float(list[10])
    nprint=int(list[11])
    
    print(mdsys.natoms)
    print(mdsys.redmass)
    print(mdsys.epsilon)
    print(mdsys.sigma)
    print(mdsys.rcut)
    print(mdsys.box)
    print(restfile)
    print(trajfile)
    print(ergfile)
    print(mdsys.nsteps)
    print(mdsys.dt)
    print(nprint)
