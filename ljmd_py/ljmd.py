#!/usr/bin/env python

from ctypes import *
from math import *
import sys

class parm(Structure):
   _mdsys_ = [ ("natoms", c_int), ("nfi", c_int), ("nsteps", c_int), ("dt", c_double), ("redmass", c_double), ("epsilon", c_double), ("sigma", c_double), ("box", c_double), ("rcut", c_double), ("ekin", c_double), ("epot", c_double), ("temp", c_double), ("rx", POINTER(c_double)), ("ry", POINTER(c_double)), ("rz", POINTER(c_double)),("vx", POINTER(c_double)),("vy", POINTER(c_double)), ("vz", POINTER(c_double)), ("fx", POINTER(c_double)),("fy", POINTER(c_double)), ("fz", POINTER(c_double))]

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
    
    sys=parm()
    mvsq2e=2390.05736153349
    
    sys.natoms=int(list[0])
    sys.redmass =float(list[1]) *mvsq2e
    sys.epsilon=float(list[2])
    sys.sigma=float(list[3])
    sys.rcut=float(list[4])
    sys.box=float(list[5])
    restfile=(list[6])
    trajfile=(list[7])
    ergfile=(list[8])
    sys.nsteps=int(list[9])
    sys.dt=float(list[10])
    nprint=int(list[11])
    
    print(sys.natoms)
    print(sys.redmass)
    print(sys.epsilon)
    print(sys.sigma)
    print(sys.rcut)
    print(sys.box)
    print(restfile)
    print(trajfile)
    print(ergfile)
    print(sys.nsteps)
    print(sys.dt)
    print(nprint)
