#!/usr/bin/env python

from ctypes import *
from math import *
from sys import *

class parm(Srtucture):
   _mdsys_ = [ ("natoms", c_int), ("nfi", c_int), ("nsteps", c_int), ("dt", c_double), ("redmass", c_double), ("epsilon", c_double), ("sigma", c_double), ("box", c_double), ("rcut", c_double), ("ekin", c_double), ("epot", c_double), ("temp", c_double), ("rx", POINTER(c_double)), ("ry", POINTER(c_double)), ("rz", POINTER(c_double)),("vx", POINTER(c_double)),("vy", POINTER(c_double)), ("vz", POINTER(c_double)), ("fx", POINTER(c_double)),("fy", POINTER(c_double)), ("fz", POINTER(c_double))]

def read():

    my_list = sys.stdin.readlines()
    x = []
    for i in my_list:
        tmp=i.split()
        dat=tmp[0]
        x.append(dat)
    
    return (x)

mvsq2e=2390.05736153349

if __name__ == "__main__":
    argv = sys.argv[:]

    dso = CDLL("./ljmd.so")
    
    sys=parm()
    
    list=read()
    
    sys.natoms=c_int(list[0])
    sys.redmass =c_double(list[1])*mvsq2e
    sys.epsilon=c_double(list[2])
    sys.sigma=c_double(list[3])
    sys.rcut=c_double(list[4])
    sys.box=c_double(list[5])
    restfile=(list[6])
    trajfile=(list[7])
    ergfile=(list[8])
    sys.nsteps=c_int(list[9])
    sys.dt=c_double(list[10])
    nprint=c_int(list[11])
    
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




















/
