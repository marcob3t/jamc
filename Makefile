# -*- Makefile -*-
SHELL=/bin/sh
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=objects_serial/%.o)
############################################

default: serial openmp mpi

serial:
	$(MAKE) $(MFLAGS) -C objects_$@

openmp:
	$(MAKE) $(MFLAGS) -C objects_$@

mpi:
	$(MAKE) $(MFLAGS) -C objects_$@

unitest:
	$(MAKE) $(MFLAGS) -C test default

clean:
	$(MAKE) $(MFLAGS) -C objects_serial clean
	$(MAKE) $(MFLAGS) -C objects_openmp clean
	$(MAKE) $(MFLAGS) -C objects_mpi clean
	$(MAKE) $(MFLAGS) -C examples clean
	$(MAKE) $(MFLAGS) -C test clean

check: serial openmp unitest mpi
	$(MAKE) $(MFLAGS) -C examples check
	#make clean

