# -*- Makefile -*-
SHELL=/bin/sh
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=objects_serial/%.o)
############################################

default: serial

serial:
	$(MAKE) $(MFLAGS) -C objects_$@

clean:
	$(MAKE) $(MFLAGS) -C objects_serial clean
	$(MAKE) $(MFLAGS) -C examples clean

check: serial
	$(MAKE) $(MFLAGS) -C examples check

test:
	$(MAKE) $(MFLAGS) -C test

