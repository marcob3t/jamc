# -*- Makefile -*-
SHELL=/bin/sh
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=Obj-serial/%.o)
############################################

default: serial
         parallel

serial:
	$(MAKE) $(MFLAGS) -C Obj-$@

#parallel:
#         $(MAKE) $(MFLAGS) -MPCC Obj-$@

clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean
	$(MAKE) $(MFLAGS) -C examples clean
        $(MAKE) $(MFLAGS) -MPCC Obj-parallel clean

check: serial
	$(MAKE) $(MFLAGS) -C examples check

#check: parallel
#       $(MAKE) $(MFLAGS) -MPCC parallel check
