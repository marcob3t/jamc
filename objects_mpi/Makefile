# -*- Makefile for unit test -*-
CC=mpic++
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -DUSE_MPI #-DTIMING
LDLIBS=-lm
INC= -I$(INCDIR) -I.

OBJDIR=.
SRCDIR=../src
INCDIR=../src
SOURCES= $(shell find $(SRCDIR) -type f -name *.c)

OBJ= $(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SOURCES:.c=.o))
HEADS=$(wildcard $(INCDIR)*.h)

default: ../ljmd-mpi.x

# linker rule
../ljmd-mpi.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compile true source files
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(HEADS)
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f ../ljmd-mpi.x $(OBJDIR)/*.o
