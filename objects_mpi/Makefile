# -*- Makefile for unit test -*-
SHELL=/bin/sh
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
../ljmd-mpi.x: $(OBJ) $(SRCDIR)/cell.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compile true source files
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(HEADS)
	$(CC) -c $(CFLAGS) $< -o $@

# special for cell
$(OBJDIR)/cell.o: $(SRCDIR)/%.cc $(HEADS)
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f ../ljmd-mpi.x *.o
