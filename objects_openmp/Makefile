# -*- Makefile -*-
SHELL=/bin/sh
CC=g++
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -fopenmp
LDLIBS=-lm

# list of source files
SRCDIR=../src
OBJDIR=./
SOURCES= $(shell find $(SRCDIR) -type f -name *.c)
OBJ= $(patsubst $(SRCDIR)/%,./%,$(SOURCES:.c=.o))

default: ../ljmd-openmp.x

clean:
	rm -f *.mod *.o ../ljmd-openmp.x

# compilation pattern rule for objects
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

# special for cell
$(OBJDIR)/cell.o: $(SRCDIR)/cell.cc $(SRCDIR)/cell.h
	$(CC) -c $(CFLAGS) -I$(SRCDIR) -I. $< -o $@

# linker rule
../ljmd-openmp.x: $(OBJ) cell.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)


