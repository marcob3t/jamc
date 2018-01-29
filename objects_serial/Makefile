# -*- Makefile -*-
SHELL=/bin/sh
CC=g++
CFLAGS=-Wall -g -fomit-frame-pointer -O3 -ffast-math -fexpensive-optimizations -DCELL #-DTIMING #-pg
LDLIBS=-lm

# list of source files
SRCDIR=../src
OBJDIR=./
SOURCES= $(shell find $(SRCDIR) -type f -name *.c)
OBJ= $(patsubst $(SRCDIR)/%,./%,$(SOURCES:.c=.o))

default: ../ljmd-serial.x

# compilation pattern rule for objects
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

# special for cell
$(OBJDIR)/cell.o: $(SRCDIR)/cell.cc
	$(CC) -c $(CFLAGS) -I$(SRCDIR) -I. $< -o $@

# linker rule
../ljmd-serial.x: $(OBJ) cell.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

clean:
	rm -f *.o ../ljmd-serial.x
