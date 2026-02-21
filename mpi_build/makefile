
GCC	=	gcc
# CC      =       /opt/openmpi-5.0.8/bin/mpicc
# CCLINK  =       /opt/openmpi-5.0.8/bin/mpicc
CC      =       mpicc
CCLINK  =       mpicc
SHELL   =       /bin/sh
CFLAGS  =       -O2

all: mpi

mpi: ../main.c tmp/nsb_stack.o tmp/node.o tmp/path.o tmp/code_path.o tmp/util.o tmp/nauty.o tmp/nautil.o tmp/naugraph.o tmp/schreier.o tmp/naurng.o tmp/gtools.o tmp/mpi_routines.o
	$(CC) $(CFLAGS) -lm -o nauty_mpi -DMPI ../main.c tmp/nsb_stack.o tmp/node.o tmp/path.o tmp/code_path.o tmp/util.o tmp/nauty.o tmp/nautil.o tmp/naugraph.o tmp/schreier.o tmp/naurng.o tmp/gtools.o tmp/mpi_routines.o -I../lib -I../src/nauty2_8_9

tmp/nsb_stack.o: ../src/nsb_stack.c ../lib/nsb_stack.h
	$(GCC) $(CFLAGS) -c ../src/nsb_stack.c -o tmp/nsb_stack.o -I../lib -I../src/nauty2_8_9

tmp/node.o: ../src/node.c ../lib/node.h
	$(GCC) $(CFLAGS) -c ../src/node.c -o tmp/node.o -I../lib -I../src/nauty2_8_9

tmp/path.o: ../src/path.c ../lib/path.h
	$(GCC) $(CFLAGS) -c ../src/path.c -o tmp/path.o -I../lib -I../src/nauty2_8_9

tmp/code_path.o: ../src/code_path.c ../lib/code_path.h
	$(GCC) $(CFLAGS) -c ../src/code_path.c -o tmp/code_path.o -I../lib -I../src/nauty2_8_9

tmp/util.o: ../src/util.c ../lib/util.h
	$(GCC) $(CFLAGS) -c ../src/util.c -o tmp/util.o -I../lib -I../src/nauty2_8_9


#Overrides
tmp/nautil.o: ../src/nautil.c ../src/nauty2_8_9/*
	$(GCC) $(CFLAGS) -c ../src/nautil.c -o tmp/nautil.o -I../lib -I../src/nauty2_8_9

tmp/naugraph.o: ../src/naugraph.c ../src/nauty2_8_9/*
	$(GCC) $(CFLAGS) -c -DMPI ../src/naugraph.c -o tmp/naugraph.o -I../lib -I../src/nauty2_8_9


# Base Nauty
tmp/schreier.o: ../src/nauty2_8_9/schreier.c ../src/nauty2_8_9/*
	$(GCC) $(CFLAGS) -c ../src/nauty2_8_9/schreier.c -o tmp/schreier.o -I../lib -I../src/nauty2_8_9

tmp/naurng.o: ../src/nauty2_8_9/naurng.c ../src/nauty2_8_9/*
	$(GCC) $(CFLAGS) -c ../src/nauty2_8_9/naurng.c -o tmp/naurng.o -I../lib -I../src/nauty2_8_9

tmp/gtools.o: ../src/nauty2_8_9/gtools.c ../src/nauty2_8_9/*
	$(GCC) $(CFLAGS) -c ../src/nauty2_8_9/gtools.c -o tmp/gtools.o -I../lib -I../src/nauty2_8_9


# MPI
tmp/nauty.o: ../src/nauty.c ../src/nauty2_8_9/*
	$(CC) $(CFLAGS) -c -lm -DMPI ../src/nauty.c -o tmp/nauty.o -I../lib -I../src/nauty2_8_9

tmp/mpi_routines.o: ../src/mpi_routines.c ../lib/mpi_routines.h ../src/nauty2_8_9/*
	$(CC) $(CFLAGS) -c -lm -DMPI ../src/mpi_routines.c -o tmp/mpi_routines.o -I../lib -I../src/nauty2_8_9
	

clean:
	rm -f nauty_mpi tmp/*.o
