include ../Makefile.inc

# OpenMP library targets 

all: lib

dot.o : dot.c
	$(CC) $(CFLAGS) -fopenmp -c dot.c

vvm.o : vvm.c
	$(CC) $(CFLAGS) -fopenmp -c vvm.c

mvv.o : mvv.c
	$(CC) $(CFLAGS) -fopenmp -c mvv.c

mmm.o : mmm.c
	$(CC) $(CFLAGS) -fopenmp -c mmm.c

dls.o : dls.c
	$(CC) $(CFLAGS) -fopenmp -c dls.c

ils.o : ils.c
	$(CC) $(CFLAGS) -fopenmp -c ils.c

lib: dot.o vvm.o mvv.o mmm.o dls.o ils.o
	ar -rc libopenmp_bblas.a  dot.o vvm.o mvv.o mmm.o dls.o ils.o
	ranlib libopenmp_bblas.a

# Default Targets for Cleaning up the Environment

clean :
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
