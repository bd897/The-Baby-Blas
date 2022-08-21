include ../Makefile.inc

# Timing Library targets 

all: lib

dot.o : dot.c
	$(CC) $(CFLAGS) -c dot.c

vvm.o : vvm.c
	$(CC) $(CFLAGS) -c vvm.c

mvv.o : mvv.c
	$(CC) $(CFLAGS) -c mvv.c

mmm.o : mmm.c
	$(CC) $(CFLAGS)  -c mmm.c

dls.o : dls.c
	$(CC) $(CFLAGS) -c dls.c

ils.o : ils.c
	$(CC) $(CFLAGS) -c ils.c

lib: dot.o vvm.o mvv.o mmm.o dls.o ils.o
	ar -rc libserial_bblas.a dot.o vvm.o mvv.o mmm.o dls.o ils.o
	ranlib libserial_bblas.a

# Default Targets for Cleaning up the Environment

clean :
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
