FC = gfortran
FFLAGS = -c -O2

CC = gcc
CFLAGS = -c -O2

all: godunov

FSOURCE = godunov.f90 
CSOURCE = riemann.c

%.o: %.f90
	$(FC) $(FFLAGS) $< 

%.o: %.c
	$(CC) $(CFLAGS) $<

OBJECTS = $(FSOURCE:.f90=.o) $(CSOURCE:.c=.o)

godunov: $(OBJECTS)
	gfortran -o godunov $(OBJECTS)


.PHONY: clean

clean:
	rm -f *.o *.mod
