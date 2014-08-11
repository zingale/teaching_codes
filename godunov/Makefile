FC = gfortran
FFLAGS = -c -O2

all: godunov

FSOURCE = godunov.f90 riemann.f90

%.o: %.f90
	$(FC) $(FFLAGS) $< 

OBJECTS = $(FSOURCE:.f90=.o)

godunov: $(OBJECTS)
	gfortran -o godunov $(OBJECTS)

.PHONY: clean

clean:
	rm -f *.o *.mod
