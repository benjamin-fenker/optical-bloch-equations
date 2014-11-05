CC=icpc
CFLAGS=-openmp -fpic -c -ansi -pedantic -Wall -O3 -xHost -W -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wwrite-strings -fshort-enums -fno-common -g
# Should really use -Werror but gsl_odeiv2 breaks when I do that
GSLFLAGS=-lgsl -lgslcblas -lm -openmp
LDFLAGS=
SOURCES=setup_optical_pumping.cc optical_pumping.cc alkali.cc eigenvector_helper.cc optical_pumping_method.cc rate_equations.cc optical_pumping_data_structures.cc density_matrix.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=opticalPumping

all: $(OBJECTS) $(EXCECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(GSLFLAGS)
.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean::
	rm -rf $(OBJECTS)
	rm -rf $(EXECUTABLE)
