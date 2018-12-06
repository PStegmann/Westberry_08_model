#
#    Copyright © 2018 Patrick Stegmann
#.   This file is part of Westberry_2008_model.
#

# This makefile works with the GNU make command, the one find on
# GNU/Linux systems and often called gmake on non-GNU systems, if you
# are using an old style make command, please see the file
# Makefile_oldstyle provided with the package.

# ======================================================================
# Declarations:
# ======================================================================

# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -O0 -g -fbounds-check
#FCFLAGS = -Ofast
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = wb08

# "make" builds all
all: $(PROGRAMS)

# dependencies:
par.o: interp.o
global.o: par.o interp.o
wb08.o: global.o par.o interp.o
wb08: global.o par.o interp.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: ../src/%.f90
	$(FC) $(FCFLAGS) -c $<

%.o: ../src/%.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
