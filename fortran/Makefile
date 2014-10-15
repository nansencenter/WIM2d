CC=gcc
FC=gfortran

FCFLAGS=-g
# FCFLAGS=-O
CFLAGS=-I. # include directory = "."
# DEPS=RTparam_hardcoded.h RTparam_fast.h # dependencies
PROG=test_advection_weno.o 

MODS=mod_gridinfo.o mod_advection_weno.o mod_file_utils.o
TARGET=test_advection_weno

# Create objects by compiling their respective .F files
%.o: %.F
	$(FC) -c -o $@ $< $(FCFLAGS) $(CFLAGS)

# Link objects to make executable $(TARGET)
default: $(PROG) $(MODS)
	$(FC) -o $(TARGET) $^ $(FCFLAGS) $(CFLAGS)

mod_advection_weno.o:	 	mod_gridinfo.o mod_file_utils.o
$(PROG): 			$(MODS)

.PHONY: clean

clean:
	rm $(TARGET) *.o *.mod

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
