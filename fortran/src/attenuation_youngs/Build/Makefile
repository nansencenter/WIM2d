# c, fortran, and fortran-to-python linker compiler
CC	:=	gcc
FC	:=	gfortran
F2PY	:=	f2py

# root directory for source code
# (relative to this one)
src	:=	../..

# keep binaries here
bin	:=	$(src)/bin

# directory where objects are kept
ODIR	:=	objs

# directory where modules are put into
MDIR	:=	mods

# header files for "include" statements
# NB relative to $(ODIR)
HDIR	:=	$(src)/../include

# set paths for compiler
# NB relative to this directory
VPATH 	:=	.:$(HDIR):$(src)/attenuation_youngs:$(src)/main


# # optimise code
# FCFLAGS=-O
# check for errors
FCFLAGS=-g

# include modules/headers from $(MDIR) or $(HDIR)
IFLAGS=-I$(MDIR) -I$(HDIR)

# Put modules in $(MDIR) when they are created
OFLAGS=-J$(MDIR)

# Flags for f2py
F2PY_FLAGS=--fcompiler=$(FC) --f77flags="$(IFLAGS) $(OFLAGS)" --include-paths $(HDIR)

# fortran program to test attenuation
PROG_F		:=	$(ODIR)/test_atten.o
TARGET_F	:=	test_atten.exec

# fortran interface program to run WIM in ipython
# - NB f2py needs the full path
# - enable access to attenuation coeffiicents in python
PROG_PY		:=	atten_young_f2py.F
TARGET_PY	:=	atten_young_f2py

OBJS		:=	$(ODIR)/mod_common_wim.o \
   			$(ODIR)/mod_atten_young.o \
			$(ODIR)/mod_atten_young_hardcoded.o \
			$(ODIR)/mod_RTparam_utilities.o

# DEPS=print_screen_info.h \
#      print_log_info.h


# Create objects by compiling their respective .F files
# (add dependency on header files also, so we recompile if headers change)
$(ODIR)/%.o:	%.F
	mkdir -p $(ODIR)
	mkdir -p $(MDIR)
	$(FC) -c -o $@ $< $(FCFLAGS) $(IFLAGS) $(OFLAGS)

# Link objects to make executable $(TARGET_F)
exec: $(PROG_F) $(OBJS)
	$(FC) -o $(TARGET_F) $^ $(FCFLAGS) $(IFLAGS)

# Link objects to make python module $(TARGET_PY).so
py:	$(PROG_PY) $(OBJS)
	$(F2PY) -c -m $(TARGET_PY) $^ $(F2PY_FLAGS)

all: exec py

# DEPENDENCIES
#
$(ODIR)/mod_atten_young.o:		$(ODIR)/mod_atten_young_hardcoded.o	$(ODIR)/mod_common_wim.o
$(ODIR)/mod_atten_young_hardcoded.o:	$(ODIR)/mod_RTparam_utilities.o
#
$(PROG_F): 		$(ODIR)/mod_atten_young.o	$(ODIR)/mod_common_wim.o
$(PROG_PY): 		$(OBJS)				$(ODIR)/mod_common_wim.o

.PHONY: clean

# clean: keep executable, but delete modules and objects
clean:
	rm -f $(ODIR)/* $(MDIR)/*

# vclean: delete executable, and delete modules and objects
vclean:
	rm -f $(TARGET_F) $(TARGET_PY).so $(ODIR)/* $(MDIR)/*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
