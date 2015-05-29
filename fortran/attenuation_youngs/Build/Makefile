# c, fortran, and fortran-to-python linker compiler
CC	:=	gcc
FC	:=	gfortran
F2PY	:=	f2py

# root directory for source code
# (relative to this one)
src	:=	../..

# keep binaries here
bin	:=	$(src)/bin

# set path for compiler
# NB relative to this directory
VPATH 	:=	.:$(src)/init:$(src)/attenuation_youngs

# directory where objects are kept
ODIR	:=	objs

# directory where modules are put into
MDIR	:=	mods

# header files for "include" statements
# NB relative to $(ODIR)
HDIR	:=	$(src)/header_files

# # optimise code
# FCFLAGS=-O
# check for errors
FCFLAGS=-g

# include modules/headers from $(MDIR) or $(HDIR)
IFLAGS=-I$(MDIR) -I$(HDIR)

# Put modules in $(MDIR) when they are created
OFLAGS=-J$(MDIR)

# Flags for f2py
F2PY_FLAGS=--fcompiler=$(FC) --f77flags="$(IFLAGS) $(OFLAGS)"

# fortran program to run WIM
# - executable called by run_WIM2d.sh
#   (cd ../run; ./run_WIM2d.sh)
PROG_F		:=	$(ODIR)/test_atten.o
TARGET_F	:=	test_atten.exe

# fortran interface program to run WIM in ipython
# - NB f2py needs the full path
# - enable access to attenuation coeffiicents in python
PROG_PY		:=	$(src)/attenuation_youngs/RTparam_outer.F
TARGET_PY	:=	RTparam_outer

OBJ_MAIN	:=	$(ODIR)/RTparam_outer.o

OBJS		:=	\
			$(ODIR)/mod_wim_prams.o \
			$(ODIR)/RTparam_fast.o \
			$(ODIR)/RTparam_hardcoded_v2.o \
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
exe: $(PROG_F) $(OBJ_MAIN) $(OBJS)
	$(FC) -o $(TARGET_F) $^ $(FCFLAGS) $(IFLAGS)
	mkdir -p $(bin)
	cp $(TARGET_F) $(bin)

# Link objects to make python module $(TARGET_PY).so
py:	$(PROG_PY) $(OBJS)
	$(F2PY) -c -m $(TARGET_PY) $^ $(F2PY_FLAGS)
	mkdir -p $(bin)
	cp $(TARGET_PY).so $(bin)

all: exe py

# DEPENDENCIES
#
$(ODIR)/RTparam_fast.o:         $(ODIR)/RTparam_hardcoded_v2.o
$(ODIR)/RTparam_hardcoded_v2.o: $(ODIR)/mod_RTparam_utilities.o
#
$(OBJ_MAIN):		$(OBJS)
$(PROG_F): 		$(OBJ_MAIN) $(ODIR)/mod_wim_prams.o
$(PROG_PY): 		$(OBJS)

.PHONY: clean

# clean: keep executable, but delete modules and objects
clean:
	rm $(ODIR)/* $(MDIR)/*

# vclean: delete executable, and delete modules and objects
vclean:
	rm $(bin)/$(TARGET_F) $(bin)/$(TARGET_PY).so $(ODIR)/* $(MDIR)/*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
