# c, fortran, and fortran-to-python linker compiler
CC	:=	gcc
FC	:=	gfortran
F2PY	:=	f2py

# root directory for source code
# (relative to this one)
src	:=	..

# keep binaries here
bin	:=	$(src)/bin
# bin	:=	test

# set path for compiler
# NB relative to this directory
VPATH 	:=	.:$(src)/misc

# directory where objects are kept
ODIR	:=	objs

# directory where modules are put into
MDIR	:=	mods

# header files for "include" statements
# NB relative to $(ODIR)
HDIR	:=	$(src)/header_files

# # optimise code
# FCFLAGS=-O
# check for errors; fPIC is to compile in linux (no effect in mac osx)
FCFLAGS=-g -fPIC

# include modules/headers from $(MDIR) or $(HDIR)
IFLAGS=-I$(MDIR) -I$(HDIR)

# Put modules in $(MDIR) when they are created
OFLAGS=-J$(MDIR)

# Flags for f2py
F2PY_FLAGS=--fcompiler=$(FC) --f77flags="$(IFLAGS) $(OFLAGS)"

# fortran program to run WIM
# - executable called by run_WIM2d.sh
#   (cd ../run; ./run_WIM2d.sh)
PROG_F		:=	$(ODIR)/p_save_grid.o
TARGET_F	:=	grid_setup.exec

# fortran interface program to run WIM in ipython
# - NB f2py needs the full path
PROG_PY		:=	save_grid_f2py.F
TARGET_PY	:=	save_grid_f2py

# OBJ_MAIN	:=	$(ODIR)/WIM2d_run.o

OBJS		:=      \
			$(ODIR)/mod_file_utils.o	\
			$(ODIR)/mod_save_grid.o

# Create objects by compiling their respective .F files
# (add dependency on header files also, so we recompile if headers change)
$(ODIR)/%.o:	%.F
	mkdir -p $(MDIR) $(ODIR)
	$(FC) -c -o $@ $< $(FCFLAGS) $(IFLAGS) $(OFLAGS)

# Link objects to make executable $(TARGET_F)
exec: $(PROG_F) $(OBJS)
	mkdir -p $(bin)
	$(FC) -o $(TARGET_F) $^ $(FCFLAGS) $(IFLAGS)
	mv $(TARGET_F) $(bin)

# Link objects to make python module $(TARGET_PY).so
py:	$(PROG_PY) $(OBJS)
	mkdir -p $(bin)
	$(F2PY) -c -m $(TARGET_PY) $^ $(F2PY_FLAGS)
	mv $(TARGET_PY).so $(bin)

all: exec py

# DEPENDENCIES
$(PROG_F): 		$(OBJS)
# $(PROG_PY): 		$(OBJS) 	$(DEPS)

.PHONY: clean vclean

# clean: keep executable, but delete modules and objects
clean:
	rm -f 	$(ODIR)/*.o \
	   	$(MDIR)/*.mod

# vclean: delete executable, and delete modules and objects
vclean:
	rm -f	$(bin)/$(TARGET_F) \
	   	$(bin)/$(TARGET_PY).so \
	   	$(ODIR)/*.o \
	   	$(MDIR)/*.mod

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
