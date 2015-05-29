# c, fortran, and fortran-to-python linker compiler
CC	:=	gcc
FC	:=	gfortran
F2PY	:=	f2py
MXC	:=	mex

# root directory for source code
# (relative to this one)
src	:=	..

# keep binaries here
bin	:=	$(src)/bin

# set path for compiler
# NB relative to this directory
VPATH 	:=	$(src)/main:$(src)/misc:$(src)/init:$(src)/numerics
VPATH 	:=	$(VPATH):$(src)/attenuation_youngs:$(src)/header_files
VPATH 	:=	$(VPATH):$(src)/interfaces

# directory where objects are kept
ODIR	:=	objs

# directory where modules are put into
MDIR	:=	mods

# header files for "include" statements
# NB needs to be in VPATH and in -I$(HDIR) in calls to compiler
HDIR	:=	$(src)/header_files

# # optimise code
# FCFLAGS=-O
# check for errors
FCFLAGS=-g -fPIC

# include modules/headers from $(MDIR) or $(HDIR)
IFLAGS=-I$(MDIR) -I$(HDIR)

# Put modules in $(MDIR) when they are created
OFLAGS=-J$(MDIR)

# Flags for f2py
F2PY_FLAGS=--fcompiler=$(FC) --f77flags="$(IFLAGS) $(OFLAGS)"

# Flags for mex
MEX_FLAGS=-largeArrayDims $(IFLAGS)

# fortran program to run WIM
# - executable called by run_WIM2d.sh
#   (cd ../run; ./run_WIM2d.sh)
PROG_F		:=	$(ODIR)/p_WIM2d.o
TARGET_F	:=	WIM2d.exe

# fortran interface program to run WIM in ipython
# - NB f2py needs the full path
# - called by run_WIM2d.py
#   (cd ../run; ipython run_WIM2d.py)
PROG_PY		:=	$(src)/interfaces/WIM2d_f2py.F
TARGET_PY	:=	WIM2d_f2py

# mex function to run interactive WIM in matlab (old version)
MX_EXT		:=	mex
MX1		:=	WIM2d_run_io_mex
MEX1		:=	$(src)/interfaces/$(MX1).F
TARGET_MEX1	:=	$(MX1).$(MX_EXT)
# mex function is $(TARGET_MEX1)maci64  for mac OSX
# mex function is $(TARGET_MEX1)a64 	for linux

# mex function to run interactive WIM in matlab
MX2		:=	WIM2d_run_io_mex_v2
MEX2		:=	$(src)/interfaces/$(MX2).F
TARGET_MEX2	:=	$(MX2).$(MX_EXT)
# mex function is $(TARGET_MEX2)maci64  for mac OSX
# mex function is $(TARGET_MEX2)a64 	for linux

OBJ_MAIN	:=	$(ODIR)/mod_WIM2d_run.o

OBJS		:=      \
			$(ODIR)/mod_file_utils.o \
                        $(ODIR)/mod_gridinfo.o \
			$(ODIR)/mod_wim_init.o \
			$(ODIR)/mod_common_wim.o \
			$(ODIR)/mod_wim_prams.o \
			$(ODIR)/RTparam_outer.o \
			$(ODIR)/RTparam_fast.o \
			$(ODIR)/RTparam_hardcoded_v2.o \
			$(ODIR)/mod_RTparam_utilities.o \
			$(ODIR)/mod_adv_atten.o \
			$(ODIR)/mod_advect.o \
			$(ODIR)/mod_wavesice.o

# Create objects by compiling their respective .F files
# (add dependency on header files also, so we recompile if headers change)
$(ODIR)/%.o:	%.F
	mkdir -p $(ODIR) $(MDIR)
	$(FC) -c -o $@ $< $(FCFLAGS) $(IFLAGS) $(OFLAGS)

# Link objects to make executable $(TARGET_F)
# (default)
exe: $(PROG_F) $(OBJ_MAIN) $(OBJS)
	mkdir -p $(bin)
	$(FC) -o $(TARGET_F) $^ $(FCFLAGS) $(IFLAGS)
	mv $(TARGET_F) $(bin)

# Link objects to make python module $(TARGET_PY).so
py:	$(PROG_PY) $(OBJ_MAIN) $(OBJS)
	mkdir -p $(bin)
	$(F2PY) -c -m $(TARGET_PY) $^ $(F2PY_FLAGS)
	mv $(TARGET_PY).so $(bin)

# Link objects to make python module $(TARGET_PY).so
mex1:	$(MEX1) $(OBJ_MAIN) $(OBJS)
	mkdir -p $(bin)
	$(MXC) $(MEX_FLAGS) $^
	mv $(TARGET_MEX1)* $(bin)

mex2:	$(MEX2) $(OBJ_MAIN) $(OBJS)
	mkdir -p $(bin)
	$(MXC) $(MEX_FLAGS) $^
	mv $(TARGET_MEX2)* $(bin)

all:	py exe mex1 mex2

# DEPENDENCIES
$(ODIR)/mod_grid_info.o:	$(ODIR)/mod_file_utils.o	grid_info.h
$(ODIR)/mod_wim_init.o:		$(ODIR)/mod_gridinfo.o  	$(ODIR)/mod_common_wim.o
$(ODIR)/mod_common_wim.o:	$(ODIR)/mod_gridinfo.o		wave_info.h
$(ODIR)/mod_advect.o:		$(ODIR)/mod_gridinfo.o  	$(ODIR)/mod_common_wim.o
$(ODIR)/mod_adv_atten.o:	$(ODIR)/mod_gridinfo.o  	$(ODIR)/mod_common_wim.o \
				$(ODIR)/mod_advect.o
$(ODIR)/mod_wavesice.o:		$(ODIR)/mod_gridinfo.o		$(ODIR)/mod_common_wim.o \
				$(ODIR)/mod_advect.o 		$(ODIR)/mod_adv_atten.o
#
$(ODIR)/RTparam_outer.o:        $(ODIR)/RTparam_fast.o
$(ODIR)/RTparam_fast.o:         $(ODIR)/RTparam_hardcoded_v2.o
$(ODIR)/RTparam_hardcoded_v2.o: $(ODIR)/mod_RTparam_utilities.o
#
$(OBJ_MAIN):			$(OBJS)

$(PROG_F): 		$(OBJ_MAIN) 	$(ODIR)/mod_gridinfo.o
$(PROG_PY): 		$(OBJS)

$(MEX1):		grid_info.h $(OBJ_MAIN)

$(MEX2):		grid_info.h $(OBJ_MAIN)

.PHONY: clean vclean mex1 mex2

# clean: keep executable, but delete modules and objects
clean:
	rm $(ODIR)/*.o $(MDIR)/*.mod

# vclean: delete executable, and delete modules and objects
vclean:
	rm $(ODIR)/*.o $(MDIR)/*.mod \
	   $(bin)/$(TARGET_F) \
	   $(bin)/$(TARGET_PY).so \
	   $(bin)/$(TARGET_MEX1)* \
	   $(bin)/$(TARGET_MEX2)*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
