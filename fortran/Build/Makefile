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
VPATH 	:=	$(VPATH):$(src)/interfaces:$(src)/extern

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
MEX_FLAGS=-largeArrayDims -g $(IFLAGS)

# fortran program to run WIM
# - executable called by run_WIM2d.sh
#   (cd ../run; ./run_WIM2d.sh)
PROG_F		:=	$(ODIR)/p_WIM2d.o
TARGET_F	:=	WIM2d.exec

# fortran interface program to run WIM in ipython
# - NB f2py needs the full path
# - called by run_WIM2d.py
#   (cd ../run; ipython run_WIM2d.py)
PROG_PY		:=	$(src)/interfaces/WIM2d_f2py.F
TARGET_PY	:=	WIM2d_f2py

# mex function to run interactive WIM in matlab
MX_EXT		:=	mex
MX2		:=	WIM2d_run_io_mex_v2
MEX2		:=	$(src)/interfaces/$(MX2).F
TARGET_MEX2	:=	$(MX2).$(MX_EXT)
# mex function is $(TARGET_MEX2)maci64  for mac OSX
# mex function is $(TARGET_MEX2)a64 	for linux

# mex function to run interactive WIM in matlab
MX3		:=	WIM2d_run_io_mex_vSdir
MEX3		:=	$(src)/interfaces/$(MX3).F
TARGET_MEX3	:=	$(MX3).$(MX_EXT)
# mex function is $(TARGET_MEX3)maci64  for mac OSX
# mex function is $(TARGET_MEX3)a64 	for linux

# mex function to run interactive WIM in matlab
MX4		:=	WIM2d_run_io_mex_vSdir_mesh
MEX4		:=	$(src)/interfaces/$(MX4).F
TARGET_MEX4	:=	$(MX4).$(MX_EXT)
# mex function is $(TARGET_MEX4)maci64  for mac OSX
# mex function is $(TARGET_MEX4)a64 	for linux

OBJ_MAIN	:=	$(ODIR)/mod_WIM2d_run.o

OBJS		:=      \
			$(ODIR)/mod_interp.o \
			$(ODIR)/mod_file_utils.o \
                        $(ODIR)/mod_gridinfo.o \
			$(ODIR)/mod_common_wim.o \
			$(ODIR)/mod_wim_init.o \
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
$(ODIR)/%.o:	%.f90
	mkdir -p $(ODIR) $(MDIR)
	$(FC) -c -o $@ $< $(FCFLAGS) $(IFLAGS) $(OFLAGS)
$(ODIR)/%.o:	%.F
	mkdir -p $(ODIR) $(MDIR)
	$(FC) -c -o $@ $< $(FCFLAGS) $(IFLAGS) $(OFLAGS)

# Link objects to make executable $(TARGET_F)
# (default)
exec: $(PROG_F) $(OBJ_MAIN) $(OBJS)
	mkdir -p $(bin)
	$(FC) -o $(TARGET_F) $^ $(FCFLAGS) $(IFLAGS)
	mv $(TARGET_F) $(bin)

# Link objects to make python module $(TARGET_PY).so
py:	$(PROG_PY) $(OBJ_MAIN) $(OBJS)
	mkdir -p $(bin)
	$(F2PY) -c -m $(TARGET_PY) $^ $(F2PY_FLAGS)
	mv $(TARGET_PY).so $(bin)

# Link objects to make python module $(TARGET_PY).so
mex2:	$(MEX2) $(OBJ_MAIN) $(OBJS)
	mkdir -p $(bin)
	$(MXC) $(MEX_FLAGS) $^
	mv $(TARGET_MEX2)* $(bin)

mex3:	$(MEX3) $(OBJ_MAIN) $(OBJS)
	mkdir -p $(bin)
	$(MXC) $(MEX_FLAGS) $^
	mv $(TARGET_MEX3)* $(bin)

mex4:	$(MEX4) $(OBJ_MAIN) $(OBJS)
	mkdir -p $(bin)
	$(MXC) $(MEX_FLAGS) $^
	mv $(TARGET_MEX4)* $(bin)

all:	py exec mex2 mex3 mex4

All:	py exec

# DEPENDENCIES
$(ODIR)/mod_grid_info.o:	$(ODIR)/mod_file_utils.o	grid_info.h
$(ODIR)/mod_wim_init.o:		$(ODIR)/mod_gridinfo.o  	$(ODIR)/mod_common_wim.o
$(ODIR)/mod_common_wim.o:	$(ODIR)/mod_gridinfo.o		wave_info.h
$(ODIR)/mod_advect.o:		$(ODIR)/mod_gridinfo.o  	$(ODIR)/mod_common_wim.o
$(ODIR)/mod_adv_atten.o:	$(ODIR)/mod_gridinfo.o  	$(ODIR)/mod_common_wim.o \
				$(ODIR)/mod_advect.o
$(ODIR)/mod_wavesice.o:		$(ODIR)/mod_gridinfo.o		$(ODIR)/mod_common_wim.o \
				$(ODIR)/mod_advect.o 		$(ODIR)/mod_adv_atten.o	 \
				$(ODIR)/mod_interp.o
# $(ODIR)/interp.o:		$(ODIR)/r8lib.o
#
$(ODIR)/RTparam_outer.o:        $(ODIR)/RTparam_fast.o
$(ODIR)/RTparam_fast.o:         $(ODIR)/RTparam_hardcoded_v2.o
$(ODIR)/RTparam_hardcoded_v2.o: $(ODIR)/mod_RTparam_utilities.o
#
$(OBJ_MAIN):			$(OBJS)

$(PROG_F): 		$(OBJ_MAIN) 	$(ODIR)/mod_gridinfo.o
$(PROG_PY): 		grid_info.h	wave_info.h		$(OBJS)

$(MEX1):		grid_info.h $(OBJ_MAIN)
$(MEX2):		grid_info.h $(OBJ_MAIN)
$(MEX3):		grid_info.h wave_info.h	$(OBJ_MAIN)
$(MEX4):		grid_info.h wave_info.h	$(OBJ_MAIN)

.PHONY: clean vclean mex2 mex3 mex4

# clean: keep executable, but delete modules and objects
clean:
	rm -f $(ODIR)/*.o $(MDIR)/*.mod

# vclean: delete executable, and delete modules and objects
vclean:
	rm -f	$(ODIR)/*.o $(MDIR)/*.mod \
	   	$(bin)/$(TARGET_F) \
	   	$(bin)/$(TARGET_PY).so \
	   	$(bin)/$(TARGET_MEX2)* \
	   	$(bin)/$(TARGET_MEX3)* \
	   	$(bin)/$(TARGET_MEX4)*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
