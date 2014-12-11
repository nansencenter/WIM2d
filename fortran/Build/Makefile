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

# mex function to run non-interactive WIM in matlab
MEX1		:=	$(src)/interfaces/WIM2d_run_mex.F
TARGET_MEX1	:=	WIM2d_run_mex.mexmaci64

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
	$(FC) -c -o $@ $< $(FCFLAGS) $(IFLAGS) $(OFLAGS)

all:	sh py mex1

# Link objects to make executable $(TARGET_F)
sh: $(PROG_F) $(OBJ_MAIN) $(OBJS)
	$(FC) -o $(TARGET_F) $^ $(FCFLAGS) $(IFLAGS)
	mv $(TARGET_F) $(bin)

# Link objects to make python module $(TARGET_PY).so
py:	$(PROG_PY) $(OBJ_MAIN) $(OBJS)
	$(F2PY) -c -m $(TARGET_PY) $^ $(F2PY_FLAGS)
	mv $(TARGET_PY).so $(bin)

# Link objects to make python module $(TARGET_PY).so
mex1:	$(MEX1) $(OBJ_MAIN) $(OBJS)
	$(MXC) $(MEX_FLAGS) $^
	mv $(TARGET_MEX1) $(bin)

# DEPENDENCIES
$(ODIR)/mod_grid_info.o:	$(ODIR)/mod_file_utils.o	grid_info.h
$(ODIR)/mod_wim_init.o:		$(ODIR)/mod_gridinfo.o  	$(ODIR)/mod_common_wim.o
$(ODIR)/mod_common_wim.o:	$(ODIR)/mod_gridinfo.o
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
$(OBJ_MAIN):			$(OBJS) \
				print_screen_info.h \
				print_log_info.h \
				WIM2d_run.h


$(PROG_F): 		$(OBJ_MAIN) 	$(ODIR)/mod_gridinfo.o
$(PROG_PY): 		$(OBJS) \
				print_screen_info.h \
				print_log_info.h \
				WIM2d_run.h

.PHONY: clean vclean mex1

# clean: keep executable, but delete modules and objects
clean:
	rm $(ODIR)/* $(MDIR)/*

# vclean: delete executable, and delete modules and objects
vclean:
	rm $(bin)/$(TARGET_F)
	rm $(bin)/$(TARGET_PY).so
	rm $(bin)/$(TARGET_MEX1)
	rm $(ODIR)/* $(MDIR)/*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
