# c, fortran, and fortran-to-python linker compiler
CC	:=	gcc
FC	:=	gfortran
F2PY	:=	f2py
MXC	:=	mex

# root directory for source code
# (relative to this one)
P       :=  $(shell pwd)
# @$(shell echo) $(P)
src	:=  $(WIM2D_PATH)/fortran/src
# @echo $(src)

# keep binaries here
bin	:=	../bin

# set path for compiler
# NB relative to this directory
VPATH 	:=  $(src)/main
VPATH   +=  $(src)/misc
VPATH   +=  $(src)/init
VPATH   +=  $(src)/numerics
VPATH 	+=  $(src)/attenuation_youngs
VPATH 	+=  $(src)/interfaces
VPATH 	+=  $(src)/extern
VPATH 	+=  .

# directory where objects are kept
# ODIR	:=	$(P)/objs
ODIR	:=	./objs

# directory where modules are put into
# MDIR	:=	$(P)/mods
MDIR	:=	./mods

# header files for "include" statements
# NB needs to be in VPATH and in -I$(HDIR) in calls to compiler
HDIR	:=	.
VPATH 	+= 	$(HDIR)

# # optimise code
# FCFLAGS=-O
# check for errors
FCFLAGS=-g -fPIC

# include modules/headers from $(MDIR) or $(HDIR)
IFLAGS=-I$(MDIR) -I$(HDIR) -I.

# Put modules in $(MDIR) when they are created
OFLAGS=-J$(MDIR)

# Flags for mex
MEX_FLAGS=-largeArrayDims -g $(IFLAGS)

# ============================================================================
# run in double precision
# - export WIM2D_TYPE=single to use single precision;
#   leave undefined or set to double (for example) to use double precision
ifneq ($(WIM2D_TYPE),single)
	FCFLAGS += -fdefault-real-8 -fdefault-double-8
endif
# ============================================================================

# Flags for f2py
F2PY_FLAGS=--fcompiler=$(FC) --f77flags="$(IFLAGS) $(OFLAGS) $(FCFLAGS)" --f90flags="$(IFLAGS) $(OFLAGS) $(FCFLAGS)"

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
MX1		:=	WIM2d_run_io_mex_v2
MEX1		:=	$(src)/interfaces/$(MX1).F
TARGET_MEX1	:=	$(MX1).mex
# mex function is $(TARGET_MEX1)maci64  for mac OSX
# mex function is $(TARGET_MEX1)a64 	for linux

# mex function to run interactive WIM in matlab
MX2		:=	WIM2d_run_io_mex_vSdir
MEX2		:=	$(src)/interfaces/$(MX2).F
TARGET_MEX2	:=	$(MX2).mex
# mex function is $(TARGET_MEX2)maci64  for mac OSX
# mex function is $(TARGET_MEX2)a64 	for linux

# mex function to run interactive WIM in matlab
MX3		:=	WIM2d_run_io_mex_vSdir_mesh
MEX3		:=	$(src)/interfaces/$(MX3).F
TARGET_MEX3	:=	$(MX3).mex
# mex function is $(TARGET_MEX3)maci64  for mac OSX
# mex function is $(TARGET_MEX3)a64 	for linux

OBJ_MAIN	:=	$(ODIR)/mod_WIM2d_run.o

OBJS		:=      \
			$(ODIR)/mod_interp.o \
			$(ODIR)/mod_year_info.o \
			$(ODIR)/mod_file_utils.o \
                        $(ODIR)/mod_gridinfo.o \
			$(ODIR)/mod_common_wim.o \
			$(ODIR)/mod_common_wim_arrays.o \
			$(ODIR)/mod_wim_init.o \
			$(ODIR)/mod_atten_young.o \
			$(ODIR)/mod_atten_young_hardcoded.o \
			$(ODIR)/mod_RTparam_utilities.o \
			$(ODIR)/mod_adv_atten.o \
			$(ODIR)/mod_common_advect.o \
			$(ODIR)/mod_advect.o \
			$(ODIR)/mod_wavesice.o

# Create objects by compiling their respective .F files
# (add dependency on header files also, so we recompile if headers change)
$(ODIR)/%.o:	%.f90
	@mkdir -p $(ODIR) $(MDIR)
	$(FC) -c -o $@ $< $(FCFLAGS) $(IFLAGS) $(OFLAGS)
$(ODIR)/%.o:	%.F
	@mkdir -p $(ODIR) $(MDIR)
	$(FC) -c -o $@ $< $(FCFLAGS) $(IFLAGS) $(OFLAGS)

# Link objects to make executable $(TARGET_F)
# (default)
lin="*********************************************************"
exec: $(PROG_F) $(OBJ_MAIN) $(OBJS)
	@echo " "
	@echo $(lin)$(lin)$(lin)
	@mkdir -p $(bin)
	$(FC) -o $(TARGET_F) $^ $(FCFLAGS) $(IFLAGS)
	mv $(TARGET_F) $(bin)
	@echo $(lin)$(lin)$(lin)
	@echo " "

# Link objects to make python module $(TARGET_PY).so
py:	$(PROG_PY) $(OBJ_MAIN) $(OBJS)
	@echo " "
	@echo $(lin)$(lin)$(lin)
	@mkdir -p $(bin)
	$(F2PY) -c -m $(TARGET_PY) $^ $(F2PY_FLAGS)
	mv $(TARGET_PY).so $(bin)
	@echo $(lin)$(lin)$(lin)
	@echo " "

# Link objects to make python module $(TARGET_PY).so
mex1:	$(MEX1) $(OBJ_MAIN) $(OBJS)
	@echo " "
	@echo $(lin)$(lin)$(lin)
	@mkdir -p $(bin)
	$(MXC) $(MEX_FLAGS) $^
	mv $(TARGET_MEX1)* $(bin)
	@echo $(lin)$(lin)$(lin)
	@echo " "

mex2:	$(MEX2) $(OBJ_MAIN) $(OBJS)
	@echo " "
	@echo $(lin)$(lin)$(lin)
	@mkdir -p $(bin)
	$(MXC) $(MEX_FLAGS) $^
	mv $(TARGET_MEX2)* $(bin)
	@echo $(lin)$(lin)$(lin)
	@echo " "

mex3:	$(MEX3) $(OBJ_MAIN) $(OBJS)
	@echo " "
	@echo $(lin)$(lin)$(lin)
	@mkdir -p $(bin)
	$(MXC) $(MEX_FLAGS) $^
	mv $(TARGET_MEX3)* $(bin)
	@echo $(lin)$(lin)$(lin)
	@echo " "

all:	py exec mex1 mex2 mex3

All:	py exec

mex:	mex1 mex2 mex3

# DEPENDENCIES
$(ODIR)/mod_grid_info.o:		$(ODIR)/mod_file_utils.o	grid_info.h
$(ODIR)/mod_wim_init.o:			$(ODIR)/mod_gridinfo.o		grid_info.h \
	                        	$(ODIR)/mod_common_wim.o	$(ODIR)/mod_common_wim_arrays.o
$(ODIR)/mod_common_wim_arrays.o:	$(ODIR)/mod_gridinfo.o		wave_info.h \
	                        	grid_info.h
$(ODIR)/mod_common_advect.o:		$(ODIR)/mod_gridinfo.o		grid_info.h
$(ODIR)/mod_advect.o:			$(ODIR)/mod_gridinfo.o 		grid_info.h \
	                        	$(ODIR)/mod_common_wim.o	$(ODIR)/mod_common_wim_arrays.o \
					$(ODIR)/mod_common_advect.o
$(ODIR)/mod_adv_atten.o:		$(ODIR)/mod_gridinfo.o  	grid_info.h \
	                        	$(ODIR)/mod_common_wim.o	$(ODIR)/mod_common_wim_arrays.o \
					$(ODIR)/mod_advect.o 
$(ODIR)/mod_wavesice.o:			$(ODIR)/mod_gridinfo.o		grid_info.h \
	                        	$(ODIR)/mod_common_wim.o	$(ODIR)/mod_common_wim_arrays.o \
					$(ODIR)/mod_advect.o 		$(ODIR)/mod_adv_atten.o	 \
					$(ODIR)/mod_interp.o            $(ODIR)/mod_atten_young.o \
				
# $(ODIR)/interp.o:		$(ODIR)/r8lib.o
#
$(ODIR)/mod_atten_young.o:		$(ODIR)/mod_atten_young_hardcoded.o	$(ODIR)/mod_common_wim.o
$(ODIR)/mod_atten_young_hardcoded.o:	$(ODIR)/mod_RTparam_utilities.o
#
$(OBJ_MAIN):			$(OBJS)

$(PROG_F):      $(OBJ_MAIN)     $(ODIR)/mod_gridinfo.o
$(PROG_PY): 	grid_info.h	wave_info.h		$(OBJS)
$(MEX1):	grid_info.h     wave_info.h     	$(OBJ_MAIN)
$(MEX2):	grid_info.h     wave_info.h		$(OBJ_MAIN)
$(MEX3):	grid_info.h     wave_info.h		$(OBJ_MAIN)

.PHONY: clean vclean mex1 mex2 mex3

# clean: keep executable, but delete modules and objects
clean:
	rm -f $(ODIR)/*.o $(MDIR)/*.mod

# vclean: delete executable, and delete modules and objects
vclean:
	rm -f	$(ODIR)/*.o $(MDIR)/*.mod \
	   	$(bin)/$(TARGET_F) \
	   	$(bin)/$(TARGET_PY).so \
	   	$(bin)/$(TARGET_MEX1)* \
	   	$(bin)/$(TARGET_MEX2)* \
	   	$(bin)/$(TARGET_MEX3)*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
