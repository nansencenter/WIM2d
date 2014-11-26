# c, fortran, and fortran-to-python linker compiler
CC	:=	gcc
FC	:=	gfortran
F2PY	:=	f2py

# root directory for source code
# (relative to this one)
src	:=	..

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

# fortran program to run WIM
# - executable called by run_WIM2d.sh
#   (cd ../run; ./run_WIM2d.sh)
PROG_F		:=	$(ODIR)/p_WIM2d.o
TARGET_F	:=	WIM2d.exe

# fortran interface program to run WIM in ipython
# - NB f2py needs the full path
# - called by run_WIM2d.py
#   (cd ../run; ipython run_WIM2d.py)
PROG_PY		:=	$(src)/main/WIM2d_run.F
TARGET_PY	:=	WIM2d_f2py

OBJ_MAIN	:=	$(ODIR)/WIM2d_run.o

OBJS		:=	$(ODIR)/mod_gridinfo.o \
			$(ODIR)/mod_wim_init.o \
			$(ODIR)/mod_common_wim.o \
			$(ODIR)/mod_file_utils.o \
			$(ODIR)/mod_wim_prams.o \
			$(ODIR)/RTparam_outer.o \
			$(ODIR)/RTparam_fast.o \
			$(ODIR)/RTparam_hardcoded_v2.o \
			$(ODIR)/RTparam_utilities.o \
			$(ODIR)/mod_adv_atten.o \
			$(ODIR)/mod_advect.o \
			$(ODIR)/mod_wavesice.o

DEPS=print_screen_info.h \
     print_log_info.h


# Create objects by compiling their respective .F files
# (add dependency on header files also, so we recompile if headers change)
$(ODIR)/%.o:	%.F
	$(FC) -c -o $@ $< $(FCFLAGS) $(IFLAGS) $(OFLAGS)

# Link objects to make executable $(TARGET_F)
default: $(PROG_F) $(OBJ_MAIN) $(OBJS)
	$(FC) -o $(TARGET_F) $^ $(FCFLAGS) $(IFLAGS)

# Link objects to make python module $(TARGET_PY).so
py:	$(PROG_PY) $(OBJS)
	$(F2PY) -c -m $(TARGET_PY) $^ $(F2PY_FLAGS)

# DEPENDENCIES
$(ODIR)/mod_wim_init.o:		$(ODIR)/mod_gridinfo.o  $(ODIR)/mod_common_wim.o
$(ODIR)/mod_common_wim.o:	$(ODIR)/mod_gridinfo.o
$(ODIR)/mod_advect.o:		$(ODIR)/mod_gridinfo.o  $(ODIR)/mod_common_wim.o
$(ODIR)/mod_adv_atten.o:	$(ODIR)/mod_gridinfo.o  $(ODIR)/mod_common_wim.o \
				$(ODIR)/mod_advect.o
$(ODIR)/mod_wavesice.o:		$(ODIR)/mod_gridinfo.o	$(ODIR)/mod_common_wim.o \
				$(ODIR)/mod_advect.o 	$(ODIR)/mod_adv_atten.o
#
$(ODIR)/RTparam_outer.o:        $(ODIR)/RTparam_fast.o
$(ODIR)/RTparam_fast.o:         $(ODIR)/RTparam_hardcoded_v2.o
$(ODIR)/RTparam_hardcoded_v2.o: $(ODIR)/RTparam_utilities.o
#
$(OBJ_MAIN):		$(OBJS) $(DEPS)
$(PROG_F): 		$(OBJ_MAIN) $(ODIR)/mod_gridinfo.o
$(PROG_PY): 		$(OBJS) $(DEPS)

.PHONY: clean

# clean: keep executable, but delete modules and objects
clean:
	rm $(ODIR)/* $(MDIR)/*

# vclean: delete executable, and delete modules and objects
vclean:
	rm $(TARGET_F) $(TARGET_PY).so $(ODIR)/* $(MDIR)/*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
