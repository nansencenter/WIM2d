# c, fortran, and fortran-to-python linker compiler
CC	:=	gcc
FC	:=	gfortran
F2PY	:=	f2py

# set path for compiler
VPATH 	:=	.:../main:../misc:../init:../numerics:../attenuation_youngs:../header_files

# directory where modules are put into
MDIR	:=	mods

# header files for "include" statements
HDIR	:=	../header_files

# # optimise code
# FCFLAGS=-O
# check for errors
FCFLAGS=-g
# include modules/headers from $(MDIR) or $(HDIR)
IFLAGS=-I$(MDIR) -I$(HDIR)
# Put modules in $(MDIR) when they are created
OFLAGS=-J$(MDIR)
F2PY_FLAGS=--fcompiler=$(FC) --f77flags="$(IFLAGS) $(OFLAGS)"

# fortran program to run WIM
# - executable called by run_WIM2d.sh
#   (cd ../run; ./run_WIM2d.sh)
PROG_F=p_WIM2d.o
TARGET_F=WIM2d.exe

# fortran interface program to run WIM in ipython
# - called by run_WIM2d.py
#   (cd ../run; ipython run_WIM2d.py)
PROG_PY=../main/WIM2d_run.F
TARGET_PY=WIM2d_f2py

OBJ_MAIN=WIM2d_run.o

OBJS=mod_gridinfo.o \
     mod_wim_init.o \
     mod_common_wim.o \
     mod_file_utils.o \
     mod_wim_prams.o \
     RTparam_outer.o \
     RTparam_fast.o \
     RTparam_hardcoded_v2.o \
     RTparam_utilities.o \
     mod_adv_atten.o \
     mod_advect.o \
     mod_wavesice.o

DEPS=print_screen_info.h \
     print_log_info.h


# Create objects by compiling their respective .F files
# (add dependency on header files also, so we recompile if headers change)
%.o:	%.F
	$(FC) -c -o $@ $< $(FCFLAGS) $(IFLAGS) $(OFLAGS)

# Link objects to make executable $(TARGET_F)
default: $(PROG_F) $(OBJ_MAIN) $(OBJS)
	$(FC) -o $(TARGET_F) $^ $(FCFLAGS) $(IFLAGS)

py:	$(PROG_PY) $(OBJS)
	$(F2PY) -c -m $(TARGET_PY) $^ $(F2PY_FLAGS)

# DEPENDENCIES
mod_wim_init.o:		mod_gridinfo.o  mod_common_wim.o
mod_common_wim.o:	mod_gridinfo.o
mod_advect.o:		mod_gridinfo.o  mod_common_wim.o
mod_adv_atten.o:	mod_gridinfo.o  mod_common_wim.o \
   			mod_advect.o
mod_wavesice.o:		mod_gridinfo.o	mod_common_wim.o \
   			mod_advect.o 	mod_adv_atten.o
#
RTparam_outer.o:        RTparam_fast.o
RTparam_fast.o:         RTparam_hardcoded_v2.o
RTparam_hardcoded_v2.o: RTparam_utilities.o
#
$(OBJ_MAIN):		$(OBJS) $(DEPS)
$(PROG_F): 		$(OBJ_MAIN) mod_gridinfo.o
$(PROG_PY): 		$(OBJS) $(DEPS)

.PHONY: clean

# clean: keep executable, but delete modules and objects
clean:
	rm *.o $(MDIR)/*

# vclean: delete executable, and delete modules and objects
vclean:
	rm $(TARGET_F) $(TARGET_PY).so *.o $(MDIR)/*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
