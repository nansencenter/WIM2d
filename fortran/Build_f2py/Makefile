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

PROG=../main_py/WIM2d_py.F 

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
     mod_advect.o

DEPS=print_screen_info.h \
     print_log_info.h

TARGET=WIM2d_py

# Create objects by compiling their respective .F files
# (add dependency on header files also, so we recompile if headers change)
%.o:	%.F
	$(FC) -c -o $@ $< $(FCFLAGS) $(IFLAGS) $(OFLAGS)

# Link objects to make executable $(TARGET).so
default: $(PROG) $(OBJS)
	$(F2PY) -c -m $(TARGET) $^ $(F2PY_FLAGS)

# DEPENDENCIES
mod_wim_init.o:		mod_gridinfo.o mod_common_wim.o
mod_common_wim.o:	mod_gridinfo.o
mod_advect.o:		mod_gridinfo.o mod_common_wim.o
mod_adv_atten.o:	mod_gridinfo.o mod_common_wim.o mod_advect.o
#
RTparam_outer.o:        RTparam_fast.o
RTparam_fast.o:         RTparam_hardcoded_v2.o
RTparam_hardcoded_v2.o: RTparam_utilities.o
#
$(PROG): 		$(OBJS) $(DEPS)

.PHONY: clean

# clean: keep executable, but delete modules and objects
clean:
	rm *.o $(MDIR)/*

# vclean: delete executable, and delete modules and objects
vclean:
	rm $(TARGET).so *.o $(MDIR)/*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
