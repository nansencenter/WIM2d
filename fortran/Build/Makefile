CC	:=	gcc
FC	:=	gfortran

# set path for compiler
VPATH 	:=	.:../main:../misc:../init:../numerics:../attenuation_youngs

# directory where modules are put into
MDIR	:=	mods

# header files for "include" statements
HDIR	:=	../header_files

# include modules from $(MDIR), and put them in $(MDIR) when they are created
# FCFLAGS=-O -I$(MDIR)
FCFLAGS=-g -I$(MDIR) -I$(HDIR)
FCFLAGS2=-J$(MDIR)

PROG=p_WIM2d.o 

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

TARGET=WIM2d.exe

# Create objects by compiling their respective .F files
%.o:	%.F
	$(FC) -c -o $@ $< $(FCFLAGS) $(FCFLAGS2)

# Link objects to make executable $(TARGET)
default: $(PROG) $(OBJS)
	$(FC) -o $(TARGET) $^ $(FCFLAGS)

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
$(PROG): 		$(OBJS)

.PHONY: clean

# clean: keep executable, but delete modules and objects
clean:
	rm *.o $(MDIR)/*

# vclean: delete executable, and delete modules and objects
vclean:
	rm $(TARGET) *.o $(MDIR)/*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
