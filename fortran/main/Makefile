CC	:=	gcc
# FC	:=	gfortran
FC	:=	gfortran-mp-4.8
VPATH 	:=	.:../misc:../init:../numerics:../attenuation_youngs

MDIR	:=	../mods

# include modules from $(MDIR), and put them in $(MDIR) when they are created
# FCFLAGS=-O -I$(MDIR)
FCFLAGS=-g -I$(MDIR)
FCFLAGS2=-J$(MDIR)

PROG=p_WIM2d.o 

OBJS=mod_gridinfo.o \
     mod_wim_init.o \
     mod_common_wim.o \
     mod_file_utils.o \
     mod_wim_prams.o

# mod_waveadv_weno.o

TARGET=WIM2d.exe

# Create objects by compiling their respective .F files
%.o:	%.F
	$(FC) -c -o $@ $< $(FCFLAGS) $(FCFLAGS2)

# Link objects to make executable $(TARGET)
default: $(PROG) $(OBJS)
	$(FC) -o $(TARGET) $^ $(FCFLAGS)

# DEPENDENCIES
# mod_waveadv_weno.o:	mod_gridinfo.o \
# 			mod_file_utils.o
mod_wim_init.o:		mod_gridinfo.o mod_common_wim.o
mod_common_wim.o:	mod_gridinfo.o
$(PROG): 		$(OBJS)

.PHONY: clean

# clean: keep executable, but delete modules and objects
clean:
	rm *.o ../mods/*

# vclean: delete executable, and delete modules and objects
vclean:
	rm $(TARGET) *.o ../mods/*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
