CC	:=	gcc
FC	:=	gfortran

# path to numerics directory
numerics	:=	../../..
src		:=	$(numerics)/..

# modules here
MDIR	:=	mods
# objects here
ODIR	:=	objs
# header files (location of grid_info.h)
# - put it in this dir for easier testing
HDIR	:=	.
# HDIR	:=	$(src)/header_files

# put executables here
bin	:=	../bin

# paths for source files
VPATH 		:=	..:$(numerics):$(src)/misc:$(src)/init:$(HDIR)

# include modules from ../mods, and put them in ../mods when they are created
# FCFLAGS=-O -I$(MDIR) -I$(HDIR)
FCFLAGS=-g -I$(MDIR) -I$(HDIR)
FCFLAGS2=-J$(MDIR)

PROG=$(ODIR)/p_test_advect.o 

OBJS=\
	$(ODIR)/mod_file_utils.o \
	$(ODIR)/mod_gridinfo.o \
	$(ODIR)/mod_common_advect.o \
	$(ODIR)/mod_remap.o \
	$(ODIR)/mod_advect_ice.o

TARGET=test_advect.exec

# Create objects by compiling their respective .F files
$(ODIR)/%.o:	%.F
	mkdir -p $(MDIR) $(ODIR)
	$(FC) -c -o $@ $< $(FCFLAGS) $(FCFLAGS2)

# Link objects to make executable $(TARGET)
# Default is 2d advection code
2d: $(PROG) $(OBJS)
	mkdir -p $(bin)
	$(FC) -o $(TARGET) $^ $(FCFLAGS)
	mv $(TARGET) $(bin)

# dependencies:
$(ODIR)/mod_gridinfo.o:		grid_info.h
$(ODIR)/mod_common_advect.o:	$(ODIR)/mod_gridinfo.o
$(ODIR)/mod_advect.o:		$(ODIR)/mod_gridinfo.o \
				$(ODIR)/mod_file_utils.o \
				$(ODIR)/mod_common_advect.o \
				$(ODIR)/mod_remap.o
$(ODIR)/mod_remap.o:		$(ODIR)/mod_gridinfo.o
$(PROG): 			$(OBJS)

.PHONY: clean vclean

clean:
	rm objs/*.o mods/*
vclean:
	rm $(bin)/$(TARGET) objs/* mods/*


# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
