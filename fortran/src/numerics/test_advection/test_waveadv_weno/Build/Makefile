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

# put executables here
bin	:=	../bin

# paths for source files
VPATH 		:=	..:$(numerics):$(src)/misc:$(src)/init:$(HDIR)

# include modules from ../mods, and put them in ../mods when they are created
# FCFLAGS=-O -I$(MDIR) -I$(HDIR)
FCFLAGS=-g -I$(MDIR) -I$(HDIR)
FCFLAGS2=-J$(MDIR)

PROG=$(ODIR)/p_test_advect.o 
PROG_1d=$(ODIR)/p_test_advect_1d.o 

OBJS=\
	$(ODIR)/mod_file_utils.o \
	$(ODIR)/mod_gridinfo.o \
	$(ODIR)/mod_common_advect.o \
	$(ODIR)/mod_advect.o

TARGET2d=test_advect.exec
TARGET1d=test_advect_1d.exec

# Create objects by compiling their respective .F files
$(ODIR)/%.o:	%.F
	mkdir -p $(MDIR) $(ODIR)
	$(FC) -c -o $@ $< $(FCFLAGS) $(FCFLAGS2)

# Link objects to make executable $(TARGET)
# Default is 2d advection code
2d: $(PROG) $(OBJS)
	mkdir -p $(bin)
	$(FC) -o $(TARGET2d) $^ $(FCFLAGS)
	mv $(TARGET2d) $(bin)

1d:	$(PROG_1d) $(OBJS)
	mkdir -p $(bin)
	$(FC) -o $(TARGET1d) $^ $(FCFLAGS)
	mv $(TARGET1d) $(bin)

all:	2d	1d

# dependencies:
$(ODIR)/mod_gridinfo.o:		grid_info.h
$(ODIR)/mod_common_advect.o:	$(ODIR)/mod_gridinfo.o
$(ODIR)/mod_advect.o:		$(ODIR)/mod_gridinfo.o \
				$(ODIR)/mod_file_utils.o \
				$(ODIR)/mod_common_advect.o
$(PROG): 			$(OBJS)
$(PROG_1d): 			$(OBJS)

.PHONY: clean vclean

clean:
	rm objs/* mods/*
vclean:
	rm $(bin)/$(TARGET) objs/* mods/*


# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
