CC	:=	gcc
FC	:=	gfortran

# path to numerics directory
numerics	:=	../..
src		:=	$(numerics)/..

# modules here
MDIR	:=	mods
# objects here
ODIR	:=	objs
# header files
HDIR	:=	$(src)/header_files
# put executables here
bin	:=	$(numerics)/test_advection/bin

# paths for source files
VPATH 		:=	..:$(numerics):$(src)/misc:$(src)/init:$(HDIR)

# include modules from ../mods, and put them in ../mods when they are created
# FCFLAGS=-O -I$(MDIR) -I$(HDIR)
FCFLAGS=-g -I$(MDIR) -I$(HDIR)
FCFLAGS2=-J$(MDIR)

PROG=$(ODIR)/test_advect.o 
PROG_1d=$(ODIR)/test_advect_1d.o 

OBJS=\
     $(ODIR)/mod_file_utils.o \
     $(ODIR)/mod_gridinfo.o \
     $(ODIR)/mod_advect.o

TARGET=test_advect.exe

# Create objects by compiling their respective .F files
$(ODIR)/%.o:	%.F
	mkdir -p $(MDIR) $(ODIR)
	$(FC) -c -o $@ $< $(FCFLAGS) $(FCFLAGS2)

# Link objects to make executable $(TARGET)
default: $(PROG) $(OBJS)
	mkdir -p $(bin)
	$(FC) -o $(TARGET) $^ $(FCFLAGS)
	mv $(TARGET) $(bin)

1d:	$(PROG_1d) $(OBJS)
	mkdir -p $(bin)
	$(FC) -o $(TARGET) $^ $(FCFLAGS)
	mv $(TARGET) $(bin)

# dependencies:
$(ODIR)/mod_gridinfo.o:	grid_info.h
$(ODIR)/mod_advect.o:	$(ODIR)/mod_gridinfo.o \
			$(ODIR)/mod_file_utils.o
$(PROG): 		$(OBJS)
$(PROG_1d): 		$(OBJS)

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
