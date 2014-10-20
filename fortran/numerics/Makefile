CC	:=	gcc
FC	:=	gfortran
VPATH 	:=	.:../misc

MDIR	:=	../mods

# include modules from ../mods, and put them in ../mods when they are created
# FCFLAGS=-O -I$(MDIR)
FCFLAGS=-g -I$(MDIR)
FCFLAGS2=-J$(MDIR)

PROG=test_waveadv_weno.o 

OBJS=mod_gridinfo.o \
     mod_waveadv_weno.o \
     mod_file_utils.o

TARGET=test_waveadv_weno.exe

# Create objects by compiling their respective .F files
%.o:	%.F
	$(FC) -c -o $@ $< $(FCFLAGS) $(FCFLAGS2)

# Link objects to make executable $(TARGET)
default: $(PROG) $(OBJS)
	$(FC) -o $(TARGET) $^ $(FCFLAGS)

mod_waveadv_weno.o:	mod_gridinfo.o \
			mod_file_utils.o
$(PROG): 		$(OBJS)

.PHONY: clean

clean:
	rm $(TARGET) *.o ../mods/*

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right of ":";
# $@ first on the left of ":";
# $< first on the right of ":";
# The .PHONY rule keeps make from doing something with a
# file named clean.
