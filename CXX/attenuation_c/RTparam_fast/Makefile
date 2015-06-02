CC=gcc
CFLAGS=-I. # include directory = "."
DEPS=RTparam_hardcoded.h RTparam_fast.h # dependencies
OBJ=test_RTparam_fast.o RTparam_fast.o RTparam_hardcoded.o 
TARGET=a.out

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

default: $(OBJ)
	$(CC) -o $(TARGET) $^ $(CFLAGS)

.PHONY: clean

clean:
	rm $(TARGET) *.o

# -I. include headers from "." directory;
# -L. include libraries from "." directory;
# $^ everything to the right;
# $@ first on the left;
# $< first on the right;
# The .PHONY rule keeps make from doing something with a
# file named clean.
