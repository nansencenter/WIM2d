CC=clang-3.6
CFLAGS=-I. # include directory = "."
CFLAGS += -I/usr/include -lm
DEPS=RTparam_hardcoded.h RTparam_fast.h # dependencies
OBJ=RTparam_outer.o RTparam_fast.o RTparam_hardcoded.o
TARGET=test_atten.app

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
