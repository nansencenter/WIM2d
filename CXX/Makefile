PROGNAME=prog_wim2d.exec

CXX = g++-4.8
CXXFLAGS += -std=c++0x -std=c++11 -pedantic  -ftemplate-depth-256 -Wno-inline -fPIC -MMD -MP -lm -v #-Wall
CXXFLAGS += -I headers -I . -I $(BOOST_DIR)/include

LDFLAGS += -Wl,-rpath,$(BOOST_DIR)/lib
LDFLAGS += -L $(BOOST_DIR)/lib -lboost_program_options -lboost_filesystem -lboost_system


CC = gcc-4.8
CFLAGS += -std=c++0x -pedantic -Wno-inline -fPIC -MMD -MP -lm -v #-Wall
CFLAGS += -I headers -I .

# SRC = \
# 	main.cpp \
# 	options.cpp


# CSRCDIR=src
# CCSRC=\
# 	$(CSRCDIR)/RTparam_outer.c \
# 	$(CSRCDIR)/RTparam_fast.c \
# 	$(CSRCDIR)/RTparam_hardcoded.c


EXEC=$(PROGNAME)
OBJECTDIR=objs/
BINARYDIR=bin/

# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
OBJS=$(CXXSRC:.cpp=.o)
DEPS=$(CXXSRC:.cpp=.d)

OBJS:=$(addprefix $(OBJECTDIR), $(OBJS))
DEPS:=$(addprefix $(OBJECTDIR), $(DEPS))

# C files
CCSRCDIR=src
CCSRC=$(wildcard $(CCSRCDIR)/*.c)
CCOBJS=$(CCSRC:.c=.o)
CCDEPS=$(CCSRC:.c=.d)

CCOBJS:=$(addprefix $(OBJECTDIR), $(CCOBJS))
CCDEPS:=$(addprefix $(OBJECTDIR), $(CCDEPS))



# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: $(EXEC)

# Delete the object files.
clean:
	$(RM) $(OBJS) $(DEPS) $(CCOBJS) $(CCDEPS)

mrproper: clean
	$(RM) $(BINARYDIR)$(PROGNAME)


default: $(CCOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJECTDIR)%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)



$(PROGNAME): $(OBJS) $(CCOBJS)
	@mkdir -p $(BINARYDIR)
	$(CXX) $(CXXLAGS) -o $(BINARYDIR)$@ $^ $(LDFLAGS)

$(OBJECTDIR)%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)


# The compilation depends on this Makefile.
$(OBJS): Makefile

-include $(DEPS)


# .PHONY: all clean mrproper

# # Delete the object files.
# clean:
# 	$(RM) -rf $(OBJECTDIR)*

# # Delete the executable
# mrproper: clean
# 	$(RM) -f $(BINARYDIR)$(EXEC)


# hello: $(DIR_OBJS)
# 	$(CXX) $(CXXLAGS) -o $@ $^ $(LDFLAGS)

# $(DIR_OBJ)%.o: %.cpp
# 	$(CXX) -o $@ -c $< $(CFLAGS)

#$(CXX) -o $@ $^ $(LDFLAGS)
