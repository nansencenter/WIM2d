PROGNAME=prog_wim2d.exec

# C++ files
CXXSRCDIR=$(WIM2D_PATH)/CXX

CXX = g++
CXXFLAGS += -std=c++0x -std=c++11 -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread -fopenmp -pedantic  -v # -MMD -MP -lm #-Wall
CXXFLAGS += -I $(CXXSRCDIR)/headers -I . -I $(BOOST_DIR)/include -I $(OPENMPI_INCLUDE_DIR)/ -I $(CXXSRCDIR)

LDFLAGS += -fopenmp
LDFLAGS += -Wl,-rpath,$(OPENMPI_LIB_DIR)/
LDFLAGS += -L $(OPENMPI_LIB_DIR)/ -lmpi_cxx -lmpi -ldl -lstdc++ -lpthread
LDFLAGS += -Wl,-rpath,$(BOOST_DIR)/lib
LDFLAGS += -L $(BOOST_DIR)/lib -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi -lboost_date_time

CC = gcc
CFLAGS += -pedantic -Wno-inline -fPIC -O3 #-v # -MMD -MP -lm #-Wall
CFLAGS += -I $(CXXSRCDIR)/headers -I . -I $(CXXSRCDIR)

EXEC=$(PROGNAME)
OBJECTDIR=objs/
DEPENDDIR=.deps/
BINARYDIR=../bin/

# C++ files
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
OBJS=$(CXXSRC:$(CXXSRCDIR)/%.cpp=$(OBJECTDIR)/%.o)
DEPS=$(CXXSRC:$(CXXSRCDIR)/%.cpp=$(DEPENDDIR)/%.d)

# C files
CCSRCDIR=$(CXXSRCDIR)/src
CCSRC=$(wildcard $(CCSRCDIR)/*.c)

CCOBJS:=$(CCSRC:$(CCSRCDIR)/%.c=$(OBJECTDIR)/%.o)
CCDEPS:=$(CCSRC:$(CCSRCDIR)/%.c=$(DEPENDDIR)/%.d)

# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: $(EXEC)

# Delete the object and dependancy files.
clean:
	$(RM) $(OBJS) $(DEPS) $(CCOBJS) $(CCDEPS)

rmproper: clean
	$(RM) $(BINARYDIR)$(PROGNAME)

$(CCOBJS): $(OBJECTDIR)%.o : $(CCSRCDIR)/%.c
	@mkdir -p $(OBJECTDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Make dependancy rules
$(CCDEPS): $(DEPENDDIR)%.d : $(CCSRCDIR)/%.c
	@mkdir -p $(OBJECTDIR)
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CC) -M $(CFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'

$(PROGNAME): $(OBJS) $(CCOBJS)
	@echo $(OBJS)
	@echo $(CCOBJS)
	@mkdir -p $(BINARYDIR)
	$(CXX) $(CXXLAGS) -o $(BINARYDIR)$@ $^ $(LDFLAGS)

$(OBJECTDIR)%.o: $(CXXSRCDIR)/%.cpp
	@mkdir -p $(OBJECTDIR)
	$(CXX) -o $@ -c $< $(CXXFLAGS)

# Make dependancy rules
$(DEPENDDIR)%.d: %.cpp
	@mkdir -p $(OBJECTDIR)
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'


# The compilation depends on this Makefile.
$(OBJS): Makefile

-include $(DEPS)

-include $(CCDEPS)

# .PHONY: all clean rmproper
