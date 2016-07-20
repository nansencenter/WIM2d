# ========================================================
# output names/locations
PROGNAME=prog_wim2d_cpp.exec
OBJECTDIR=objs/
DEPENDDIR=.deps/
BINARYDIR=../bin/
# ========================================================


# ========================================================
# C++ files
CXXSRCDIR=$(WIM2D_PATH)/CXX/src/cpp
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
OBJS=$(CXXSRC:$(CXXSRCDIR)/%.cpp=$(OBJECTDIR)/%.o)
DEPS=$(CXXSRC:$(CXXSRCDIR)/%.cpp=$(DEPENDDIR)/%.d)

CXX = g++
CXXFLAGS += -std=c++0x -std=c++11 -ftemplate-depth-256 -Wno-inline
CXXFLAGS += -fPIC -O3 -pthread -fopenmp -pedantic  -v # -MMD -MP -lm #-Wall
CXXFLAGS += -I . -I $(BOOST_DIR)/include -I $(OPENMPI_INCLUDE_DIR)/ -I $(CXXSRCDIR) -I $(CCSRCDIR)
# ========================================================


# ========================================================
# C files
CCSRCDIR=$(WIM2D_PATH)/CXX/src/c
CCSRC=$(wildcard $(CCSRCDIR)/*.c)
CCOBJS:=$(CCSRC:$(CCSRCDIR)/%.c=$(OBJECTDIR)/%.o)
CCDEPS:=$(CCSRC:$(CCSRCDIR)/%.c=$(DEPENDDIR)/%.d)

CC = gcc
CFLAGS += -pedantic -Wno-inline -fPIC -O3 #-v # -MMD -MP -lm #-Wall
CFLAGS += -I $(CCSRCDIR) -I . -I $(CXXSRCDIR)
# ========================================================


# ========================================================
# Libraries
LDFLAGS += -fopenmp
LDFLAGS += -Wl,-rpath,$(OPENMPI_LIB_DIR)/
LDFLAGS += -L $(OPENMPI_LIB_DIR)/ -lmpi_cxx -lmpi -ldl -lstdc++ -lpthread
LDFLAGS += -Wl,-rpath,$(BOOST_DIR)/lib
LDFLAGS += -L $(BOOST_DIR)/lib -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi -lboost_date_time
# ========================================================


# ========================================================
# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: $(PROGNAME)

# Delete the object and dependancy files.
clean:
	$(RM) $(OBJS) $(DEPS) $(CCOBJS) $(CCDEPS)

# Delete the object and dependancy files, and the executable
rmproper: clean
	$(RM) $(BINARYDIR)$(PROGNAME)
# ========================================================


# ================================================================
# Rules for .c files
$(CCOBJS): $(OBJECTDIR)%.o : $(CCSRCDIR)/%.c
	@mkdir -p $(OBJECTDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Make dependancy rules for .c files
$(CCDEPS): $(DEPENDDIR)%.d : $(CCSRCDIR)/%.c
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CC) -M $(CFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'
# ================================================================


# ================================================================
# Rules for .cpp files
$(OBJECTDIR)%.o: $(CXXSRCDIR)/%.cpp
	@mkdir -p $(OBJECTDIR)
	$(CXX) -o $@ -c $< $(CXXFLAGS)

# Make dependancy rules for .cpp files
$(DEPENDDIR)%.d: $(CXXSRCDIR)/%.cpp
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'
# ================================================================


# ========================================================
# make executable
$(PROGNAME): $(OBJS) $(CCOBJS)
	@mkdir -p $(BINARYDIR)
	$(CXX) $(CXXLAGS) -o $(BINARYDIR)$@ $^ $(LDFLAGS)
# ========================================================


# The compilation depends on this Makefile.
$(OBJS): Makefile

-include $(DEPS)

-include $(CCDEPS)
