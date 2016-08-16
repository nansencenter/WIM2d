# ========================================================
# kernel type
KERNEL := $(shell uname -s)

ifeq ($(KERNEL),Linux)
	LIBRARYNAME=libwim.so.1.0
else
	LIBRARYNAME=libwim.dylib
endif

# ========================================================

# ========================================================
# output names/locations
PROGNAME=prog_wim2d_cpp.exec
OBJECTDIR=objs/
DEPENDDIR=.deps/
BINARYDIR=../bin/
# ========================================================

WIM2D_CXX=$(WIM2D_PATH)/CXX

# ========================================================
# C++ files
CXXSRCDIR=$(WIM2D_CXX)/src
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
OBJS=$(CXXSRC:$(CXXSRCDIR)/%.cpp=$(OBJECTDIR)/%.o)
DEPS=$(CXXSRC:$(CXXSRCDIR)/%.cpp=$(DEPENDDIR)/%.d)

CXX = g++
#CXXFLAGS += -std=c++0x -std=c++11 -ftemplate-depth-256 -Wno-inline
CXXFLAGS += -std=c++11 -pedantic -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread -v
CXXFLAGS += -fPIC -O3 -pthread -fopenmp -pedantic  -v # -MMD -MP -lm #-Wall
CXXFLAGS += -I . -I $(BOOST_DIR)/include -I $(OPENMPI_INCLUDE_DIR)/ -I $(WIM2D_CXX)/include
# ========================================================


# ========================================================
# C files
CCSRCDIR=$(WIM2D_PATH)/CXX/src
CCSRC=$(wildcard $(CCSRCDIR)/*.c)
CCOBJS:=$(CCSRC:$(CCSRCDIR)/%.c=$(OBJECTDIR)/%.o)
CCDEPS:=$(CCSRC:$(CCSRCDIR)/%.c=$(DEPENDDIR)/%.d)

CC = gcc
CFLAGS += -pedantic -Wno-inline -fPIC -O3 #-v # -MMD -MP -lm #-Wall
CFLAGS += -I $(CCSRCDIR) -I . -I $(WIM2D_CXX)/include
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
all: exec #$(PROGNAME)

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
exec: $(OBJS) $(CCOBJS)
	@mkdir -p $(BINARYDIR)
	$(CXX) $(CXXLAGS) -o $(BINARYDIR)$(PROGNAME) $^ $(LDFLAGS)
# ========================================================

# ========================================================
# make wim library

LIBRARYDIR=$(WIM2D_PATH)/CXX/lib/

LDFLAGSLIB = $(LDFLAGS)
ifeq ($(KERNEL),Linux)
	LDFLAGSLIB += -shared -Wl,-soname,libwim.so.1
	#LDFLAGSLIB += -lm -ldl #-lhwloc
else
	LDFLAGSLIB += -dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-current_version,1.0,-install_name,$(LIBRARYDIR)$(LIBRARYNAME)
	#LDFLAGSLIB += -ldl -lstdc++
endif

LOBJS := $(filter-out $(OBJECTDIR)/main.o, $(OBJS))

ifeq ($(KERNEL),Linux)
lib: $(LOBJS) $(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGSLIB)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libwim.so
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libwim.so.1
else
lib: $(LOBJS) $(CCOBJS)
	@mkdir -p $(LIBRARYDIR)
	@echo "============================================"
	@echo Creating shared library $(LIBRARYDIR)$(LIBRARYNAME)
	@echo "============================================"
	$(CXX) $(CXXFLAGS) -o $(LIBRARYDIR)$(LIBRARYNAME) $^ $(LDFLAGSLIB)
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libwim.dylib.0
	@ln -sf $(LIBRARYNAME) $(LIBRARYDIR)libwim.dylib.1
endif
# ========================================================

# The compilation depends on this Makefile.
$(OBJS): Makefile

-include $(DEPS)
-include $(CCDEPS)
