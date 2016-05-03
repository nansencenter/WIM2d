PROGNAME=prog_wim2d.exec

CXX = g++
CXXFLAGS += -std=c++0x -std=c++11 -ftemplate-depth-256 -Wno-inline -fPIC -O3 -pthread -fopenmp -pedantic  -v # -MMD -MP -lm #-Wall
CXXFLAGS += -I headers -I . -I $(BOOST_DIR)/include

LDFLAGS += -fopenmp
#LDFLAGS += -Wl,-rpath,$(OPENMPI_LIB_DIR)/
#LDFLAGS += -L $(OPENMPI_LIB_DIR)/ -lmpi_cxx -lmpi -ldl -lstdc++ -lpthread

LDFLAGS += -Wl,-rpath,$(BOOST_DIR)/lib
LDFLAGS += -L $(BOOST_DIR)/lib -lboost_program_options -lboost_filesystem -lboost_system -lboost_mpi -lboost_date_time

CC = gcc
CFLAGS += -pedantic -Wno-inline -fPIC -O3 #-v # -MMD -MP -lm #-Wall
CFLAGS += -I headers -I .

EXEC=$(PROGNAME)
OBJECTDIR=objs/
DEPENDDIR=.deps/
BINARYDIR=bin/

# C++ files
CXXSRCDIR=.
CXXSRC=$(wildcard $(CXXSRCDIR)/*.cpp)
OBJS=$(CXXSRC:.cpp=.o)
DEPS=$(CXXSRC:.cpp=.d)

OBJS:=$(addprefix $(OBJECTDIR), $(OBJS))
DEPS:=$(addprefix $(DEPENDDIR), $(DEPS))

# C files
CCSRCDIR=src
CCSRC=$(wildcard $(CCSRCDIR)/*.c)

CCOBJS:=$(CCSRC:$(CCSRCDIR)/%.c=$(OBJECTDIR)%.o)
CCDEPS:=$(CCSRC:$(CCSRCDIR)/%.c=$(DEPENDDIR)%.d)

# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: $(EXEC)

# Delete the object and dependancy files.
clean:
	$(RM) $(OBJS) $(DEPS) $(CCOBJS) $(CCDEPS)

mrproper: clean
	$(RM) $(BINARYDIR)$(PROGNAME)

$(CCOBJS): $(OBJECTDIR)%.o : $(CCSRCDIR)/%.c
	@mkdir -p $(OBJECTDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Make dependancy rules
$(CCDEPS): $(DEPENDDIR)%.d : $(CCSRCDIR)/%.c
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CC) -M $(CFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'

$(PROGNAME): $(OBJS) $(CCOBJS)
	@mkdir -p $(BINARYDIR)
	$(CXX) $(CXXLAGS) -o $(BINARYDIR)$@ $^ $(LDFLAGS)

$(OBJECTDIR)%.o: %.cpp
	@mkdir -p $(OBJECTDIR)
	$(CXX) -o $@ -c $< $(CXXFLAGS)

# Make dependancy rules
$(DEPENDDIR)%.d: %.cpp
	@mkdir -p $(DEPENDDIR)
	@$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< | sed "s^$*.o^$(OBJECTDIR)$*.o^" > $@'


# The compilation depends on this Makefile.
$(OBJS): Makefile

-include $(DEPS)

-include $(CCDEPS)


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


# OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

# $(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
#     $(CC) $(CFLAGS) -c $< -o $@
#     @echo "Compiled "$<" successfully!"

# $(BINDIR)/$(TARGET): $(OBJECTS)
#     $(LINKER) $@ $(LFLAGS) $(OBJECTS)
#     @echo "Linking complete!"
