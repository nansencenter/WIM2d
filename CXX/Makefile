PRG_NAME=prog_wim2d.app
SOURCES=\
	main.cpp \
	options.cpp


CXX = clang++-3.6
CXXFLAGS += -std=c++0x -std=c++11 -pedantic  -ftemplate-depth-256 -Wno-inline -Wall -fPIC -stdlib=libstdc++ -MMD -MP -v #-lboost_program_options #-v
# CXXFLAGS += -I headers
LDFLAGS += -lboost_program_options -lboost_system -lboost_filesystem -I/opt/local/boost-1.58/include/boost # -L/usr/lib -lboost_program_options


# CXX = g++
# CXXFLAGS += -std=c++0x -std=c++11 -pedantic -static -pthread -lboost_thread -ftemplate-depth-256 -Wno-inline -Wall -fPIC  #-lboost_program_options #-v
# LDFLAGS = -lboost_program_options -lboost_filesystem -lboost_system
# LDLIBS += -lboost_program_options
#LDFLAGS += -I/opt/local/boost-1.58/include -L/opt/local/boost-1.58/lib  # -L/usr/lib -lboost_program_options

# -include $(all_objs:.o=.d)

DEPS=$(patsubst %.cpp,%.d,$(SOURCES))
OBJECTS=$(patsubst %.cpp,%.o,$(SOURCES))

# Rules to always execute.
.PHONY: all clean mrproper

# Default action.
all: $(PRG_NAME)

# Delete the object files.
clean:
	$(RM) $(OBJECTS) $(DEPS)

# Delete the object files and the compiled program.
mrproper: clean
	$(RM) $(PRG_NAME)

# Linkage.
$(PRG_NAME): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

# The compilation depends on this Makefile.
$(OBJECTS): Makefile

-include $(DEPS)


#g++ main.cpp -o first -lboost_program_options -v
