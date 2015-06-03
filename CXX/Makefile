PRG_NAME=prog_wim2d.app
SOURCES=\
	main.cpp \
	options.cpp



# CXX = clang++-3.6
# CXXFLAGS += -std=c++0x -std=c++11 -pedantic  -ftemplate-depth-256 -Wno-inline -Wall -fPIC -stdlib=libstdc++ -MMD -MP -lm -v #-lboost_program_options #-v
# # CXXFLAGS += -I headers
# LDFLAGS += -lboost_program_options -lboost_system -lboost_filesystem -I/opt/local/boost-1.58/include/boost # -L/usr/lib -lboost_program_options


CXX = g++-4.8
CXXFLAGS += -std=c++0x -std=c++11 -pedantic  -ftemplate-depth-256 -Wno-inline -fPIC -MMD -MP -lm -v #-Wall
LDFLAGS += -lboost_program_options -lboost_filesystem -lboost_system

#CFLAGS = -O3 -Wall
#LDFLAGS += -lboost_program_options -lboost_system -lboost_filesystem
#LDFLAGS += -I/opt/local/boost-1.58/include -L/opt/local/boost-1.58/lib  # -L/usr/lib -lboost_program_options
#LDLIBS = lboost_program_options -lboost_filesystem -lboost_system
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
# $(PRG_NAME): $(OBJECTS)
# 	$(CXX) $(LDFLAGS) -o $@ $^


$(PRG_NAME) : $(OBJECTS)
	$(CXX) $(CXXLAGS) -o $@ $^ $(LDFLAGS)


# The compilation depends on this Makefile.
$(OBJECTS): Makefile

-include $(DEPS)


#g++ -std=c++11 main.cpp options.cpp -o test -lboost_program_options -lboost_system -lboost_filesystem -v
#g++ main.cpp -o first -lboost_program_options -v
