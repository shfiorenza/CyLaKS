CC=mpicc
CXX=mpicxx

SRCDIR = src
OBJDIR = obj
BINDIR = bin
SRCEXT = cpp

COMPILE_FLAGS = -std=c++14 -fopenmp
RCOMPILE_FLAGS = -D NDEBUG -O2 -march=native
DCOMPILE_FLAGS = -D DEBUG -Og -g
LINK_FLAGS = -fopenmp 
GSLINCS = `pkg-config --cflags gsl`
GSLLIBS = `pkg-config --libs gsl`
YAMLINCS = -I./libs/yaml-cpp/include
YAMLLIBS = -Wl,-rpath=./libs/yaml-cpp/static -L./libs/yaml-cpp/static -lyaml-cpp
INCLUDES = $(GSLINCS) $(YAMLINCS)
LIBS = $(GSLLIBS) $(YAMLLIBS)

# If we're running on Summit, add appropriate flags
ifeq ($(LOC), Summit)
	COMPILE_FLAGS += -xCORE-AVX2
	LINK_FLAGS += -lm
endif

# Combine compiler and linker flags
ifeq ($(CFG),release)
	export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(RCOMPILE_FLAGS)
	export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS)
else
	export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(DCOMPILE_FLAGS)
	export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS)
endif

SOURCES = $(shell find $(SRCDIR) -name '*.$(SRCEXT)' -printf '%T@\t%p\n' \
			  									| sort -k 1nr | cut -f2-)
SIM_SRC = $(SRCDIR)/sim.cpp
MAIN_SOURCES = $(SIM_SRC)

# These are the common sources
SRCS = $(filter-out $(MAIN_SOURCES) $(EXCLUDE_SOURCES), $(SOURCES))
OBJECTS = $(SRCS:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.o)
SIM_OBJ = $(SIM_SRC:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.o)

.PHONY: dirs
dirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

.PHONY: clean
clean:
	rm -f -r $(OBJDIR)
	rm -f -r $(BINDIR)
	rm -f sim

.PHONY : clean-output
clean-output:
	rm -f *.file

sim: dirs $(BINDIR)/sim; cp $(BINDIR)/sim sim

$(BINDIR)/sim: $(OBJECTS) $(SIM_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS) $(LIBS)

# source file rules
$(OBJDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
