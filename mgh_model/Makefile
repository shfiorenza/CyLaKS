# <variant> library not yet supported on intel compiler. Use gcc/8.2.0 for now.
# Ideally, <variant> should be depricated by polymorphism in the future. 
ifeq ($(LOC), summit)
    CC=gcc # icc
    CXX=g++ # icpc
else
    CC=gcc
    CXX=g++
endif

SRCDIR = src
OBJDIR = obj
BINDIR = bin
SRCEXT = cpp

COMPILE_FLAGS = -std=c++17 
RCOMPILE_FLAGS = -D NDEBUG -O2 -march=native
DCOMPILE_FLAGS = -D DEBUG -O0 -g
LINK_FLAGS =
YAMLINCS = -I./libs/yaml-cpp/include
YAMLLIBS = -Wl,-rpath=./libs/yaml-cpp/static -L./libs/yaml-cpp/static -lyaml-cpp
# If we're running on Summit, add appropriate flags
ifeq ($(LOC), summit)
    GSLINCS = -I/projects/shfi4480/overlap_analysis/mgh_model/libs/gsl/include
    GSLLIBS = -L/projects/shfi4480/overlap_analysis/mgh_model/libs/gsl/lib -lgsl -lgslcblas -lm
#    COMPILE_FLAGS += -xCORE-AVX2
#    LINK_FLAGS += -lm
else
    GSLINCS = -I./libs/gsl/include
    GSLLIBS = -Wl,-rpath=./libs/gsl/static -L./libs/gsl/static -lgsl
    LINK_FLAGS += -static
endif
INCLUDES = $(GSLINCS) $(YAMLINCS) 
LIBS = $(GSLLIBS) $(YAMLLIBS)

# Combine compiler and linker flags
ifeq ($(CFG),release)
    export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(RCOMPILE_FLAGS)
    export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS)
else
    export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(DCOMPILE_FLAGS)
    export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS)
endif

SOURCES = $(shell find $(SRCDIR) -name '*.$(SRCEXT)' -printf '%T@\t%p\n' | sort -k 1nr | cut -f2-)
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
	rm -f *.file *.log

sim: dirs $(BINDIR)/sim; cp $(BINDIR)/sim sim

$(BINDIR)/sim: $(OBJECTS) $(SIM_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS) $(LIBS)

# source file rules
$(OBJDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
