SRCDIR = src
OBJDIR = obj
BINDIR = bin
SRCEXT = cpp

COMPILE_FLAGS = -std=c++11
RCOMPILE_FLAGS = -D NDEBUG -O2 -march=native -g
DCOMPILE_FLAGS = -D DEBUG -O0 -g
LINK_FLAGS = 
RLINK_FLAGS = 
DLINK_FLAGS =

ifeq ($(THREADING),omp)
	COMPILE_FLAGS += -fopenmp
	LINK_FLAGS += -fopenmp
endif

ifeq ($(THREADING),iomp)
	COMPILE_FLAGS += -openmp -parallel-source-info=2
	LINK_FLAGS += -openmp -parallel-source-info=2
endif

ifeq ($(NOGRAPH),true)
	COMPILE_FLAGS += -D NOGRAPH
	LINK_FLAGS += -D NOGRAPH
endif

UNAME_S:=$(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	CC = gcc-5
	CXX = g++-5
	GSLINCS = `pkg-config --cflags gsl`
	GSLLIBS = `pkg-config --libs gsl`
	FFTWINCS = `pkg-config --cflags fftw3`
	GLFW3INCS = -I/opt/X11/include -I/usr/X11R6/include
	GLFW3LIBS = -L/opt/X11/lib -L/usr/local/Cellar/glfw/3.2.1/lib -lglfw3 -framework OpenGL -lglew 
	YAMLINCS = -I./Lib/yaml-cpp/include
	YAMLLIBS = -Wl,-rpath,./Lib/yaml-cpp/osxlib -L./Lib/yaml-cpp/osxlib/ -lyaml-cpp
	PNGLIBS = `/usr/local/bin/libpng-config --ldflags`
	INCLUDES = $(GLFW3INCS) $(YAMLINCS) $(GSLINCS) $(FFTWINCS)
	LIBS = $(PNGLIBS) $(GLFW3LIBS) $(YAMLLIBS) $(GSLLIBS)
else ifeq ($(LOC),Janus)
	GSLINCS = `pkg-config --cflags gsl`
	GSLLIBS = `pkg-config --libs gsl`
	GLFW3INCS = 
	GLFW3LIBS = 
	FFTWINCS = `pkg-config --cflags fftw3`
	YAMLINCS = -I./januslib/yaml-cpp/include
	YAMLLIBS = -Wl,-rpath,./januslib/yaml-cpp/static -L./januslib/yaml-cpp/static -lyaml-cpp
	PNGLIBS = `/usr/bin/libpng-config --ldflags`
	INCLUDES = $(GLFW3INCS) $(YAMLINCS) $(GSLINCS) $(FFTWINCS)
	LIBS = $(GLFW3LIBS) $(YAMLLIBS) $(GSLLIBS) $(PNGLIBS)
else ifeq ($(LOC),Summit)
	GSLINCS = `pkg-config --cflags gsl`
	GSLLIBS = `pkg-config --libs gsl`
	GLFW3INCS = 
	GLFW3LIBS = 
	FFTWINCS = `pkg-config --cflags fftw3`
	YAMLINCS = -I./Lib/yaml-cpp/include
	YAMLLIBS = -Wl,-rpath,./Lib/yaml-cpp/summitlib -L./Lib/yaml-cpp/summitlib/ -lyaml-cpp
	PNGLIBS = `/usr/bin/libpng-config --ldflags`
	INCLUDES = $(GLFW3INCS) $(YAMLINCS) $(GSLINCS) $(FFTWINCS)
	LIBS = $(GLFW3LIBS) $(YAMLLIBS) $(GSLLIBS) $(PNGLIBS)
else ifeq ($(LOC),Pando)
	GSLINCS = `pkg-config --cflags gsl`
	GSLLIBS = `pkg-config --libs gsl`
	GLFW3INCS = 
	GLFW3LIBS = 
	FFTWINCS = `pkg-config --cflags fftw3`
	YAMLINCS = -I./pandolib/yaml-cpp/include
	YAMLLIBS = -Wl,-rpath,./pandolib/yaml-cpp/static -L./pandolib/yaml-cpp/static -lyaml-cpp
	PNGLIBS = `/usr/bin/libpng-config --ldflags`
	INCLUDES = $(GLFW3INCS) $(YAMLINCS) $(GSLINCS) $(FFTWINCS)
	LIBS = $(GLFW3LIBS) $(YAMLLIBS) $(GSLLIBS) $(PNGLIBS)
else
	GSLINCS = -I/usr/local/include
	GSLLIBS = -L/usr/local/libs -lgsl -lgslcblas -lm
	YAMLINCS = -I./libs/yaml-cpp/include
	YAMLLIBS = -Wl,-rpath,./libs/yaml-cpp/static -L./libs/yaml-cpp/static -lyaml-cpp
	PNGLIBS = `/usr/bin/libpng-config --ldflags`
	INCLUDES = $(GSLINCS) $(YAMLINCS)
	LIBS = $(GSLLIBS) $(YAMLLIBS) 
endif

print-%: ; @echo $*=$($*)

SHELL = /bin/bash

.SUFFIXES:

ifneq ($(LIBS),)
	COMPILE_FLAGS += $(shell pkg-config --cflags $(LIBS))
	LINK_FLAGS += $(shell pkg-config --libs $(LIBS))
endif

# Special stuff for intel compiler
CC=$(CXX)
ifeq ($(CC),icpc)
	COMPILE_FLAGS += -Wno-deprecated
	#RCOMPILE_FLAGS += -openmp -DBOB_OMP
else
	COMPILE_FLAGS += -Wno-deprecated-declarations -Wno-deprecated
	#RCOMPILE_FLAGS += -fopenmp -DBOB_OMP
endif

# Combine compiler and linker flags
ifeq ($(CFG),release)
	export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(RCOMPILE_FLAGS)
	export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS) $(RLINK_FLAGS)
else
	export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(DCOMPILE_FLAGS)
	export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS) $(DLINK_FLAGS)
endif

# build information on all sources
ifeq ($(UNAME_S),Darwin)
	SOURCES = $(shell find $(SRCDIR) -name '*.$(SRCEXT)' | sort -k 1nr | cut -f2-)
else
	SOURCES = $(shell find $(SRCDIR) -name '*.$(SRCEXT)' -printf '%T@\t%p\n' \
			  										| sort -k 1nr | cut -f2-)
endif

# fallback case
rwildcard = $(foreach d, $(wildcard $1*), $(call rwildcard,$d/,$2) \
									$(filter $(subst *,%,$2), $d))
ifeq ($(SOURCES),)
	SOURCES := $(call rwildcard, $(SRCDIR), *.$(SRCEXT))
endif

# Now we have to figure out which we are building of the program list, since that matters
# for things like not compiling more than one main
SIM_SRC = $(SRCDIR)/sim.cpp

MAIN_SOURCES = $(SIM_SRC)

# These are the common sources
SRCS = $(filter-out $(MAIN_SOURCES) $(EXCLUDE_SOURCES), $(SOURCES))

OBJECTS = $(SRCS:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.o)
SIM_OBJ = $(SIM_SRC:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.o)
DEPS = $(OBJECTS:.o=.d)

.PHONY: dirs
dirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

.PHONY: clean
clean:
	rm -f -r $(OBJDIR)
	rm -f -r $(BINDIR)

.PHONY : clean-output
clean-output :
	rm -f *.file sim

sim: dirs $(BINDIR)/sim; cp $(BINDIR)/sim sim

$(BINDIR)/sim: $(OBJECTS) $(SIM_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS) $(LIBS)

# add dependencies
-include $(DEPS)

# source file rules
$(OBJDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
