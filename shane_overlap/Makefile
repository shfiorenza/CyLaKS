CC =		g++
CFLAGS =	-O2 -g
LNKFLAGS =
GLXFLAGS =      -L/usr/X11R6/lib/ -lXext -lXmu -lX11 -lGL -lGLU
GSLFLAGS =      -L/usr/local/lib/ -I/usr/local/include/gsl/ -lgsl -lgslcblas
 
FFTFLAGS =      -L/sw/lib/

%o:%cpp
		$(CC) $(CFLAGS) -cpp $<


.PHONEY : clean
clean :
	rm -f *.o

OVERLAP_OBJS =	main.o gfopen.o\
		kMC.o ran3.o\

sim:		$(OVERLAP_OBJS)
		$(CC) $(LNKFLAGS) -o sim $(OVERLAP_OBJS) $(GSLFLAGS) -lm

$(OVERLAP_OBJS):master_header.h constants.h macros.h parameters.h prototypes.h structures.h

BUILD_OBJS =	

build_mts:	$(BUILD_OBJS) build_microtubules.o
		$(CC) $(LNKFLAGS) -o build_mts build_microtubules.o $(BUILD_OBJS) -lm

$(BUILD_OBJS):	master_header.h constants.h macros.h parameters.h prototypes.h structures.h
