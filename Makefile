CC =		gcc
CFLAGS =	-O2 -g
LNKFLAGS = 	
GLXFLAGS =      -L/usr/X11R6/lib/ -lXext -lXmu -lX11 -lGL -lGLU
GSLFLAGS =      -L/usr/local/lib/ -I/usr/local/include/gsl/ -lgsl -lgslcblas
FFTFLAGS =      -L/sw/lib/

%o:%c
		$(CC) $(CFLAGS) -c $<

.PHONEY : clean
clean :
	rm -f *.o

OVERLAP_OBJS =	main.o parse_parameters.o parse_tokens.o\
		overlap_growth.o output.o\
		gfopen.o allocate.o ran3.o error_exit.o

overlap:	$(OVERLAP_OBJS)
		$(CC) $(LNKFLAGS) -o overlap $(OVERLAP_OBJS) $(GSLFLAGS) -lm

$(OVERLAP_OBJS):hskuan.h constants.h macros.h parameters.h prototypes.h

BUILD_OBJS =	gfopen.o allocate.o error_exit.o

build_micro_hskuan:	$(BUILD_OBJS) build_microtubule.o
		$(CC) $(LNKFLAGS) -o build_micro_hskuan build_microtubule.o $(BUILD_OBJS) -lm

$(BUILD_OBJS):	hskuan.h constants.h macros.h parameters.h prototypes.h

ANA_OBJS =	gfopen.o allocate.o error_exit.o

analysis:	$(ANA_OBJS) analysis.o
		$(CC) $(LNKFLAGS) -o analysis analysis.o $(ANA_OBJS) -lm
