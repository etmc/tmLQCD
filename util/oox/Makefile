CC=gcc
CXX=g++
CFLAGS=-O2 -fexpensive-optimizations -fomit-frame-pointer # -mfpmath=sse -msse2 
LIBS=-lm
OBJECTS_OOX=oox.o
INCLUDE=-I./


# variables for oox_ga executable
# if you want to compile with ga lib support
# please adjust the GALIBPATH variable
# to the toplevel dir of galib
# it is assumed that you compiled the library
# such that a libga.a file is present in the 
# ./ga subdir of galib
GALIBPATH=/usr1/scratch/annube/galib247
LIBS_GA=${LIBS} -L${GALIBPATH}/ga -lga
CFLAGS_GA=${CFLAGS} -DWITHGALIB
INCLUDE_GA=${INCLUDE} -I${GALIBPATH}
OBJECTS_OOX_GA=oox_ga.o oox_gawrapper.o


all: oox oox_ga

oox: ${OBJECTS_OOX} Makefile
	${CXX} ${OBJECTS_OOX} -o $@ ${CFLAGS} ${LIBS}

oox_ga: ${OBJECTS_OOX_GA} Makefile
	${CXX} ${OBJECTS_OOX_GA} -o $@ ${CFLAGS_GA} ${LIBS_GA}

oox_gawrapper.o: oox_gawrapper.cxx
	${CXX} ${CFLAGS_GA} -o $@ -c $< ${INCLUDE_GA}

oox_ga.o: oox.c
	${CC} ${CFLAGS_GA} -o $@ -c $< ${INCLUDE_GA}

clean:
	rm oox oox_ga *.o

.SUFFIXES:

%.o: %.c
	${CC} ${CFLAGS}	-o $@ -c $< ${INCLUDE}

%.o: %.cxx
	${CXX} ${CFLAGS} -o $@ -c $< ${INCLUDE}
