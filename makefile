#CC = mpicc.lam
#CFLAGS = --ansi -pedantic -Wall -O3 -malign-double -fomit-frame-pointer \
#	-ffloat-store -DSSE -fasm
CC = mpcc_r
#CFLAGS = -q64 -qipa -qsrcmsg -O3 -qhot -DMPI -bmaxdata:0x70000000
CFLAGS = -q64 -DMPI -DPARALLELT -bmaxdata:0x70000000
COMPILE = ${CC} ${CFLAGS} 

TARGETS = clover_eo hybrid linsolve xchange \
	expo hybrid_update observables start \
	geometry_eo linalg_eo ranlxs sw io

all: hybrid

${addsuffix .o, ${TARGETS}}: %.o: %.c %.h makefile global.h
	${COMPILE} -c $< -o $@

hybrid: ${addsuffix .o, ${TARGETS}} global.h
	 ${COMPILE} ${addsuffix .o, ${TARGETS}} -o hybrid -lm

clean:
	rm -f *.o hybrid
