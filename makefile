CC = mpicc.lam
CFLAGS = --ansi -pedantic -Wall -O3 -malign-double -fomit-frame-pointer \
	-ffloat-store -DSSE -fasm
COMPILE = ${CC} ${CFLAGS} 

TARGETS = clover_eo hybrid linsolve xchange \
	expo hybrid_update observables start \
	geometry_eo linalg_eo ranlxs sw

all: hybrid

${addsuffix .o, ${TARGETS}}: %.o: %.c %.h Makefile global.h
	${COMPILE} -c $< -o $@

hybrid: ${addsuffix .o, ${TARGETS}} global.h
	 ${COMPILE} ${addsuffix .o, ${TARGETS}} -o hybrid -lm

clean:
	rm -f *.o hybrid
