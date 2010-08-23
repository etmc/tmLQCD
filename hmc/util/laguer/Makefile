CXX=g++
CXXFLAGS=-g -O3
CLNDIR=${HOME}/daten/workdir/cln/

chebyRoot: chebyRoot.C Makefile chebyRoot.H
	${CXX} $< -g -o $@ -I${CLNDIR}/include/ -L${CLNDIR}/lib -lcln -lm

clean:
	rm -f *.o chebyRoot *.dat *.log *~
