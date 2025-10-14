TARGETS = scalar_prod_r_test

USESF = yes

OS = -os3

# gcc shouldn't see this options, that's why we don't use CGLAGS here
NLCCFLAGS = -D_STD_C99_COMPLEX_CHECKED -D_STD_C99_COMPLEX -Dapenext
INCLUDES = -I../
# workaround to let nlcc not see the non-standard complex.h
NLCCINCLUDES = -I${NROOT}/include/nlibc/ ${INCLUDES}

NLCCOPTS = -gp ${NLCCFLAGS} ${NLCCINCLUDES}
ifdef USESF
  MPPOPTS = -sf -v
  SHAKEROPTS = -n -z 
else 
  MPPOPTS = -v
  SHAKEROPTS = +a -z 
endif
SOFANOPTS = --rr

# needed due to a bug in nlcc
NLCCOS = -OS3

NLCC = nlcc-0.5.2
MPP = mpp
SOFAN = sofan
SHAKER = shaker
M4 = m4
CCDEP = gcc
DEPFLAGS = -MM -MQ $*.sasm ${CFLAGS} ${INCLUDES}

DEPFILES = $(addsuffix .d, ${TARGETS})
MEMFILES = $(addsuffix .mem, ${TARGETS}) $(addsuffix -sofan.mem, ${TARGETS}) \
	   $(addsuffix .no, ${TARGETS}) $(addsuffix -sofan.no, ${TARGETS})
ASMFILES = $(addsuffix .sasm, ${TARGETS}) $(addsuffix .masm, ${TARGETS}) $(addsuffix -sofan.masm, ${TARGETS})
NCDFILES = $(addsuffix .ncd, ${TARGETS}) $(addsuffix -sofan.ncd, ${TARGETS})
SFOUTFILES = $(addsuffix .svn-out, ${TARGETS}) $(addsuffix .svn-out%, ${TARGETS}) \
             $(addsuffix .sf_log, ${TARGETS}) $(addsuffix .sf_log%, ${TARGETS}) \
             $(addsuffix .sf_log0, ${TARGETS}) $(addsuffix .sf_log0%, ${TARGETS}) \
             $(addsuffix .err-sf, ${TARGETS}) $(addsuffix .svn-out, ${TARGETS}) \
             $(addsuffix .dmo, ${TARGETS}) \
	     $(addsuffix -sofan.svn-out, ${TARGETS}) $(addsuffix -sofan.svn-out%, ${TARGETS}) \
             $(addsuffix -sofan.sf_log, ${TARGETS}) $(addsuffix -sofan.sf_log%, ${TARGETS}) \
             $(addsuffix -sofan.sf_log0, ${TARGETS}) $(addsuffix -sofan.sf_log0%, ${TARGETS}) \
             $(addsuffix -sofan.err-sf, ${TARGETS}) $(addsuffix -sofan.svn-out, ${TARGETS}) \
             $(addsuffix -sofan.dmo, ${TARGETS})
GCCBINARIES = $(addsuffix .gccbin, ${TARGETS})

all: $(addsuffix -sofan.mem, ${TARGETS})
allgcc:  $(addsuffix .gccbin, ${TARGETS})

-include $(DEPFILES)

%.mem: %.masm
	${SHAKER} ${SHAKEROPTS} $<

%.masm: %.sasm
	${MPP} ${OS} ${MPPOPTS} $<

%-sofan.masm: %.masm
	${SOFAN} ${SOFANOPTS} $< $@

%.sasm: %.c Makefile
	${NLCC} ${NLCCOPTS} ${NLCCOS} -S $<

%.ncd: %.mem
	dispminit $< > $@

%-sofan.perf: %-sofan.ncd
	nperf -asm=$*.sasm -c -l -a $< > $@ || (rm -f $@; exit 1)

# beware, this is not very general
%.gccbin: %.c
	gcc -I../ $< -o $@

$(DEPFILES): %.d: %.c Makefile
	$(CCDEP) ${DEPFLAGS} ${INCLUDES} $< > $@

clean:
	rm -f ${ASMFILES} ${MEMFILES} ${NCDFILES} ${GCCBINARIES}

distclean: clean
	rm -f ${DEPFILES} ${SFOUTFILES}

.SECONDARY:
.DELETE_ON_ERROR:
