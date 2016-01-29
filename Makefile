# Please create Makefile.local according to Makefile.local.example
CONFIGFILE ?= Makefile.local
include $(CONFIGFILE)
include Makefile.rules

.PHONY: all clean

all:
	@(cd lua; CC="$(CC)" make generic)
	@(cd depends; make "CC=$(CC)" -f Makefile.pfft)
	@(cd bigfile; make "CC=$(CC)" "MPICC=$(CC)")
	@(cd libfastpm; make)
	@(cd src; make)

clean:
	@(cd libfastpm; make clean)
	@(cd src; make clean)

