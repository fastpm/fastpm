# Please create Makefile.local according to Makefile.local.example
CONFIGFILE ?= Makefile.local
include $(CONFIGFILE)
include Makefile.rules

.PHONY: all clean

all:
	@(cd lua; CC="$(CC)" $(MAKE) generic)
	@(cd depends; $(MAKE) "CC=$(CC)" -f Makefile.pfft)
	@(cd bigfile; $(MAKE) "CC=$(CC)" "MPICC=$(CC)")
	@(cd libfastpm; $(MAKE))
	@(cd src; $(MAKE))

clean:
	@(cd libfastpm; $(MAKE) clean)
	@(cd src; $(MAKE) clean)

deep-clean: clean
	(cd lua; make clean)
	(cd depends; rm -rf double single install src)
	(cd bigfile; make clean)
