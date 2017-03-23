# Please create Makefile.local according to Makefile.local.example
CONFIGFILE ?= Makefile.local
include $(CONFIGFILE)
include Makefile.rules

.PHONY: all clean

all:
	@(cd lua; CC="$(CC)" $(MAKE) generic)
	@(cd depends; $(MAKE))
	@(cd libfastpm; $(MAKE))
	@(cd libfastpmio; $(MAKE))
	@(cd src; $(MAKE))
	@(cd tests; $(MAKE))

clean:
	@(cd libfastpm; $(MAKE) clean)
	@(cd libfastpmio; $(MAKE) clean)
	@(cd depends; $(MAKE) clean)
	@(cd src; $(MAKE) clean)
	@(cd tests; $(MAKE) clean)

deep-clean: clean
	(cd lua; make clean)
	(cd depends; rm -rf double single install src)
	(cd depends/bigfile; make clean)
	(cd depends/mpsort; make clean)
