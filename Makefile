include Makefile.local
include Makefile.rules

.PHONY: all clean

all:
	@(cd lua; CC="$(CC)" make generic)
	@(make -f Makefile.pfft)
	@(cd libfastpm; make)
	@(cd src; make)

clean:
	@(cd libfastpm; make clean)
	@(cd src; make clean)

