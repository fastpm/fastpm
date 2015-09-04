CPPFLAGS += -I depends/install/include
LDFLAGS += -L depends/install/lib

CPPFLAGS += -I lua
LDFLAGS += -L lua

SOURCES = fastpm.c pmpfft.c pmghosts.c pmpaint.c pmstore.c pm2lpt.c \
		msg.c


fastpm: $(SOURCES:%.c=.objs/%.o)
	mpicc -g -O0 -o fastpm $(SOURCES:%.c=.objs/%.o) \
			$(LDFLAGS) -llua -lgsl -lpfft -lfftw3_mpi -lfftw3 -lm 

include $(SOURCES:%.c=.deps/%.d)

.objs/%.o : %.c
	@if ! [ -d .objs ]; then mkdir .objs; fi
	mpicc -g -O0 -c $(CPPFLAGS) -o $@ $<

.deps/%.d : %.c
	@if ! [ -d .deps ]; then mkdir .deps; fi
	echo -n .objs/ > $@
	mpicc -M $(CPPFLAGS) -o - $< >> $@
