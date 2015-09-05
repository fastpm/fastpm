CPPFLAGS += -I depends/install/include
LDFLAGS += -L depends/install/lib

CPPFLAGS += -I lua
LDFLAGS += -L lua

SOURCES = fastpm.c pmpfft.c pmghosts.c pmpaint.c pmstore.c pm2lpt.c \
		readparams.c msg.c power.c


fastpm: $(SOURCES:%.c=.objs/%.o)
	mpicc -g -O0 -o fastpm $(SOURCES:%.c=.objs/%.o) \
			$(LDFLAGS) -llua -lgsl \
			-lpfft -lfftw3_mpi -lfftw3 \
			-lpfftf -lfftw3f_mpi -lfftw3f \
			-lm 

-include $(SOURCES:%.c=.deps/%.d)

.objs/%.o : %.c
	@if ! [ -d .objs ]; then mkdir .objs; fi
	mpicc -g -O0 -c $(CPPFLAGS) -o $@ $<

.deps/%.d : %.c
	@if ! [ -d .deps ]; then mkdir .deps; fi
	@if mpicc -M $(CPPFLAGS) -o - $< > $@.tmp ; then \
		(echo -n .objs/; cat $@.tmp ;) > $@ ; \
	fi;
	@rm $@.tmp
