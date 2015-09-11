
CC=mpicc
DEPCMD = gcc -MG -MT .objs/$(<:%.c=%.o) -MM $(CPPFLAGS)

OPTIMIZE = -O0 -g

CPPFLAGS += -I lua
LDFLAGS += -L lua

GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH = $(GSL_DIR) depends/install

CPPFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LDFLAGS += $(foreach dir, $(DIR_PATH), -L$(dir)/lib) 

SOURCES = fastpm.c pmpfft.c pmghosts.c pmpaint.c pmstore.c pm2lpt.c \
		readparams.c msg.c power.c pmsteps.c pmtimer.c pmio-runpb.c

PFFTLIB = depends/install/lib/libpfft.a
LUALIB = lua/liblua.a

fastpm: $(PFFTLIB) $(LUALIB) $(SOURCES:%.c=.objs/%.o) 
	$(CC) $(OPTIMIZE) -o fastpm $(SOURCES:%.c=.objs/%.o) \
			$(LDFLAGS) -llua -lgsl \
			-lpfft -lfftw3_mpi -lfftw3 \
			-lpfftf -lfftw3f_mpi -lfftw3f \
			-lm 

$(LUALIB): lua/Makefile
	(cd lua; CC="$(CC)" make generic)

$(PFFTLIB): depends/install_pfft.sh
	# FIXME: some configure flags may not work. 
	# shall we adopt autotools?
	@if ! [ -f $(PFFTLIB) ]; then \
		MPICC=$(CC) sh depends/install_pfft.sh $(PWD)/depends/install | tail;\
	else \
		touch $(PFFTLIB); \
	fi;
	
-include $(SOURCES:%.c=.deps/%.d)

.objs/%.o : %.c
	@if ! [ -d .objs ]; then mkdir .objs; fi
	$(CC) $(OPTIMIZE) -c $(CPPFLAGS) -o $@ $<

.deps/%.d : %.c
	@if ! [ -d .deps ]; then mkdir .deps; fi
	@if ! $(DEPCMD) -o $@ $< ; then \
		rm $@; \
	fi;

clean:
	rm -rf .objs
	rm -rf .deps
