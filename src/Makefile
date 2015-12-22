include ../Makefile.local
include ../Makefile.rules

GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH = $(GSL_DIR) ../depends/install

CPPFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LDFLAGS += $(foreach dir, $(DIR_PATH), -L$(dir)/lib) 

CPPFLAGS += -I../api/ -I../lua/


SOURCES = main.c power.c readparams.c runpb.c readgrafic.c 
#
# the following files are derived and will be purged by make.
SOURCES += lua-preface.c

.PHONY: all

all: fastpm testlib 

fastpm: $(LUALIB) $(SOURCES:%.c=.objs/%.o) ../libfastpm/libfastpm.a
	$(CC) $(OPTIMIZE) $(OPENMP) -o fastpm $(SOURCES:%.c=.objs/%.o) ../libfastpm/libfastpm.a \
			 ../lua/liblua.a $(LDFLAGS) -lgsl -lgslcblas \
			-lpfft_omp -lfftw3_mpi -lfftw3_omp -lfftw3 \
			-lpfftf_omp -lfftw3f_mpi -lfftw3f_omp -lfftw3f \
			-lm 

testlib : testlib.o
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ testlib.o \
			../libfastpm/libfastpm.a \
			$(LDFLAGS) \
			-lgsl -lgslcblas \
			-lpfft_omp -lfftw3_mpi -lfftw3_omp -lfftw3 \
			-lpfftf_omp -lfftw3f_mpi -lfftw3f_omp -lfftw3f \
			-lm 

%.c : %.lua
	xxd -i $^ > $@

clean:
	rm -rf .objs
	rm -rf .deps

