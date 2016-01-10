include ../Makefile.local
include ../Makefile.rules

GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl

DIR_PATH = $(GSL_DIR) ../depends/install

CPPFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LDFLAGS += $(foreach dir, $(DIR_PATH), -L$(dir)/lib) 

CPPFLAGS += -I../api/ -I../lua/


FASTPM_SOURCES = main.c power.c readparams.c runpb.c readgrafic.c 
TEST_SOURCES = test2lpt.c testpm.c

SOURCES = $(FASTPM_SOURCES) $(TEST_SOURCES)
.PHONY: all

all: fastpm test2lpt testpm

fastpm: $(LUALIB) $(FASTPM_SOURCES:%.c=.objs/%.o) .objs/lua-dump.o .objs/lua-preface.o ../libfastpm/libfastpm.a
	$(CC) $(OPTIMIZE) $(OPENMP) -o fastpm $(FASTPM_SOURCES:%.c=.objs/%.o) .objs/lua-dump.o .objs/lua-preface.o ../libfastpm/libfastpm.a \
			 ../lua/liblua.a $(LDFLAGS) -lgsl -lgslcblas \
			-lpfft_omp -lfftw3_mpi -lfftw3_omp -lfftw3 \
			-lpfftf_omp -lfftw3f_mpi -lfftw3f_omp -lfftw3f \
			-lm 

test2lpt : test2lpt.o
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ test2lpt.o \
			../libfastpm/libfastpm.a \
			$(LDFLAGS) \
			-lgsl -lgslcblas \
			-lpfft_omp -lfftw3_mpi -lfftw3_omp -lfftw3 \
			-lpfftf_omp -lfftw3f_mpi -lfftw3f_omp -lfftw3f \
			-lm 

testpm : testpm.o
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ testpm.o \
			../libfastpm/libfastpm.a \
			$(LDFLAGS) \
			-lgsl -lgslcblas \
			-lpfft_omp -lfftw3_mpi -lfftw3_omp -lfftw3 \
			-lpfftf_omp -lfftw3f_mpi -lfftw3f_omp -lfftw3f \
			-lm 

-include $(SOURCES:%.c=.deps/%.d)

%.c : %.lua
	xxd -i $^ > $@

clean:
	rm -rf .objs
	rm -rf .deps
	rm -rf fastpm testpm test2lpt

