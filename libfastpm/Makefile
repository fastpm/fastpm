-include ../Makefile.local
-include ../Makefile.rules

CPPFLAGS += -I ../api

DIR_PATH = $(GSL_DIR) ../depends/install

CPPFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LDFLAGS += $(foreach dir, $(DIR_PATH), -L$(dir)/lib) 

LIBSOURCES = libfastpm.c \
    solver-2lpt.c \
    solver-pm.c \
    solver-pm-steps.c \
    solver-pm-grav.c \
    solver-pm-model.c \
    solver-pm-model-pt.c \
    solver-pm-model-linear.c \
    solver-pm-model-pm.c \
    transfer.c \
    utils.c \
    utils-ps.c \
    utils-ic.c \
    cosmology.c \
    logging.c \
    prof.c \
    vpm.c  \
    pmpfft.c  \
    pmghosts.c  \
    pmpaint.c \
    pmstore.c \
    pm2lpt.c \
    pmapi.c \
	hmc.c

PFFTLIBS0 = libpfft_omp.a libpfftf_omp.a
PFFTLIBS0 += libfftw3_omp.a libfftw3_mpi.a libfftw3.a
PFFTLIBS0 += libfftw3f_omp.a libfftw3f_mpi.a libfftw3f.a
# prefix the path
PFFTLIBS = $(PFFTLIBS0:%.a=../depends/install/lib/%.a)

all: libfastpm.a 

libfastpm.a : $(LIBSOURCES:%.c=.objs/%.o)
	$(AR) rc $@ $^
	../ar-join $@ $(PFFTLIBS)

#libfastpm.so : libfastpm.a $(PFFTLIBS)
#	$(CC) -shared -fvisibility=hidden -o $@ $^ 

-include $(LIBSOURCES:%.c=.deps/%.d)

clean:
	rm -rf .objs
	rm -rf .deps
	rm -rf libfastpm.a

