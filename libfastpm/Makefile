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
    pmapi.c 

libfastpm.a : $(LIBSOURCES:%.c=.objs/%.o)
	$(AR) rcs $@ $(LIBSOURCES:%.c=.objs/%.o) 

-include $(LIBSOURCES:%.c=.deps/%.d)

clean:
	rm -rf .objs
	rm -rf .deps
	rm -rf libfastpm.a

