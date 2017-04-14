include ../Makefile.local
include ../Makefile.rules

GSL_DIR   ?= 

DIR_PATH = $(GSL_DIR)

CPPFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LDFLAGS += $(foreach dir, $(DIR_PATH), -L$(dir)/lib) 

CPPFLAGS += -I../api/ -I../lua/ -I../bigfile -I../depends/install -I../mpsort

TEST_SOURCES = testpm.c \
               testconstrained.c \
               testlightcone.c \
               testrecorder.c \
               testgradient-c2r.c \
               testgradient-paint.c \
               testgradient-readout.c \
               testgradient-r2c.c \
               testgradient-force.c \
               testgradient-force-simple.c \

TESTS = $(TEST_SOURCES:%.c=%)

LIBFASTPM_LIBS = ../libfastpmio/libfastpm-io.a ../libfastpm/libfastpm.a ../depends/libfastpm-dep.a

SOURCES = $(TEST_SOURCES)
.PHONY: all

all: $(TESTS)

testpm : .objs/testpm.o $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

testconstrained : .objs/testconstrained.o $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

testlightcone : .objs/testlightcone.o $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

testrecorder: .objs/testrecorder.o $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

testgradient-c2r: .objs/testgradient-c2r.o $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

testgradient-r2c: .objs/testgradient-r2c.o $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

testgradient-paint: .objs/testgradient-paint.o $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

testgradient-readout: .objs/testgradient-readout.o $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

testgradient-force: .objs/testgradient-force.o $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

testgradient-force-simple: .objs/testgradient-force-simple.o $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

-include $(SOURCES:%.c=.deps/%.d)

clean:
	rm -rf .objs
	rm -rf .deps
	rm -rf $(TESTS)

