include ../Makefile.local
include ../Makefile.rules

GSL_DIR   ?= 

DIR_PATH = $(GSL_DIR)

CPPFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LDFLAGS += $(foreach dir, $(DIR_PATH), -L$(dir)/lib) 

CPPFLAGS += -I../api/ -I../lua/ -I../bigfile -I../depends/install 


FASTPM_SOURCES = fastpm.c power.c readparams.c runpb.c readgrafic.c io.c
TEST_SOURCES = test2lpt.c testpm.c

SOURCES = $(FASTPM_SOURCES) $(TEST_SOURCES)
.PHONY: all
FASTPM_LUASOURCES = lua-dump.lua lua-preface.lua

all: fastpm-lua fastpm test2lpt testpm

fastpm-lua: fastpm-lua.o readparams.o $(FASTPM_LUASOURCES:%.lua=.objs/%.o) \
     ../lua/liblua.a ../libfastpm/libfastpm.a
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
		$(LDFLAGS) -lm 

fastpm: $(FASTPM_SOURCES:%.c=.objs/%.o) $(FASTPM_LUASOURCES:%.lua=.objs/%.o) \
     ../lua/liblua.a ../libfastpm/libfastpm.a ../bigfile/libbigfile.a ../bigfile/libbigfile-mpi.a
	$(CC) $(OPTIMIZE) $(OPENMP) -o fastpm $^ \
		$(LDFLAGS) -lgsl -lgslcblas -lm

test2lpt : test2lpt.o ../libfastpm/libfastpm.a
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

testpm : testpm.o ../libfastpm/libfastpm.a
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

-include $(SOURCES:%.c=.deps/%.d)

%.c : %.lua
	xxd -i $^ > $@

clean:
	rm -rf .objs
	rm -rf .deps
	rm -rf fastpm fastpm-lua testpm test2lpt

