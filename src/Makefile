include ../Makefile.local
include ../Makefile.rules

GSL_DIR   ?= 

DIR_PATH = $(GSL_DIR)

CPPFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LDFLAGS += $(foreach dir, $(DIR_PATH), -L$(dir)/lib) 

CPPFLAGS += -I../api/ -I../lua/ -I../bigfile -I../depends/install -I../mpsort
LUA ?= ../lua/lua

FASTPM_SOURCES = fastpm.c lua-main.c lua-runtime.c runpb.c readgrafic.c io.c lua-config.c
FASTPM_LUA_SOURCES = fastpm-lua.c lua-main.c lua-runtime.c lua-config.c
TEST_SOURCES = testpm.c testlightcone.c

SOURCES = $(FASTPM_SOURCES) $(TEST_SOURCES)
.PHONY: all

LUASOURCES = \
	lua-runtime-dump.lua \
	lua-runtime-config.lua \
	lua-runtime-fastpm.lua \

all: fastpm-lua fastpm testpm testlightcone

fastpm-lua: $(FASTPM_LUA_SOURCES:%.c=.objs/%.o) \
     ../lua/liblua.a ../libfastpm/libfastpm.a
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
		$(LDFLAGS) -lm 

fastpm: $(FASTPM_SOURCES:%.c=.objs/%.o)\
     ../lua/liblua.a ../libfastpm/libfastpm.a ../bigfile/libbigfile.a ../bigfile/libbigfile-mpi.a
	$(CC) $(OPTIMIZE) $(OPENMP) -o fastpm $^ \
		$(LDFLAGS) -lgsl -lgslcblas -lm

testpm : .objs/testpm.o ../libfastpm/libfastpm.a
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

testlightcone : .objs/testlightcone.o .objs/io.o ../libfastpm/libfastpm.a ../bigfile/libbigfile.a ../bigfile/libbigfile-mpi.a
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ \
	    $(LDFLAGS) -lgsl -lgslcblas -lm

-include $(SOURCES:%.c=.deps/%.d)

lua-runtime.c : $(LUASOURCES) make-lua-runtime.sh
	bash make-lua-runtime.sh $(LUASOURCES) > $@

lua-config.c lua-config.h : $(LUASOURCES) make-lua-config.lua
	$(LUA) make-lua-config.lua lua-config

clean:
	rm -rf .objs
	rm -rf .deps
	rm -rf fastpm fastpm-lua testpm testlightcone test2lpt

