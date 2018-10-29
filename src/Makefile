include ../Makefile.local
include ../Makefile.rules

DIR_PATH = 

CPPFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LDFLAGS += $(foreach dir, $(DIR_PATH), -L$(dir)/lib) 

CPPFLAGS += -I../api/ -I../lua/ -I../depends/install/include
LUA ?= ../lua/lua

FASTPM_SOURCES = fastpm.c param-mpi.c param.c lua-main.c lua-runtime.c runpb.c readgrafic.c lua-config.c stacktrace.c
FASTPM_LUA_SOURCES = fastpm-lua.c param.c lua-main.c lua-runtime.c lua-config.c

FASTPM_FOF_SOURCES = fastpm-fof.c param-mpi.c param.c lua-main.c lua-runtime.c lua-config.c stacktrace.c

SOURCES = $(FASTPM_SOURCES) $(TEST_SOURCES)
.PHONY: all

LUASOURCES = \
	lua-runtime-dump.lua \
	lua-runtime-config.lua \
	lua-runtime-fastpm.lua \

LIBFASTPM_LIBS = ../libfastpmio/libfastpm-io.a ../libfastpm/libfastpm.a ../depends/libfastpm-dep.a

all: fastpm-lua \
     fastpm-fof \
     fastpm

fastpm-lua: $(FASTPM_LUA_SOURCES:%.c=.objs/%.o) ../lua/liblua.a $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ $(LDFLAGS) -lm 

fastpm: $(FASTPM_SOURCES:%.c=.objs/%.o) ../lua/liblua.a $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ $(LDFLAGS) $(GSL_LIBS) -lm

fastpm-fof: $(FASTPM_FOF_SOURCES:%.c=.objs/%.o) ../lua/liblua.a $(LIBFASTPM_LIBS)
	$(CC) $(OPTIMIZE) $(OPENMP) -o $@ $^ $(LDFLAGS) $(GSL_LIBS) -lm

-include $(SOURCES:%.c=.deps/%.d)

lua-runtime.c : $(LUASOURCES) make-lua-runtime.sh
	bash make-lua-runtime.sh $(LUASOURCES) > $@

lua-config.c lua-config.h : $(LUASOURCES) make-lua-config.lua
	$(LUA) make-lua-config.lua lua-config

clean:
	rm -rf .objs
	rm -rf .deps
	rm -rf fastpm fastpm-lua $(TESTS)

