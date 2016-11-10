-include ../Makefile.local
-include ../Makefile.rules

CPPFLAGS += -I ../api -I ../mpsort/

DIR_PATH = $(GSL_DIR) ../depends/install

CPPFLAGS += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LDFLAGS += $(foreach dir, $(DIR_PATH), -L$(dir)/lib) 
CPPFLAGS += -I../api/ -I../depends/install

FILES = $(shell git ls-files)
FILES += $(shell git ls-files ../api)

LIBSOURCES = io.c \

all: libfastpm-io.a 

libfastpm-io.a : $(LIBSOURCES:%.c=.objs/%.o)
	$(AR) rc $@ $^

version.c : $(FILES)
	GIT=`git describe --always --dirty --abbrev=10`; \
	echo 'const char * LIBFASTPM_IO_VERSION = "1.0.'$${GIT}'";' > $@

-include $(LIBSOURCES:%.c=.deps/%.d)

clean:
	rm -rf .objs
	rm -rf .deps
	rm -rf libfastpm.a

