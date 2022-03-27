.PHONY: all

CC?=cc
PREFIX?=/usr
PIC?=-fPIC

all: libchealpix.a

install: libchealpix.a
	install -d $(PREFIX)
	install -d $(PREFIX)/lib
	install -d $(PREFIX)/include
	install -d $(PREFIX)/include/chealpix
	install libchealpix.a $(PREFIX)/lib/libchealpix.a
	install chealpix.h $(PREFIX)/include/chealpix/chealpix.h
	
libchealpix.a: chealpix.c chealpix.h
	$(CC) $(CFLAGS) $(PIC) -c -o chealpix.o chealpix.c
	ar r $@ chealpix.o
	ranlib $@

clean:
	rm -f libchealpix.a *.o
