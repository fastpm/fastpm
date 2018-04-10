.PHONY: all
CC?=cc
PREFIX?=/usr
PIC?=-fPIC

all: libkdcount.a

install: libkdcount.a 
	install -d $(PREFIX)
	install -d $(PREFIX)/lib
	install -d $(PREFIX)/include
	install -d $(PREFIX)/include/kdcount
	install libkdcount.a $(PREFIX)/lib/libkdcount.a
	install kdtree.h $(PREFIX)/include/kdcount/kdtree.h

kd_enum.o: kd_enum.c kdtree.h
	$(CC) $(CFLAGS) $(PIC) -o $@ -c kd_enum.c

kd_fof.o: kd_fof.c kdtree.h
	$(CC) $(CFLAGS) $(PIC) -o $@ -c kd_fof.c

kdtree.o: kdtree.c kdtree.h
	$(CC) $(CFLAGS) $(PIC) -o $@ -c kdtree.c

libkdcount.a: kdtree.o kd_enum.o kd_fof.o
	ar r $@ $^
	ranlib $@

clean:
	rm -f *.a *.o
