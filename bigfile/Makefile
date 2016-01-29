CC=mpiicc
MPICC=mpicc
PREFIX=/usr
all: libbigfile.a libbigfile-mpi.a

install: libbigfile.a libbigfile-mpi.a
	install libbigfile.a $(PREFIX)/lib/libbigfile.a
	install libbigfile-mpi.a $(PREFIX)/lib/libbigfile-mpi.a
	install bigfile-mpi.h $(PREFIX)/include/bigfile-mpi.h
	install bigfile.h $(PREFIX)/include/bigfile-mpi.h

bigfile.o: bigfile.c bigfile.h
	$(CC) $(CFLAGS) -o $@ -c bigfile.c
libbigfile.a: bigfile.o
	ar r $@ bigfile.o
	ranlib $@
bigfile-mpi.o: bigfile-mpi.c bigfile-mpi.h
	$(MPICC) $(CFLAGS) -o $@ -c bigfile-mpi.c
libbigfile-mpi.a: bigfile-mpi.o
	ar r $@ bigfile-mpi.o
	ranlib $@
clean: 
	rm *.a *.o
