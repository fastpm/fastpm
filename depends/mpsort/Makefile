CC?=cc
MPICC?=mpicc
PREFIX=/usr

all: libradixsort.a libmpsort-mpi.a

install: libradixsort.a libmpsort-mpi.a
	install -d $(PREFIX)/lib
	install -d $(PREFIX)/include
	install libradixsort.a $(PREFIX)/lib/libradixsort.a
	install libmpsort-mpi.a $(PREFIX)/lib/libmpsort-mpi.a
	install mpsort.h $(PREFIX)/include/mpsort.h

clean:
	rm -f *.o *.a	
tests: main main-mpi bench-mpi

main: main.c libmpsort-omp.a libradixsort.a
	$(CC) -o main $^
main-mpi: main-mpi.c libmpsort-mpi.a libradixsort.a
	$(MPICC) -o main-mpi $^
bench-mpi: bench-mpi.c libmpsort-mpi.a libradixsort.a
	$(MPICC) -o bench-mpi $^

libradixsort.a: radixsort.c
	$(CC) -c -o radixsort.o radixsort.c
	ar r libradixsort.a radixsort.o
	ranlib libradixsort.a

libmpsort-omp.a: mpsort-omp.c
	$(CC) -c -o mpsort-omp.o mpsort-omp.c
	ar r libmpsort-omp.a mpsort-omp.o
	ranlib libmpsort-omp.a

libmpsort-mpi.a: mpsort-mpi.c
	$(MPICC) -c -o mpsort-mpi.o mpsort-mpi.c
	ar r libmpsort-mpi.a mpsort-mpi.o
	ranlib libmpsort-mpi.a

