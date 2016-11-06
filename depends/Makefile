-include ../Makefile.local
-include ../Makefile.rules

.PHONY: all clean
BIGFILE_LIBS = libbigfile.a libbigfile-mpi.a
MPSORT_LIBS = libradixsort.a libmpsort-mpi.a
PFFT_LIBS = libfftw3.a \
			libfftw3_mpi.a \
			libfftw3_omp.a \
			libfftw3f.a \
			libfftw3f_mpi.a \
			libfftw3f_omp.a \
			libpfft_omp.a \
			libpfftf_omp.a \

.PHONY: all

all : libfastpm-dep.a

libfastpm-dep.a: $(BIGFILE_LIBS:%.a=install/lib/%.a) \
	  $(MPSORT_LIBS:%.a=install/lib/%.a) \
	  $(PFFT_LIBS:%.a=install/lib/%.a)
	../ar-join $@ $^

$(BIGFILE_LIBS:%.a=install/lib/%.a):
	(cd bigfile; $(MAKE) install "PREFIX=$(PWD)/install" "CC=$(CC)" "MPICC=$(CC)")

$(MPSORT_LIBS:%.a=install/lib/%.a):
	(cd mpsort; $(MAKE) install "PREFIX=$(PWD)/install" "CC=$(CC)" "MPICC=$(CC)")

$(PFFT_LIBS:%.a=install/lib/%.a):
	($(MAKE) "CC=$(CC)" -f Makefile.pfft)

