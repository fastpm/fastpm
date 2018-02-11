-include ../Makefile.local
-include ../Makefile.rules

.PHONY: all clean
BIGFILE_LIBS = libbigfile.a libbigfile-mpi.a
MPSORT_LIBS = libradixsort.a libmpsort-mpi.a
ifdef OPENMP
PFFT_LIBS = libfftw3.a \
			libfftw3_mpi.a \
			libfftw3_omp.a \
			libfftw3f.a \
			libfftw3f_mpi.a \
			libfftw3f_omp.a \
			libpfft_omp.a \
			libpfftf_omp.a \

else
PFFT_LIBS = libfftw3.a \
			libfftw3_mpi.a \
			libfftw3f.a \
			libfftw3f_mpi.a \
			libpfft.a \
			libpfftf.a \

endif

.PHONY: all

all : libfastpm-dep.a

libfastpm-dep.a: $(BIGFILE_LIBS:%.a=install/lib/%.a) \
	  $(MPSORT_LIBS:%.a=install/lib/%.a) \
	  $(PFFT_LIBS:%.a=install/lib/%.a)
	../ar-join $@ $^

.INTERMEDIATE: .pfft .bigfile .mpsort

$(BIGFILE_LIBS:%.a=install/lib/%.a): .bigfile

.bigfile:
	(cd bigfile; $(MAKE) install "PREFIX=$(PWD)/install" "CC=$(CC) $(OPTIMIZE) $(CFLAGS)" "MPICC=$(CC) $(OPTIMIZE) $(CFLAGS)")

$(MPSORT_LIBS:%.a=install/lib/%.a): .mpsort

.mpsort:
	(cd mpsort; $(MAKE) install "PREFIX=$(PWD)/install" "CC=$(CC) $(OPTIMIZE) $(CFLAGS)" "MPICC=$(CC) $(OPTIMIZE) $(CFLAGS)")

$(PFFT_LIBS:%.a=install/lib/%.a): .pfft

.pfft:
	($(MAKE) "OPENMP=$(OPENMP)" "CC=$(CC)" "CFLAGS=$(OPTIMIZE) $(CFLAGS)" -f Makefile.pfft "PFFT_CONFIGURE_FLAGS=$(PFFT_CONFIGURE_FLAGS)")

clean:
	rm -f libfastpm-dep.a
