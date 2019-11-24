void radix_sort(void * base, size_t nmemb, size_t size,
        void (*radix)(const void * ptr, void * radix, void * arg),
        size_t rsize,
        void * arg);

#ifdef _OPENMP
/* openmp support */
void mpsort_omp(void * base, size_t nmemb, size_t size,
        void (*radix)(const void * ptr, void * radix, void * arg), 
        size_t rsize, 
        void * arg);

#endif

#ifdef MPI_VERSION
/* MPI support */
#define MPSORT_DISABLE_SPARSE_ALLTOALLV (1 << 1)
#define MPSORT_DISABLE_GATHER_SORT (1 << 3)
#define MPSORT_REQUIRE_GATHER_SORT (1 << 4)
#define MPSORT_REQUIRE_SPARSE_ALLTOALLV (1 << 6)

void mpsort_mpi_set_options(int options);
int mpsort_mpi_has_options(int options);
void mpsort_mpi_unset_options(int options);

void mpsort_mpi_impl(void * base, size_t nmemb, size_t elsize,
        void (*radix)(const void * ptr, void * radix, void * arg),
        size_t rsize,
        void * arg, MPI_Comm comm,
        const int line, const char * file);

#define mpsort_mpi(base, nmemb, elsize, radix, rsize, arg, comm) \
    mpsort_mpi_impl(base, nmemb, elsize, radix, rsize, arg, comm, \
    __LINE__, __FILE__)

void mpsort_mpi_newarray_impl(void * base, size_t nmemb,
        void * out, size_t outnmemb,
        size_t elsize,
        void (*radix)(const void * ptr, void * radix, void * arg),
        size_t rsize,
        void * arg, MPI_Comm comm,
        const int line, const char * file);

#define mpsort_mpi_newarray(base, nmemb, out, outnmemb, elsize, \
    radix, rsize, arg, comm) \
    mpsort_mpi_newarray_impl(base, nmemb, out, outnmemb, elsize, \
    radix, rsize, arg, comm, __LINE__, __FILE__)

void mpsort_mpi_report_last_run();

#ifdef __INTEL_COMPILER
#warning MPSORT: detected an Intel Compiler.
#warning MPSORT: As of Oct 27 2019, icc frequently produces buggier code than gcc when interfacing with MPI and multithreading.
#warning MPSORT: If your application misbehaves, consider recompiling fresh with gcc to rule out compiler bugs.
#endif
#endif
