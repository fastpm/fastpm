void radix_sort(void * base, size_t nmemb, size_t size,
        void (*radix)(const void * ptr, void * radix, void * arg),
        size_t rsize,
        void * arg);

#ifdef _OPENMP
void mpsort_omp(void * base, size_t nmemb, size_t size,
        void (*radix)(const void * ptr, void * radix, void * arg), 
        size_t rsize, 
        void * arg);

#endif
#ifdef MPI_COMM_WORLD
void mpsort_mpi(void * base, size_t nmemb, size_t elsize,
        void (*radix)(const void * ptr, void * radix, void * arg), 
        size_t rsize, 
        void * arg, MPI_Comm comm);
void mpsort_mpi_newarray(void * base, size_t nmemb, 
        void * out, size_t outnmemb,
        size_t size,
        void (*radix)(const void * ptr, void * radix, void * arg), 
        size_t rsize, 
        void * arg, MPI_Comm comm);
#endif
