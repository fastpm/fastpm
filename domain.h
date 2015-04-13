
void domain_init(int Nmesh, double BoxSize);
void domain_decompose(Particles* const particles, void * scratch, size_t scratch_bytes);
void domain_finalize();
