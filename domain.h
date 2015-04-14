
void domain_init(int Nmesh, double BoxSize);
void domain_decompose(Particles* const particles, void * scratch, size_t scratch_bytes);
void domain_decompose_min(Snapshot * const particles, void * scratch, size_t scratch_bytes);
void domain_finalize();
int domain_create_ghosts(Particles* const particles, double eps, void * scratch, size_t scratch_bytes);
void domain_annihilate_ghosts(Particles* const particles, int nghosts, float (* f3)[3], 
        void * scratch, size_t scratch_bytes);
