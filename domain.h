
void domain_init(int Nmesh, double BoxSize);
void domain_decompose(Particles * particles);
void domain_finalize();
void domain_wrap(Particles * particles);
int domain_create_ghosts(Particles * particles, double eps);
void domain_annihilate_ghosts(Particles * particles, int nghosts, float (* f3)[3]);
