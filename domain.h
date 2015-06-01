
void domain_init(double boxsize, int nc);
void domain_set_size(int pm_nc_factor);
void domain_decompose(Particles * particles);
void domain_free();
void domain_wrap(Particles * particles);
int domain_create_ghosts(Particles * particles, double eps);
void domain_annihilate_ghosts(Particles * particles, int nghosts, float (* f3)[3]);
