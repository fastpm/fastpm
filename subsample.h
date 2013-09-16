#ifndef SUBSAMPLE_H
#define SUBSAMPLE_H 1

//void write_subsample(const char filename[], const int fac, Snapshot const * const snapshot, void* const mem, size_t mem_size);

void subsample_init(const double subsample_factor, const unsigned int seed);

void write_random_sabsample(const char filename[], Snapshot const * const snapshot, void* const mem, size_t mem_size);
#endif
