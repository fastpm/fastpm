#ifndef MOVE_H
#define MOVE_H 1

#include <stdlib.h>
#include "particle.h"

void move_particles2(Particles* const particles, const float BoxSize, void* const buf, const size_t size);

#endif
