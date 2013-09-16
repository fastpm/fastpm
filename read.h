#ifndef READ_H
#define READ_H 1

#include "particle.h"

int read_snapshot(const char filename[], Snapshot* snapshot, void* buf, size_t size);

#endif
