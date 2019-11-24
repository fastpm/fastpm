#ifndef INTERNAL_H
#define INTERNAL_H
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


typedef int (*_compar_fn_t)(const void * r1, const void * r2, size_t rsize);
typedef void (*_bisect_fn_t)(void * r, const void * r1, const void * r2, size_t rsize);

struct crstruct {
    void * base;
    size_t nmemb;
    size_t size;
    size_t rsize;
    void * arg;
    void (*radix)(const void * ptr, void * radix, void * arg);
    _compar_fn_t compar;
    _bisect_fn_t bisect;
};

void _setup_radix_sort(
        struct crstruct *d,
        void * base,
        size_t nmemb,
        size_t size,
        void (*radix)(const void * ptr, void * radix, void * arg),
        size_t rsize,
        void * arg);
#endif
