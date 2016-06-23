#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "internal.h"

/*****
 * msort.c is stripped from glibc git. 
 * qsort_r is not available till 2.8 but
 * most systems has older glibc 
 * ****/
#include "stdlib/msort.c"


/****
 * sort by radix;
 * internally this uses qsort_r of glibc.
 *
 **** */

void radix_sort(void * base, size_t nmemb, size_t size, 
        void (*radix)(const void * ptr, void * radix, void * arg), 
        size_t rsize, 
        void * arg) {
    struct crstruct d;
    _setup_radix_sort(&d, size, radix, rsize, arg);
    mpsort_qsort_r(base, nmemb, size, _compute_and_compar_radix, &d);
}


/* implementation ; internal */
int _compute_and_compar_radix(const void * p1, const void * p2, void * arg) {
    struct crstruct * d = arg;
    char r1[d->rsize], r2[d->rsize];
    d->radix(p1, r1, d->arg);
    d->radix(p2, r2, d->arg);
    int c1 = d->compar(r1, r2, d->rsize);
    return c1;
}


#define DEFTYPE(type) \
static int _compar_radix_ ## type ( \
        const type * u1,  \
        const type * u2,  \
        void * junk) { \
    return (signed) (*u1 > *u2) - (signed) (*u1 < *u2); \
} \
static void _bisect_radix_ ## type ( \
        type * u, \
        const type * u1,  \
        const type * u2,  \
        void * junk) { \
    *u = *u1 + ((*u2 - *u1) >> 1); \
}
DEFTYPE(uint16_t)
DEFTYPE(uint32_t)
DEFTYPE(uint64_t)

static int _compar_radix(const void * r1, const void * r2, size_t rsize, int dir) {
    int i;
    /* from most significant */
    const unsigned char * u1 = r1;
    const unsigned char * u2 = r2;
    if(dir < 0) {
        u1 += rsize - 1;
        u2 += rsize - 1;;
    }
    for(i = 0; i < rsize; i ++) {
        if(*u1 < *u2) return -1;
        if(*u1 > *u2) return 1;
        u1 += dir;
        u2 += dir;
    }
    return 0;
}
static int _compar_radix_u8(const void * r1, const void * r2, size_t rsize, int dir) {
    int i;
    /* from most significant */
    const uint64_t * u1 = r1;
    const uint64_t * u2 = r2;
    if(dir < 0) {
        u1 = (const uint64_t *) ((const char*) u1 + rsize - 8);
        u2 = (const uint64_t *) ((const char*) u2 + rsize - 8);
    }
    for(i = 0; i < rsize; i += 8) {
        if(*u1 < *u2) return -1;
        if(*u1 > *u2) return 1;
        u1 += dir;
        u2 += dir;
    }
    return 0;
}
static int _compar_radix_le(const void * r1, const void * r2, size_t rsize) {
    return _compar_radix(r1, r2, rsize, -1);
}
static int _compar_radix_be(const void * r1, const void * r2, size_t rsize) {
    return _compar_radix(r1, r2, rsize, +1);
}
static int _compar_radix_le_u8(const void * r1, const void * r2, size_t rsize) {
    return _compar_radix_u8(r1, r2, rsize, -1);
}
static int _compar_radix_be_u8(const void * r1, const void * r2, size_t rsize) {
    return _compar_radix_u8(r1, r2, rsize, +1);
}
static void _bisect_radix(void * r, const void * r1, const void * r2, size_t rsize, int dir) {
    int i;
    const unsigned char * u1 = r1;
    const unsigned char * u2 = r2;
    unsigned char * u = r;
    unsigned int carry = 0;
    /* from most least significant */
    for(i = 0; i < rsize; i ++) {
        unsigned int tmp = (unsigned int) *u2 + *u1 + carry;
        if(tmp >= 256) carry = 1;
        else carry = 0;
        *u = tmp;
        u -= dir;
        u1 -= dir;
        u2 -= dir;
    }
    u += dir;
    for(i = 0; i < rsize; i ++) {
        unsigned int tmp = *u + carry * 256;
        carry = tmp & 1;
        *u = (tmp >> 1) ;
        u += dir;
    }
}
static void _bisect_radix_le(void * r, const void * r1, const void * r2, size_t rsize) {
    _bisect_radix(r, r1, r2, rsize, -1);
}
static void _bisect_radix_be(void * r, const void * r1, const void * r2, size_t rsize) {
    _bisect_radix(r, r1, r2, rsize, +1);
}

void _setup_radix_sort(
        struct crstruct *d,
        size_t size,
        void (*radix)(const void * ptr, void * radix, void * arg), 
        size_t rsize, 
        void * arg) {
    const char deadbeef[] = "deadbeef";
    const uint32_t * ideadbeef = (uint32_t *) deadbeef;
    d->rsize = rsize;
    d->arg = arg;
    d->radix = radix;
    d->size = size;
    switch(rsize) {
        case 2:
            d->compar = (_compar_fn_t) _compar_radix_uint16_t;
            d->bisect = (_bisect_fn_t) _bisect_radix_uint16_t;
            break;
        case 4:
            d->compar = (_compar_fn_t) _compar_radix_uint32_t;
            d->bisect = (_bisect_fn_t) _bisect_radix_uint32_t;
            break;
        case 8:
            d->compar = (_compar_fn_t) _compar_radix_uint64_t;
            d->bisect = (_bisect_fn_t) _bisect_radix_uint64_t;
            break;
        default:
            if(ideadbeef[0] != 0xdeadbeef) {
                if(rsize % 8 == 0) {
                    d->compar = _compar_radix_le_u8;
                } else{
                    d->compar = _compar_radix_le;
                }
                d->bisect = _bisect_radix_le;
            } else {
                if(rsize % 8 == 0) {
                    d->compar = _compar_radix_be_u8;
                } else{
                    d->compar = _compar_radix_be;
                }
                d->bisect = _bisect_radix_be;
            }
    }
}

