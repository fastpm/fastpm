/*
 * returns index of the last item satisfying
 * [item] < P,
 *
 * returns -1 if [all] < P
 * */

static ptrdiff_t _bsearch_last_lt(void * P,
    void * base, size_t nmemb,
    struct crstruct * d) {

    if (nmemb == 0) return -1;

    unsigned char tmpradix[d->rsize];
    ptrdiff_t left = 0;
    ptrdiff_t right = nmemb - 1;

    d->radix((char*) base, tmpradix, d->arg);
    if(d->compar(tmpradix, P, d->rsize) >= 0) {
        return - 1;
    }
    d->radix((char*) base + right * d->size, tmpradix, d->arg);
    if(d->compar(tmpradix, P, d->rsize) < 0) {
        return nmemb - 1;
    }

    /* left <= i <= right*/
    /* [left] < P <= [right] */
    while(right > left + 1) {
        ptrdiff_t mid = ((right - left + 1) >> 1) + left;
        d->radix((char*) base + mid * d->size, tmpradix, d->arg);
        /* if [mid] < P , move left to mid */
        /* if [mid] >= P , move right to mid */
        int c1 = d->compar(tmpradix, P, d->rsize);
        if(c1 < 0) {
            left = mid;
        } else {
            right = mid;
        }
    }
    return left;
}

/*
 * returns index of the last item satisfying
 * [item] <= P,
 *
 * */
static ptrdiff_t _bsearch_last_le(void * P,
    void * base, size_t nmemb,
    struct crstruct * d) {

    if (nmemb == 0) return -1;

    unsigned char tmpradix[d->rsize];
    ptrdiff_t left = 0;
    ptrdiff_t right = nmemb - 1;

    d->radix((char*) base, tmpradix, d->arg);
    if(d->compar(tmpradix, P, d->rsize) > 0) {
        return -1;
    }
    d->radix((char*) base + right * d->size, tmpradix, d->arg);
    if(d->compar(tmpradix, P, d->rsize) <= 0) {
        return nmemb - 1;
    }

    /* left <= i <= right*/
    /* [left] <= P < [right] */
    while(right > left + 1) {
        ptrdiff_t mid = ((right - left + 1) >> 1) + left;
        d->radix((char*) base + mid * d->size, tmpradix, d->arg);
        /* if [mid] <= P , move left to mid */
        /* if [mid] > P , move right to mid*/
        int c1 = d->compar(tmpradix, P, d->rsize);
        if(c1 <= 0) {
            left = mid;
        } else {
            right = mid;
        }
    }
    return left;
}

/*
 * do a histogram of mybase, based on bins defined in P.
 * P is an array of radix of length Plength,
 * myCLT, myCLE are of length Plength + 2
 *
 * myCLT[i + 1] is the count of items less than P[i]
 * myCLE[i + 1] is the count of items less than or equal to P[i]
 *
 * myCLT[0] is always 0
 * myCLT[Plength + 1] is always mynmemb
 *
 * */
static void _histogram(unsigned char * P, int Plength, void * mybase, size_t mynmemb,
        ptrdiff_t * myCLT, ptrdiff_t * myCLE,
        struct crstruct * d) {
    int it;

    if(myCLT) {
        myCLT[0] = 0;
        for(it = 0; it < Plength; it ++) {
            /* No need to start from the beginging of mybase, since myubase and P are both sorted */
            ptrdiff_t offset = myCLT[it];
            myCLT[it + 1] = _bsearch_last_lt(P + it * d->rsize,
                            ((char*) mybase) + offset * d->size,
                            mynmemb - offset, d)
                            + 1 + offset;
        }
        myCLT[it + 1] = mynmemb;
    }
    if(myCLE) {
        myCLE[0] = 0;
        for(it = 0; it < Plength; it ++) {
            /* No need to start from the beginging of mybase, since myubase and P are both sorted */
            ptrdiff_t offset = myCLE[it];
            myCLE[it + 1] = _bsearch_last_le(P + it * d->rsize,
                            ((char*) mybase) + offset * d->size,
                            mynmemb - offset, d)
                            + 1 + offset;
        }
        myCLE[it + 1] = mynmemb;
    }
}

struct piter {
    int * stable;
    int * narrow;
    int Plength;
    unsigned char * Pleft;
    unsigned char * Pright;
    struct crstruct * d;
};
static void piter_init(struct piter * pi,
        unsigned char * Pmin, unsigned char * Pmax, int Plength,
        struct crstruct * d) {
    pi->stable = malloc(Plength * sizeof(int));
    pi->narrow = malloc(Plength * sizeof(int));
    pi->d = d;
    pi->Pleft = malloc(d->rsize * Plength);
    pi->Pright = malloc(d->rsize * Plength);
    pi->Plength = Plength;

    int i;
    for(i = 0; i < pi->Plength; i ++) {
        pi->stable[i] = 0;
        pi->narrow[i] = 0;
        memcpy(&pi->Pleft[i * d->rsize], Pmin, d->rsize);
        memcpy(&pi->Pright[i * d->rsize], Pmax, d->rsize);
    }
}
static void piter_destroy(struct piter * pi) {
    free(pi->stable);
    free(pi->narrow);
    free(pi->Pleft);
    free(pi->Pright);
}

/*
 * this will bisect the left / right in piter.
 * note that piter goes [left, right], thus we need
 * to maintain an internal status to make sure we go over
 * the additional 'right]'. (usual bisect range is
 * '[left, right)' )
 * */
static void piter_bisect(struct piter * pi, unsigned char * P) {
    struct crstruct * d = pi->d;
    int i;
    for(i = 0; i < pi->Plength; i ++) {
        if(pi->stable[i]) continue;
        if(pi->narrow[i]) {
            /* The last iteration, test Pright directly */
            memcpy(&P[i * d->rsize],
                &pi->Pright[i * d->rsize],
                d->rsize);
            pi->stable[i] = 1;
        } else {
            /* ordinary iteration */
            d->bisect(&P[i * d->rsize],
                    &pi->Pleft[i * d->rsize],
                    &pi->Pright[i * d->rsize], d->rsize);
            /* in case the bisect can't move P beyond left,
             * the range is too small, so we set flag narrow,
             * and next iteration we will directly test Pright */
            if(d->compar(&P[i * d->rsize],
                &pi->Pleft[i * d->rsize], d->rsize) <= 0) {
                pi->narrow[i] = 1;
            }
        }
#if 0
        printf("bisect %d %u %u %u\n", i, *(int*) &P[i * d->rsize],
                *(int*) &pi->Pleft[i * d->rsize],
                *(int*) &pi->Pright[i * d->rsize]);
#endif
    }
}
static int piter_all_done(struct piter * pi) {
    int i;
    int done = 1;
#if 0
#pragma omp single
    for(i = 0; i < pi->Plength; i ++) {
        printf("P %d stable %d narrow %d\n",
            i, pi->stable[i], pi->narrow[i]);
    }
#endif
    for(i = 0; i < pi->Plength; i ++) {
        if(!pi->stable[i]) {
            done = 0;
            break;
        }
    }
    return done;
}

/*
 * bisection acceptance test.
 *
 * test if the counts satisfies CLT < C <= CLE.
 * move Pleft / Pright accordingly.
 * */
static void piter_accept(struct piter * pi, unsigned char * P,
        ptrdiff_t * C, ptrdiff_t * CLT, ptrdiff_t * CLE) {
    struct crstruct * d = pi->d;
    int i;
#if 0
    for(i = 0; i < pi->Plength + 1; i ++) {
        printf("counts %d LT %ld C %ld LE %ld\n",
                i, CLT[i], C[i], CLE[i]);
    }
#endif
    for(i = 0; i < pi->Plength; i ++) {
        if( CLT[i + 1] < C[i + 1] && C[i + 1] <= CLE[i + 1]) {
            pi->stable[i] = 1;
            continue;
        } else {
            if(CLT[i + 1] >= C[i + 1]) {
                /* P[i] is too big */
                memcpy(&pi->Pright[i * d->rsize], &P[i * d->rsize], d->rsize);
            } else {
                /* P[i] is too small */
                memcpy(&pi->Pleft[i * d->rsize], &P[i * d->rsize], d->rsize);
            }
        }
    }
}
