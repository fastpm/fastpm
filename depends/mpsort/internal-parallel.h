
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

    char tmpradix[d->rsize];
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

    char tmpradix[d->rsize];
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
static void _histogram(char * P, int Plength, void * mybase, size_t mynmemb, 
        ptrdiff_t * myCLT, ptrdiff_t * myCLE,
        struct crstruct * d) {
    int it;

    myCLT[0] = 0;
    myCLE[0] = 0;
    for(it = 0; it < Plength; it ++) {
        myCLT[it + 1] = _bsearch_last_lt(P + it * d->rsize, mybase, mynmemb, d) + 1;
        myCLE[it + 1] = _bsearch_last_le(P + it * d->rsize, mybase, mynmemb, d) + 1;
    }
    myCLT[it + 1] = mynmemb;
    myCLE[it + 1] = mynmemb;
}

#if 0
/* 
 * solve for the communication layout based on
 *
 * C: the desired number of items per task
 * GL_CLT[t,i+1]: the offset of lt P[i] in task t 
 * GL_CLE[t,i+1]: the offset of le P[i] in task t
 * 
 * the result is saved in
 *
 * GL_C[t, i]: the offset of sending to task i in task t.
 *
 * this routine requires GL_ to scale with NTask * NTask;
 * won't work with 1,000 + ranks.
 * */
static void _solve_for_layout (
        int NTask, 
        ptrdiff_t * C,
        ptrdiff_t * GL_CLT, 
        ptrdiff_t * GL_CLE, 
        ptrdiff_t * GL_C) {
    int NTask1 = NTask + 1;
    int i, j;
    /* first assume we just send according to GL_CLT */
    for(i = 0; i < NTask + 1; i ++) {
        for(j = 0; j < NTask; j ++) {
            GL_C[j * NTask1 + i] = GL_CLT[j * NTask1 + i];
        }
    }

    /* Solve for each receiving task i 
     *
     * this solves for GL_C[..., i + 1], which depends on GL_C[..., i]
     *
     * and we have GL_C[..., 0] == 0 by definition.
     *
     * this cannot be done in parallel wrt i because of the dependency. 
     *
     *  a solution is guaranteed because GL_CLE and GL_CLT
     *  brackes the total counts C (we've found it with the
     *  iterative counting.
     *
     * */

    for(i = 0; i < NTask; i ++) {
        ptrdiff_t sure = 0;

        /* how many will I surely receive? */
        for(j = 0; j < NTask; j ++) {
            ptrdiff_t sendcount = GL_C[j * NTask1 + i + 1] - GL_C[j * NTask1 + i];
            sure += sendcount;
        }
        /* let's see if we have enough */
        ptrdiff_t deficit = C[i + 1] - C[i] - sure;

        for(j = 0; j < NTask; j ++) {
            /* deficit solved */
            if(deficit == 0) break;
            if(deficit < 0) {
                fprintf(stderr, "serious bug: more items than there should be: deficit=%ld\n", deficit);
                abort();
            }
            /* how much task j can supply ? */
            ptrdiff_t supply = GL_CLE[j * NTask1 + i + 1] - GL_C[j * NTask1 + i + 1];
            if(supply < 0) {
                fprintf(stderr, "serious bug: less items than there should be: supply =%ld\n", supply);
                abort();
            }
            if(supply <= deficit) {
                GL_C[j * NTask1 + i + 1] += supply;
                deficit -= supply;
            } else {
                GL_C[j * NTask1 + i + 1] += deficit;
                deficit = 0;
            }
        }
    }

#if 0
    for(i = 0; i < NTask; i ++) {
        for(j = 0; j < NTask + 1; j ++) {
            printf("%d %d %d, ", 
                    GL_CLT[i * NTask1 + j], 
                    GL_C[i * NTask1 + j], 
                    GL_CLE[i * NTask1 + j]);
        }
        printf("\n");
    }
#endif

}
#endif

struct piter {
    int * stable;
    int * narrow;
    int Plength;
    char * Pleft;
    char * Pright;
    struct crstruct * d;
};
static void piter_init(struct piter * pi, 
        char * Pmin, char * Pmax, int Plength,
        struct crstruct * d) {
    pi->stable = calloc(Plength, sizeof(int));
    pi->narrow = calloc(Plength, sizeof(int));
    pi->d = d;
    pi->Pleft = calloc(Plength, d->rsize);
    pi->Pright = calloc(Plength, d->rsize);
    pi->Plength = Plength;
    
    int i;
    for(i = 0; i < pi->Plength; i ++) {
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
static void piter_bisect(struct piter * pi, char * P) {
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
                &pi->Pleft[i * d->rsize], d->rsize) == 0) {
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
static void piter_accept(struct piter * pi, char * P, 
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
