#ifdef __TEST_PERMUTE__
#include <stdio.h>
#include <stdlib.h>
#endif
/*
   These functions permute and ipermute array given by integer
    permute array ind.

   A scratch space of size np / 8 bytes is used (allocated on the
    stack)
    
   permute:    OUT[i]       = IN[ind[i]]     i = 0 .. N-1
   invpermute: OUT[ind[i]] = IN[i]           i = 0 .. N-1
*/
static void permute(void * data, int np, size_t elsize, int * ind){
    unsigned char * done = malloc(np / 8 + 1);
    if(!done) {
        fastpm_raise(-1, "no memory for 'done'\n");
    }
    char * q = data;
    char * temp1 = alloca(elsize);
    int j, i;
    memset(done, 0, np / 8 + 1);
    for(i = 0; i < np ; i ++){
        /* if this item is already permuted by a previous ring.  */
        if(done[i >> 3] & (1 << (i & 7))) continue;

        /* shuffle items on the permuting-ring starting from i */
        /* this works too when the ring is of length 1. */

        /* mov the old item at i to temp1 to bootstrap the shuffling */
        if (i >= np || i < 0) {
            abort();
        }
        memcpy(temp1, &q[i * elsize], elsize);

        /* loop till we are back to the head of the ring */

        int ii = i;
        for(j = ind[ii]; j != i; ii = j, j = ind[j]) {
            if (j >= np || j < 0) {
                abort();
            }
            if (ii >= np || ii < 0) {
                abort();
            }
            memcpy(&q[ii * elsize], &q[j * elsize], elsize);
            /* now j contains the correct item */;
            done[ii >> 3] |= 1 << (ii & 7);
        }
        if (ii >= np || ii < 0) {
            abort();
        }
        /* now move the saved item to the end of the ring */
        memcpy(&q[ii * elsize], temp1, elsize);
        done[ii >> 3] |= 1 << (ii & 7);
    }
    free(done);
}

static void ipermute(void * data, int np, size_t elsize, int * ind){
    unsigned char * done = malloc(np / 8 + 1);
    if(!done) {
        fastpm_raise(-1, "no memory for 'done'\n");
    }
    char * q = data;
    char * temp1 = alloca(elsize);
    char * temp2 = alloca(elsize);
    int j, i;
    memset(done, 0, np / 8 + 1);
    for(i = 0; i < np ; i ++){
        /* if this item is already permuted by a previous ring.  */
        if(done[i >> 3] & (1 << (i & 7))) continue;

        /* shuffle items on the permuting-ring starting from i */
        /* this works too when the ring is of length 1. */

        /* mov the old item at i to temp1 to bootstrap the shuffling */
        memcpy(temp1, &q[i * elsize], elsize);

        /* temp1 always contains the item being moved to the correct location */
        /* loop till we are back to the head of the ring */

        for(j = ind[i]; j != i; j = ind[j]) {
            /* swap temp1 and j */
            memcpy(temp2, &q[j * elsize], elsize);
            memcpy(&q[j * elsize], temp1, elsize);
            memcpy(temp1, temp2, elsize);
            /* now j contains the correct item */;
            done[j >> 3] |= 1 << (j & 7);
        }
        /* no need of a full swap to finish the head of the ring */
        memcpy(&q[j * elsize], temp1, elsize);
        done[j >> 3] |= 1 << (j & 7);
    }
    free(done);
}

#ifdef __TEST_PERMUTE__
void test(int * ind) {
    int i;
    int * data = malloc(sizeof(int) * 4);
    for(i = 0; i < 4; i ++) {
        data[i] = i;
    }
    permute(data, 4, sizeof(data[0]), ind);
    for(i = 0; i < 4; i ++){
        if(data[i] != ind[i]) {
            abort();
        }
    }
}
void main() {
    int ind[4];
    int i, j;
    for(ind[0] = 0; ind[0] < 4; ind[0]++)
    for(ind[1] = 0; ind[1] < 4; ind[1]++)
    for(ind[2] = 0; ind[2] < 4; ind[2]++)
    for(ind[3] = 0; ind[3] < 4; ind[3]++) {
        for(i = 0; i < 4; i ++) {
            for(j = 0; j < 4; j ++) {
                if(i == j) continue;
                if(ind[i] == ind[j]) goto next;
            }
        }
        for(i = 0; i < 4; i ++) {
            printf("%d ", ind[i]);
        }
        printf("\n");
        test(ind);
        next:
        continue;
    }
}

#endif
