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
    unsigned char * done = alloca(np / 8) + 1;
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
        memcpy(temp1, &q[i * elsize], elsize);

        /* loop till we are back to the head of the ring */

        int i0 = i;
        for(j = ind[i]; j != i0; i = j, j = ind[j]) {
            memcpy(&q[i * elsize], &q[j * elsize], elsize);
            /* now j contains the correct item */;
            done[i >> 3] |= 1 << (i & 7);
        }
        /* now move the saved item to the end of the ring */
        memcpy(&q[i * elsize], temp1, elsize);
        done[i >> 3] |= 1 << (i & 7);
    }
}

static void ipermute(void * data, int np, size_t elsize, int * ind){
    unsigned char * done = alloca(np / 8) + 1;
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
}

#ifdef __TEST_PERMUTE__
void main() {
    char data[5] = {0, 1, 2, 3, 4};
    char odata[5] = {0, 1, 2, 3, 4};
    int ind[5] = {1, 2, 3, 4, 0};
    memcpy(odata, data, 5);
    permute(data, 5, sizeof(char), ind);
    int i;
    for(i = 0; i < 5; i ++){
        printf("data[%d] <= data[%d], %d\n", i, ind[i], data[i]);
    }
    memcpy(data, odata, 5);
    ipermute(data, 5, sizeof(char), ind);
    for(i = 0; i < 5; i ++){
        printf("data[%d] <= data[%d], %d\n", ind[i], i, data[i]);
    }
}
#endif
