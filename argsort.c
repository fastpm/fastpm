/*this is not tested. Finish it up! */
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>

typedef int (*__compar_d_fn_t) (const void *, const void *, void *);

void
argsort(void *b, size_t n, size_t s, __compar_d_fn_t cmp, void *arg, 
        ptrdiff_t * argind, void * scratch, size_t scratch_bytes)
{
    struct msort_param p;

    size = n * sizeof (void *);
    if (size > scratch_bytes) 
    {
        abort();
        return;
    } else {
        p.t = scratch;
    }

    p.s = s;
    p.var = 3;
    p.cmp = cmp;
    p.arg = arg;

    /* Indirect sorting.  */
    char *ip = (char *) b;
    void **tp = (void **) argind;
    void **t = tp;

    while ((void *) t < (void *) (tp + n))
    {
        *t++ = ip;
        ip += s;
    }
    p.s = sizeof (void *);
    msort_with_tmp (&p, p.t, n);

    ptrdiff_t *tind = argind;

    while ((void*) tind < (void *) (tp + n))
    {
        *tind -= ip;
        tind ++;
    }

}

static void
msort_with_tmp (const struct msort_param *p, void *b, size_t n)
{
    char *b1, *b2;
    size_t n1, n2;

    if (n <= 1)
        return;

    n1 = n / 2;
    n2 = n - n1;
    b1 = b;
    b2 = (char *) b + (n1 * p->s);

    msort_with_tmp (p, b1, n1);
    msort_with_tmp (p, b2, n2);

    char *tmp = p->t;
    const size_t s = p->s;
    __compar_d_fn_t cmp = p->cmp;
    void *arg = p->arg;
    switch (p->var)
    {
        case 3:
            while (n1 > 0 && n2 > 0)
            {
                if ((*cmp) (*(const void **) b1, *(const void **) b2, arg) <= 0)
                {
                    *(void **) tmp = *(void **) b1;
                    b1 += sizeof (void *);
                    --n1;
                }
                else
                {
                    *(void **) tmp = *(void **) b2;
                    b2 += sizeof (void *);
                    --n2;
                }
                tmp += sizeof (void *);
            }
            break;
    }

    if (n1 > 0)
        memcpy (tmp, b1, n1 * s);
    memcpy (b, p->t, (n - n2) * s);
}
