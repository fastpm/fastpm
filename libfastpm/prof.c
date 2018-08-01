#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include <fastpm/libfastpm.h>
#include <fastpm/prof.h>
#include <fastpm/logging.h>

struct FastPMClock {
    double tcum;
    double twait;
    double t0;
    char file[128];
    char func[128];
    char name[128];
    struct FastPMClock * next;
};

static FastPMClock * head = NULL;

FastPMClock * 
fastpm_clock_create(const char * file, const char * func, const char * name) 
{
    FastPMClock * p = malloc(sizeof(FastPMClock));
    p->t0 = 0;
    p->tcum = 0;
    p->twait = 0;
    strncpy(p->file, file, 120);
    strncpy(p->func, func, 120);
    strncpy(p->name, name, 120);
    p->next = NULL;
    return p;
}

FastPMClock * 
fastpm_clock_find(const char * file, const char * func, const char * name) 
{
    FastPMClock * p;
    if(head == NULL) {
        goto notfound;
    }
    for(p = head; p; p=p->next) {
        if( (0 == strcmp(p->file, file))
         && (0 == strcmp(p->func, func))
         && (0 == strcmp(p->name, name))) {
            return p;
        }
    }
notfound:
    p = fastpm_clock_create(file, func, name);
    p->next = head;
    head = p;
    return p;
}

void fastpm_clock_in(FastPMClock * clock) 
{
    clock->t0 = MPI_Wtime();
}

void fastpm_clock_out(FastPMClock * clock) 
{
    double t1 = MPI_Wtime();
    clock->tcum += (t1 - clock->t0);
}

void fastpm_clock_out_barrier(FastPMClock * clock, MPI_Comm comm) 
{
    double t1 = MPI_Wtime();
    clock->tcum += (t1 - clock->t0);

    MPI_Barrier(comm);
    double t2 = MPI_Wtime();
    clock->twait += t2 - t1;
}

static void fastpm_clock_stat_one(FastPMClock * clock, MPI_Comm comm)
{
    double min = 0;
    double max = 0;
    double mean = 0;
    int N;
    MPI_Comm_size(comm, &N);
    MPI_Allreduce(&clock->tcum, &min, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&clock->tcum, &max, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&clock->tcum, &mean, 1, MPI_DOUBLE, MPI_SUM, comm);
    mean /= N;
    fastpm_info("%8.2f %8.2f %8.2f  %12s : %s : %s \n", 
            min, max, mean, clock->name, clock->func, clock->file);

}
static int
fastpm_clock_cmp(FastPMClock * c1, FastPMClock * c2)
{
    int i;
    i = strcmp(c1->file, c2->file);
    if(i) return i;
    i = strcmp(c1->func, c2->func);
    if(i) return i;
    i = strcmp(c1->name, c2->name);
    if(i) return i;
    return 0;
}

static void fastpm_clock_sort() 
{
    FastPMClock * h = NULL;
    FastPMClock * q = NULL;
    FastPMClock * p = head;
    FastPMClock * pn;
    FastPMClock * qn;
 
    h = head;
    if(!head) return;

    head = head->next;
    h->next = NULL;
 
    for(p = head; p; p=pn) {
        pn = p->next;
        if(fastpm_clock_cmp(p, h) < 0) {
            p->next = h;
            h = p;            
            continue;
        } 
        for(q = h; q->next; q=qn) {
            qn = q->next; 
            if(fastpm_clock_cmp(p, q->next) < 0) {
                p->next = q->next;
                q->next = p;
                goto next;
            }
        }
        /* no inserting, append */
        q->next = p;
        p->next = NULL;
next:
        continue;
    }

    head = h;
}

void fastpm_clock_stat(MPI_Comm comm)
{
    FastPMClock * p;
    int ThisTask;
    MPI_Comm_rank(comm, &ThisTask);
    FastPMClock foo = {0};

    fastpm_info("%8s %8s %8s \n", "min", "max", "mean");

    fastpm_clock_sort();

    const char * lastfunc = "";
    if(ThisTask == 0) {
        for(p = head; p ; p=p->next) {
            MPI_Bcast(p, sizeof(foo), MPI_BYTE, 0, comm); 
            if(strcmp(lastfunc, p->func) != 0) {
                fastpm_ilog(INFO, "-------------------------------------"
                            "-------------------------------------\n");
            }
            lastfunc = p->func;
            fastpm_clock_stat_one(p, comm);
        }
        foo.file[0] = 0;
        MPI_Bcast(&foo, sizeof(foo), MPI_BYTE, 0, comm); 
    } else {
        while(1) {
            MPI_Bcast(&foo, sizeof(foo), MPI_BYTE, 0, comm); 
            if(foo.file[0] == 0) {
                break;
            }
            p = fastpm_clock_find(foo.file, foo.func, foo.name);
            fastpm_clock_stat_one(p, comm);
        }
    }
}

