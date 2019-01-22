typedef void (*reduce_func)(
    FastPMStore * src,
    ptrdiff_t isrc,
    FastPMStore * dest,
    ptrdiff_t idest,
    int ci,
    void * userdata);

typedef struct PMGhostData {
    PM * pm;
    FastPMStore * source; /* original particles */
    FastPMStore * p;  /* ghost particles */

    double Below[3];
    double Above[3];

    /* private members */
    int * Nsend;
    int * Osend;
    int * Nrecv;
    int * Orecv;
    void * send_buffer;
    void * recv_buffer;

    /* iterator status */
    ptrdiff_t ipar;
    int * ighost_to_ipar;
    int rank;
    int * reason; /* relative offset causing the ghost */
} PMGhostData;

PMGhostData * 
pm_ghosts_create(PM * pm, FastPMStore * p, FastPMColumnTags attributes, int support);

PMGhostData * 
pm_ghosts_create_full(PM * pm, FastPMStore * p,
        FastPMColumnTags attributes,
        double below[],
        double above[]
        );
void
pm_ghosts_send(PMGhostData * pgd, FastPMColumnTags attributes);

void pm_ghosts_reduce(PMGhostData * pgd, FastPMColumnTags attribute, reduce_func reduce, void * userdata);

void pm_ghosts_free(PMGhostData * pgd);

void
pm_ghosts_has_ghosts(PMGhostData * pgd, uint8_t * has_ghosts);
