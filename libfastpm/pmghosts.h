typedef struct PMGhostData {
    PM * pm;
    FastPMStore * p;
    size_t np;
    size_t np_upper;
    size_t nghosts;
    enum FastPMPackFields attributes;
    fastpm_posfunc get_position;

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
    ptrdiff_t * reason; /* relative offset causing the ghost */
    size_t elsize;
} PMGhostData;

typedef void (*pm_iter_ghosts_func)(PM * pm, PMGhostData * ppd);

PMGhostData * 
pm_ghosts_create(PM * pm, FastPMStore * p, enum FastPMPackFields attributes, fastpm_posfunc get_position);

void pm_ghosts_reduce(PMGhostData * pgd, enum FastPMPackFields attributes);

typedef
void (*pm_ghosts_reduce_func)(PMGhostData * pgd,
            enum FastPMPackFields attributes,
            ptrdiff_t index,
            void * buf, void * userdata);

void
pm_ghosts_reduce_any(PMGhostData * pgd, enum FastPMPackFields attributes,
        pm_ghosts_reduce_func func, void * userdata);

void pm_ghosts_free(PMGhostData * pgd);

