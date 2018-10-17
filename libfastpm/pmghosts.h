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
pm_ghosts_create(PM * pm, FastPMStore * p, enum FastPMPackFields attributes);

PMGhostData * 
pm_ghosts_create_full(PM * pm, FastPMStore * p,
        enum FastPMPackFields attributes,
        double below[],
        double above[]
        );
void
pm_ghosts_send(PMGhostData * pgd, enum FastPMPackFields attributes);

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

