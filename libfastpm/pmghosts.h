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
pm_ghosts_create(PM * pm, FastPMStore * p, FastPMColumnTags attributes);

PMGhostData * 
pm_ghosts_create_full(PM * pm, FastPMStore * p,
        FastPMColumnTags attributes,
        double below[],
        double above[]
        );
void
pm_ghosts_send(PMGhostData * pgd, FastPMColumnTags attributes);

void pm_ghosts_reduce(PMGhostData * pgd, FastPMFieldDescr field);

void pm_ghosts_free(PMGhostData * pgd);

