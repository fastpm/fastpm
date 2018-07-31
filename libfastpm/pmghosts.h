typedef struct PMGhostData {
    PM * pm;
    FastPMStore * p;
    size_t np;
    size_t np_upper;
    size_t nghosts;

    /* the get_position member is used to determine the target rank */
    fastpm_posfunc get_position;
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
    size_t elsize;
} PMGhostData;

PMGhostData * 
pm_ghosts_create(PM * pm, FastPMStore * p, enum FastPMPackFields attributes, fastpm_posfunc get_position);

PMGhostData * 
pm_ghosts_create_full(PM * pm, FastPMStore * p,
        enum FastPMPackFields attributes,
        fastpm_posfunc get_position,
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

