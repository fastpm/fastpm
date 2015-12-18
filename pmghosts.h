typedef struct PMGhostData {
    PM * pm;
    void * pdata;
    size_t np;
    size_t np_upper;
    size_t nghosts;
    int attributes;
    void   (*get_position)(void * pdata, ptrdiff_t index, double pos[3]);

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
    int ReductionAttributes;
    size_t elsize;
} PMGhostData;

typedef void (*pm_iter_ghosts_func)(PM * pm, PMGhostData * ppd);

PMGhostData * pm_ghosts_create(PM * pm, PMStore * p, int attributes, 
        void (*getpos)(void * pdata, ptrdiff_t index, double pos[3]));

void pm_ghosts_reduce(PMGhostData * pgd, int attributes);
void pm_ghosts_free(PMGhostData * pgd);

