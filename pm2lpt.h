
typedef double (*pkfunc)(double k, void * data);

void 
pm_2lpt_init(PM * pm, PMStore * p, int Ngrid, double BoxSize, MPI_Comm comm);

void 
pm_2lpt_main(PM * pm, PMStore * p, MPI_Comm comm);

/* The following functions fill the gaussian field*/
void 
pm_2lpt_fill_gaussian_gadget(PM * pm, int seed, pkfunc pk, void * pkdata);
void 
pm_2lpt_fill_gaussian_fast(PM * pm, int seed, pkfunc pk, void * pkdata);
void 
pm_2lpt_fill_gaussian_slow(PM * pm, int seed, pkfunc pk, void * pkdata);
