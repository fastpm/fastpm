typedef double (*pkfunc)(double k, void * data);

/* The following functions fill the gaussian field*/
void 
pm_ic_fill_gaussian_gadget(PM * pm, int seed, pkfunc pk, void * pkdata);
void 
pm_ic_fill_gaussian_fast(PM * pm, int seed, pkfunc pk, void * pkdata);
void 
pm_ic_fill_gaussian_slow(PM * pm, int seed, pkfunc pk, void * pkdata);
void 
pm_ic_read_gaussian(PM * pm, char * filename, pkfunc pk, void * pkdata);
