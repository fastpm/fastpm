typedef double (*pkfunc)(double k, void * data);

/* The following functions fill the gaussian field*/
void 
pmic_fill_gaussian_gadget(PM * pm, float_t * delta_k, int seed, pkfunc pk, void * pkdata);
void 
pmic_fill_gaussian_fast(PM * pm, float_t * delta_k, int seed, pkfunc pk, void * pkdata);
void 
pmic_fill_gaussian_slow(PM * pm, float_t * delta_k, int seed, pkfunc pk, void * pkdata);
void 
pmic_read_gaussian(PM * pm, float_t * delta_k, char * filename, pkfunc pk, void * pkdata);
