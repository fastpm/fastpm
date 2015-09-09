
typedef double (*pkfunc)(double k, void * data);

void 
pm_2lpt_main(PMStore * p, int Ngrid, double BoxSize, 
    pkfunc pk, int seed, void * pkdata, MPI_Comm comm);

