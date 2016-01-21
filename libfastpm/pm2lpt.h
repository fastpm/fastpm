
void 
pm_2lpt_init(PM * pm, PMStore * p, int Ngrid, double BoxSize, MPI_Comm comm);

void 
pm_2lpt_solve(PM * pm, FastPMFloat * delta_k, PMStore * p, double shift[3]);

void 
pm_2lpt_evolve(double aout, PMStore * p, double Omega, int zaonly);
