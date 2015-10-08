
void 
pm_2lpt_init(PM * pm, PMStore * p, int Ngrid, double BoxSize, MPI_Comm comm);

void 
pm_2lpt_main(PM * pm, PMStore * p, MPI_Comm comm);

void 
pm2lpt_set_initial(double aout, PMStore * p, double shift[3], double Omega);
