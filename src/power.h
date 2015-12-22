#ifndef POWER_H
#define POWER_H 1

double PowerSpec(const double k);
double PowerSpecWithData(const double k, void * data);
void power_init(const char filename[], const double a_init, const double sigma8, const double omega_m, const double omega_lambda, MPI_Comm comm);
#endif
