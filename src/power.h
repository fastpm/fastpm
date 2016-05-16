#ifndef POWER_H
#define POWER_H 1

double PowerSpec(const double k);
double PowerSpecWithData(const double k, void * data);
void power_init(const char filename[], const double sigma8, MPI_Comm comm);
#endif
