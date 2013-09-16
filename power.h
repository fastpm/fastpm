#ifndef POWER_H
#define POWER_H 1

double PowerSpec(const double k);
void power_init(const char filename[], const double a_init, const double sigma8, const double omega_m, const double omega_lambda);
double GrowthFactor(double astart, double aend);

#endif
