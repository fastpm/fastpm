
double 
GrowthFactor(double astart, double aend, double Omega, double OmegaLambda);

static inline double Qfactor(double a, double Omega, double OmegaLambda) { // Q\equiv a^3 H(a)/H0.
    return sqrt(Omega/(a*a*a)+OmegaLambda)*a*a*a;
}

