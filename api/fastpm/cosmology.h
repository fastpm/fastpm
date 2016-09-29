#ifdef __cplusplus
extern "C" {
#endif
extern double HubbleConstant;
extern double HubbleDistance;

typedef struct {
    double OmegaM;
    double OmegaLambda;
} FastPMCosmology;

double GrowthFactor(double a, FastPMCosmology c);
double GrowthFactor2(double a, FastPMCosmology c);

double DLogGrowthFactor(double a, FastPMCosmology c);
double DLogGrowthFactor2(double a, FastPMCosmology c);

double HubbleEa(double a, FastPMCosmology c);
double DHubbleEaDa(double a, FastPMCosmology c);
double D2HubbleEaDa2(double a, FastPMCosmology c);
double DGrowthFactorDa(double a, FastPMCosmology c);
double D2GrowthFactorDa2(double a, FastPMCosmology c);

double ComovingDistance(double a, FastPMCosmology c);
double OmegaA(double a, FastPMCosmology c);

#ifdef __cplusplus
}
#endif
