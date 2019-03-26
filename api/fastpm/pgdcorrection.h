
FASTPM_BEGIN_DECLS

typedef struct {
    FastPMPainterType PainterType;
    int PainterSupport;
    double alpha0;
    double A;
    double B;
    double kl;
    double ks;
} FastPMPGDCorrection;

double
fastpm_pgdc_get_alpha(FastPMPGDCorrection * pgdc, double a);

double
fastpm_pgdc_get_kl(FastPMPGDCorrection * pgdc, double a);

double
fastpm_pgdc_get_ks(FastPMPGDCorrection * pgdc, double a);

void
fastpm_pgdc_calculate(FastPMPGDCorrection * pgdc,
    PM * pm,
    FastPMStore * p,
    FastPMFloat * delta_k, double a, double fac);

FASTPM_END_DECLS
