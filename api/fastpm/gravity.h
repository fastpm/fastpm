FASTPM_BEGIN_DECLS

#define FASTPM_CRITICAL_DENSITY 27.7455 /* 1e10 Msun /h*/

typedef struct {
    FastPMKernelType KernelType;
    FastPMDealiasingType DealiasingType;
    FastPMPainterType PainterType;
    int PainterSupport;
} FastPMGravity;

void
fastpm_gravity_calculate(FastPMGravity * gravity, PM * pm, FastPMStore * p, FastPMFloat * delta_k);

void
gravity_apply_kernel_transfer(FastPMGravity * gravity, PM * pm, FastPMFloat * delta_k, FastPMFloat * canvas, FastPMFieldDescr field);

FASTPM_END_DECLS
