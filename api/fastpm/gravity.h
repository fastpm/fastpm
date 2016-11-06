FASTPM_BEGIN_DECLS
typedef struct {
    FastPMKernelType KernelType;
    FastPMDealiasingType DealiasingType;
    FastPMPainterType PainterType;
    int PainterSupport;
    int ComputePotential;
} FastPMGravity;

void
fastpm_gravity_calculate(FastPMGravity * gravity, PM * pm, FastPMStore * p, FastPMFloat * delta_k);

void
fastpm_gravity_calculate_gradient(FastPMGravity * gravity,
    PM * pm,
    FastPMStore * grad_acc,
    FastPMStore * p,
    FastPMStore * grad_pos
);

FASTPM_END_DECLS
