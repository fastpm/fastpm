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
gravity_apply_kernel_transfer(FastPMGravity * gravity, PM * pm, FastPMFloat * delta_k, FastPMFloat * canvas, int d);

FASTPM_END_DECLS
