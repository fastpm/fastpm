FASTPM_BEGIN_DECLS

#define FASTPM_CRITICAL_DENSITY 27.7455 /* 1e10 Msun /h*/

void
fastpm_solver_compute_force(FastPMSolver * fastpm,
    FastPMPainter * painter,
    FastPMDealiasingType dealias,
    FastPMKernelType kernel,
    FastPMFloat * delta_k);

void
gravity_apply_kernel_transfer(FastPMKernelType kernel, PM * pm, FastPMFloat * delta_k, FastPMFloat * canvas, FastPMFieldDescr field);

FASTPM_END_DECLS
