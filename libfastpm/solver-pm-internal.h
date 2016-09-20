void fastpm_calculate_forces(FastPMSolver * fastpm, FastPMFloat * delta_k);

void 
fastpm_kick_store(FastPMSolver * fastpm, 
              FastPMStore * pi, FastPMStore * po,
              double af);

void 
fastpm_drift_store(FastPMSolver * fastpm,
               FastPMStore * pi, FastPMStore * po,
               double af);

void 
fastpm_set_snapshot(FastPMSolver * fastpm,
                FastPMStore * p, FastPMStore * po,
                double aout);
