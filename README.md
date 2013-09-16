cola_halo
=========

Parallel COLA cosmological simulation + 2LPT initial condition generator + FoF halo finder

This is a cosmological N-body simulation code doing the following on the fly, designed to generate hundreds of realizations:
1. Generates random Gaussian initial condition with 2LPT -- second-order Lagrangian Perturbation Theory (Scoccimarro 1998, Crocce, Pueblas, & Scoccimarro 2006).
2. Time evolve N-body particles with the COmoving Lagrangian Acceleration (COLA) method (Tassev, Zaldarriaga, & Eisenstein 2013)
3. Finds dark-matter haloes with the Friends-of-Friends algorithm (Davis et al. 1985).

# References 
Crocce, M., Pueblas, S., & Scoccimarro, R. 2006, MNRAS, 373, 369
Davis, M., Efstathiou, G., Frenk, C.~S., & White, S.D.M. 1985, ApJ, 292, 371 
Scoccimarro, R. 1998, MNRAS, 299, 1097
Tassev, S., Zaldarriaga, M., & Eisenstein, D.~J. 2013, JCAP, 6, 36 
