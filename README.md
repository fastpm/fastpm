cola_halo
=========

Parallel COLA cosmological simulation + 2LPT initial condition
generator + FoF halo finder

This is a parallel cosmological N-body simulation code running the
following on the fly, designed to generate hundreds of realizations:

1. Generates random Gaussian initial condition with 2LPT --
second-order Lagrangian Perturbation Theory (Scoccimarro 1998; Crocce,
Pueblas, & Scoccimarro 2006,[[1]]).
2. Time evolve N-body particles with the COmoving Lagrangian
Acceleration (COLA) method (Tassev, Zaldarriaga, & Eisenstein
2013,[[2]])
3. Finds dark-matter haloes with the Friends-of-Friends algorithm
(Davis et al. 1985, [Washington University N-body shop][3]).


This code is based on publicly avaiable codes ([1], [2], [3]),
modified and redistributed under GPLv3.

# Authors
This code is assembled and parallelized (COLA & FOF) by Jun Koda.

# References 
1. Crocce, M., Pueblas, S., & Scoccimarro, R. 2006, MNRAS, 373, 369
2. Davis, M., Efstathiou, G., Frenk, C.~S., & White, S.D.M. 1985, ApJ, 292, 371 
3. Scoccimarro, R. 1998, MNRAS, 299, 1097
4. Tassev, S., Zaldarriaga, M., & Eisenstein, D.~J. 2013, JCAP, 6, 36 

[1]: http://cosmo.nyu.edu/roman/2LPT/
[2]: https://bitbucket.org/tassev/colacode/
[3]: http://www-hpcc.astro.washington.edu/tools/fof.html

