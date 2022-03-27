/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2012 Krzysztof M. Gorski, Eric Hivon, Martin Reinecke,
 *                          Benjamin D. Wandelt, Anthony J. Banday,
 *                          Matthias Bartelmann,
 *                          Reza Ansari & Kenneth M. Ganga
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------*/
/*
 * chealpix.h
 */

#ifndef CHEALPIX_H
#define CHEALPIX_H

#ifdef __cplusplus
extern "C" {
#endif

/* -------------------- */
/* Constant Definitions */
/* -------------------- */

#ifndef HEALPIX_NULLVAL
#define HEALPIX_NULLVAL (-1.6375e30)
#endif /* HEALPIX_NULLVAL */

/* pixel operations */
/* ---------------- */
void ang2pix_nest(long nside, double theta, double phi, long *ipix);
void ang2pix_ring(long nside, double theta, double phi, long *ipix);

void pix2ang_nest(long nside, long ipix, double *theta, double *phi);
void pix2ang_ring(long nside, long ipix, double *theta, double *phi);

void nest2ring(long nside, long ipnest, long *ipring);
void ring2nest(long nside, long ipring, long *ipnest);

long nside2npix(long nside);
long npix2nside(long npix);

void ang2vec(double theta, double phi,   double *vec);
void vec2ang(const double *vec, double *theta, double *phi);

void vec2pix_nest(long nside, const double *vec, long *ipix);
void vec2pix_ring(long nside, const double *vec, long *ipix);

void pix2vec_nest(long nside, long ipix, double *vec);
void pix2vec_ring(long nside, long ipix, double *vec);

/* operations on Nside values up to 2^29 */

typedef int64_t hpint64;

void ang2pix_nest64(hpint64 nside, double theta, double phi, hpint64 *ipix);
void ang2pix_ring64(hpint64 nside, double theta, double phi, hpint64 *ipix);

void pix2ang_nest64(hpint64 nside, hpint64 ipix, double *theta, double *phi);
void pix2ang_ring64(hpint64 nside, hpint64 ipix, double *theta, double *phi);

void nest2ring64(hpint64 nside, hpint64 ipnest, hpint64 *ipring);
void ring2nest64(hpint64 nside, hpint64 ipring, hpint64 *ipnest);

hpint64 nside2npix64(hpint64 nside);
long npix2nside64(hpint64 npix);

void vec2pix_nest64(hpint64 nside, const double *vec, hpint64 *ipix);
void vec2pix_ring64(hpint64 nside, const double *vec, hpint64 *ipix);

void pix2vec_nest64(hpint64 nside, hpint64 ipix, double *vec);
void pix2vec_ring64(hpint64 nside, hpint64 ipix, double *vec);

/* FITS operations */
/* --------------- */

float *read_healpix_map (const char *infile, long *nside, char *coordsys,
  char *ordering);

void write_healpix_map (const float *signal, long nside, const char *filename,
  char nest, const char *coordsys);

long get_fits_size(const char *filename, long *nside, char *ordering);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CHEALPIX_H */
