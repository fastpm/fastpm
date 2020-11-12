#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

#define LENGTH_FERMI_DIRAC_TABLE 4000 //sets the length of the table on which CDF
                                      //will be evaluated
#define MAX_FERMI_DIRAC          20.0 //Upper limit of F-D distribution in units of
                                      // p/T

#define NEXT(n, i)  (((n) + (i)/(n)) >> 1) // Needed for Healpix routine

//needed for healpix
unsigned int isqrt(int number) {
    unsigned int n  = 1;
    unsigned int n1 = NEXT(n, number);

    while(abs(n1 - n) > 1) {
        n  = n1;
        n1 = NEXT(n, number);
    }
    while(n1*n1 > number)
        n1--;
    return n1;
}


//Converts from pixel number to unit vector (needed for healpix)
void pix2vec (int pix, double *vec, int n_side)
{
    double z, phi;
    int nside_ = n_side;
    long ncap_=nside_*(nside_-1)*2;
    long npix_=12*nside_*nside_;
    double fact2_  = 4./npix_;
    if (pix<ncap_) /* North Polar cap */
    {
        int iring = (int)(0.5*(1+isqrt(1+2*pix))); /* counted from North pole */
        int iphi  = (pix+1) - 2*iring*(iring-1);

        z = 1.0 - (iring*iring)*fact2_;
        phi = (iphi-0.5) * 0.5*M_PI/iring;
    }
    else if (pix<(npix_-ncap_)) /* Equatorial region */
    {
        double fact1_  = (nside_<<1)*fact2_;
        int ip  = pix - ncap_;
        int iring = ip/(4*nside_) + nside_; /* counted from North pole */
        int iphi  = ip%(4*nside_) + 1;
        /* 1 if iring+nside is odd, 1/2 otherwise */
        double fodd = ((iring+nside_)&1) ? 1 : 0.5;

        int nl2 = 2*nside_;
        z = (nl2-iring)*fact1_;
        phi = (iphi-fodd) * M_PI/nl2;
    }
    else /* South Polar cap */
    {
        int ip = npix_ - pix;
        int iring = (int)(0.5*(1+isqrt(2*ip-1))); /* counted from South pole */
        int iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

        z = -1.0 + (iring*iring)*fact2_;
        phi = (iphi-0.5) * 0.5*M_PI/iring;
    }

    double v[3];
    v[0] = sin(acos(z))*cos(phi);
    v[1] = sin(acos(z))*sin(phi);
    v[2] = z;

    /* rotate the vector to break degeneracy with grid axes */
    vec[0] = 0.5*v[0] -0.5*v[1] +0.70710678*v[2];
    vec[1] = 0.85355339*v[0] +0.14644661*v[1] -0.5*v[2];
    vec[2] = 0.14644661*v[0] +0.85355339*v[1] +0.5*v[2];
}

/* make params struct for functions we will integrate with gsl */
typedef struct kernel_params {
    double m[3];
    int n;
} kernel_params;

double fermi_dirac_kernel_vol_1(double x)
{
    /* Functional form of FD distribution for 1 particle
       over the full volume (i.e. no *x^2) */
    return 1. / (exp(x) + 1);
}

double fermi_dirac_kernel_vol(double x, void * p)
{
    /* FD kernel in the case of multiple particles
       over the full volume (i.e. no *x^2) */
    struct kernel_params * params = (struct kernel_params *)p;

    double result = fermi_dirac_kernel_vol_1(x);
    // only loop from i=1. i=0 term has been included in result already (r=1)
    for(int i = 1; i < params->n; i ++) {
        double r = params->m[i] / params->m[0];
        result += (r*r*r*r) * fermi_dirac_kernel_vol_1(x * r);
    }
    return result;
}

/* now define the 1D kernel and moments */
double fermi_dirac_kernel(double x, void * p)
{
    return x*x * fermi_dirac_kernel_vol(x, p);
}

/* Samples the low velocity end more effectively */
double low_vel_kernel(double x, void * p)
{
    return x * fermi_dirac_kernel_vol(x, p);
}

/* Needed for calculatin the velocity dispersion */
double fermi_dirac_dispersion(double x, void * p)
{
     return x*x*x*x * fermi_dirac_kernel_vol(x, p);
}


void 
divide_fd(double *vel_table, double *mass, 
          kernel_params * params, 
          int n_shells, 
          int lvk)
{
    double fermi_dirac_vel_ncdm[LENGTH_FERMI_DIRAC_TABLE];
    double fermi_dirac_cdf_ncdm[LENGTH_FERMI_DIRAC_TABLE]; //stores CDF

    int i,j;
    double v_bin,u,bin_average,bin_mass;
    double vel_edges[n_shells];
    double result,error;
    gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (1000);

    gsl_function F,G,H;
    if(lvk)
        F.function = &low_vel_kernel;
    else
        F.function = &fermi_dirac_kernel;
    G.function = &fermi_dirac_kernel;
    H.function = &fermi_dirac_dispersion;
    
    F.params = (void*) params;
    G.params = (void*) params;
    H.params = (void*) params;

    // Initialize the CDF table
    for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++){
        fermi_dirac_vel_ncdm[i] = MAX_FERMI_DIRAC * i / (LENGTH_FERMI_DIRAC_TABLE - 1.0);
        gsl_integration_qags (&F, 0, fermi_dirac_vel_ncdm[i], 0,1e-7 , 1000,
			  w, &result, &error);
        fermi_dirac_cdf_ncdm[i] = result;
    }

    //Normalize to 1
    for(i = 0; i < LENGTH_FERMI_DIRAC_TABLE; i++)
        fermi_dirac_cdf_ncdm[i] /= fermi_dirac_cdf_ncdm[LENGTH_FERMI_DIRAC_TABLE - 1];


    //Define the edges of the velocity bins (lower edge held to 0)
    for(i=0;i<n_shells;i++){
        v_bin = (i+1)/(n_shells*1.0);
        j=0;
        while(j < LENGTH_FERMI_DIRAC_TABLE - 2)
            if(v_bin > fermi_dirac_cdf_ncdm[j + 1])
                j++;
            else
                break;

        u = (v_bin - fermi_dirac_cdf_ncdm[j]) / (fermi_dirac_cdf_ncdm[j + 1] - fermi_dirac_cdf_ncdm[j]);

        vel_edges[i] = fermi_dirac_vel_ncdm[j] * (1 - u) + fermi_dirac_vel_ncdm[j + 1] * u;
    }

    //Get the bin velocities and bin masses
    double total_mass;
    gsl_integration_qags (&G, 0, fermi_dirac_vel_ncdm[LENGTH_FERMI_DIRAC_TABLE-1], 0, 1e-7, 1000,
			  w, &result, &error);
    total_mass = result;
    for(i=0;i<n_shells;i++){
        //Special case for lowest bin - left edge has to be zero
        if (i==0){
            gsl_integration_qags (&H, 0,vel_edges[i], 0, 1e-7, 1000,
			  w, &result, &error);
            bin_average = result;
            gsl_integration_qags (&G, 0, vel_edges[i], 0, 1e-7, 1000,
			  w, &result, &error);
            bin_average /= result;

            bin_mass = result/total_mass;
        }
        else{
            gsl_integration_qags (&H, vel_edges[i-1], vel_edges[i], 0, 1e-7, 1000,
			  w, &result, &error);
            bin_average = result;
            gsl_integration_qags (&G, vel_edges[i-1], vel_edges[i], 0, 1e-7, 1000,
			  w, &result, &error);
            bin_average /= result;

            bin_mass = result/total_mass;
        }
        vel_table[i] = sqrt(bin_average);
        mass[i] = bin_mass;
    }

}


void divide_sphere_healpix(double *vec_table, int n_side)  //so is the direction totally fixed? could we rotate the divisions relative to north pole (i.e. relative to the box)?
{
    int i,k;
    double v_sq[3];
    for(k=0;k<3;k++)
        v_sq[k] = 0.;
    
    // Get unit vectors for all pixels
    for(i=0;i<12*n_side*n_side;i++){
        pix2vec(i,&vec_table[i*3],n_side);
        for(k=0;k<3;k++)
            v_sq[k] += vec_table[i*3+k]*vec_table[i*3+k];
    }
    // Isotropize the velocity dispersion - important for low n_side    [should we remove this if for high n_side???]
    for(k=0;k<3;k++){
        v_sq[k] /= 12*n_side*n_side;
        v_sq[k] /= 1./3.; // Set each direction to have 1/3 of total dispersion
    }
    for(i=0;i<12*n_side*n_side;i++)
        for(k=0;k<3;k++)
            vec_table[i*3+k] /= sqrt(v_sq[k]);
}

void divide_sphere_fibonacci(double *vec_table, int n_side)
{
    int n_fibonacci = n_side;
    double lat, lon;
    int i, j;
    //int N_tot = 2 * n_fibonacci + 1;

    for (i=-n_fibonacci;i<n_fibonacci+1;i++){
        lat = asin(2.0*i/(2.0*n_fibonacci+1));
        lon = 2.0*M_PI*i*2.0/(1.0+sqrt(5.0));
        j = i + n_fibonacci;
        vec_table[j*3+0] = cos(lat)*cos(lon);
        vec_table[j*3+1] = cos(lat)*sin(lon);
        vec_table[j*3+2] = sin(lat);
    }
}

void divide_sphere(FastPMncdmSphereScheme ncdm_sphere_scheme,
                   double *vec_table, int n_side)
{
    switch (ncdm_sphere_scheme){
        case FASTPM_NCDM_SPHERE_HEALPIX:
            divide_sphere_healpix(vec_table, n_side);
        break;
        case FASTPM_NCDM_SPHERE_FIBONACCI:
            divide_sphere_fibonacci(vec_table, n_side);
        break;
        default:
            fastpm_raise(-1, "Wrong ncdm sphere scheme.\n");
    }
}

static void _fastpm_ncdm_init_fill(FastPMncdmInitData* nid);

//create and fill the table
FastPMncdmInitData* 
fastpm_ncdm_init_create(
    double BoxSize,
    FastPMCosmology * c,
    double z,
    int n_shells,
    int n_side,
    int lvk,
    FastPMncdmSphereScheme ncdm_sphere_scheme)
{
    FastPMncdmInitData* nid = malloc(sizeof(nid[0]));

    nid->BoxSize = BoxSize;
    nid->m_ncdm_sum = 0;
    for(int i = 0; i < c->N_ncdm; i ++) {
        nid->m_ncdm[i] = c->m_ncdm[i];
        nid->m_ncdm_sum += c->m_ncdm[i];
    }
    nid->n_ncdm = c->N_ncdm;
    nid->Omega_ncdm = c->Omega_ncdm;
    nid->z = z;
    fastpm_info("ncdm reference redshift = %g\n", z);
    nid->n_shells = n_shells;
    nid->n_side = n_side;
    nid->lvk = lvk;
    nid->ncdm_sphere_scheme = ncdm_sphere_scheme;
    
    /* this is the total number of velocity vectors produced */
    switch (ncdm_sphere_scheme){
        case FASTPM_NCDM_SPHERE_HEALPIX:
            nid->n_sphere = 12 * n_side * n_side;
        break;
        case FASTPM_NCDM_SPHERE_FIBONACCI:
            nid->n_sphere = 2 * n_side + 1;
        break;
        default:
            fastpm_raise(-1, "Wrong ncdm sphere scheme.\n");
    }
    
    nid->n_split = n_shells * nid->n_sphere;
    
    /* recall vel is a pointer to a 3 element array
    so lets make into multidim array of dim (n_split x 3) */
    //for now store all (3) neutrinos in the same nid table.
    nid->vel = calloc(nid->n_split, sizeof(double[3]));
    nid->mass = calloc(nid->n_split, sizeof(double));   //a 1d array to store all the masses
    
    _fastpm_ncdm_init_fill(nid);
    
    return nid;
}

//free the table
void 
fastpm_ncdm_init_free(FastPMncdmInitData* nid)
{
    free(nid->vel);
    free(nid->mass);
    free(nid);
}

//append the table
static void
_fastpm_ncdm_init_fill(FastPMncdmInitData* nid)
{   
    int n_ncdm = nid->n_ncdm;
    int n_shells = nid->n_shells;
    int n_side = nid->n_side;
    int n_sphere = nid->n_sphere;
    int lvk = nid->lvk;
    
    double *vel_table, *vec_table, *masstab;
    vel_table = malloc(sizeof(double)*n_shells);
    /* masstab is the distribution integrated per shell, sums to 1 */
    masstab = calloc(n_shells,sizeof(double));
    vec_table = malloc(sizeof(double)*n_sphere*3);
    
    /* define the kernel params for dividing fd */
    kernel_params* params = malloc(sizeof(params[0]));
    for(int i = 0; i < n_ncdm; i ++) {
        params->m[i] = nid->m_ncdm[i];
    }
    params->n = n_ncdm;
    
    divide_fd(vel_table, masstab, params, n_shells, lvk);
    divide_sphere(nid->ncdm_sphere_scheme, vec_table, n_side);
    
    /* velocity_conversion_factor converts the dimless exponent
       in the FD distribution to velocity in FastPM units:
       conjugate momentum (a^2 xdot, where x is comoving dist).
       kTc = 50.3 eV/c^2 km/s */
    double m0 = nid->m_ncdm[0];
    double velocity_conversion_factor = 50.3 / m0 / HubbleConstant;

    int i, j, d;
    int s = 0;                          //row num  (s for split index)
    for(i = 0; i < n_sphere; i ++){
        for(j = 0; j < n_shells; j ++){
            /* define mass st sum over split gives 1 */
            nid->mass[s] = masstab[j] / n_sphere;
            for(d = 0; d < 3; d ++){
                nid->vel[s][d] = vel_table[j] * vec_table[i*3 + d] * velocity_conversion_factor;
            }
            s++;
        }
    }
    free(vel_table); 
    free(masstab);
    free(vec_table);
}

void
fastpm_split_ncdm(FastPMncdmInitData * nid,
                  FastPMStore * src,
                  FastPMStore * dest,
                  MPI_Comm comm)
{
    /*
    Takes store src and splits according to the ncdm 
    init data, and uses this to populate an dest store.
    src is fully populated after calling the setup_lpt 
    func, but dest has only been init'd.
    */
    
    dest->np = src->np * nid->n_split;    // FIXME: ? Is proc imbalance an issue that needs to be accounted for?
    if(dest->np > dest->np_upper) {
        fastpm_raise(-1, "exceeding limit on the number limit of particles for the ncdm species \n");
    }

    size_t np_total = fastpm_store_get_np_total(src, comm);     //np_total is num of UNsplit ncdms

    /* avg mass of UNsplit ncdm site */
    double M0 = nid->Omega_ncdm * FASTPM_CRITICAL_DENSITY * pow(nid->BoxSize, 3) / np_total;
    double Mestate = M0 / nid->n_ncdm;
    fastpm_info("average mass of a ncdm particle is %g\n", Mestate);

    /* copy and amend meta-data */
    memmove(&dest->meta, &src->meta, sizeof(src->meta));

    dest->meta.M0 = 0.; /* will have a mass per particle */

    /* Compute the factor used to displace ncdm so that the expanded spheres almost touch.
       This assumes the largest speed particle is in the final entry of vel table */
    int n_shells = nid->n_shells;
    int n_sphere = nid->n_sphere;
    double vthm_max = 0;
    for (int d = 0; d < 3; d ++)
        vthm_max += nid->vel[n_sphere*n_shells-1][d]*nid->vel[n_sphere*n_shells-1][d];
    vthm_max = sqrt(vthm_max);
    double disp_factor = 0.5 * nid->BoxSize / nid->n_ncdm / vthm_max * (n_shells - 1) / n_shells;

    ptrdiff_t i, s, d;
    ptrdiff_t r = 0;    //row index of dest
    int c;
    double vthm;
    for(i = 0; i < src->np; i ++) {    //loop thru each src particle
        for(s = 0; s < nid->n_split; s ++) {  //loop thru split velocities
            /* copy all cols el by el */
            for(c = 0; c < 32; c ++) {
                if (!dest->columns[c] || !src->columns[c]) continue;

                size_t elsize = dest->_column_info[c].elsize;
                memcpy(dest->columns[c] + r * elsize,
                       src->columns[c] + i * elsize,
                       elsize);
            }

            /* give id, mass and add thm vel */
            dest->id[r] = s * src->meta._q_size + src->id[i];
            dest->mass[r] = nid->mass[s] * M0;    // ensures sum over split gives M0

            for(d = 0; d < 3; d ++){
                vthm = nid->vel[s][d];
                dest->v[r][d] = vthm;
                dest->x[r][d] += vthm * disp_factor;
                if (dest->q) dest->q[r][d] += vthm * disp_factor;
            }
            r ++;
        }
    }
}

