#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include <fastpm/thermalvelocity.h>

#define LENGTH_FERMI_DIRAC_TABLE 4000 //sets the length of the table on which CDF
                                      //will be evaluated
#define MAX_FERMI_DIRAC          20.0 //Upper limit of F-D distribution in units of
                                      // p/T

#define NEXT(n, i)  (((n) + (i)/(n)) >> 1) // Needed for Healpix routine

//should i initialize all funcs here?????????????????????

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

  vec[0] = sin(acos(z))*cos(phi);
  vec[1] = sin(acos(z))*sin(phi);
  vec[2] = z;
}


// Functional form of F-D distribution
double fermi_dirac_kernel(double x)
{
  return x * x / (exp(x) + 1);
}

// Samples the low velocity end more effectively
double low_vel_kernel(double x)
{
  return x / (exp(x) + 1);
}

// Needed for calculatin the velocity dispersion
double fermi_dirac_dispersion(double x)
{
  return x * x * x * x/ (exp(x) + 1);
}




void divide_fd(double *vel_table, double *mass, int n_shells)
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
  F.function = &low_vel_kernel;
  G.function = &fermi_dirac_kernel;
  H.function = &fermi_dirac_dispersion;

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
  gsl_integration_qags (&G, 0, fermi_dirac_vel_ncdm[LENGTH_FERMI_DIRAC_TABLE-1], 0,1e-7 , 1000,
			  w, &result, &error);
  total_mass = result;
  for(i=0;i<n_shells;i++){
    //Special case for lowest bin - left edge has to be zero
     if (i==0){
       gsl_integration_qags (&H, 0,vel_edges[i], 0,1e-7 , 1000,
			  w, &result, &error);
       bin_average = result;
       gsl_integration_qags (&G, 0, vel_edges[i], 0,1e-7 , 1000,
			  w, &result, &error);
       bin_average /= result;

       bin_mass = result/total_mass;
     }
     else{
       gsl_integration_qags (&H, vel_edges[i-1],vel_edges[i], 0,1e-7 , 1000,
			  w, &result, &error);
       bin_average = result;
       gsl_integration_qags (&G, vel_edges[i-1], vel_edges[i], 0,1e-7 , 1000,
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
  // Isotropize the velocity dispersion - important for low n_side
  for(k=0;k<3;k++){
    v_sq[k] /= 12*n_side*n_side;
    v_sq[k] /= 1./3.; // Set each direction to have 1/3 of total dispersion
  }
  for(i=0;i<12*n_side*n_side;i++)
    for(k=0;k<3;k++)
      vec_table[i*3+k] /= sqrt(v_sq[k]);
}

void divide_sphere_fibonacci(double *vec_table, int n_fibonacci)
{
  double lat, lon;
  int i,j,k;
  int N_tot = 2*n_fibonacci+1;


  for (i=-n_fibonacci;i<n_fibonacci+1;i++){
    lat = asin(2.0*i/(2.0*n_fibonacci+1));
    lon = 2.0*M_PI*i*2.0/(1.0+sqrt(5.0));
    j = i + n_fibonacci;
    vec_table[j*3+0] = cos(lat)*cos(lon);
    vec_table[j*3+1] = cos(lat)*sin(lon);
    vec_table[j*3+2] = sin(lat);

  }

}


//all below is for healpix... can add fibonacci later.!!!!!!!!!!!!

static void _fastpm_ncdm_init_fill(FastPMncdmInitData* nid);

//create and fill the table
FastPMncdmInitData* 
fastpm_ncdm_init_create(double m_ncdm, double z, int n_shells, int n_side)
{
    FastPMncdmInitData* nid = malloc(sizeof(nid[0]));

    nid->m_ncdm = m_ncdm;
    nid->z = z;
    nid->n_shells = n_shells;    //is this ok or should i point to these for now?
    nid->n_side = n_side;
    
    int n_split = 12 * n_shells * n_side*n_side;     //this is the total number of velocity vectors produced
    /* recall vel is a pointer to a 3 element array
    so lets make into multidim array of dim (n_split x 3) */
    nid->vel = calloc(n_split, sizeof(double[3]));
    nid->mass = calloc(n_split, sizeof(double));   //a 1d array to store all the masses.

    _fastpm_ncdm_init_fill(nid);
    
    return nid;
}

//free the table
void 
fastpm_ncdm_init_free(FastPMncdmInitData* nid)
{   //hmmmm seems overkill, but maybe useful to give this a func name for when we use in another module
    free(nid->vel);
    free(nid->mass);
    free(nid);
}

//append the table
static void
_fastpm_ncdm_init_fill(FastPMncdmInitData* nid)    ///call in create.  no need for it otherwise.
{   
    int n_shells = nid->n_shells;
    int n_side = nid->n_side;
    
    double *vel_table, *vec_table, *masstab;
    vel_table = malloc(sizeof(double)*n_shells);
    masstab = calloc(n_shells,sizeof(double));
    vec_table = malloc(sizeof(double)*12*n_side*n_side*3);
    
    divide_fd(vel_table,masstab,n_shells);

    divide_sphere_healpix(vec_table,n_side);

    double velocity_conversion_factor = 50.3*(1.+ nid->z)*(1./nid->m_ncdm); //In km/s in Gadget internal units
    //velocity_conversion_factor *= sqrt(1. + redshift); // Now in Gadget I/O units

    int i, j, k;
    int r=0;                          //row num
    for(i=0;i<12*n_side*n_side;i++){
        for(j=0;j<n_shells;j++){
            nid->mass[r] = masstab[j]/(12.*n_side*n_side);
            //printf("%g\n",masstab[j]/(12.*n_side*n_side));
            //printf("hiii%g\n",nid->mass[1]);
            for(k=0;k<3;k++){
                nid->vel[r][k] = vel_table[j]*vec_table[i*3+k]*velocity_conversion_factor;    // *(nid->vel+r)[k]
                //if (j == n_shells-1) printf("%d byyyyy%g\n",k, nid->mass[1]);
            }
            r++;
        }
    }
    
    //free all above table memories?
    free(vel_table); 
    free(masstab);
    free(vec_table);
}

//nid input store, output

//void
//fastpm_split_ncdm(FastPMncdmInitData* nid, FastPMStore * p, FastPMStore * q)
//{
    /*Takes store p, splits all particles according 
    to the ncdm init data, giving store q.*/
    

//}