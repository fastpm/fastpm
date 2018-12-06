#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>

#include <fastpm/libfastpm.h>

//include Fermi-Dirac integration table for neutrinos
#include <fastpm/Ftable.h>

#define STEF_BOLT 2.851e-48  // in units: [ h * (10^10Msun/h) * s^-3 * K^-4 ]
#define rho_crit 27.7455     //rho_crit0 in mass/length (not energy/length)
#define LIGHT 9.722e-15      // in units: [ h * (Mpc/h) * s^-1 ]      //need to update all these units and change Omega_g

#define kB 8.617330350e-5    //boltzman in eV/K   i th

double HubbleDistance = 2997.92458; /* Mpc/h */   /*this c*1e5 in SI units*/
double HubbleConstant = 100.0; /* Mpc/h / km/s*/     //OTHER WAY ROUND!

//move struct to h
//typedef struct {
//    double h;
//    double Omega_cdm;
//    double Omega_Lambda;
//    double T_cmb;    /*related to omegaR*/
//    double N_eff;  //N_ur;      /*this is N_eff*/ //actually i might just make this number o fmassless neutrinos, just adds to rad
//    double M_nu; ///    m_ncdm[3]; for now assume 3 nus of same mass
//    int N_nu;  //N_ncdm;
//} FastPMCosmology;

double interpolate(const double xa[], const double ya[], size_t size, double xi)
{
    /*
    Linearly interpolate at xi from data xa[] and ya[]
    */
    
    gsl_interp* interp 
        = gsl_interp_alloc(gsl_interp_linear, size);  //linear
    gsl_interp_accel* acc 
        = gsl_interp_accel_alloc();
    
    gsl_interp_init(interp, xa, ya, size);
    
    double yi = gsl_interp_eval(interp, xa, ya, xi, acc);
    
    gsl_interp_free(interp);
    gsl_interp_accel_free (acc);
    
    return yi;
}

double Omega_g(FastPMCosmology * c)   
{
    /* This is really Omega_g0. All Omegas in this code mean 0, exceot nu!!, I would change notation, but anyway 
    Omega_g0 = 4 \sigma_B T_{CMB}^4 8 \pi G / (3 c^3 H0^2)*/
    return 4 * STEF_BOLT * pow(c->T_cmb, 4) / pow(LIGHT, 3) / rho_crit / pow(c->h, 2);
}

double getFtable(int F_id, double y)
{
    //Not using size as an arg, it's globally defined
    //Gets the interpolated value, unless y=0 for which result is hard coded
    //F_id: 1 for F, 2 for F', 3 for F''
    //F at y means interp at x ... kind of confusing mixing x and y, but whatever.
    
    if (y == 0.) {               //does this work for double?
        if (F_id == 1) {
            return 5.6822;    //check
        }else{
            return 0.;
        }
    }else{
        return interpolate(Ftable[0], Ftable[F_id], Fsize, y);
    }
}

double Fconst(FastPMCosmology * c)
{
    /*
    This is a cosmology dependent constant which is the argument divided by a of F, DF, DDF used repeatedly in code
    To be clear, evaluate F at Fconst*a
    */
    
    if (c->T_cmb == 0.){
        return 0.;                 //really get infinity, but no gamma means no nu, so effectively F_const=0. Exact value doesnt matter, because the Omega_g factor in Omega_nu will make it 0 regarldess of  what F gives, just need to avoif nan.
    }else{
        double Gamma_nu = pow(c->N_eff / c->N_nu, 1./4.) * pow( 4./11., 1./3.);   //nu to photon temp ratio today
        double T_nu = Gamma_nu * c->T_cmb;
        return c->M_nu / c->N_nu / (kB * T_nu);
    }
}


double OmegaNuTimesHubbleEaSq(double a, FastPMCosmology * c)   
{
    /* Use interpolation to find Omega_nu(a) * E(a)^2 */
    
    double Fc = Fconst(c);
    double F = getFtable(1, Fc*a);    //row 1 for F
    double A = 15. / pow(M_PI, 4) * pow(4./11., 4./3.) * c->N_eff * Omega_g(c);
    
    return A / (a*a*a*a) * F;
}

//lots of repn in the below funcs, should we define things globally (seems messy), or object orient?
double DOmegaNuTimesHubbleEaSqDa(double a, FastPMCosmology * c)
{
    double Fc = Fconst(c);
    double DF = getFtable(2, Fc*a);   //row 2 for F'
    double A = 15. / pow(M_PI, 4) * pow(4./11., 4./3.) * c->N_eff * Omega_g(c);
    
    double OnuESq = OmegaNuTimesHubbleEaSq(a,c);
    
    return -4. / a * OnuESq + A / (a*a*a*a) * Fc * DF;
}

double D2OmegaNuTimesHubbleEaSqDa2(double a, FastPMCosmology * c)
{
    double Fc = Fconst(c);
    double DDF = getFtable(3, Fc*a);   //row 3 for F''
    double A = 15. / pow(M_PI, 4) * pow(4./11., 4./3.) * c->N_eff * Omega_g(c);
    
    double OnuESq = OmegaNuTimesHubbleEaSq(a,c);
    double DOnuESqDa = DOmegaNuTimesHubbleEaSqDa(a,c);
    
    return -12. / (a*a) * OnuESq - 8. / a * DOnuESqDa + A / (a*a*a*a) * (Fc*Fc) * DDF;
}

double HubbleEa(double a, FastPMCosmology * c)
{
    /* H(a) / H0 */
    return sqrt(Omega_g(c) / (a*a*a*a) + c->Omega_cdm / (a*a*a) + OmegaNuTimesHubbleEaSq(a, c) + c->Omega_Lambda);
}

double Omega_cdm_a(double a, FastPMCosmology * c)
{
    //as a func of a. I dont like this. I would use 0s in struct instead. ask yu.
    double E = HubbleEa(a, c);
    return c->Omega_cdm / (a*a*a) / (E*E);
}

double OmegaA(double a, FastPMCosmology * c) {
    //this is what I called Omega_cdm_a above. choose which to use later.
    return Omega_cdm_a(a, c);
}

double DHubbleEaDa(double a, FastPMCosmology * c)
{
    /* d E / d a*/
    double E = HubbleEa(a, c);
    double DOnuESqDa = DOmegaNuTimesHubbleEaSqDa(a,c);
    
    return 0.5 / E * ( - 4 * Omega_g(c) / pow(a,5) - 3 * c->Omega_cdm / pow(a,4) + DOnuESqDa );
}

double D2HubbleEaDa2(double a, FastPMCosmology * c)
{
    double E = HubbleEa(a,c);
    double dEda = DHubbleEaDa(a,c);
    double D2OnuESqDa2 = D2OmegaNuTimesHubbleEaSqDa2(a,c);

    return 0.5 / E * ( 20 * Omega_g(c) / pow(a,6) + 12 * c->Omega_cdm / pow(a,5) + D2OnuESqDa2 - 2 * pow(dEda,2) );
}

double OmegaSum(double a, FastPMCosmology* c)
{
    //should always equal 1. good for testing.
    double sum = Omega_g(c) / pow(a, 4);
    sum += c->Omega_cdm / pow(a, 3);
    sum += OmegaNuTimesHubbleEaSq(a, c);
    sum += c->Omega_Lambda;
    return sum / pow(HubbleEa(a, c), 2);
}

static int growth_ode(double a, const double y[], double dyda[], void *params)
{
    //is yy needed?
    FastPMCosmology* c = (FastPMCosmology * ) params;
    
    const double E = HubbleEa(a, c);
    const double dEda = DHubbleEaDa(a, c);
    
    double dydlna[4];
    dydlna[0] = y[1];
    dydlna[1] = - (2. + a / E * dEda) * y[1] + 1.5 * Omega_cdm_a(a, c) * y[0];
    dydlna[2] = y[3];
    dydlna[3] = - (2. + a / E * dEda) * y[3] + 1.5 * Omega_cdm_a(a, c) * (y[2] - y[0]*y[0]);
    
    //divide by  a to get dyda
    for (int i=0; i<4; i++){
        dyda[i] = dydlna[i] / a;
    }
    
    return GSL_SUCCESS;
}

//moved to .h
//typedef struct{
//    double y0;
//    double y1;
//    double y2;
//    double y3;
//} ode_soln;

static ode_soln growth_ode_solve(double a, FastPMCosmology * c)
{
    /* NOTE that the analytic COLA growthDtemp() is 6 * pow(1 - c.OmegaM, 1.5) times growth() */
    /* This returns an array of {d1, F1, d2, F2}
        Is there a nicer way to do this within one func and no object?*/

    gsl_odeiv2_system F;
    F.function = &growth_ode;
    F.jacobian = NULL;
    F.dimension = 4;
    F.params = (void*) c;
    
    gsl_odeiv2_driver * drive 
        = gsl_odeiv2_driver_alloc_standard_new(&F,
                                               gsl_odeiv2_step_rkf45, 
                                               1e-5, 
                                               1e-8,
                                               1e-8,
                                               1,
                                               1);
    
     /* We start early to avoid lambda. ??????????????*/
    double aini = 1e-5;
    
    //Note the normalisation of D is arbitrary as we will only use it to calcualte growth fractor.
    
    //MD initial conditions for now.
    double yini[4] = {aini, aini, -3./7.*aini*aini, -6./7.*aini*aini};
    //double yini[4] = {1e7*aini, 1e52*aini, -3./7.*aini*aini*1e97, -6./7.*aini*aini/1e50};
   
    //int stat = 
    gsl_odeiv2_driver_apply(drive, &aini, a, yini);
    //if (stat != GSL_SUCCESS) {     //If succesful, stat will = GSL_SUCCESS and yinit will become the final values...
        //endrun(1,"gsl_odeiv in growth: %d. Result at %g is %g %g\n",stat, curtime, yinit[0], yinit[1]); //need to define endrun
    //    printf(stat);    //quick fix for now 
    //}
    
    gsl_odeiv2_driver_free(drive);
    /*Store derivative of D if needed.*/
    //if(dDda) {
    //    *dDda = yinit[1]/pow(a,3)/(hubble_function(a)/CP->Hubble);
    //}
    
    ode_soln soln;
    soln.y0 = yini[0];
    soln.y1 = yini[1];
    soln.y2 = yini[2];
    soln.y3 = yini[3];
    
    return soln;
}

double growth(double a, FastPMCosmology * c) {
    //d1
    return growth_ode_solve(a, c).y0;
}

double DgrowthDlna(double a, FastPMCosmology * c) {
    //F1
    return growth_ode_solve(a, c).y1;
}

double growth2(double a, FastPMCosmology * c) {
    //d2
    return growth_ode_solve(a, c).y2;
}

double Dgrowth2Dlna(double a, FastPMCosmology * c) {
    //F2
    return growth_ode_solve(a, c).y3;
}

double GrowthFactor(double a, FastPMCosmology * c) { // growth factor for LCDM
    return growth(a, c) / growth(1., c);           //[this is D(a)/D_today for LCDM]
}

double DLogGrowthFactor(double a, FastPMCosmology * c) {
    /* dlnD1/dlna */
    double d1 = growth(a, c);
    double F1 = DgrowthDlna(a, c);
    
    return F1 / d1;
}

double GrowthFactor2(double a, FastPMCosmology * c) {
    /* Normalised D2. Is this correct???*/
    double d0 = growth(1., c);  //0 for today
    return growth2(a, c) / growth2(1., c); //(d0*d0); // ??????????????????;
}

double DLogGrowthFactor2(double a, FastPMCosmology * c) {
    /* dlnD2/dlna */
    double d2 = growth2(a, c);
    double F2 = Dgrowth2Dlna(a, c);
    
    return F2 / d2;
}

double DGrowthFactorDa(double a, FastPMCosmology * c) {
    double d0 = growth(1., c);
    double F1 = DgrowthDlna(a, c);
    return F1 / a / d0;
}

double D2GrowthFactorDa2(double a, FastPMCosmology * c) {
    double E = HubbleEa(a, c);
    double dEda = DHubbleEaDa(a, c);
    
    double d1 = growth(a, c);
    double F1 = DgrowthDlna(a, c);
    double d0 = growth(1., c);  //0 for today
    
    double ans = 0.;
    ans -= (3. + a / E * dEda) * F1;
    ans += 1.5 * Omega_cdm_a(a, c) * d1;
    ans /= d0 * a*a;
    return ans;
}

static double
comoving_distance_int(double a, void * params)
{
    FastPMCosmology * c = (FastPMCosmology * ) params;
    return 1. / (a * a * HubbleEa(a, c));
}

double ComovingDistance(double a, FastPMCosmology * c) {

    /* We tested using ln_a doesn't seem to improve accuracy */

    int WORKSIZE = 100000;

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;

    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = &comoving_distance_int;
    F.params = (void*) c;

    gsl_integration_qag(&F, a, 1., 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41,
            workspace, &result, &abserr); 
            //lowered tol by /10 to avoid error from round off (maybe I need to make my paras more accurate)

    gsl_integration_workspace_free(workspace);

    return result;
}

#ifdef TEST_COSMOLOGY
int main() {
    /* the old COLA growthDtemp is 6 * pow(1 - c.OmegaM, 1.5) times growth */
    double a;
    FastPMCosmology c[1] = {{
        .OmegaM = 0.3,
        .OmegaLambda = 0.7
    }};

    printf("OmegaM D dD/da d2D/da2 D2 E dE/dA d2E/da2 \n");
    for(c->OmegaM = 0.1; c->OmegaM < 0.6; c->OmegaM += 0.1) {
        double f = 6 * pow(1 - c->OmegaM, 1.5);
        c->OmegaLambda = 1 - c->OmegaM;
        double a = 0.8;
        printf("%g %g %g %g %g %g %g %g %g\n",
            c->OmegaM, 
            ComovingDistance(a, c),
            GrowthFactor(a, c),
            DGrowthFactorDa(a, c),
            D2GrowthFactorDa2(a, c),

            GrowthFactor2(a, c),
            HubbleEa(a, c),
            DHubbleEaDa(a, c),
            D2HubbleEaDa2(a, c)
            );
    }
}


#endif
