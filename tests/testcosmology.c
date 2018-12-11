#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>

//include Fermi-Dirac integration table for neutrinos
#include <fastpm/Ftable.h>  //""

//now that I moved this code over to cosmology.c, need to include cosmology.h here in testcosmology.c in order to avoid onfilciting type errors (because there would otherwise be multiple defns of the same thing!)
//#include <fastpm/cosmology.h>

//double HubbleDistance = 2997.92458; /* Mpc/h */   /*this c*1e5 in SI units*/
//double HubbleConstant = 100.0; /* Mpc/h / km/s*/     //OTHER WAY ROUND!

#define STEF_BOLT 2.851e-48  // in units: [ h * (10^10Msun/h) * s^-3 * K^-4 ]
#define rho_crit 27.7455     //rho_crit0 in mass/length (not energy/length)
#define LIGHT 9.722e-15      // in units: [ h * (Mpc/h) * s^-1 ]      //need to update all these units and change Omega_g

#define kB 8.617330350e-5    //boltzman in eV/K   i think these units work best as we'll define neutrino mass in eV

/*Define FastPMCosmology Class for Nu code.
Had to change the syntax comapred to cosmology.h version use typedef and put name at end.
This removes error
*/

typedef FastPMCosmology FastPMCosmologyNu;


double HubbleEaNu(double a, FastPMCosmologyNu * c)
{
    /* H(a) / H0 */
    return sqrt(Omega_g(c) / (a*a*a*a) + c->Omega_cdm / (a*a*a) + Omega_ncdmTimesHubbleEaSq(a, c) + c->Omega_Lambda);
}

double DHubbleEaDaNu(double a, FastPMCosmologyNu * c)
{
    /* d E / d a*/
    double E = HubbleEaNu(a, c);
    double DOnuESqDa = DOmega_ncdmTimesHubbleEaSqDa(a,c);
    
    return 0.5 / E * ( - 4 * Omega_g(c) / pow(a,5) - 3 * c->Omega_cdm / pow(a,4) + DOnuESqDa );
}

double D2HubbleEaDa2Nu(double a, FastPMCosmologyNu * c)
{
    double E = HubbleEaNu(a,c);
    double dEda = DHubbleEaDaNu(a,c);
    double D2OnuESqDa2 = D2Omega_ncdmTimesHubbleEaSqDa2(a,c);

    return 0.5 / E * ( 20 * Omega_g(c) / pow(a,6) + 12 * c->Omega_cdm / pow(a,5) + D2OnuESqDa2 - 2 * pow(dEda,2) );
}


static int growth_odeNu(double a, const double y[], double dyda[], void *params)
{
    //is yy needed?
    FastPMCosmologyNu* c = (FastPMCosmologyNu * ) params;
    
    const double E = HubbleEaNu(a, c);
    const double dEda = DHubbleEaDaNu(a, c);
    
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

//typedef struct{
//    double y0;
//    double y1;
//    double y2;
//    double y3;
//} ode_soln;

static ode_soln growth_ode_solveNu(double a, FastPMCosmologyNu * c)
{
    /* NOTE that the analytic COLA growthDtemp() is 6 * pow(1 - c.OmegaM, 1.5) times growth() */
    /* This returns an array of {d1, F1, d2, F2}
        Is there a nicer way to do this within one func and no object?*/

    gsl_odeiv2_system FF;
    FF.function = &growth_odeNu;
    FF.jacobian = NULL;
    FF.dimension = 4;
    FF.params = (void*) c;
    
    gsl_odeiv2_driver * drive 
        = gsl_odeiv2_driver_alloc_standard_new(&FF,
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

double growthNu(double a, FastPMCosmologyNu * c) {
    //d1
    return growth_ode_solveNu(a, c).y0;
}

double DgrowthDlnaNu(double a, FastPMCosmologyNu * c) {
    //F1
    return growth_ode_solveNu(a, c).y1;
}

double growth2Nu(double a, FastPMCosmologyNu * c) {
    //d2
    return growth_ode_solveNu(a, c).y2;
}

double Dgrowth2DlnaNu(double a, FastPMCosmologyNu * c) {
    //F2
    return growth_ode_solveNu(a, c).y3;
}

double OmegaANu(double a, FastPMCosmologyNu * c) {
    //this is what I called Omega_cdm_a above. choose which to use later.
    return Omega_cdm_a(a, c);
}

double GrowthFactorNu(double a, FastPMCosmologyNu * c) { // growth factor for LCDM
    return growthNu(a, c) / growthNu(1., c);           //[this is D(a)/D_today for LCDM]
}

double DLogGrowthFactorNu(double a, FastPMCosmologyNu * c) {
    /* dlnD1/dlna */
    double d1 = growthNu(a, c);
    double F1 = DgrowthDlnaNu(a, c);
    
    return F1 / d1;
}

double GrowthFactor2Nu(double a, FastPMCosmologyNu * c) {
    /* Normalised D2. Is this correct???*/
    // double d0 = growthNu(1., c);  //0 for today
    return growth2Nu(a, c) / growth2Nu(1., c); //(d0*d0); // ??????????????????;
}

double DLogGrowthFactor2Nu(double a, FastPMCosmologyNu * c) {
    /* dlnD2/dlna */
    double d2 = growth2Nu(a, c);
    double F2 = Dgrowth2DlnaNu(a, c);
    
    return F2 / d2;
}

double DGrowthFactorDaNu(double a, FastPMCosmologyNu * c) {
    double d0 = growthNu(1., c);
    double F1 = DgrowthDlnaNu(a, c);
    return F1 / a / d0;
}

double D2GrowthFactorDa2Nu(double a, FastPMCosmologyNu * c) {
    double E = HubbleEaNu(a, c);
    double dEda = DHubbleEaDaNu(a, c);
    
    double d1 = growthNu(a, c);
    double F1 = DgrowthDlnaNu(a, c);
    double d0 = growthNu(1., c);  //0 for today
    
    double ans = 0.;
    ans -= (3. + a / E * dEda) * F1;
    ans += 1.5 * Omega_cdm_a(a, c) * d1;
    ans /= d0 * a*a;
    return ans;
}

static double
comoving_distance_intNu(double a, void * params)
{
    FastPMCosmologyNu * c = (FastPMCosmologyNu * ) params;
    return 1. / (a * a * HubbleEaNu(a, c));
}

double ComovingDistanceNu(double a, FastPMCosmologyNu * c) {

    /* We tested using ln_a doesn't seem to improve accuracy */

    int WORKSIZE = 100000;

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;

    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = &comoving_distance_intNu;
    F.params = (void*) c;

    gsl_integration_qag(&F, a, 1., 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41,
            workspace, &result, &abserr); 
            //lowered tol by /10 to avoid error from round off (maybe I need to make my paras more accurate)

    gsl_integration_workspace_free(workspace);

    return result;
}


int main(int argc, char * argv[]) {

    MPI_Init(&argc, &argv);

    libfastpm_init();

    MPI_Comm comm = MPI_COMM_WORLD;

    fastpm_set_msg_handler(fastpm_default_msg_handler, comm, NULL);

    /* do no rad first */
    FastPMCosmology c[1] = {{
        .h=0.6772,
        .Omega_cdm=0.3,
        .Omega_Lambda=0.7,
        .T_cmb=0.,
        .N_eff=3.046,
        .m_ncdm={0, 0, 0},               //(assuming 3 nus of mass 1ev, this is the sum of their masses)
        .N_nu=0,
        .N_ncdm=0
    }};

    printf("OmegaM X D dD/da d2D/da2 D2 E dE/dA d2E/da2 f1 f2 \n");
    for(c->Omega_cdm = 0.1; c->Omega_cdm < 0.6; c->Omega_cdm += 0.1) {
        c->Omega_Lambda = 1 - c->Omega_cdm;
        double a = 0.8;
        
        printf("%g %g %g %g %g %g %g %g %g %g %g\n",
            c->Omega_cdm, 
            ComovingDistance(a, c),
            //growth(a, c),
            GrowthFactor(a, c),
            DGrowthFactorDa(a, c),
            D2GrowthFactorDa2(a, c),

            GrowthFactor2(a, c),
            HubbleEa(a, c),
            DHubbleEaDa(a, c),
            D2HubbleEaDa2(a, c),
               
            DLogGrowthFactor2(a, c),
            DLogGrowthFactor(a, c)
            );
    }
    
    FastPMCosmology cNu[1] = {{     //Nu labels that this is a test with neutrinos in the simulation NOT using the FastPMCosmology object.
        .h=0.6772,
        .Omega_cdm=0.3,
        .Omega_Lambda=0.7,
        .T_cmb=2.73,
        .N_eff=3.046,
        .m_ncdm={1., 0, 0},               //(assuming 3 nus of mass 1ev, this is the sum of their masses)
        .N_nu=3,
        .N_ncdm=1
    }};

    printf("OmegaM X D dD/da d2D/da2 D2 E dE/dA d2E/da2 f1 f2 \n");
    for(cNu->Omega_cdm = 0.1; cNu->Omega_cdm < 0.6; cNu->Omega_cdm += 0.1) {
        cNu->Omega_Lambda = 1 - cNu->Omega_cdm - Omega_r(cNu) - Omega_ncdmTimesHubbleEaSq(1.,cNu);
        double a = 0.8;
        
        printf("%g %g %g %g %g %g %g %g %g %g %g\n",
            cNu->Omega_cdm, 
            ComovingDistance(a, cNu),
            //growth(a, c),
            GrowthFactor(a, cNu),
            DGrowthFactorDa(a, cNu),
            D2GrowthFactorDa2(a, cNu),

            GrowthFactor2(a, cNu),
            HubbleEa(a, cNu),
            DHubbleEaDa(a, cNu),
            D2HubbleEaDa2(a, cNu),
               
            DLogGrowthFactor2(a, cNu),
            DLogGrowthFactor(a, cNu)
            );
    }

    /*
    double xa[100];
    double ya[100];
    for (int i=0; i<100; i+=1){
        xa[i]=i;
        ya[i]=i;
    }
    printf("%g",interpolate(xa,ya,4.3));
    */
    libfastpm_cleanup();
    MPI_Finalize();
    return 0;
}
