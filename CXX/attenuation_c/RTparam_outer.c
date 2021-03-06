//Standard files;
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//My header files;
#include <RTparam_fast.h>
#include <RTparam_hardcoded.h>
#include <RTparam_outer.h>

#define PI 3.14159265358979
#define EPS 1.0e-12
#define MAXIT 100


#if 0

/*construct main function*/
int main() {

  //inputs to RTparam_outer;
  double h        = 2.0;   // thickness [m];
  double period   = 10.0;  // wave period [s];
  double visc_rp  = 13.0;  // damping parameter [Pa.(m.s^{-1})^{-1}];
  //
  double E     = 5.45e9;// Young's modulus [Pa];
  double g     = 9.81;  // gravity [ms^{-2}];
  double rhow  = 1025;  // density of water [kg.m^{-3}];
  double rho   = 0.9;   // = (rhoi/rhow), rhoi = density of ice;
  double nu    = 0.3;   // Poisson's ratio;
  //
  double rhoi,om,guess;
  double *params;

  //outputs from RTparam_outer;
  double ac,modT,argR,argT;
  double damping,kice,kwtr,int_adm;

  om     = 2*PI/period;
  rhoi   = rho*rhow;
  guess  = om*om/g;

  //allocate memory for pointer params, and assign values;
  params    = malloc( 5*sizeof(double) );
  params[0] = E;
  params[1] = g;
  params[2] = rhow;
  params[3] = rhoi;
  params[4] = nu;

  //call RTparam;
  RTparam_outer(&damping,&kice,&kwtr,&int_adm,
                &ac,&modT,&argR,&argT,
                h,om,visc_rp,guess,params);

  //free memory from params;
  free(params);

  //check results;
  printf("\noutput: damping   = %0.9e\n", damping);
  printf("\noutput: kice      = %0.9e\n", kice);
  printf("\noutput: kwtr      = %0.9e\n", kwtr);
  printf("\noutput: int_adm   = %0.9e\n", int_adm);
  printf(" ");
  printf("\noutput: ac        = %0.9e\n", ac);
  printf("\noutput: |T|       = %0.9e\n", modT);
  printf("\noutput: Arg[R]    = %0.9e\n", argR);
  printf("\noutput: Arg[T]    = %0.9e\n", argT);
}

#endif

//function [damping,kice,kwtr,int_adm,...
//            alp_scat,modT,argR,argT] =...
//               RT_param_outer(h,om,E,visc_rp,guess)
int RTparam_outer(double *damping,double *kice,double *kwtr,double *int_adm,
                  double *ac,double *modT,double *argR,double *argT,
                  double h,double om,double visc_rp,double guess,double *params) {
//int RTparam_outer(double *,double *,double *,double *,
//                  double *,double *,double *,double *,
//                  double,double,double,double);

   double alp_nd,h_nd,zeta_nd;
   double Hw_nd,varpi,L;
   double damping_nd,visc_rp_nd;
   double H_nd = 4.0;//infinite depth used in scattering calculation;
   //
   double ki,kw,avc,BG1,BG2;
   double E,g,rhow,rhoi,nu,rho;
   double D;

   E     = params[0];
   g     = params[1];
   rhow  = params[2];
   rhoi  = params[3];
   nu    = params[4];
   rho   = rhoi/rhow;
   //
   D        = E*(h*h*h)/12/(1-nu*nu);
   L        = exp( 0.2*log(D/rhow/om/om) );
   alp_nd   = om*om/g*L;
   h_nd     = h/L;
   zeta_nd  = rho*h_nd;

   //get wavenumber for ice;
   varpi = 1/alp_nd-zeta_nd;
   //[ki,BG2,avc]   = gen_root_ice(varpi,H_nd,guess*L);
   gen_root_ice(&ki,&BG2,&avc,varpi,H_nd,guess*L);
   *kice = ki/L;

   //get wavenumber for water;
   varpi = 1/alp_nd;
   Hw_nd = H_nd+zeta_nd;
   //[kw,BG1] = gen_root_wtr(varpi,Hw_nd,alp_nd);
   gen_root_wtr(&kw,&BG1,varpi,Hw_nd,alp_nd);
   *kwtr  = kw/L;

   //get intrinsic admittance;
   //|R|^2+int_adm*|T|^2=1
   printf("\noutput: BG1   = %0.9e\n", BG1);
   printf("\noutput: BG2   = %0.9e\n", BG2);
   *int_adm  = BG1/BG2;

   //get viscous attenuation;
   visc_rp_nd  = visc_rp/rhow/om/L;
   damping_nd  = avc*visc_rp_nd;
   *damping    = damping_nd/L;

   //get interpolated attenuation coefficient;
   RTparam_fast(ac,modT,argR,argT,
                alp_nd,h_nd,*int_adm);
}

/*******************************************************************/
//%%function [ki,BG2,avc]=gen_root_ice(del,H,guess)
//finds the root of the ice dispersion relation nearest to 'guess'.
//Newton-Rhapson method solver
//of f  = Lam*k*sinh(k*H)-cosh(k*H) = 0, Lam=k^4+del
int gen_root_ice(double *ki2, double *BG2,double *avc,
                 double del, double H,double guess) {
//int gen_root_ice(double *, double *,double *,
//                 double, double,double);

   double fac,k0,dk,ki;
   double res,Lam,Lampr,denom;

   fac   = 1.0;
   k0    = guess;

   //Call dispersion relation function;
   //dk    = NR_corr_term(k0,del,H,fac);
   NR_corr_term(&dk,&Lam,&Lampr,k0,del,H,fac);
   ki = k0-dk;

   while(fabs(dk) > EPS) {
     k0 = ki;

     //Call dispersion relation function;
     //dk = NR_corr_term(k0,del,H,fac);
     NR_corr_term(&dk,&Lam,&Lampr,k0,del,H,fac);
     ki = k0-dk;
   }

   //Call dispersion relation function;
   //[dk,Lam,Lampr] = NR_corr_term(ki,del,H,fac);
   NR_corr_term(&dk,&Lam,&Lampr,ki,del,H,fac);
   //
   denom = H*(Lam*Lam*ki*ki-1)+Lampr;
   res   = -ki/denom;

   //Outputs;
   *ki2  = ki;
   *BG2  = Lam*Lam*res;
   *avc  = ki/denom;// -1i*dk/d(visc_rp) - non-dimensional
}

/*******************************************************************/
//%%function [kw,BG1]=gen_root_wtr(del,H,guess)
//finds the root of the water dispersion relation nearest to 'guess'.
//Newton-Rhapson method solver
//of f  = del*k*sinh(k*H)-cosh(k*H) = 0
int gen_root_wtr(double *kw2, double *BG1,
                 double del, double H,double guess) {

   double fac,k0,dk,kw;
   double res,Lam,Lampr,denom;

   fac   = 0.0;
   k0    = guess;

   //Call dispersion relation function;
   //dk    = NR_corr_term(k0,del,H,fac);
   NR_corr_term(&dk,&Lam,&Lampr,k0,del,H,fac);
   kw = k0-dk;

   while(fabs(dk) > EPS) {
      k0       = kw;

      //Call dispersion relation function;
      //dk       = NR_corr_term(k0,del,H,fac);
      NR_corr_term(&dk,&Lam,&Lampr,k0,del,H,fac);
      kw       = k0-dk;
   }

   //Call dispersion relation function;
   //[dk,Lam,Lampr] = NR_corr_term(kw,del,H,fac);
   NR_corr_term(&dk,&Lam,&Lampr,kw,del,H,fac);
   denom          = H*(Lam*Lam*kw*kw-1)+Lampr;
   res            = -kw/denom;
   //printf("(denom,res,Bg2)=(%.9e,%.9e,%.9e)",res,denom,Lam*Lam*res);

   //Outputs;
   *kw2  = kw;
   *BG1  = Lam*Lam*res;
}

// function [dk,Lam,Lampr] = NR_corr_term(k,del,H,fac)
// %% dk=f/f_k, where f has the same zeros as of the dispersion function,
// %% is the correction term in the Newton-Rhapson method for finding zeros in f.
int NR_corr_term(double *dk,double *Lam2,double *Lampr2,
                 double k, double del, double H,double fac) {
//int NR_corr_term(double *,double *,double *,
//                 double, double, double,double);

   double f,df;
   double Lam,Lampr,x,k4;

   k4    = k*k*k*k;
   Lam   = fac*k4+del;
   Lampr = 5*fac*k4+del;
   x     = 7.5;

   if(fabs(k*H)<=x) {
      f  = Lam*k*sinh(k*H)-cosh(k*H);
      df = Lam*(k*H)*cosh(k*H)+(Lampr-H)*sinh(k*H);
   }
   else {
      f  = Lam*k*tanh(k*H)-1;
      df = Lam*k*H+(Lampr-H)*tanh(k*H);
   }

   //Outputs;
   *dk      = f/df;
   *Lam2    = Lam;
   *Lampr2  = Lampr;
}
