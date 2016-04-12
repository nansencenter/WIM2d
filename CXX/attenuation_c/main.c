#include <RTparam_outer.h>

#define PI 3.14159265358979
#define EPS 1.0e-12
#define MAXIT 100

int main(int argc, char **argv)
{
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
