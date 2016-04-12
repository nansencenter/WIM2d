//Standard files;
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//My header files;
#include <RTparam_fast.h>
#include <RTparam_hardcoded.h>

#define PI 3.14159265358979
#define EPS 1.0e-8
#define MAXIT 100

/*construct main function*/
int main() {
  double int_adm  = 1.33;
  double rho      = 0.9;
  //int    OPT      = 2;
  //int    LOW      = 1;
  double alp_nd,hnd;//inputs to RTparam;
  double ac,modT,argR,argT;//outputs from RTparam;
  int OPT,LOW;

  for(LOW=0;LOW<=0;LOW++) {
  for(OPT=1;OPT<=5;OPT++) {

  if(LOW==1) {
     hnd  = 0.1;
     if(OPT==1) {
        alp_nd  = 1.0e-4;
     }
     else if(OPT==2) {
        alp_nd  = 1.0e-1;
     }
     else if(OPT==3) {
        alp_nd   = 1.2;
     }
     else if(OPT==4) {
        alp_nd  = 8.918121856041489;
     }
     else if(OPT==5) {
        alp_nd  = 72.701156608510175;
     }
  }
  else if(LOW==0) {
     hnd  = .23;
     if(OPT==1) {
        alp_nd  = 1.0e-4;
     }
     else if(OPT==2) {
        alp_nd  = 1.0e-1;
     }
     else if(OPT==3) {
        alp_nd  = 1.2;
     }
     else if(OPT==4) {
        alp_nd  = 6.703563588309956;
     }
     else if(OPT==5) {
        alp_nd  = 38.239307354260035;
     }
  }
  //call RTparam;
  RTparam_fast(&ac,&modT,&argR,&argT,alp_nd,hnd,int_adm);

  //check results;
  printf("output: ac     = %0.9e\n", ac);
  printf("output: |T|    = %0.9e\n", modT);
  printf("output: Arg[R] = %0.9e\n", argR);
  printf("output: Arg[T] = %0.9e\n", argT);
  printf("************************************\n");

  }
  }
}
