/* To Compile:
 * 'mex imag_root_ice.c' (on a machine with gsl installed) */
/*CALL: x=imag_root_ice(del,H,i0,i1)
	->finds the imaginary root of the dispersion relation for ice,
	  f=del*w*sin(w)+H^5*cos(w), inside the interval [i0,i1]*/
#include <math.h>
//#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <RTparam.h>
#include <RTparam_hardcoded.h>

#define PI 3.14159265358979
#define EPS 1.0e-8
#define MAXIT 100




/*construct main function*/
int main() {
  double int_adm  = 1.33;
  int    OPT      = 5;
  int    LOW      = 0;
  double ac,modT,argR,argT;//outputs from RTparam_inner;
  double tx,ty;            //inputs to RTparam_inner;

  if(LOW==1) {
     ty  = 0.537243574;
     if(OPT==1) {
        tx  = 0.081382113;
     }
     else if(OPT==2) {
        tx  = 0.463351327;
     }
     else if(OPT==3) {
        tx  = 0.722706232;
     }
     else if(OPT==4) {
        tx  = 0.6;
     }
     else if(OPT==5) {
        tx  = 0.512044103;
     }
  }
  else if(LOW==0) {
     ty  = -0.700000000;
     if(OPT==1) {
        tx  = 0.081382113;
     }
     else if(OPT==2) {
        tx  = 0.463351327;
     }
     else if(OPT==3) {
        tx  = 0.722706232;
     }
     else if(OPT==4) {
        tx  = 0.6;
     }
     else if(OPT==5) {
        tx  = 0.457159276;
     }
  }

  RTparam_inner(&ac,&modT,&argR,&argT,
                tx,ty,OPT,LOW,int_adm);

  if(1) {
     //test inputs;
     printf("\ninputs: (tx,ty)=(%f,%f)\n",tx,ty);

     ////test chebys vector;
     //printf("\nvector of chebyshev coeffients:\n");
     //for(r=0;r<Ntot;r++) {
     //   printf("%f\n", chebys[r]);
     //}

     //test outputs;
     //printf("\noutput: z      = %f\n", z);
     printf("\noutput: ac     = %0.9e\n", ac);
     printf("\noutput: |T|    = %0.9e\n", modT);
     printf("\noutput: Arg[R] = %0.9e\n", argR);
     printf("\noutput: Arg[T] = %0.9e\n", argT);
  }
}

int RTparam_inner(double *ac2,double *modT2,double *argR2,double *argT2,
                  double tx,double ty,int OPT,int LOW,double int_adm) {

  /***********************************************************/
  //declare "input/output variables";
  double z;
  double ac,modR,modT,argR,argT;
  double Rr,Ri,Tr,Ti;

  //auxiliary variables;
  double *chebys;//use malloc/free later
  double *Zthings;//use malloc/free later
  int    Ncx,Ncy,ncol;
  int    i,j,r;
  int    Nthings;
  int    INTERP_MODE[5] = {1,1,3,2,1};
  int    IM;

  //Order of chebyshev expansion for lower thicknesses;
  int Ncx_L[5] = {10,10,10,10,3}; //for alp_nd;
  int Ncy_L[5] = {10,10,10,10,10};//for h_nd;

  //Order of chebyshev expansion for higher thicknesses;
  int Ncx_H[5] = {10,10,10,10,4};//for alp_nd;
  int Ncy_H[5] = {10,10,10,10,10};//for h_nd;
  /***********************************************************/

  //Mode of interpolation;
  if(INTERP_MODE[OPT-1]==1) {
     //Interpolate log(ac), Arg[R,T];
     Nthings   = 3;
  }
  else if(INTERP_MODE[OPT-1]==2) {
     //Interpolate (ac), Arg[R,T];
     Nthings   = 3;
  }
  else if(INTERP_MODE[OPT-1]==3) {
     //Interpolate Re[R,T], Im[R,T];
     Nthings   = 4;
  }

  if (LOW==1) {
     //Lower non-dimensional thicknesses;

     /*******************************************************/
     if(OPT==1) {
        Ncx    = Ncx_L[OPT-1];
        Ncy    = Ncy_L[OPT-1];

        /**************************************************/
        //Actual interpolation;

        //get vector of chebyshev coefficients into;
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        //get interpolated quantities;
        Zthings   = malloc( Nthings*sizeof(double) );

        printf("OPT=%d\n",OPT);
        printf("IM=%d\n",INTERP_MODE[OPT-1]);
        printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy);
        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_L1(chebys,ncol);
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        //free memory from chebys;
        free(chebys);
        /**************************************************/

        /**************************************************/
        //interpret interpolated quantities;
        if(INTERP_MODE[OPT-1]==1) {
           ac     = exp(Zthings[0]);
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==2) {
           ac     = Zthings[0];
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==3) {
           Rr     = Zthings[0];
           Ri     = Zthings[1];
           Tr     = Zthings[2];
           Ti     = Zthings[3];
           //
           argR   = atan2(Ri,Rr);
           argT   = atan2(Ti,Tr);
           modR   = sqrt(Rr*Rr+Ri*Ri);
           modT   = sqrt(Tr*Tr+Ti*Ti);
           ac     = -2*log(1-modR*modR);
        }
        /**************************************************/

        //free memory from Zthings;
        free(Zthings);

     }
     /*******************************************************/
     else if(OPT==2) {
        Ncx    = Ncx_L[OPT-1];
        Ncy    = Ncy_L[OPT-1];

        /**************************************************/
        //Actual interpolation;

        //get vector of chebyshev coefficients into;
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        //get interpolated quantities;
        Zthings   = malloc( Nthings*sizeof(double) );

        printf("OPT=%d\n",OPT);
        printf("IM=%d\n",INTERP_MODE[OPT-1]);
        printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy);
        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_L2(chebys,ncol);
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        //free memory from chebys;
        free(chebys);
        /**************************************************/

        /**************************************************/
        //interpret interpolated quantities;
        if(INTERP_MODE[OPT-1]==1) {
           ac     = exp(Zthings[0]);
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==2) {
           ac     = Zthings[0];
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==3) {
           Rr     = Zthings[0];
           Ri     = Zthings[1];
           Tr     = Zthings[2];
           Ti     = Zthings[3];
           //
           argR   = atan2(Ri,Rr);
           argT   = atan2(Ti,Tr);
           modR   = sqrt(Rr*Rr+Ri*Ri);
           modT   = sqrt(Tr*Tr+Ti*Ti);
           ac     = -2*log(1-modR*modR);
        }
        /**************************************************/

        //free memory from Zthings;
        free(Zthings);

     }
     /*******************************************************/
     else if(OPT==3) {
        Ncx    = Ncx_L[OPT-1];
        Ncy    = Ncy_L[OPT-1];

        /**************************************************/
        //Actual interpolation;

        //get vector of chebyshev coefficients into;
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        //get interpolated quantities;
        Zthings   = malloc( Nthings*sizeof(double) );

        printf("OPT=%d\n",OPT);
        printf("IM=%d\n",INTERP_MODE[OPT-1]);
        printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy);
        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_L3(chebys,ncol);
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        //free memory from chebys;
        free(chebys);
        /**************************************************/

        /**************************************************/
        //interpret interpolated quantities;
        if(INTERP_MODE[OPT-1]==1) {
           ac     = exp(Zthings[0]);
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==2) {
           ac     = Zthings[0];
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==3) {
           Rr     = Zthings[0];
           Ri     = Zthings[1];
           Tr     = Zthings[2];
           Ti     = Zthings[3];
           //
           argR   = atan2(Ri,Rr);
           argT   = atan2(Ti,Tr);
           modR   = sqrt(Rr*Rr+Ri*Ri);
           modT   = sqrt(Tr*Tr+Ti*Ti);
           ac     = -2*log(1-modR*modR);
        }
        /**************************************************/

        //free memory from Zthings;
        free(Zthings);

     }
     /*******************************************************/
     else if(OPT==4) {
        Ncx    = Ncx_L[OPT-1];
        Ncy    = Ncy_L[OPT-1];

        /**************************************************/
        //Actual interpolation;

        //get vector of chebyshev coefficients into;
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        //get interpolated quantities;
        Zthings   = malloc( Nthings*sizeof(double) );

        printf("OPT=%d\n",OPT);
        printf("IM=%d\n",INTERP_MODE[OPT-1]);
        printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy);
        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_L4(chebys,ncol);
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        //free memory from chebys;
        free(chebys);
        /**************************************************/

        /**************************************************/
        //interpret interpolated quantities;
        if(INTERP_MODE[OPT-1]==1) {
           ac     = exp(Zthings[0]);
           argR   = Zthings[1];
           argT   = Zthings[2];
        }
        else if(INTERP_MODE[OPT-1]==2) {
           ac     = Zthings[0];
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==3) {
           Rr     = Zthings[0];
           Ri     = Zthings[1];
           Tr     = Zthings[2];
           Ti     = Zthings[3];
           //
           argR   = atan2(Ri,Rr);
           argT   = atan2(Ti,Tr);
           modR   = sqrt(Rr*Rr+Ri*Ri);
           modT   = sqrt(Tr*Tr+Ti*Ti);
           ac     = -2*log(1-modR*modR);
        }
        /**************************************************/

        //free memory from Zthings;
        free(Zthings);
     }

     /*******************************************************/
     else if(OPT==5) {
        Ncx    = Ncx_L[OPT-1];
        Ncy    = Ncy_L[OPT-1];

        /**************************************************/
        //Actual interpolation;

        //get vector of chebyshev coefficients into;
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        //get interpolated quantities;
        Zthings   = malloc( Nthings*sizeof(double) );

        printf("OPT=%d\n",OPT);
        printf("IM=%d\n",INTERP_MODE[OPT-1]);
        printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy);
        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_H5(chebys,ncol);
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        //free memory from chebys;
        free(chebys);
        /**************************************************/

        /**************************************************/
        //interpret interpolated quantities;
        if(INTERP_MODE[OPT-1]==1) {
           ac     = exp(Zthings[0]);
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==2) {
           ac     = Zthings[0];
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==3) {
           Rr     = Zthings[0];
           Ri     = Zthings[1];
           Tr     = Zthings[2];
           Ti     = Zthings[3];
           //
           argR   = atan2(Ri,Rr);
           argT   = atan2(Ti,Tr);
           modR   = sqrt(Rr*Rr+Ri*Ri);
           modT   = sqrt(Tr*Tr+Ti*Ti);
           ac     = -2*log(1-modR*modR);
        }
        /**************************************************/

        //free memory from Zthings;
        free(Zthings);

     }
  } else {
     //Higher non-dimensional thicknesses;

     /*******************************************************/
     if(OPT==1) {
        Ncx    = Ncx_H[OPT-1];
        Ncy    = Ncy_H[OPT-1];

        /**************************************************/
        //Actual interpolation;

        //get vector of chebyshev coefficients into;
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        //get interpolated quantities;
        Zthings   = malloc( Nthings*sizeof(double) );

        printf("OPT=%d\n",OPT);
        printf("IM=%d\n",INTERP_MODE[OPT-1]);
        printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy);
        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_H1(chebys,ncol);
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        //free memory from chebys;
        free(chebys);
        /**************************************************/

        /**************************************************/
        //interpret interpolated quantities;
        if(INTERP_MODE[OPT-1]==1) {
           ac     = exp(Zthings[0]);
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==2) {
           ac     = Zthings[0];
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==3) {
           Rr     = Zthings[0];
           Ri     = Zthings[1];
           Tr     = Zthings[2];
           Ti     = Zthings[3];
           //
           argR   = atan2(Ri,Rr);
           argT   = atan2(Ti,Tr);
           modR   = sqrt(Rr*Rr+Ri*Ri);
           modT   = sqrt(Tr*Tr+Ti*Ti);
           ac     = -2*log(1-modR*modR);
        }
        /**************************************************/

        //free memory from Zthings;
        free(Zthings);

     }
     /*******************************************************/
     else if(OPT==2) {
        Ncx    = Ncx_H[OPT-1];
        Ncy    = Ncy_H[OPT-1];

        /**************************************************/
        //Actual interpolation;

        //get vector of chebyshev coefficients into;
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        //get interpolated quantities;
        Zthings   = malloc( Nthings*sizeof(double) );

        printf("OPT=%d\n",OPT);
        printf("IM=%d\n",INTERP_MODE[OPT-1]);
        printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy);
        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_H2(chebys,ncol);
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        //free memory from chebys;
        free(chebys);
        /**************************************************/

        /**************************************************/
        //interpret interpolated quantities;
        if(INTERP_MODE[OPT-1]==1) {
           ac     = exp(Zthings[0]);
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==2) {
           ac     = Zthings[0];
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==3) {
           Rr     = Zthings[0];
           Ri     = Zthings[1];
           Tr     = Zthings[2];
           Ti     = Zthings[3];
           //
           argR   = atan2(Ri,Rr);
           argT   = atan2(Ti,Tr);
           modR   = sqrt(Rr*Rr+Ri*Ri);
           modT   = sqrt(Tr*Tr+Ti*Ti);
           ac     = -2*log(1-modR*modR);
        }
        /**************************************************/

        //free memory from Zthings;
        free(Zthings);

     }
     /*******************************************************/
     else if(OPT==3) {
        Ncx    = Ncx_H[OPT-1];
        Ncy    = Ncy_H[OPT-1];

        /**************************************************/
        //Actual interpolation;

        //get vector of chebyshev coefficients into;
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        //get interpolated quantities;
        Zthings   = malloc( Nthings*sizeof(double) );

        printf("OPT=%d\n",OPT);
        printf("IM=%d\n",INTERP_MODE[OPT-1]);
        printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy);
        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_H3(chebys,ncol);
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        //free memory from chebys;
        free(chebys);
        /**************************************************/

        /**************************************************/
        //interpret interpolated quantities;
        if(INTERP_MODE[OPT-1]==1) {
           ac     = exp(Zthings[0]);
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==2) {
           ac     = Zthings[0];
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==3) {
           Rr     = Zthings[0];
           Ri     = Zthings[1];
           Tr     = Zthings[2];
           Ti     = Zthings[3];
           //
           argR   = atan2(Ri,Rr);
           argT   = atan2(Ti,Tr);
           modR   = sqrt(Rr*Rr+Ri*Ri);
           modT   = sqrt(Tr*Tr+Ti*Ti);
           ac     = -2*log(1-modR*modR);
        }
        /**************************************************/

        //free memory from Zthings;
        free(Zthings);

     }
     /*******************************************************/
     else if(OPT==4) {
        Ncx    = Ncx_H[OPT-1];
        Ncy    = Ncy_H[OPT-1];

        /**************************************************/
        //Actual interpolation;

        //get vector of chebyshev coefficients into;
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        //get interpolated quantities;
        Zthings   = malloc( Nthings*sizeof(double) );

        printf("OPT=%d\n",OPT);
        printf("IM=%d\n",INTERP_MODE[OPT-1]);
        printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy);
        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_H4(chebys,ncol);
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        //free memory from chebys;
        free(chebys);
        /**************************************************/

        /**************************************************/
        //interpret interpolated quantities;
        if(INTERP_MODE[OPT-1]==1) {
           ac     = exp(Zthings[0]);
           argR   = Zthings[1];
           argT   = Zthings[2];
        }
        else if(INTERP_MODE[OPT-1]==2) {
           ac     = Zthings[0];
           argR   = Zthings[1];
           argT   = Zthings[2];
           modT   = sqrt(exp(-ac/2.0)/int_adm);
        }
        else if(INTERP_MODE[OPT-1]==3) {
           Rr     = Zthings[0];
           Ri     = Zthings[1];
           Tr     = Zthings[2];
           Ti     = Zthings[3];
           //
           argR   = atan2(Ri,Rr);
           argT   = atan2(Ti,Tr);
           modR   = sqrt(Rr*Rr+Ri*Ri);
           modT   = sqrt(Tr*Tr+Ti*Ti);
           ac     = -2*log(1-modR*modR);
        }
        /**************************************************/

        //free memory from Zthings;
        free(Zthings);
     }

     /*******************************************************/
     else if(OPT==5) {
        Ncx = Ncx_H[OPT-1];
        Ncy = Ncy_H[OPT-1];
        IM  = INTERP_MODE[OPT-1];

        /**************************************************/
        //Actual interpolation;

        //get vector of chebyshev coefficients into;
        chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );

        //get interpolated quantities;
        Zthings   = malloc( Nthings*sizeof(double) );

        printf("OPT=%d\n",OPT);
        printf("IM=%d\n",IM);
        printf("(Ncx,Ncy)=(%d,%d)\n",Ncx,Ncy);
        for(ncol=0;ncol<Nthings;ncol++) {
           Amn_fxn_H5(chebys,ncol);
           OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
           Zthings[ncol]   = z;
        }

        //free memory from chebys;
        free(chebys);
        /**************************************************/

        /**************************************************/
        //interpret interpolated quantities;
        RTparam_get_ac(ac2,modT2,argR2,argT2,Zthings,IM,int_adm);
        //if(INTERP_MODE[OPT-1]==1) {
        //   ac     = exp(Zthings[0]);
        //   argR   = Zthings[1];
        //   argT   = Zthings[2];
        //   modT   = sqrt(exp(-ac/2.0)/int_adm);
        //}
        //else if(INTERP_MODE[OPT-1]==2) {
        //   ac     = Zthings[0];
        //   argR   = Zthings[1];
        //   argT   = Zthings[2];
        //   modT   = sqrt(exp(-ac/2.0)/int_adm);
        //}
        //else if(INTERP_MODE[OPT-1]==3) {
        //   Rr     = Zthings[0];
        //   Ri     = Zthings[1];
        //   Tr     = Zthings[2];
        //   Ti     = Zthings[3];
        //   //
        //   argR   = atan2(Ri,Rr);
        //   argT   = atan2(Ti,Tr);
        //   modR   = sqrt(Rr*Rr+Ri*Ri);
        //   modT   = sqrt(Tr*Tr+Ti*Ti);
        //   ac     = -2*log(1-modR*modR);
        //}
        /**************************************************/

        //free memory from Zthings;
        free(Zthings);

     }
  }

  //change values of outputs;
  //*ac2   = ac;
  //*modT2 = modT;
  //*argR2 = argR;
  //*argT2 = argT;

  if(1) {
     //test inputs;
     printf("\ninputs: (tx,ty)=(%f,%f)\n",tx,ty);

     ////test chebys vector;
     //printf("\nvector of chebyshev coeffients:\n");
     //for(r=0;r<Ntot;r++) {
     //   printf("%f\n", chebys[r]);
     //}

     //test outputs;
     //printf("\noutput: z      = %f\n", z);
     printf("\noutput: ac     = %f\n", ac);
     printf("\noutput: |T|    = %f\n", modT);
     printf("\noutput: Arg[R] = %f\n", argR);
     printf("\noutput: Arg[T] = %f\n", argT);
  }


}

/*********************************************************
*********************************************************/
//int OP_chebinterp2d(double *z,double tx,double ty,
//                     double chebys[],int Ncx,int Ncy) {
int OP_chebinterp2d(double *z,double tx,double ty,
                     double *chebys,int Ncx,int Ncy) {

  double Amn, am,z0;
  double Tm0_x, Tm1_x, Tm_x;
  double Tn0_y, Tn1_y, Tn_y;
  int    s,nx,ny;

  //printf("\n(in OP_chebinterp2d)(tx,ty)=(%f,%f)\n",tx,ty);
  z0     = 0.0;
  Tm0_x  = 1.0;
  Tm1_x  = tx;
  //printf("\ntx,Tm1_x,Tm0_x=%f,%f,%f\n",tx,Tm1_x,Tm0_x);

  /**************************************************
  * Need to calculate z=\sum_{m=0}^Ncx a_m*T_m(tx),
  *  where a_m=\sum_{n=0}^Ncy A_{mn}*T_n(ty);
  **************************************************/

  //m=0 term;
  nx    = 0;
  Tm_x  = Tm0_x;
  am    = 0.0;
  //
  ny    = 0;
  Tn0_y = 1.0;
  s     = nx + ny*(Ncx+1);//NB s starts from 0;
  Amn   = chebys[s];
  am    = am+Amn*Tn0_y;
  //printf("\n{nx,ny,s}=%d,%d,%d\n", nx,ny,s);
  //
  //printf("\ns=%d, Amn=%f,Tn=%f,am=%f\n", s,Amn,Tn0_y,am);
  //
  ny    = 1;
  Tn1_y = ty;
  s     = nx + ny*(Ncx+1);//NB s starts from 0;
  Amn   = chebys[s];
  am    = am+Amn*Tn1_y;
  //printf("\n{nx,ny,s}=%d,%d,%d\n", nx,ny,s);
  //
  //printf("\ns=%d, Amn=%f,Tn=%f,am=%f\n", s,Amn,Tn1_y,am);
  //
  for(ny=2;ny<=Ncy;ny++) {
     Tn_y  = 2*ty*Tn1_y-Tn0_y;
     //printf("\nny,ty,Tn_y,Tn1_y,Tn0_y=%d,%f,%f,%f,%f\n",ny,ty,Tn_y,Tn1_y,Tn0_y);
     //
     Tn0_y = Tn1_y;
     Tn1_y = Tn_y;
     //
     s     = nx + ny*(Ncx+1);//NB s starts from 0;
     //printf("\n{nx,ny,s}=%d,%d,%d\n", nx,ny,s);
     //
     Amn   = chebys[s];
     am    = am+Amn*Tn_y;
     //printf("\ns=%d, Amn=%f,Tn=%f,am=%f\n", s,Amn,Tn_y,am);
  }

  //update z;
  z0 = z0+am*Tm_x;
  //printf("z (m=0): %f\n",z0);

  /*******************************************************
  * finished m=0 term;
  * now do m=1 term;
  *******************************************************/

  nx    = 1;
  Tm_x  = Tm1_x;
  am    = 0.0;
  //
  ny    = 0;
  Tn0_y = 1.0;
  s     = nx + ny*(Ncx+1);
  Amn   = chebys[s];
  am    = am+Amn*Tn0_y;
  //
  ny    = 1;
  Tn1_y = ty;
  s     = nx + ny*(Ncx+1);
  Amn   = chebys[s];
  am    = am+Amn*Tn1_y;

  for(ny=2;ny<=Ncy;ny++) {
     Tn_y  = 2*ty*Tn1_y-Tn0_y;
     Tn0_y = Tn1_y;
     Tn1_y = Tn_y;
     //
     s     = nx + ny*(Ncx+1);//NB s starts from 0;
     Amn   = chebys[s];
     am    = am+Amn*Tn1_y;
  }

  //update z;
  z0 = z0+am*Tm_x;
  //printf("z (m=1): %f\n",z0);

  /*******************************************************
  * Finished m=1 term;
  * Now sum over rest of terms (m=2 to m=Ncx);
  *******************************************************/
  for(nx=2;nx<=Ncx;nx++) {
     Tm_x  = 2*tx*Tm1_x-Tm0_x;
     Tm0_x = Tm1_x;
     Tm1_x = Tm_x;
     am    = 0.0;
     //
     ny    = 0;
     Tn0_y = 1.0;
     s     = nx + ny*(Ncx+1);
     Amn   = chebys[s];
     am    = am+Amn*Tn0_y;
     //
     ny    = 1;
     Tn1_y = ty;
     s     = nx + ny*(Ncx+1);
     Amn   = chebys[s];
     am    = am+Amn*Tn1_y;
     //
     for(ny=2;ny<=Ncy;ny++) {
        Tn_y  = 2*ty*Tn1_y-Tn0_y;
        Tn0_y = Tn1_y;
        Tn1_y = Tn_y;
        //
        s     = nx + ny*(Ncx+1);//NB s starts from 0;
        Amn   = chebys[s];
        am    = am+Amn*Tn1_y;
     }

     //update z;
     z0 = z0+am*Tm_x;
     //printf("z (m=%d): %f\n",nx,z0);
   }

   *z = z0;
   //printf("\n(in OP_chebinterp2d) output (z) = %f\n",*z);

   return;
}

/***********************************************/
int RTparam_get_ac(double *ac2,double *modT2,double *argR2,double *argT2,
      double *Zthings,int IM,double int_adm) {
//int RTparam_get_ac(double *,double *,double *,double *,
//    double *,int,double);

  double ac,modR,modT,argR,argT;
  double Rr,Ri,Tr,Ti;

  if(IM==1) {
     ac     = exp(Zthings[0]);
     argR   = Zthings[1];
     argT   = Zthings[2];
     modT   = sqrt(exp(-ac/2.0)/int_adm);
  }
  else if(IM==2) {
     ac     = Zthings[0];
     argR   = Zthings[1];
     argT   = Zthings[2];
     modT   = sqrt(exp(-ac/2.0)/int_adm);
  }
  else if(IM==3) {
     Rr     = Zthings[0];
     Ri     = Zthings[1];
     Tr     = Zthings[2];
     Ti     = Zthings[3];
     //
     argR   = atan2(Ri,Rr);
     argT   = atan2(Ti,Tr);
     modR   = sqrt(Rr*Rr+Ri*Ri);
     modT   = sqrt(Tr*Tr+Ti*Ti);
     ac     = -2*log(1-modR*modR);
  }

  *ac2   = ac;
  *modT2 = modT;
  *argR2 = argR;
  *argT2 = argT;

}
/***********************************************/
