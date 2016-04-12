/* To Compile:
 * 'mex imag_root_ice.c' (on a machine with gsl installed) */
/*CALL: x=imag_root_ice(del,H,i0,i1)
	->finds the imaginary root of the dispersion relation for ice,
	  f=del*w*sin(w)+H^5*cos(w), inside the interval [i0,i1]*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "mex.h"

//my header files;
#include <RTparam.h>
#include <RTparam_hardcoded.h>

#define PI 3.14159265358979
#define EPS 1.0e-8
#define MAXIT 100




/*construct main function*/
int main() {

  //declare "input/output variables";
  double tx,ty,z;
  int    Ncx,Ncy,ncol;
  double *chebys;//use malloc/free later

  //auxiliary variables;
  int i,j,r,Ntot;

  //initialise other variables;
  tx     = -0.945346551;
  ty     = 0.708408642;
  Ncx    = 10;
  Ncy    = 10;
  Ntot   = (Ncx+1)*(Ncy+1);

  //get vector of chebyshev coefficients into;
  chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );
  ncol   = 0;//get column 0,1,2
  Amn_fxn_L4(chebys,ncol);
  //Amn_fxn(chebys);

  //do the interpolation;
  OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);

  if(1) {
     //test inputs;
     printf("\ninputs: (tx,ty)=(%f,%f)\n",tx,ty);

     //test chebys vector;
     printf("\nvector of chebyshev coeffients:\n");
     for(r=0;r<Ntot;r++) {
        printf("%f\n", chebys[r]);
     }

     //test outputs;
     if(ncol==0) {
        printf("\noutput: ac = %f\n", z);
     }
     else if(ncol==1) {
        printf("\noutput: Arg[R] = %f\n", z);
     }
     else if(ncol==2) {
        printf("\noutput: Arg[T] = %f\n", z);
     }
  }

  //free memory from chebys;
  free(chebys);

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
int Amn_fxn(double *chebys) {

  double Amn[4][3] = {
     {1.1,2.1,3.1},
     {4.1,5.1,6.1},
     {7.1,8.1,9.1},
     {10.1,11.1,12.1}
  };

  int M  = 4;
  int N  = 3;
  int i,j,r;

  for(i=0;i<M;i++) {
     for(j=0;j<N;j++) {
        r         = i+M*j;
        chebys[r] = Amn[i][j];
     }
  }
}
/***********************************************/
