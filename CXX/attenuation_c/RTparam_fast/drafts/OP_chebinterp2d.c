/* To Compile:
 * 'mex imag_root_ice.c' (on a machine with gsl installed) */
/*CALL: x=imag_root_ice(del,H,i0,i1)
	->finds the imaginary root of the dispersion relation for ice,
	  f=del*w*sin(w)+H^5*cos(w), inside the interval [i0,i1]*/
#include <math.h>
//#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979
#define EPS 1.0e-8
#define MAXIT 100


/*************************************************/
/*function declarations TODO: put in header file?*/
//int OP_chebinterp2d(double *,double,double,
//                     double[],int,int);
int OP_chebinterp2d(double *,double,double,
                     double *,int,int);
/*************************************************/

/*construct main function*/
int main() {

  //declare "input/output variables";
  double tx,ty,z;
  int    Ncx,Ncy;
  //double chebys[(Ncx+1)*(Ncy+1)];
  //double chebys[12];
  double *chebys;//use malloc/free later

  //auxiliary variables;
  int i,j,r,Ntot;

  //const double Amn[Ncx+1][Ncy+1] = {
  double Amn[4][3] = {
     {1.1,2.1,3.1},
     {4.1,5.1,6.1},
     {7.1,8.1,9.1},
     {10.1,11.1,12.1}
  };

  //initialise other variables;
  tx     = 0.2;
  ty     = -0.12;
  printf("\nInputs: (tx,ty)=(%f,%f)\n",tx,ty);
  Ncx    = 3;
  Ncy    = 2;
  Ntot   = (Ncx+1)*(Ncy+1);

  //Arrange vector of chebyshev coefficients into;
  chebys = malloc( (Ncx+1)*(Ncy+1)*sizeof(double) );
  for(i=0;i<(Ncx+1);i++) {
     for(j=0;j<(Ncy+1);j++) {

        r         = i+(Ncx+1)*j;
        chebys[r] = Amn[i][j];
        //printf("(%d,%d)(rows=%d)->r=%d\n", i,j,Ncx+1,r);
        //printf("%f can be written %e\n", chebys[r], chebys[r]);
        //printf("(%d,%d)->%f\n", i,j,Amn[i][j]);
        //printf("%f\n", chebys[r]);

     }
  }

  //test chebys vector;
  printf("\nVector of chebyshev coeffients:\n");
  for(r=0;r<Ntot;r++) {
     printf("%f\n", chebys[r]);
  }

  //CALL OP_chebinterp2d;
  OP_chebinterp2d(&z,tx,ty,chebys,Ncx,Ncy);
  printf("\nOutput: z = %f\n", z);
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

