#include "mex.h"
#include "math.h"
#include "string.h"

/*#define TEST*/
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _b : _a; })
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

void weno3pdV2(double* gin, double* u, double* scuy,
                       double* scp2i, double* scp2, double* saoout,
                       double dt, int nx, int nbdy)
{

    double cq00=-1./2 ,cq01=3./2, cq10=1./2, cq11=1./2, ca0=1./3, ca1=2./3, eps=1e-12;
    double q0, q1, a0, a1, q;
    int im1, im2, ip1;
    int nxext   =  nx+2*nbdy;
    double *ful = calloc(nxext,sizeof(double));
    double *fuh = calloc(nxext,sizeof(double));
    double *gt  = calloc(nxext,sizeof(double));


    /*fluxes in x direction*/
    for (int i = 2; i < nxext-1; i++)
    {
        im1 = i-1;

        if (u[i] > 0.)
        {
            /*coefficents to calc higher-order fluxes*/
            im2 = im1-1;
            q0 = cq00*gin[im2]+cq01*gin[im1];
            q1 = cq10*gin[im1]+cq11*gin[i];
            a0 = ca0;
            a1 = ca1*(fabs(gin[im2]-gin[im1])+eps)
               /(fabs(gin[im1]-gin[i])+eps);

            /*lower-order fluxes*/
            ful[i] = u[i]*gin[im1]*scuy[i];
        }
        else
        {
            /*coefficents to calc higher-order fluxes*/
            ip1 = i+1;
            q0 = cq11*gin[im1]+cq10*gin[i];
            q1 = cq01*gin[i]+cq00*gin[ip1];
            a0 = ca1;
            a1 = ca0*(fabs(gin[im1]-gin[i])+eps)/(fabs(gin[i]-gin[ip1])+eps);

            // lower-order fluxes
            ful[i] = u[i]*gin[i]*scuy[i];
        }

        // higher-order fluxes
        fuh[i] = (u[i]*(a0*q0+a1*q1)*scuy[i]/(a0+a1))-ful[i];
    }


    /* update field with low order fluxes*/
    for (int i = 0; i < nxext-1; i++)
        gt[i] = gin[i]-dt*(ful[(i+1)]-ful[i])*scp2i[i];

    q = 0.25/dt;

    /* obtain fluxes in x direction with limited high order correction fluxes*/
    for (int i = 1; i < nxext; i++)
        fuh[i] = ful[i]+
           max(-q*gt[i]*scp2[i],
                    min(q*gt[(i-1)]*scp2[(i-1)],fuh[i]));

    /* compute the spatial advective operator*/
    for (int i = 0; i < nxext-1; i++)
        saoout[i] = -(fuh[(i+1)]-fuh[i])*scp2i[i];

    /* deallocate memory for temp arrays*/
    free(ful);
    free(fuh);
    free(gt);
}


void padVar(double* u, double* upad, int advopt_, int nx, int nbdy)
{

   int nxext   = nx+2*nbdy;

    for (int i = 0; i < nxext; i++)
    {
        if ((nbdy-1 < i) && (i < nx+nbdy))
            upad[i] = u[(i-nbdy)];

        if (advopt_ == 1)/*x-periodic*/
        {
            /*NB change from matlab to C ordering*/
            /* make periodic in i
               - far-left cells */
            if ((i < nbdy) )
                upad[i] = u[(nx-nbdy+i)];

            /* - far-right cells */
            if ((nx+nbdy-1 < i))
                upad[i] = u[(i-nx-nbdy)];
        }/*"x-periodic" */
    }/*i*/
}


#if !defined(TEST)
void waveAdvWeno(double *h, double *u,double *LANDMASK,
      double *scp2, double *scp2i, double *scuy,
      double dt, int nx, int nbdy,int advopt)
#else
void waveAdvWeno(double *h, double *u,double *LANDMASK,
      double *scp2, double *scp2i, double *scuy,
      double dt, int nx, int nbdy,int advopt,double *test_array)
#endif
{
    int nxext = nx+2*nbdy;

    /*allocate memory of temp arrays*/
    /*-declare pointers*/
    double *sao,*u_pad,*h_pad,*scp2_pad,*scp2i_pad,*scuy_pad,*scvx_pad,*hp;
    sao         = malloc(nxext*sizeof(double));/*result of malloc is a pointer*/
    u_pad       = malloc(nxext*sizeof(double));
    h_pad       = malloc(nxext*sizeof(double));
    scp2_pad    = malloc(nxext*sizeof(double));
    scp2i_pad   = malloc(nxext*sizeof(double));
    scuy_pad    = malloc(nxext*sizeof(double));
    hp          = malloc(nxext*sizeof(double));


    /*boundary conditions*/
    int advopt_uv    = 1;/*x-periodic*/
    int advopt_grid  = 1;/*x-periodic*/
    padVar(u, u_pad, advopt_uv,nx,nbdy);
    padVar(scp2, scp2_pad,advopt_grid,nx,nbdy);
    padVar(scp2i, scp2i_pad,advopt_grid,nx,nbdy);
    padVar(scuy, scuy_pad,advopt_grid,nx,nbdy);
    padVar(h, h_pad,advopt,nx,nbdy);/*this field uses the variable boundary condition advopt*/
#if defined(TEST)
    /*get array to test, but take transpose*/
    for (int i = 0; i < nxext; i++)
        test_array[i]   = h_pad[i];
#endif

    /* prediction step*/
    weno3pdV2(h_pad, u_pad, scuy_pad, scp2i_pad, scp2_pad,
          sao,dt,nx,nbdy);

    if (nbdy>3)
    {
       /* need to loop over full padded domain*/
       for (int i = 0; i < nxext; i++)
           hp[i] = h_pad[i]+dt*sao[i];
    }
    else
    {
       printf("%s","\nAdvection (WENO): 'nbdy' should be >=4\n");
       abort();
    }

    /* correction step*/
    weno3pdV2(hp, u_pad, scuy_pad, scp2i_pad, scp2_pad,
          sao,dt,nx,nbdy);


    for (int i = 0; i < nx; i++)
    {
            double tmp = 0.5*(h_pad[(i+nbdy)]
                              +hp[(i+nbdy)]+dt*sao[(i+nbdy)]);

            /*mask land cells*/
            h[i]   = tmp*(1-LANDMASK[i]);
    }



    /*deallocate memory of temp arrays*/
    free(sao);
    free(u_pad);
    free(h_pad);
    free(scp2_pad);
    free(scp2i_pad);
    free(scuy_pad);
    free(hp);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/*
   main call is to 
      waveAdvWeno(double* h, double* u,int* LANDMASK,
         double* scp2, double* scp2i, double* scuy,
         double dt, int nx, int nbdy,int adv_opt)

    matlab call is:
    h_new = waveadv_weno_mex(nx,dt,advopt,h, u, LANDMASK, scp2, scp2i, scuy)
*/

    /*---------------------------------------------------------------------------------------------*/
    /*check no of arguments*/
#if !defined (TEST)
    if ( nlhs != 1 )
    {
        mexErrMsgIdAndTxt("waveadv_weno_mex:Outputs",
              "no of output arguments should be 1.\nexample matlab call:\nh_new = waveadv_weno_1d_mex(nx,dt,advopt,\n   h, u,LANDMASK, scp2, scp2i, scuy)\n");
    }
#else/*add a test output for debugging*/
    if ( nlhs != 2 )
    {
        mexErrMsgIdAndTxt("waveadv_weno_mex:Outputs",
              "no of output arguments should be 2.\nexample matlab call:\n[h_new,test] = waveadv_weno_1d_mex(nx,dt,advopt,\n   h, u, LANDMASK, scp2, scp2i, scuy)\n");
    }
#endif

    if ( nrhs != 9 )
    {
        mexErrMsgIdAndTxt("waveadv_weno_mex:Inputs",
              "no of input arguments should be 9.\nexample matlab call:\n[h_new,test] = waveadv_weno_1d_mex(nx,dt,advopt,\n   h, u, LANDMASK, scp2, scp2i, scuy)\n");
    }
    /*---------------------------------------------------------------------------------------------*/

    int nbdy = 4;  /* number of ghost cells to use*/

    /*---------- Input ----------*/
    /*scalar inputs*/
    mwSize nx        = (mwSize)   mxGetScalar(prhs[0]);   /* grid size in x dirn*/
    double dt        = (double)   mxGetScalar(prhs[1]);   /* timestep*/
    int    advopt    = (int)      mxGetScalar(prhs[2]);  /* periodicity option:
                                                             0=not periodic; 1=x-periodic*/

    /*check size of input arays*/
    for (int j=3;j<9;j++)
    {
        mwSize M   = mxGetM(prhs[j]);
        mwSize N   = mxGetN(prhs[j]);
        if ( (M!=nx) || (N!=1) )
        {
            mexErrMsgIdAndTxt("waveadv_weno_mex:Inputs",
              "incorrect size of input argument");
        }
    }


    /*pointers to input arrays*/
    /*NB all in matlab ordering*/
    double *h        = (double *) mxGetPr(prhs[3]);
    double *u        = (double *) mxGetPr(prhs[4]);
    double *LANDMASK = (double *) mxGetPr(prhs[5]);
    double *scp2     = (double *) mxGetPr(prhs[6]);
    double *scp2i    = (double *) mxGetPr(prhs[7]);
    double *scuy     = (double *) mxGetPr(prhs[8]);


    /*define output array*/
    plhs[0] = mxCreateDoubleMatrix(nx, 1, mxREAL);/*pointer to output array*/
    double *h_new   = (double *) mxGetPr(plhs[0]);
    for (mwIndex i=0;i<nx;i++)
        h_new[i]    = h[i];/*initialise to h*/
#if defined(TEST)
    plhs[1] = mxCreateDoubleMatrix(nx+2*nbdy, 1, mxREAL);/*pointer to output array*/
    double *test_array   = (double *) mxGetPr(plhs[1]);
#endif

    /*do advection
     * h_new is a pointer that is modified by waveAdvWeno
     * - this is also the output */
#if !defined(TEST)
    waveAdvWeno(h_new, u,LANDMASK, scp2, scp2i, scuy, 
        dt, (int) nx, nbdy,advopt);
#else
    waveAdvWeno(h_new, u,LANDMASK, scp2, scp2i, scuy,
        dt, (int) nx, nbdy,advopt,test_array);
#endif

    return;
}
