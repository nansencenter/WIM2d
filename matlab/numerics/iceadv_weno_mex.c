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

void weno3pd(double *gin, double *u, double *v, double *scuy,
             double *scvx,double *scp2i,double *scp2,
             int *pmask,int *umask,int *vmask,int *isp,int *isu,int *isv,
             int *ifp,int *ifu,int *ifv,int *ilp,int *ilu,int *ilv,
             double *saoout, double dt, int nx, int ny, int nbdy)
{

    double cq00=-1./2 ,cq01=3./2, cq10=1./2, cq11=1./2, ca0=1./3, ca1=2./3, eps=1e-12;
    double q0, q1, a0, a1, q;
    int im1, im2, ip1, jm1, jm2, jp1;
    int nxext,nyext;
    
    nxext = nx+2*nbdy;
    nyext = ny+2*nbdy;
    double *ful   = calloc(nxext*nyext,sizeof(double));
    double *fuh   = calloc(nxext*nyext,sizeof(double));
    double *fvl   = calloc(nxext*nyext,sizeof(double));
    double *fvh   = calloc(nxext*nyext,sizeof(double));
    double *gt    = calloc(nxext*nyext,sizeof(double));


    /*fluxes in x direction*/
    for (int j = 0; j < nyext; j++)
        for (int l = 0; l < isu[j]; l++)
            for (int i = max(2,ifu[j+l*nyext]+nbdy-1); i < max(nxext-1,ilu[j+l*nyext]+nbdy-1); i++)
            {
                im1 = i-1;

                if (u[i*nyext+j] > 0.)
                {
                    /*coefficents to calc higher-order fluxes*/
                    im2 = im1-1;
                    q0 = cq00*gin[im2*nyext+j]+cq01*gin[im1*nyext+j];
                    q1 = cq10*gin[im1*nyext+j]+cq11*gin[i*nyext+j];
                    a0 = ca0;
                    a1 = ca1*(fabs(gin[im2*nyext+j]-gin[im1*nyext+j])+eps)
                       /(fabs(gin[im1*nyext+j]-gin[i*nyext+j])+eps);

                    /*lower-order fluxes*/
                    ful[i*nyext+j] = u[i*nyext+j]*gin[im1*nyext+j]*scuy[i*nyext+j];
                }
                else
                {
                    /*coefficents to calc higher-order fluxes*/
                    ip1 = i+1;
                    q0 = cq11*gin[im1*nyext+j]+cq10*gin[i*nyext+j];
                    q1 = cq01*gin[i*nyext+j]+cq00*gin[ip1*nyext+j];
                    a0 = ca1;
                    a1 = ca0*(fabs(gin[im1*nyext+j]-gin[i*nyext+j])+eps)/(fabs(gin[i*nyext+j]-gin[ip1*nyext+j])+eps);

                    // lower-order fluxes
                    ful[i*nyext+j] = u[i*nyext+j]*gin[i*nyext+j]*scuy[i*nyext+j];
                }

                // higher-order fluxes
                fuh[i*nyext+j] = (u[i*nyext+j]*(a0*q0+a1*q1)*scuy[i*nyext+j]/(a0+a1))-ful[i*nyext+j];
            }

    /*fluxes in y direction*/
    for (int j = 2; j < nyext-1; j++)
        for (int l = 0; l < isv[j]; l++)
            for (int i = max(0,ifv[j+l*nyext]+nbdy-1); i < max(nxext,ilv[j+l*nyext]+nbdy-1); i++)
            {
                jm1 = j-1;

                if (v[i*nyext+j] > 0.)
                {
                    jm2 = jm1-1;
                    q0 = cq00*gin[i*nyext+jm2]+cq01*gin[i*nyext+jm1];
                    q1 = cq10*gin[i*nyext+jm1]+cq11*gin[i*nyext+j];
                    a0 = ca0;
                    a1 = ca1*(fabs(gin[i*nyext+jm2]-gin[i*nyext+jm1])+eps)
                       /(fabs(gin[i*nyext+jm1]-gin[i*nyext+j])+eps);
                    fvl[i*nyext+j] = v[i*nyext+j]*gin[i*nyext+jm1]*scvx[i*nyext+j];
                }
                else
                {
                    jp1 = j+1;
                    q0 = cq11*gin[i*nyext+jm1]+cq10*gin[i*nyext+j];
                    q1 = cq01*gin[i*nyext+j]+cq00*gin[i*nyext+jp1];
                    a0 = ca1;
                    a1 = ca0*(fabs(gin[i*nyext+jm1]-gin[i*nyext+j])+eps)
                       /(fabs(gin[i*nyext+j]-gin[i*nyext+jp1])+eps);
                    fvl[i*nyext+j] = v[i*nyext+j]*gin[i*nyext+j]*scvx[i*nyext+j];
                }

                fvh[i*nyext+j] = (v[i*nyext+j]*(a0*q0+a1*q1)*scvx[i*nyext+j]/(a0+a1))-fvl[i*nyext+j];
            }


    /* update field with low order fluxes*/
    for (int j = 0; j < nyext; j++)
        for (int l = 0; l < isp[j]; l++)
            for (int i = max(0,ifp[j+l*nyext]+nbdy-1); i < max(nxext-1,ilp[j+l*nyext]+nbdy-1); i++)
                gt[i*nyext+j] = gin[i*nyext+j]-dt*(ful[(i+1)*nyext+j]
                      -ful[i*nyext+j]+fvl[i*nyext+j+1]-fvl[i*nyext+j])*scp2i[i*nyext+j];

    q = 0.25/dt;

    /* obtain fluxes in x direction with limited high order correction fluxes*/
    for (int j = 0; j < nyext; j++)
        for (int l = 0; l < isu[j]; l++)
            for (int i = max(1,ifu[j+l*nyext]+nbdy-1); i < max(nxext,ilu[j+l*nyext]+nbdy-1); i++)
                fuh[i*nyext+j] = ful[i*nyext+j]+
                   max(-q*gt[i*nyext+j]*scp2[i*nyext+j],
                            min(q*gt[(i-1)*nyext+j]*scp2[(i-1)*nyext+j],fuh[i*nyext+j]));

    /* obtain fluxes in y direction with limited high order correction fluxes*/
    for (int j = 1; j < nyext; j++)
        for (int l = 0; l < isv[j]; l++)
            for (int i = max(0,ifv[j+l*nyext]+nbdy-1); i < max(nxext,ilv[j+l*nyext]+nbdy-1); i++)
                fvh[i*nyext+j]=fvl[i*nyext+j]+
                   max(-q*gt[i*nyext+j]*scp2[i*nyext+j],
                            min(q*gt[i*nyext+j-1]*scp2[i*nyext+j-1],fvh[i*nyext+j]));


    /* compute the spatial advective operator*/
    for (int j = 0; j < nyext; j++)
        for (int l = 0; l < isp[j]; l++)
            for (int i = max(0,ifp[j+l*nyext]+nbdy-1); i < max(nxext-1,ilp[j+l*nyext]+nbdy-1); i++)
                saoout[i*nyext+j] = -(fuh[(i+1)*nyext+j]-fuh[i*nyext+j]+fvh[i*nyext+j+1]
                    -fvh[i*nyext+j])*scp2i[i*nyext+j];

    /* deallocate memory for temp arrays*/
    free(ful);
    free(fuh);
    free(fvl);
    free(fvh);
    free(gt);
}


void padVar(double* u, double* upad, int advopt_, int nx, int ny, int nbdy)
{

   int nxext   = nx+2*nbdy;
   int nyext   = ny+2*nbdy;

    for (int i = 0; i < nxext; i++)
    {
        for (int j = 0; j < nyext; j++)
        {

            if ((nbdy-1 < i) && (i < nx+nbdy) && (nbdy-1 < j) && (j < ny+nbdy))
            {
                /*NB change from matlab to C ordering*/
                upad[i*nyext+j] = u[(i-nbdy)+nx*(j-nbdy)];
            }

            if (advopt_ != 0)/*either y-periodic or xy-periodic*/
            {
                /*NB change from matlab to C ordering*/
                /* make periodic in j
                   - lower cells*/
                bool i_inner = ((nbdy-1 < i) && (i < nx+nbdy));
                if ((j < nbdy) && i_inner)
                    upad[i*nyext+j] = u[(i-nbdy)+nx*(ny-nbdy+j)];

                /* - upper cells */
                if ((ny+nbdy-1 < j) && i_inner)
                    upad[i*nyext+j] = u[(i-nbdy)+nx*(j-ny-nbdy)];
            }

            if (advopt_ == 1)/*xy-periodic*/
            {
                /*NB change from matlab to C ordering*/
                /* make periodic in i
                   - far-left cells */
                bool j_inner = ((nbdy-1 < j) && (j < ny+nbdy));
                if ((i < nbdy) && j_inner )
                    upad[i*nyext+j] = u[(nx-nbdy+i)+nx*(j-nbdy)];

                /* - far-right cells */
                if ((nx+nbdy-1 < i) && j_inner )
                    upad[i*nyext+j] = u[(i-nx-nbdy)+nx*(j-nbdy)];

                /* TR */
                if ((nx+nbdy-1 < i) && (ny+nbdy-1 < j))
                    upad[i*nyext+j] = u[(i-nx-nbdy)+nx*(j-ny-nbdy)];

                /* BL */
                if ((i < nbdy) && (j < nbdy))
                    upad[i*nyext+j] = u[(i+nx-nbdy)+nx*j];

                /* BR */
                if ((nx+nbdy-1 < i) && (j < nbdy))
                    upad[i*nyext+j] = u[(i-nx-nbdy)+nx*(ny-nbdy+j)];

                /* TL */
                if ((i < nbdy) && (ny+nbdy-1 < j))
                    upad[i*nyext+j] = u[(i+nx-nbdy)+nx*(j-ny-nbdy)];
            }/*"xy-periodic" */
        }/*j*/
    }/*i*/
}


#if !defined(TEST)
void iceAdvWeno(double *h, double *u, double *v,
      double *scp2, double *scp2i, double *scuy, double *scvx,
      int *pmask,int *umask,int *vmask,int *isp,int *isu,int *isv,
      int *ifp,int *ifu,int *ifv,int *ilp,int *ilu,int *ilv,
      double dt, int nx, int ny, int nbdy)
#else
void iceAdvWeno(double *h, double *u, double *v,
      double *scp2, double *scp2i, double *scuy, double *scvx,
      int *pmask,int *umask,int *vmask,int *isp,int *isu,int *isv,
      int *ifp,int *ifu,int *ifv,int *ilp,int *ilu,int *ilv,
      double dt, int nx, int ny, int nbdy, double *test_array)
#endif
{
    int nxext = nx+2*nbdy;
    int nyext = ny+2*nbdy;

    /*allocate memory of temp arrays (also initialise to 0 everywhere)*/
    double *sao       = calloc(nxext*nyext,sizeof(double));/*result of calloc (sao) is a pointer*/
    double *u_pad     = calloc(nxext*nyext,sizeof(double));
    double *v_pad     = calloc(nxext*nyext,sizeof(double));
    double *h_pad     = calloc(nxext*nyext,sizeof(double));
    double *scp2_pad  = calloc(nxext*nyext,sizeof(double));
    double *scp2i_pad = calloc(nxext*nyext,sizeof(double));
    double *scuy_pad  = calloc(nxext*nyext,sizeof(double));
    double *scvx_pad  = calloc(nxext*nyext,sizeof(double));
    double *hp        = calloc(nxext*nyext,sizeof(double));


    /*boundary conditions*/
    /*NB padVar changes from matlab to C ordering
     * (takes the transpose)*/
    int advopt       = 0;/*not-periodic*/
    int advopt_uv    = 0;/*not-periodic*/
    int advopt_grid  = 1;/*xy-periodic*/
    padVar(u, u_pad, advopt_uv,nx,ny,nbdy);
    padVar(v, v_pad,advopt_uv,nx,ny,nbdy);
    padVar(scp2, scp2_pad,advopt_grid,nx,ny,nbdy);
    padVar(scp2i, scp2i_pad,advopt_grid,nx,ny,nbdy);
    padVar(scuy, scuy_pad,advopt_grid,nx,ny,nbdy);
    padVar(scvx, scvx_pad,advopt_grid,nx,ny,nbdy);
    padVar(h, h_pad,advopt,nx,ny,nbdy);/*this field uses the variable boundary condition advopt*/
#if defined(TEST)
    /*get array to test, but take transpose*/
    for (int i = 0; i < nxext; i++)
        for (int j = 0; j < nyext; j++)
            test_array[i+nxext*j]   = h_pad[i*nyext+j];
#endif

    /* prediction step*/
    weno3pd(h_pad, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad,
            pmask,umask,vmask,isp,isu,isv,
            ifp,ifu,ifv,ilp,ilu,ilv,
            sao,dt,nx,ny,nbdy);

    if (nbdy>3)
    {
        /* need to loop over full padded domain*/
        /*NB i,j are for h_pad, which is padded*/
        /*ifp is relative: ifp=1-nbdy->i=0 or 1st padded row,
         * so compare i,ifp+nbdy-1 */
        /* NB ifp,ilp have matlab ordering*/
        /* NB hp,h_pad have C ordering*/
        for (int j = 0; j < nyext; j++)
            for (int l = 0; l < isp[j]; l++)
                for (int i = max(0,ifp[j+l*nyext]+nbdy-1); i < max(nxext,ilp[j+l*nyext]+nbdy-1); i++)
                    hp[i*nyext+j] = h_pad[i*nyext+j]+dt*sao[i*nyext+j];
    }
    else
    {
       printf("%s","\nAdvection (WENO): 'nbdy' should be >=3\n");
       abort();
    }

    /* correction step*/
    weno3pd(hp, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad,
            pmask,umask,vmask,isp,isu,isv,
            ifp,ifu,ifv,ilp,ilu,ilv,
            sao,dt,nx,ny,nbdy);


    /*NB i,j are for h, which is not padded*/
    /*ifp is relative: ifp=1-nbdy->i=-nbdy or 1st padded row,
     * so i=0 corresponds to ifp=1,
     * so compare i,ifp-1 */
    for (int j = 0; j < ny; j++)
        for (int l = 0; l < isp[j]; l++)
            for (int i = max(0,ifp[j+l*nyext]-1); i < max(nx,ilp[j+l*nyext]-1); i++)
                h[i+nx*j]   = 0.5*(h_pad[(i+nbdy)*nyext+j+nbdy]
                                  +hp[(i+nbdy)*nyext+j+nbdy]+dt*sao[(i+nbdy)*nyext+j+nbdy]);


    /*deallocate memory of temp arrays*/
    free(sao);
    free(u_pad);
    free(v_pad);
    free(h_pad);
    free(scp2_pad);
    free(scp2i_pad);
    free(scuy_pad);
    free(scvx_pad);
    free(hp);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/*
    main call is to 
    iceAdvWeno(double *h, double *u, double *v,
          double *scp2, double *scp2i, double *scuy, double *scvx,
          int *pmask,int *umask,int *vmask,int *isp,int *isu,int *isv,
          int *ifp,int *ifu,int *ifv,int *ilp,int *ilu,int *ilv,
          double dt, int nx, int ny, int nbdy)

    matlab call is:
    h_new = iceadv_weno_mex(nx,ny,dt,advopt,nbdy,h, u, v, scp2, scp2i, scuy, scvx,...
                            pmask,umask,vmask,isp,isu,isv,ifp,ifu,ifv,ilp,ilu,ilv)
*/

    /*---------------------------------------------------------------------------------------------*/
    /*check no of arguments*/
#if !defined (TEST)
    if ( nlhs != 1 )
    {
        mexErrMsgIdAndTxt("iceadv_weno_mex:Outputs",
              "no of output arguments should be 1.\nexample matlab call:\nh_new = iceadv_weno_mex(nx,ny,dt,nbdy,\n   h, u, v,LANDMASK, scp2, scp2i, scuy, scvx)\n");
    }
#else/*add a test output for debugging*/
    if ( nlhs != 2 )
    {
        mexErrMsgIdAndTxt("iceadv_weno_mex:Outputs",
              "no of output arguments should be 2.\nexample matlab call:\n[h_new,test] = iceadv_weno_mex(nx,ny,dt,nbdy,\n   h, u, v,LANDMASK, scp2, scp2i, scuy, scvx)\n");
    }
#endif

    if ( nrhs != 23 )
    {
        mexErrMsgIdAndTxt("iceadv_weno_mex:Inputs",
              "no of input arguments should be 23.\nexample matlab call:\nh_new = iceadv_weno_mex(nx,ny,dt,nbdy,\n   h, u, v,LANDMASK, scp2, scp2i, scuy, scvx)\n");
    }
    /*---------------------------------------------------------------------------------------------*/

    /*---------- Input ----------*/
    /*scalar inputs*/
    mwSize nx       = (mwSize)   mxGetScalar(prhs[0]);  /* grid size in x dirn*/
    mwSize ny       = (mwSize)   mxGetScalar(prhs[1]);  /* grid size in y dirn*/
    double dt       = (double)   mxGetScalar(prhs[2]);  /* timestep*/
    mwSize nbdy     = (mwSize)   mxGetScalar(prhs[3]);  /* no of ghost cells (>=4)*/
    mwSize nxext    = nx+2*nbdy;
    mwSize nyext    = ny+2*nbdy;
    mwSize max_no_sections_p,max_no_sections_u,max_no_sections_v;


    /*check size of input arays*/
    for (int j=4;j<11;j++)
    {
        /*h, u, v, scp2, scp2i, scuy, scvx*/
        mwSize M   = mxGetM(prhs[j]);
        mwSize N   = mxGetN(prhs[j]);
        if ( (M!=nx) | (N!=ny) )
        {
            mexErrMsgIdAndTxt("iceadv_weno_mex:Inputs",
              "incorrect size of input argument (5-11)");
        }
    }
    for (int j=11;j<14;j++)
    {
        /*pmask,umask,vmask*/
        mwSize M   = mxGetM(prhs[j]);
        mwSize N   = mxGetN(prhs[j]);
        if ( (M!=nxext) | (N!=nyext) )
        {
            mexErrMsgIdAndTxt("iceadv_weno_mex:Inputs",
              "incorrect size of input argument (12-14)");
        }
    }
    for (int j=14;j<17;j++)
    {
        /*isp,isu,isv*/
        mwSize M   = mxGetM(prhs[j]);
        mwSize N   = mxGetN(prhs[j]);
        if ( (M!=nyext) | (N!=1) )
        {
            mexErrMsgIdAndTxt("iceadv_weno_mex:Inputs",
              "incorrect size of input argument (15-20)");
        }
    }
    for (int j=17;j<23;j++)
    {
        /*ifp,ifu,ifv,ilp,ilu,ilv*/
        mwSize M   = mxGetM(prhs[j]);
        mwSize N   = mxGetN(prhs[j]);
        if (j==17)
        {
            max_no_sections_p = N;
            mwSize N2 = mxGetN(prhs[j+3]);
            if ( N != N2 )
            {
                mexErrMsgIdAndTxt("iceadv_weno_mex:Inputs",
                  "ifp,ilp should have the same no of columns");
            }
        }
        if (j==18)
        {
            max_no_sections_u = N;
            mwSize N2 = mxGetN(prhs[j+3]);
            if ( N != N2 )
            {
                mexErrMsgIdAndTxt("iceadv_weno_mex:Inputs",
                  "ifu,ilu should have the same no of columns");
            }
        }
        if (j==19)
        {
            max_no_sections_v = N;
            mwSize N2 = mxGetN(prhs[j+3]);
            if ( N != N2 )
            {
                mexErrMsgIdAndTxt("iceadv_weno_mex:Inputs",
                  "ifv,ilv should have the same no of columns");
            }
        }
        if ( M != nyext )
            mexErrMsgIdAndTxt("iceadv_weno_mex:Inputs",
              "incorrect size of input argument (18-23)");
    }


    /*pointers to input arrays*/
    /*NB all in matlab ordering*/
    double *h        = (double *) mxGetPr(prhs[4]);
    double *u        = (double *) mxGetPr(prhs[5]);
    double *v        = (double *) mxGetPr(prhs[6]);
    double *scp2     = (double *) mxGetPr(prhs[7]);
    double *scp2i    = (double *) mxGetPr(prhs[8]);
    double *scuy     = (double *) mxGetPr(prhs[9]);
    double *scvx     = (double *) mxGetPr(prhs[10]);
    int    *pmask    = (int    *) mxGetPr(prhs[11]);
    int    *umask    = (int    *) mxGetPr(prhs[12]);
    int    *vmask    = (int    *) mxGetPr(prhs[13]);
    int    *isp      = (int    *) mxGetPr(prhs[14]);
    int    *isu      = (int    *) mxGetPr(prhs[15]);
    int    *isv      = (int    *) mxGetPr(prhs[16]);
    int    *ifp      = (int    *) mxGetPr(prhs[17]);
    int    *ifu      = (int    *) mxGetPr(prhs[18]);
    int    *ifv      = (int    *) mxGetPr(prhs[19]);
    int    *ilp      = (int    *) mxGetPr(prhs[20]);
    int    *ilu      = (int    *) mxGetPr(prhs[21]);
    int    *ilv      = (int    *) mxGetPr(prhs[22]);


    /*define output array*/
    plhs[0] = mxCreateDoubleMatrix(nx, ny, mxREAL);/*pointer to output array*/
    double *h_new   = (double *) mxGetPr(plhs[0]);
    for (mwIndex i=0;i<nx*ny;i++)
        h_new[i]    = h[i];/*initialise to h*/
#if defined(TEST)
    plhs[1] = mxCreateDoubleMatrix(nxext, nyext, mxREAL);/*pointer to output array*/
    double *test_array   = (double *) mxGetPr(plhs[1]);
#endif

    /*do advection
     * h_new is a pointer that is modified by iceAdvWeno
     * - this is also the output */
#if !defined(TEST)
    iceAdvWeno(h_new, u, v, scp2, scp2i, scuy, scvx,
        pmask,umask,vmask,isp,isu,isv,ifp,ifu,ifv,ilp,ilu,ilv,
        dt, (int) nx, (int) ny, (int) nbdy);
#else
    iceAdvWeno(h_new, u, v, scp2, scp2i, scuy, scvx,
        pmask,umask,vmask,isp,isu,isv,ifp,ifu,ifv,ilp,ilu,ilv,
        dt, (int) nx, (int) ny, (int) nbdy, test_array);
#endif

    return;
}
