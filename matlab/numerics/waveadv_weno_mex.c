#include "mex.h"
#include "math.h"
#include "string.h"

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _b : _a; })
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

void weno3pdV2(double* gin, double* u, double* v, double* scuy,
                       double* scvx, double* scp2i, double* scp2, double* saoout,
                       double dt, int nx, int ny, int nbdy)
{

    double cq00=-1./2 ,cq01=3./2, cq10=1./2, cq11=1./2, ca0=1./3, ca1=2./3, eps=1e-12;
    double q0, q1, a0, a1, q;
    int im1, im2, ip1, jm1, jm2, jp1, ymargin;
    int nxext,nyext;
    
    nxext = nx+nbdy;
    nyext = ny+nbdy;
    double *ful   = calloc(nxext*nyext,sizeof(double));
    double *fuh   = calloc(nxext*nyext,sizeof(double));
    double *fvl   = calloc(nxext*nyext,sizeof(double));
    double *fvh   = calloc(nxext*nyext,sizeof(double));
    double *gt    = calloc(nxext*nyext,sizeof(double));


    /*fluxes in x direction*/
    for (int i = 2; i < nxext-1; i++)
    {
        for (int j = 0; j < nyext; j++)
        {
            double q0, q1, a0, a1, q;
            int im1, im2, ip1, jm1, jm2, jp1, ymargin;

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
    }

    /*fluxes in y direction*/
    for (int i = 0; i < nxext; i++)
    {
        for (int j = 2; j < nyext-1; j++)
        {
            double q0, q1, a0, a1, q;
            int im1, im2, ip1, jm1, jm2, jp1;

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
    }


    /* update field with low order fluxes*/
    for (int i = 0; i < nxext-1; i++)
    {
        for (int j = 0; j < nyext; j++)
        {
            gt[i*nyext+j] = gin[i*nyext+j]-dt*(ful[(i+1)*nyext+j]
                  -ful[i*nyext+j]+fvl[i*nyext+j+1]-fvl[i*nyext+j])*scp2i[i*nyext+j];
        }
    }

    q = 0.25/dt;

    /* obtain fluxes with limited high order correction fluxes*/
    for (int i = 1; i < nxext; i++)
    {
        for (int j = 0; j < nyext; j++)
        {
            fuh[i*nyext+j] = ful[i*nyext+j]+
               max(-q*gt[i*nyext+j]*scp2[i*nyext+j],
                        min(q*gt[(i-1)*nyext+j]*scp2[(i-1)*nyext+j],fuh[i*nyext+j]));
        }
    }

    /* obtain fluxes with limited high order correction fluxes*/
    for (int i = 0; i < nxext; i++)
    {
        for (int j = 1; j < nyext; j++)
        {
            fvh[i*nyext+j]=fvl[i*nyext+j]+
               max(-q*gt[i*nyext+j]*scp2[i*nyext+j],
                        min(q*gt[i*nyext+j-1]*scp2[i*nyext+j-1],fvh[i*nyext+j]));
        }
    }

    /* compute the spatial advective operator*/
    for (int i = 0; i < nxext-1; i++)
    {
        for (int j = 0; j < nyext; j++)
        {
            saoout[i*nyext+j] = -(fuh[(i+1)*nyext+j]-fuh[i*nyext+j]+fvh[i*nyext+j+1]
                -fvh[i*nyext+j])*scp2i[i*nyext+j];
        }
    }

    /* deallocate memory for temp arrays*/
    free(ful);
    free(fuh);
    free(fvl);
    free(fvh);
    free(gt);
}


void padVar(double* u, double* upad, int advopt_, int nx, int ny, int nbdy)
{

   int nxext   = nx+nbdy;
   int nyext   = ny+nbdy;

    for (int i = 0; i < nxext; i++)
    {
        for (int j = 0; j < nyext; j++)
        {

            if ((nbdy-1 < i) && (i < nx+nbdy) && (nbdy-1 < j) && (j < ny+nbdy))
            {
                upad[i*nyext+j] = u[(i-nbdy)*ny+j-nbdy];
            }

            if (advopt_ != 0)/*either y-periodic or xy-periodic*/
            {
                /* make periodic in j
                   - lower cells*/
                bool i_inner = ((nbdy-1 < i) && (i < nx+nbdy));
                if ((j < nbdy) && i_inner)
                    upad[i*nyext+j] = u[(i-nbdy)*ny+ny-nbdy+j];

                /* - upper cells */
                if ((ny+nbdy-1 < j) && i_inner)
                    upad[i*nyext+j] = u[(i-nbdy)*ny+j-ny-nbdy];
            }

            if (advopt_ == 1)/*xy-periodic*/
            {
                /* make periodic in i
                   - far-left cells */
                bool j_inner = ((nbdy-1 < j) && (j < ny+nbdy));
                if ((i < nbdy) && j_inner )
                    upad[i*nyext+j] = u[(nx-nbdy+i)*ny+j-nbdy];

                /* - far-right cells */
                if ((nx+nbdy-1 < i) && j_inner )
                    upad[i*nyext+j] = u[(i-nx-nbdy)*ny+j-nbdy];

                /* TR */
                if ((nx+nbdy-1 < i) && (ny+nbdy-1 < j))
                    upad[i*nyext+j] = u[(i-nx-nbdy)*ny+j-ny-nbdy];

                /* BL */
                if ((i < nbdy) && (j < nbdy))
                    upad[i*nyext+j] = u[(i+nx-nbdy)*ny+j];

                /* BR */
                if ((nx+nbdy-1 < i) && (j < nbdy))
                    upad[i*nyext+j] = u[(i-nx-nbdy)*ny+ny-nbdy+j];

                /* TL */
                if ((i < nbdy) && (ny+nbdy-1 < j))
                    upad[i*nyext+j] = u[(i+nx-nbdy)*ny+j-ny-nbdy];
            }/*"xy-periodic" */
        }/*j*/
    }/*i*/
}


void waveAdvWeno(double* h, double* u, double* v,double* LANDMASK,
      double* scp2, double* scp2i, double* scuy, double* scvx,
      double dt, int nx, int ny, int nbdy,int advopt)
{
    int nxext = nx+nbdy;
    int nyext = ny+nbdy;

    /*allocate memory of temp arrays*/
    double *sao        = calloc(nxext*nyext,sizeof(double));
    double *u_pad      = calloc(nxext*nyext,sizeof(double));
    double *v_pad      = calloc(nxext*nyext,sizeof(double));
    double *h_pad      = calloc(nxext*nyext,sizeof(double));
    double *scp2_pad   = calloc(nxext*nyext,sizeof(double));
    double *scp2i_pad  = calloc(nxext*nyext,sizeof(double));
    double *scuy_pad   = calloc(nxext*nyext,sizeof(double));
    double *scvx_pad   = calloc(nxext*nyext,sizeof(double));
    double *hp         = calloc(nxext*nyext,sizeof(double));


    /*boundary conditions*/
    int advopt_uv    = 1;/*xy-periodic*/
    int advopt_grid  = 1;/*xy-periodic*/
    padVar(u, u_pad, advopt_uv,nx,ny,nbdy);
    padVar(v, v_pad,advopt_uv,nx,ny,nbdy);
    padVar(scp2, scp2_pad,advopt_grid,nx,ny,nbdy);
    padVar(scp2i, scp2i_pad,advopt_grid,nx,ny,nbdy);
    padVar(scuy, scuy_pad,advopt_grid,nx,ny,nbdy);
    padVar(scvx, scvx_pad,advopt_grid,nx,ny,nbdy);
    padVar(h, h_pad,advopt,nx,ny,nbdy);/*this field uses the variable boundary condition advopt*/

    /* prediction step*/
    weno3pdV2(h_pad, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad,
          sao,dt,nx,ny,nbdy);

    if (nbdy>3)
    {
       /* if using nghost>3, don't need to apply boundary conditions
           before correction step,
           but need to loop over full padded domain*/

       for (int i = 0; i < nxext; i++)
       {
           for (int j = 0; j < nyext; j++)
           {
               hp[i*nyext+j] = h_pad[i*nyext+j]+dt*sao[i*nyext+j];
           }
       }
    }
    else
    {
       printf("%s","\nAdvection (WENO): 'nbdy' should be >=3\n");
       abort();
    }

    /* correction step*/
    weno3pdV2(hp, u_pad, v_pad, scuy_pad, scvx_pad, scp2i_pad, scp2_pad,
          sao,dt,nx,ny,nbdy);


    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            h[i*ny+j] = 0.5*(h_pad[(i+nbdy)*nyext+j+nbdy]+hp[(i+nbdy)*nyext+j+nbdy]+dt*sao[(i+nbdy)*nyext+j+nbdy]);

            /*mask land cells*/
            h[i*ny+j] *= 1-LANDMASK[i*ny+j];
        }
    }

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
      waveAdvWeno(double* h, double* u, double* v,int* LANDMASK,
         double* scp2, double* scp2i, double* scuy, double* scvx,
         double dt, int nx, int ny, int nbdy,int adv_opt)

    matlab call is:
    h_new = waveAdvWeno(h, u, v,LANDMASK, scp2, scp2i, scuy, scvx,
               dt, nx, ny, nbdy,advopt)
*/

    /*---------- Input ----------*/
    /* pointers to input arrays - want to change to vectors*/
    double **prhs0 = (double **) mxGetPr(prhs[0]);
    double **prhs1 = (double **) mxGetPr(prhs[1]);
    double **prhs2 = (double **) mxGetPr(prhs[2]);
    double **prhs3 = (double **) mxGetPr(prhs[3]);
    double **prhs4 = (double **) mxGetPr(prhs[4]);
    double **prhs5 = (double **) mxGetPr(prhs[5]);
    double **prhs6 = (double **) mxGetPr(prhs[6]);
    double **prhs7 = (double **) mxGetPr(prhs[7]);

    /*scalar inputs*/
    double dt        = (double)   mxGetScalar(prhs[8]);   /* timestep*/
    mwSize nx        = (mwSize)   mxGetScalar(prhs[9]);   /* grid size in x dirn*/
    mwSize ny        = (mwSize)   mxGetScalar(prhs[10]);  /* grid size in y dirn*/
    mwSize nbdy      = (mwSize)   mxGetScalar(prhs[11]);  /* number of ghost cells to use*/
    int    advopt    = (int)      mxGetScalar(prhs[12]);  /* periodicity option:
                                                             0=not periodic; 1=x-y-periodic; 2=y-periodic*/

    /*convert to vectors*/
    double *h        = calloc(nx*ny,sizeof(double));
    double *u        = calloc(nx*ny,sizeof(double));
    double *v        = calloc(nx*ny,sizeof(double));
    double *LANDMASK = calloc(nx*ny,sizeof(double));
    double *scp2     = calloc(nx*ny,sizeof(double));
    double *scp2i    = calloc(nx*ny,sizeof(double));
    double *scuy     = calloc(nx*ny,sizeof(double));
    double *scvx     = calloc(nx*ny,sizeof(double));

    mwIndex i,j;
    for (i=0;i<nx;i++)
    {
        for (j=0;j<ny;j++)
        {
            /*TODO check ordering - matlab vs c*/
            h        [i*ny+j] = prhs0[i][j]; /* field to be advected*/
            u        [i*ny+j] = prhs1[i][j]; /* u velocity field*/
            v        [i*ny+j] = prhs2[i][j]; /* v velocity field*/
            LANDMASK [i*ny+j] = prhs3[i][j]; /* 1 on land, 0 else*/
            scp2     [i*ny+j] = prhs4[i][j]; /* area of grid cell (scpx*scpy)*/
            scp2i    [i*ny+j] = prhs5[i][j]; /* 1/scp2*/
            scuy     [i*ny+j] = prhs6[i][j]; /* dy at u points*/
            scvx     [i*ny+j] = prhs7[i][j]; /* dx at v points*/
        }
    }


    /*do advection
     * h is a pointer that is modified by waveAdvWeno
     * - this is now the output */
    waveAdvWeno(h, u, v,LANDMASK, scp2, scp2i, scuy, scvx,
      dt, (int) nx, (int) ny, (int) nbdy,advopt);

    /*define and set output array*/
    plhs[0]          = mxCreateDoubleMatrix(nx, ny, mxREAL);
    double **h_new   = (double **) mxGetPr(plhs[0]);
    for (i=0;i<nx;i++)
       for (j=0;j<ny;j++)
          h_new[i][j]   = h[i*ny+j];
        
    /*deallocate memory in arrays*/
    free(h);
    free(u);
    free(v);
    free(LANDMASK);
    free(scp2);
    free(scp2i);
    free(scuy);
    free(scvx);
    return;
}
