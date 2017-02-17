function h  = iceadv_weno(h,u,v,grid_prams,masks,dt,nbdy)
%% advection_weno.m
%% Author: Timothy Williams
%% Date:   20140821, 05:40:03 CEST
%%adapted from HYCOM advection code in Fortran 77;

%
% --- ------------------------------------------------------------------
% --- Advection is done with flux limited 3rd order WENO in space and
% --- 2nd order Runge-Kutta in time
% --- ------------------------------------------------------------------
%

%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% modified [TW 28.1.2014] for wave advection
%% - iceadv builds up ice (& thus waves) at coast;
%% - iceadv_v2 lets waves go into/out of land,
%%   then cancels the waves on land at end of routine;
%% - just adds the wave cancellation to Luke's original fix
%%   for waves on boundaries;
%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% ARGUMENTS:
%% in/out:  h is qty to be advected;
%%
%% in:      u,v,dt are speed components and time step;
%% in:      scuy,scvx is mesh size at u (v) points in y (x) direction
%%           - see common_blocks.h;
%% in:      scp2, scp2i are grid box area at p points, and its inverse;
%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if sum(abs(h(:)))==0
   %% h = 0 everywhere
   %% - nothing to do
   return;
end

ii       = grid_prams.nx;
jj       = grid_prams.ny;
scuy     = grid_prams.scuy;
scvx     = grid_prams.scvx;
scp2i    = grid_prams.scp2i;
scp2     = grid_prams.scp2;
LANDMASK = grid_prams.LANDMASK;
%clear grid_prams;

idm   = ii;
jdm   = jj;

%% size of boundary
%% - need >=4 ghost cells
%%    - 3rd order in space
%%    - 2nd order in time (prediction + correction steps)
%% - if we use 3, we need to apply the boundary conditions between
%%    the prediction & correction steps
nbdy        = 4;
ireal       = nbdy+(1:idm)';  %%non-ghost i indices
jreal       = nbdy+(1:jdm)';  %%non-ghost j indices  
TEST_PLOT   = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%assign values in ghost cells
%%-do this beforehand in case want to parallelise??

adv_opt_uv     = 0;%not periodic
adv_opt_grid   = 1;%x-y periodic

u     = pad_var(u     ,adv_opt_uv   ,nbdy);
v     = pad_var(v     ,adv_opt_uv   ,nbdy);
scvx  = pad_var(scvx  ,adv_opt_grid ,nbdy);
scuy  = pad_var(scuy  ,adv_opt_grid ,nbdy);
scp2  = pad_var(scp2  ,adv_opt_grid ,nbdy);
scp2i = pad_var(scp2i ,adv_opt_grid ,nbdy);

%%ADV_OPT determines how h is extended to ghost cells:
%% 0: zeros in ghost cells
%% 1: periodic in i,j
%% 2: periodic in j only (zeros in ghost cells i<0,i>ii)
ADV_OPT  = 0;
h        = pad_var(h,ADV_OPT,nbdy);
%jtst  = (ii+2*nbdy)+(-12:0);
%tst2d = h(jtst,4),pause
%pcolor(h),colorbar,GEN_pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Prediction step
sao   = weno3pd_v2(h,u,v,scuy,scvx,scp2i,scp2,dt,nbdy,masks);
%tst2d = sao(jtst,4),pause

if nbdy>=4
   %% no need to enforce periodicity between prediction and correction steps
   %% - loop over all cells (including ghost cells)
   hp = 0*h;
   for j_ = 1-nbdy:jj+nbdy
      j  = j_+nbdy;%%1-nbdy->1
      for l=1:masks.isp(j)
         for i_=max(1-nbdy,masks.ifp(j,l)):min(ii+nbdy,masks.ilp(j,l))
            i  = i_+nbdy;%%1-nbdy->1
            hp(i,j)  = h(i,j)+dt*sao(i,j);
         end%i - rows
      end%l - sections
   end%j - columns
%elseif nbdy==3
%   %% enforce periodicity between prediction and correction steps
%   %% - only loop over non-ghost cells
%   hp = zeros(ii,jj);
%   for i_ = 1:ii
%   for j_ = 1:jj
%      i           = i_+nbdy;%%1-nbdy->1
%      j           = j_+nbdy;%%1-nbdy->1
%      hp(i_,j_)   = h(i,j)+dt*sao(i,j);
%   end%j
%   end%i
%   hp = pad_var(hp,ADV_OPT,nbdy);
else
   error('nbdy should be >=4');
end

if TEST_PLOT
   imshow_array(h,nbdy,'h (initial,padded)');pause;
   imshow_array(sao,nbdy,'sao');pause;
   imshow_array(hp,nbdy,'h (prediction)');pause;
end

% --- Correction step
sao   = weno3pd_v2(hp,u,v,scuy,scvx,scp2i,scp2,dt,nbdy,masks);

%%set waves on land to 0
if TEST_PLOT
   %larger loop for demo
   for i_ = 1-nbdy:ii+nbdy
   for j_ = 1-nbdy:jj+nbdy
      i        = i_+nbdy;%%1-nbdy->1
      j        = j_+nbdy;%%1-nbdy->1
      h(i,j)   = .5*(h(i,j)+hp(i,j)+dt*sao(i,j));
   end%j
   end%i
   imshow_array(h,nbdy,'h (corrected)');pause;
else
   for j_ = 1:jj
      j  = j_+nbdy;%%1-nbdy->1
      for l=1:masks.isp(j)
         for i_=max(1-nbdy,masks.ifp(j,l)):min(ii+nbdy,masks.ilp(j,l))
            i        = i_+nbdy;%%1-nbdy->1
            h(i,j)   = .5*(h(i,j)+hp(i,j)+dt*sao(i,j));
         end%i-rows
      end%l-sections
   end%j-columns
end
h  = h(ireal,jreal).*(1-LANDMASK);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_pad = pad_var(u,OPT,nbdy)

[ii,jj]  = size(u);%ii,jj,nbdy
u_pad    = zeros(ii+2*nbdy,jj+2*nbdy);
ivec     = 1:ii;
jvec     = 1:jj;
%%
u_pad(ivec+nbdy,jvec+nbdy)  = u;

if (OPT==1)
   %%make things periodic

   %%make it periodic in i
   ivec1 = 1-nbdy:0;
   ivec2 = ii+1:ii+nbdy;
   %%
   u_pad(ivec1+nbdy,jvec+nbdy) = u(ii-nbdy:ii-1,:);
   u_pad(ivec2+nbdy,jvec+nbdy) = u(1:nbdy,:);

   %%make it periodic in j
   jvec1 = 1-nbdy:0;
   jvec2 = jj+1:jj+nbdy;
   %%
   u_pad(ivec+nbdy,jvec1+nbdy) = u(1:ii,jj-nbdy:jj-1);
   u_pad(ivec+nbdy,jvec2+nbdy) = u(1:ii,1:nbdy);

   %%BR,TL
   ivec1 = ii+1:ii+nbdy;
   ivec2 = 1-nbdy:0;
   jvec1 = jj+1:jj+nbdy;
   jvec2 = 1-nbdy:0;
   %%
   u_pad(ivec1+nbdy,jvec1+nbdy)  = u(1:nbdy,1:nbdy);
   u_pad(ivec2+nbdy,jvec2+nbdy)  = u(ii-nbdy:ii-1,jj-nbdy:jj-1);

   %%BL,TR
   ivec1 = ii+1:ii+nbdy;
   ivec2 = 1-nbdy:0;
   jvec1 = 1-nbdy:0;
   jvec2 = jj+1:jj+nbdy;
   %%
   u_pad(ivec1+nbdy,jvec1+nbdy)  = u(1:nbdy,jj-nbdy:jj-1);
   u_pad(ivec2+nbdy,jvec2+nbdy)  = u(ii-nbdy:ii-1,1:nbdy);

elseif (OPT==2)
   %%make it periodic in j (y) only
   jvec1 = 1-nbdy:0;
   jvec2 = jj+1:jj+nbdy;
   %%
   u_pad(ivec+nbdy,jvec1+nbdy) = u(1:ii,jj-nbdy:jj-1);
   u_pad(ivec+nbdy,jvec2+nbdy) = u(1:ii,1:nbdy);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sao = weno3pd_v2(g,u,v,scuy,scvx,scp2i,scp2,dt,nbdy,masks)
%
% --- ------------------------------------------------------------------
% --- By a weighted essentially non-oscillatory scheme with up to 3th
% --- order accuracy, obtain the spatial advective operator of a
% --- 2-dimensional field defined at the scalar points of a C-grid. The
% --- fluxes are limited to make the scheme positive definite.
% --- Advective velocities in the i- and j-direction are defined at u-
% --- and v-points, respectively.
% --- ------------------------------------------------------------------
 
ii    = size(g,1)-2*nbdy; 
jj    = size(g,2)-2*nbdy; 
idm   = ii;
jdm   = jj;
sao   = 0*g;

cq00  = -1./2.;
cq01  = 3./2.;
cq10  = 1./2.;
cq11  = 1./2.;
ca0   = 1./3.;
ca1   = 2./3.;
eps   = 1.e-12;

ful   = zeros(idm+2*nbdy,jdm+2*nbdy);
fuh   = zeros(idm+2*nbdy,jdm+2*nbdy);
fvl   = zeros(idm+2*nbdy,jdm+2*nbdy);
fvh   = zeros(idm+2*nbdy,jdm+2*nbdy);
gt    = zeros(idm+2*nbdy,jdm+2*nbdy);
jtst  = idm+2*nbdy+(-12:0)';

%
% --- Compute grid cell boundary fluxes. Split in a low order flux
% --- (donor cell) and a high order correction flux.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fluxes in x dirn
for j_ = 1-nbdy:jj+nbdy
   j  = j_+nbdy;
   for l=1:masks.isu(j)
      for i_ = max(1-nbdy+2,masks.ifu(j,l)):min(ii+nbdy-1,masks.ilu(j,l))
         i     = i_+nbdy;
         im1   = i-1;
       
         if (u(i,j)>0.)
            im2   = im1-masks.umask(im1,j);%i-2
       
            q0 = cq00*g(im2,j)+cq01*g(im1,j);
            q1 = cq10*g(im1,j)+cq11*g(i  ,j);

            a0 = ca0;
            a1 = ca1*(abs(g(im2,j)-g(im1,j))+eps)/(abs(g(im1,j)-g(i  ,j))+eps);

            ful(i,j) = u(i,j)*g(im1,j)*scuy(i,j);

          else
             ip1  = i+masks.umask(i+1,j);
       
             q0   = cq11*g(im1,j)+cq10*g(i  ,j);
             q1   = cq01*g(i  ,j)+cq00*g(ip1,j);
       
             a0   = ca1;
             a1   = ca0*(abs(g(im1,j)-g(i  ,j))+eps)/(abs(g(i  ,j)-g(ip1,j))+eps);
       
             ful(i,j)   = u(i,j)*g(i  ,j)*scuy(i,j);
          end

          fuh(i,j)   = u(i,j)*(a0*q0+a1*q1)/(a0+a1)*scuy(i,j)-ful(i,j);

      end%i-rows
   end%l-sections
end%j-columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%fluxes in y dirn:
for j_ = 1-nbdy+2:jj+nbdy-1
   j  = j_+nbdy;
   for l=1:masks.isv(j);
      for i_ = max(1-nbdy,masks.ifv(j,l)):min(ii+nbdy,masks.ilv(j,l))
         i     = i_+nbdy;
         jm1   = j-1;

         if (v(i,j)>0.)
            jm2   = jm1-masks.vmask(i,jm1);%j-2

            q0 = cq00*g(i,jm2)+cq01*g(i,jm1);
            q1 = cq10*g(i,jm1)+cq11*g(i,j  );

            a0 = ca0;
            a1 = ca1*(abs(g(i,jm2)-g(i,jm1))+eps)/(abs(g(i,jm1)-g(i,j  ))+eps);

            fvl(i,j) = v(i,j)*g(i,jm1)*scvx(i,j);

         else
            jp1   = jm1+masks.vmask(i,j+1);%j+1

            q0 = cq11*g(i,jm1)+cq10*g(i,j  );
            q1 = cq01*g(i,j  )+cq00*g(i,jp1);

            a0 = ca1;
            a1 = ca0*(abs(g(i,jm1)-g(i,j  ))+eps)/(abs(g(i,j  )-g(i,jp1))+eps);

            fvl(i,j) = v(i,j)*g(i,j  )*scvx(i,j);
         end

         fvh(i,j) = v(i,j)*(a0*q0+a1*q1)/(a0+a1)*scvx(i,j)-fvl(i,j);
      end%i-rows
   end%l-sections
end%j-columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Update field with low order fluxes.
for j_ = 1-nbdy:jj+nbdy-1
   j  = j_+nbdy;
   for l=1:masks.isp(j)
      for i_ = max(1-nbdy,masks.ifp(j,l)):min(ii+nbdy-1,masks.ilp(j,l))
         i        = i_+nbdy;
         gt(i,j)  = g(i,j)-dt*(ful(i+1,j)-ful(i,j) +fvl(i,j+1)-fvl(i,j))*scp2i(i,j);
      end%i-rows
   end%l-sections
end%j-columns
%tst2d = gt(jtst,4),pause

% --- Obtain fluxes with limited high order correction fluxes.
q  = .25/dt;
for j_ = 1-nbdy:jj+nbdy
   j  = j_+nbdy;
   for l=1:masks.isu(j)
      for i_ = max(1-nbdy+1,masks.ifu(j,l)):min(ii+nbdy,masks.ilu(j,l))
         i        = i_+nbdy;
         fuh(i,j) = ful(i,j)+max(   -q*gt(i  ,j)*scp2(i  ,j),...
                                    min( q*gt(i-1,j)*scp2(i-1,j),...
                                          fuh(i,j))   );
      end%i-rows
   end%l-sections
end%j-columns

for j_ = 1-nbdy+1:jj+nbdy
   j  = j_+nbdy;
   for l=1:masks.isv(j)
      for i_ = max(1-nbdy,masks.ifv(j,l)):min(ii+nbdy,masks.ilv(j,l))
         i        = i_+nbdy;
         fvh(i,j) = fvl(i,j)+max(   -q*gt(i,j  )*scp2(i,j  ),...
                                    min( q*gt(i,j-1)*scp2(i,j-1),...
                                         fvh(i,j)) );
      end%i-rows
   end%l-sections
end%j-columns

% --- Compute the spatial advective operator.
for j_ = 1-nbdy:jj+nbdy-1
   j  = j_+nbdy;
   for l=1:masks.isp(j)
      for i_ = max(1-nbdy,masks.ifp(j,l)):min(ii+nbdy-1,masks.ilp(j,l))
         i        = i_+nbdy;
         sao(i,j) = -(fuh(i+1,j)-fuh(i,j)+fvh(i,j+1)-fvh(i,j))*scp2i(i,j);
      end%i-rows
   end%l-sections
end%j-columns

return

function imshow_array(arr,nbdy,ttl)

[nx,ny]  = size(arr);

%make (0,0) a corner
yvec     = (1:ny+1)-nbdy-1;
xvec     = (1:nx+1)-nbdy-1;
[Y,X]    = meshgrid(yvec,xvec);

%pad array, since pcolor in matlab ignores last row & column
clf;
arr2              = zeros(nx+1,ny+1);
arr2(1:nx,1:ny)   = arr;
P                 = pcolor(X,Y,arr2);

%plot boundary of 'real' domain
hold on;
x1 = nx+1-2*nbdy;
y1 = ny+1-2*nbdy;
cc = 'k';
xx = [0,x1]; plot(xx,0*xx,cc,'linewidth',2); plot(xx,y1+0*xx,cc,'linewidth',2);
yy = [0,y1]; plot(0*yy,yy,cc,'linewidth',2); plot(x1+0*yy,yy,cc,'linewidth',2); 

% colorbar and title
set(P,'edgecolor','none');
colorbar;
title(ttl);
fn_fullscreen;
hold off;
