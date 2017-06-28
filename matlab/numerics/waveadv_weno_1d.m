function h  = waveadv_weno_1d(h,u,grid_prams,dt,adv_options)
%% waveadve_weno_1d.m
%% Author: Timothy Williams
%% Date:   20150405
%% adapted from HYCOM advection code in Fortran 77;

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
%% in:      u,dt are speed and time step;
%% in:      adv_options.ADV_OPT determines periodicity: 0, not; 1, is
%%
%% in:      grid_prams is structure containing grid information:
%%                    nx: 49
%%                    ny: 1 (NOT USED)
%%                  scuy: [49x1 double]   - mesh size at u points in y direction
%%                  scp2: [49x1 double]   - grid box area at p points
%%                 scp2i: [49x1 double]   - inverse of scp2
%%              LANDMASK: [49x1 double]   - 1 on land, 0 on water
%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if sum(abs(h(:)))==0
   %% h=0 everywhere
   %% - nothing to do
   return;
end

ADV_OPT   = adv_options.ADV_OPT;
if ADV_OPT==2
   ADV_OPT   = 0;%1d, so periodicity in y is N/A
end

ii       = grid_prams.nx;
scuy     = grid_prams.scuy(:,1);
scp2i    = grid_prams.scp2i(:,1);
scp2     = grid_prams.scp2(:,1);
LANDMASK = grid_prams.LANDMASK(:,1);
clear grid_prams;

%% size of boundary
%% - need >=4 ghost cells
%%    - 3rd order in space
%%    - 2nd order in time (prediction + correction steps)
%% - if use 3 ghost cells, need to apply boundary conditions between the prediction and correction steps
nbdy        = 4;
idm         = ii;
ireal       = nbdy+(1:idm)';  %%non-ghost i indices
TEST_PLOT   = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%assign values in ghost cells
%%-do this beforehand in case want to parallelise??

%%make all these periodic in i,j (x,y)
advopt_u    = 1;%periodic in x
advopt_grid = 1;%periodic in x
u           = pad_var_1d(u    ,advopt_u   ,nbdy);
scuy        = pad_var_1d(scuy ,advopt_grid,nbdy);%dy: only needed for 1/dx=scp2/scuy
scp2        = pad_var_1d(scp2 ,advopt_grid,nbdy);
scp2i       = pad_var_1d(scp2i,advopt_grid,nbdy);
%%
CFL   = abs(u.*(scuy./scp2)*dt);
if max(CFL)>1
   error(['CFL too large (',num2str(max(CFL)),') - reduce time step (or speed)']);
end

%%ADV_OPT determines how h is extended to ghost cells:
%% 0: zeros in ghost cells
%% 1: periodic in i,j
%% 2: periodic in j only (zeros in ghost cells i<0,i>ii)
h     = pad_var_1d(h,ADV_OPT,nbdy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Prediction step
sao   = weno3pd_v2_1d(h,u,scuy,scp2i,scp2,dt,nbdy);

hp = 0*h;
if nbdy>=4
   %% loop over all cells
   %% - no need to enforce periodicity between prediction and correction steps
   for i_ = 1-nbdy:ii+nbdy
      i     = i_+nbdy;%%1-nbdy->1
      hp(i) = h(i)+dt*sao(i);
   end%i
elseif nbdy==3
   %% enforce periodicity between prediction and correction steps
   for i_ = 1:ii
      i     = i_+nbdy;%%1-nbdy->1
      hp(i) = h(i)+dt*sao(i);
   end%i
   hp = pad_var_1d(hp(ireal),ADV_OPT,nbdy);
else
   error('nbdy should be >=3');
end

if TEST_PLOT
   test_plot(h,nbdy,'h (initial,padded)'); pause;
   test_plot(sao,nbdy,'sao'); pause;
   test_plot(hp,nbdy,'h (predicted)'); pause;
end

% --- Correction step
sao   = weno3pd_v2_1d(hp,u,scuy,scp2i,scp2,dt,nbdy);
if TEST_PLOT
   %% large loop for demo
   for i_ = 1-nbdy:ii+nbdy
      i     = i_+nbdy;
      h(i)  = .5*(h(i)+hp(i)+dt*sao(i));
   end%i
   test_plot(h,nbdy,'h (corrected)'); pause;
else
   for i_ = 1:ii
      i     = i_+nbdy;
      h(i)  = .5*(h(i)+hp(i)+dt*sao(i));
   end%i
end

%%set waves on land to 0
h  = h(ireal).*(1-LANDMASK);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_pad = pad_var_1d(u,ADV_OPT,nbdy)

ii       = length(u);
u_pad    = zeros(ii+2*nbdy,1);
ivec     = 1:ii;
%%
u_pad(ivec+nbdy)  = u;
%%
if (ADV_OPT==1)
   %%make things periodic

   %%make it periodic in i
   ivec1 = 1-nbdy:0;
   ivec2 = ii+1:ii+nbdy;
   %%
   u_pad(ivec1+nbdy) = u(ii-nbdy:ii-1);
   u_pad(ivec2+nbdy) = u(1:nbdy);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sao = weno3pd_v2_1d(g,u,scuy,scp2i,scp2,dt,nbdy)
%sao      = weno3pd_v2_1d(h,u,scuy,scp2i,scp2,dt,nbdy);
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
idm   = ii;
sao   = 0*g;

cq00  = -1./2.;
cq01  = 3./2.;
cq10  = 1./2.;
cq11  = 1./2.;
ca0   = 1./3.;
ca1   = 2./3.;
eps   = 1.e-12;

ful   = zeros(idm+2*nbdy,1);
fuh   = zeros(idm+2*nbdy,1);
gt    = zeros(idm+2*nbdy,1);

%
% --- Compute grid cell boundary fluxes. Split in a low order flux
% --- (donor cell) and a high order correction flux.
%
for i_ = 1-nbdy+2:ii+nbdy-1
   i     = i_+nbdy;%%1-nbdy->1
   im1   = i-1;
 
   if (u(i)>0.)
      im2   = im1-1;%i-2
 
      q0 = cq00*g(im2)+cq01*g(im1);
      q1 = cq10*g(im1)+cq11*g(i  );

      a0 = ca0;
      a1 = ca1*(abs(g(im2)-g(im1))+eps)/(abs(g(im1)-g(i))+eps);

      ful(i) = u(i)*g(im1)*scuy(i);

    else
       ip1  = i+1;
 
       q0   = cq11*g(im1)+cq10*g(i  );
       q1   = cq01*g(i  )+cq00*g(ip1);
 
       a0   = ca1;
       a1   = ca0*(abs(g(im1)-g(i))+eps)/(abs(g(i)-g(ip1))+eps);
 
       ful(i)   = u(i)*g(i)*scuy(i);
    end

    fuh(i)   = u(i)*(a0*q0+a1*q1)/(a0+a1)*scuy(i)-ful(i);

end%i

% --- Update field with low order fluxes.
for i_ = 1-nbdy:ii+nbdy-1
   i     = i_+nbdy;%%1-nbdy->1
   %%
   gt(i) = g(i)-dt*(ful(i+1)-ful(i))*scp2i(i);%scuy is already inside ful
end%i

% --- Obtain fluxes with limited high order correction fluxes.
q  = .25/dt;
for i_ = 1-nbdy+1:ii+nbdy
   i     = i_+nbdy;%%1-nbdy->1

   fuh(i) = ful(i)+max( -q*gt(i)*scp2(i),...
                        min(  q*gt(i-1)*scp2(i-1),...
                              fuh(i)));
end%i

% --- Compute the spatial advective operator.
for i_ = 1-nbdy:ii+nbdy-1
   i        = i_+nbdy;%%1-nbdy->1
   sao(i)   = -(fuh(i+1)-fuh(i))*scp2i(i);
end%i

return

function test_plot(h,nbdy,ttl)
nx    = length(h);
xvec  = (0:nx)'-nbdy;%corners

%%plot steps
X           = xvec(1)+zeros(2*nx,1);%[0 1 1 2 2 ... nx+1]
X(2:2:end)  = xvec(2:end);
X(3:2:end)  = xvec(2:end-1);

H           = zeros(2*nx,1);%[0 1 1 2 2 ... nx+1]
H(1:2:end)  = h;
H(2:2:end)  = h;
plot(X,H,'k');
xlim(X([1,end]));
Y     = get(gca,'ylim');
yav   = mean(Y);
dy    = Y(2)-yav;
Y     = yav+1.25*dy*[-1,1];
set(gca,'ylim',Y);

hold on;
plot(0*Y,Y,':r');
x1 = nx-2*nbdy;
plot(x1+0*Y,Y,':r');
title(ttl);
hold off;
