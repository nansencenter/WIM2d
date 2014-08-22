function h  = advection_weno(h,u,v,scuy,scvx,scp2i,scp2,dt,LANDMASK)
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

ii    = size(h,1);
jj    = size(h,1);
idm   = ii;
jdm   = jj;
%%
nbdy  = 3;              %%size of boundary - need 3 ghost cells since 3rd order in time
ireal = nbdy+(1:idm)';  %%non-ghost i indices
jreal = nbdy+(1:jdm)';  %%non-ghost j indices  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%assign values in ghost cells
%%-do this beforehand in case want to parallelise??
%     real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
%    &  h,u,v,scuy,scvx,scp2i,scp2
%     real dt
tmp            = h;
h              = zeros(idm+2*nbdy,jdm+2*nbdy);
h(ireal,jreal) = tmp;
%%
tmp            = u;
u              = zeros(idm+2*nbdy,jdm+2*nbdy);
u(ireal,jreal) = tmp;
%%
tmp            = v;
v              = zeros(idm+2*nbdy,jdm+2*nbdy);
v(ireal,jreal) = tmp;
%%
dx                = max(max(scvx));
tmp               = scvx;
scvx              = dx+zeros(idm+2*nbdy,jdm+2*nbdy);
scvx(ireal,jreal) = tmp;
%%
dy                = max(max(scuy));
tmp               = scuy;
scuy              = dy+zeros(idm+2*nbdy,jdm+2*nbdy);
scuy(ireal,jreal) = tmp;
%%
Area              = max(max(scp2));
tmp               = scp2;
scp2              = Area+zeros(idm+2*nbdy,jdm+2*nbdy);
scp2(ireal,jreal) = tmp;
%%
Areai                = max(max(scp2i));
tmp                  = scp2i;
scp2i                = Areai+zeros(idm+2*nbdy,jdm+2*nbdy);
scp2i(ireal,jreal)   = tmp;
clear tmp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: sao,hp
%     real dtm
%     integer i,j,l


%sao(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)=0

%%Comment in fortran - don't think valid anymore
% --- Use a modified time step since velocities are in m/s while scale
% --- factors are in cm
%      dtm=dt*1.e2
dtm   = dt;

% --- Prediction step
sao      = weno3pd_v2(h,u,v,scuy,scvx,scp2i,scp2,dtm);
margin   = nbdy;
%     do j=1-margin,jj+margin
%       do i=1-margin,ii+margin!![TW 28.1.2014]
for i_ = 1-margin:ii+margin
for j_ = 1-margin:jj+margin
   i  = i_+nbdy;%%1-nbdy->1
   j  = j_+nbdy;%%1-nbdy->1

   hp(i,j)  = h(i,j)+dtm*sao(i,j);
end%j
end%i

% --- Correction step
sao   = weno3pd_v2(hp,u,v,scuy,scvx,scp2i,scp2,dtm);
%     do j=1-margin,jj+margin
%       do i=1-margin,ii+margin!![TW 28.1.2014]
for i_ = 1-margin:ii+margin
for j_ = 1-margin:jj+margin
   i  = i_+nbdy;%%1-nbdy->1
   j  = j_+nbdy;%%1-nbdy->1

   %%if not on land keep as it is;
   h(i,j)  = .5*(h(i,j)+hp(i,j)+dtm*sao(i,j));
end%j
end%i

%%set waves on land to 0
h  = h(ireal,jreal).*(1-LANDMASK);

return

function sao = weno3pd_v2(g,u,v,scuy,scvx,scp2i,scp2,dt)
%
% --- ------------------------------------------------------------------
% --- By a weighted essentially non-oscillatory scheme with up to 3th
% --- order accuracy, obtain the spatial advective operator of a
% --- 2-dimensional field defined at the scalar points of a C-grid. The
% --- fluxes are limited to make the scheme positive definite.
% --- Advective velocities in the i- and j-direction are defined at u-
% --- and v-points, respectively.
% --- ------------------------------------------------------------------
 
%     real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
%    &  g,sao,u,v,scuy,scvx,scp2i,scp2
%     real dt

%     real cq00,cq01,cq10,cq11,ca0,ca1,eps
%     parameter (cq00=-1./2.,cq01= 3./2.,                               &
%    &           cq10= 1./2.,cq11= 1./2.,                               &
%    &           ca0=1./3.,ca1=2./3.,                                   &
%    &           eps=1.e-12)
nbdy  = 3;
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

%     real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::               &
%    &  ful,fuh,fvl,fvh,gt
%     real q0,q1,a0,a1,q
%     integer i,j,l,im1,im2,ip1,jm1,jm2,jp1
ful   = zeros(idm+2*nbdy,jdm+2*nbdy);
fuh   = zeros(idm+2*nbdy,jdm+2*nbdy);
fvl   = zeros(idm+2*nbdy,jdm+2*nbdy);
fvh   = zeros(idm+2*nbdy,jdm+2*nbdy);
gt    = zeros(idm+2*nbdy,jdm+2*nbdy);

%
% --- Compute grid cell boundary fluxes. Split in a low order flux
% --- (donor cell) and a high order correction flux.
%

%     do j=0,jj+2
%       do i=0,ii+2
for i_ = 0:ii+2
for j_ = 0:jj+2
   i  = i_+nbdy;%%1-nbdy->1
   j  = j_+nbdy;%%1-nbdy->1
   %%
   ful(i,j) = 0;
   fuh(i,j) = 0;
   fvl(i,j) = 0;
   fvh(i,j) = 0;
end%j
end%i

%     do j=0,jj+1
%       do i=0,ii+2%%[TW 28.1.2014]
for i_ = 0:ii+2
for j_ = 0:jj+1
   i  = i_+nbdy;%%1-nbdy->1
   j  = j_+nbdy;%%1-nbdy->1
   %%
   im1   = i-1;
 
   if (u(i,j)>0.)
      %!iu is a water mask (water at point and to left of point);
      %!for waves we make it 1 everywhere and mask waves that go on land later;
      %!im2   = im1-iu(im1,j);
      im2   = im1-1;%i-2
 
      q0 = cq00*g(im2,j)+cq01*g(im1,j);
      q1 = cq10*g(im1,j)+cq11*g(i  ,j);

      a0 = ca0;
      a1 = ca1*(abs(g(im2,j)-g(im1,j))+eps)/(abs(g(im1,j)-g(i  ,j))+eps);

      ful(i,j) = u(i,j)*g(im1,j)*scuy(i,j);

    else
       %!ip1  = i+iu(i+1,j);
       ip1  = i+1;
 
       q0   = cq11*g(im1,j)+cq10*g(i  ,j);
       q1   = cq01*g(i  ,j)+cq00*g(ip1,j);
 
       a0   = ca1;
       a1   = ca0*(abs(g(im1,j)-g(i  ,j))+eps)/(abs(g(i  ,j)-g(ip1,j))+eps);
 
       ful(i,j)   = u(i,j)*g(i  ,j)*scuy(i,j);
    end

    fuh(i,j)   = u(i,j)*(a0*q0+a1*q1)/(a0+a1)*scuy(i,j)-ful(i,j);

end%j
end%i

%     do j=0,jj+2
%       jm1=j-1
%       do i=0,ii+1!![TW 28.1.2014]
for i_ = 0:ii+1
for j_ = 0:jj+2
   i     = i_+nbdy;%%1-nbdy->1
   j     = j_+nbdy;%%1-nbdy->1
   jm1   = j-1;

   if (v(i,j)>0.)
      %!iv is a water mask (water at point and beneath point);
      %!for waves we make it 1 everywhere and mask waves that go on land later;
      %!jm2   = jm1-iv(i,jm1);
      jm2   = jm1-1;

      q0 = cq00*g(i,jm2)+cq01*g(i,jm1);
      q1 = cq10*g(i,jm1)+cq11*g(i,j  );

      a0 = ca0;
      a1 = ca1*(abs(g(i,jm2)-g(i,jm1))+eps)/(abs(g(i,jm1)-g(i,j  ))+eps);

      fvl(i,j) = v(i,j)*g(i,jm1)*scvx(i,j);

   else
      %!jp1   = j+iv(i,j+1);
      jp1   = j+1;

      q0 = cq11*g(i,jm1)+cq10*g(i,j  );
      q1 = cq01*g(i,j  )+cq00*g(i,jp1);

      a0 = ca1;
      a1 = ca0*(abs(g(i,jm1)-g(i,j  ))+eps)/(abs(g(i,j  )-g(i,jp1))+eps);

      fvl(i,j) = v(i,j)*g(i,j  )*scvx(i,j);
   end

   fvh(i,j) = v(i,j)*(a0*q0+a1*q1)/(a0+a1)*scvx(i,j)-fvl(i,j);
end%j
end%i

% --- Update field with low order fluxes.
%     do j=0,jj+1
%       do i=0,ii+1!![TW 28.1.2014]
for i_ = 0:ii+1
for j_ = 0:jj+1
   i     = i_+nbdy;%%1-nbdy->1
   j     = j_+nbdy;%%1-nbdy->1
   %%
   gt(i,j)  = g(i,j)-dt*(ful(i+1,j)-ful(i,j) +fvl(i,j+1)-fvl(i,j))*scp2i(i,j);
end%j
end%i

% --- Obtain fluxes with limited high order correction fluxes.
q  = .25/dt;
%     do j=1,jj
%       do i=1,ii+1!![TW 28.1.2014]
for i_ = 0:ii+1
for j_ = 0:jj
   i     = i_+nbdy;%%1-nbdy->1
   j     = j_+nbdy;%%1-nbdy->1

   fuh(i,j) = ful(i,j)+max(-q*gt(i  ,j)*scp2(i  ,j), min( q*gt(i-1,j)*scp2(i-1,j),fuh(i,j)));
end%j
end%i

%     do j=1,jj+1
%       do i=1,ii!![TW 28.1.2014]
for i_ = 0:ii
for j_ = 0:jj+1
   i     = i_+nbdy;%%1-nbdy->1
   j     = j_+nbdy;%%1-nbdy->1

   fvh(i,j) = fvl(i,j)+max(-q*gt(i,j  )*scp2(i,j  ), min( q*gt(i,j-1)*scp2(i,j-1),fvh(i,j)));
end%j
end%i

% --- Compute the spatial advective operator.
%     do j=1,jj
%       do i=1,ii!![TW 28.1.2014]
for i_ = 0:ii
for j_ = 0:jj+1
   i     = i_+nbdy;%%1-nbdy->1
   j     = j_+nbdy;%%1-nbdy->1

   sao(i,j) = -(fuh(i+1,j)-fuh(i,j)+fvh(i,j+1)-fvh(i,j))*scp2i(i,j);
end%j
end%i

return
