function [Hs,Tp,mwd] = fn_spectral_integrals(om_vec,wavdir,Sdir)
%% function for testing:
%% split energy between fwd and back directions
%% also test integration schemes (const panel vs Simpson's)

nw    = length(om_vec);
ndir  = length(wavdir);
%%
sz = size(Sdir);
nx = sz(1);
ny = sz(2);

if nw==1
   wt_om = 1;
else
   %%simpson's rule 
   %%NB needs odd number of points
   %%NB om_vec needs to be equally spaced
   dom    = abs(om_vec(2)-om_vec(1));
   wt_om  = wt_simpsons(nw,dom);
end

%% integrate over frequencies:
mom0  = zeros(nx,ny,ndir);
mom2  = zeros(nx,ny,ndir);
for w=1:nw
   om    = om_vec(w);
   mom0  = mom0+wt_om(w)*abs(Sdir(:,:,:,w));
   mom2  = mom2+wt_om(w)*om^2*abs(Sdir(:,:,:,w));
      %%take abs so Sdir is always positive (can become negative due to advection)
end


%% integrate over directions (splitting fwd/back):
th_vec   = -pi/180*(wavdir+90);%%-90 deg (from W) -> 0 radians; 0 deg (from N) -> -pi/2 radians
dth      = th_vec(2)-th_vec(1);
%%
Mom0  = zeros(nx,ny);
Tp    = zeros(nx,ny);
mwd   = zeros(nx,ny);
for wth=1:ndir
   Mom0  = Mom0+mom0(:,:,wth)*th_vec(wth)*dth;
end
%%
mom0  = dth*sum(mom0,3);
mom2  = dth*sum(mom2,3);
%%
Hs       = 4*sqrt(mom0);
jp       = find(Hs>0);
Tp(jp)   = 2*pi*sqrt(mom0(jp)./mom2(jp));
mwd(jp)  = -90 -180/pi*sqrt(Mom0(jp)./mom0(jp));
