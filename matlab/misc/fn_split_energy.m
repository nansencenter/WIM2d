function [Ep,Em,Et1,Et2] = fn_split_energy(om_vec,wavdir,Sdir)
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
Edir  = zeros(nx,ny,ndir);
for w=1:nw
   Edir  = Edir+wt_om*abs(Sdir(:,:,:,w));
      %%take abs so Sdir is always positive (can become negative due to advection)
end

%% integrate over directions (splitting fwd/back):
th_vec   = -pi/180*(wavdir+90);%%-90 deg (from W) -> 0 radians; 0 deg (from N) -> -pi/2 radians
dth      = th_vec(2)-th_vec(1);
jp       = find(cos(th_vec)>=0);
jm       = find(cos(th_vec)<=0);
%%
nd0   = length(jp);
wt_th = wt_simpsons(nd0,dth);
%%
Ep = zeros(nx,ny);
for jj=1:length(jp)
   j  = jp(jj);
   Ep = Ep+wt_th(jj)*Edir(:,:,j);
end
%%
Em = zeros(nx,ny);
for jj=1:length(jm)
   j  = jm(jj);
   Em = Em+wt_th(jj)*Edir(:,:,j);
end
%%
Et1      = zeros(nx,ny);
Et2      = zeros(nx,ny);
wt2      = wt_simpsons(ndir+1,dth); % Simpson's rule
wt2(1)   = 2*wt2(1);                % double 1st element since repeated
wt2(end) = [];                      % don't need last element since it's same as 1st
for j=1:ndir
   Et1 = Et1+dth*Edir(:,:,j);    % standard integration for periodic functions
   Et2 = Et2+wt2(j)*Edir(:,:,j); % Simpson's rule
end
