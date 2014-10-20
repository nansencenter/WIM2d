%% adv_atten_timestep_simple.m
%% Author: Timothy Williams
%% Date: 20141018, 18:04:46 CEST

function S = adv_atten_timestep_simple(grid_prams,ice_prams,s1,dt)

ndir        = s1.ndir;
wavdir      = s1.wavdir;
S           = s1.Sdir;%% size(S)=[nx,ny,ndir] - do 1 freq at a time;
ag_eff      = s1.ag_eff;
atten_dim   = s1.atten_dim;
ICE_MASK    = s1.ICE_MASK;
clear s1;

theta = pi/180*(90-wavdir);%%waves-to, anti-clockwise, radians

%%advection;
for jth  = 1:ndir
   u  = ag_eff*cos(theta(jth));
   v  = ag_eff*sin(theta(jth));
   S(:,:,jth)  = waveadv_weno(S(:,:,jth),u,v,grid_prams,dt);
end

%%attenuation
nx = grid_prams.nx;
ny = grid_prams.ny;
clear grid_prams;

for i = 1:nx
for j = 1:ny
   %% atten_dim = ENERGY attenuation coeff [m^{-1}]
   if ICE_MASK(i,j)>0
      S(i,j,:) = S(i,j,:)*...
         exp(-atten_dim(i,j)*ag_eff(i,j)*dt);
   end
end
end
