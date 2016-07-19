%% atten_simple.m
%% Author: Timothy Williams
%% Date: 20141018, 18:04:46 CEST

function [S,S_freq,tau_x,tau_y] = ...
   atten_simple(grid_prams,ice_prams,s1,dt)

nx = grid_prams.nx;
ny = grid_prams.ny;
clear grid_prams;

ndir        = s1.ndir;
wavdir      = s1.wavdir;
S           = s1.Sdir;%% size(S)=[nx,ny,ndir] - do 1 freq at a time;
ag_eff      = s1.ag_eff;
atten_dim   = s1.atten_dim + s1.damp_dim;%%treat scattered E and damped E in same way
ICE_MASK    = s1.ICE_MASK;
clear s1;

theta = -pi/180*(90+wavdir);%%waves-to, anti-clockwise, radians

%% weights for integral over directions
%% NB using radians for mwd;
if ndir>1
   wt_theta = ones(ndir,1)*(2*pi/ndir);
else
   wt_theta = 1;
end
S_freq   = zeros(nx,ny);
tau_x    = zeros(nx,ny);
tau_y    = zeros(nx,ny);

chk   = 1;
for i = 1:nx
for j = 1:ny
   %% atten_dim = ENERGY attenuation coeff [m^{-1}]
   if ICE_MASK(i,j)>0
      S_th        = squeeze(S(i,j,:));
      source      = -atten_dim(i,j)*ag_eff(i,j)*S_th;%% m^{-1}*[m/s]*[m^2s] = m^2
      tau_x(i,j)  = -(cos(theta).*wt_theta)'*source;
      tau_y(i,j)  = -(sin(theta).*wt_theta)'*source;
         %% tau_x,tau_y need to be multiplied by rho_wtr*g/phase_vel
         %%  and integrated over frequency as well;
         %% units: [m^2]*[kg/m^3]*[m/s^2]*[s/m]*s^{-1}
         %%         = kg/m/s^2 = Pa

         %% NB take '-' because here we have calc'd the stress
         %% ON the waves (cf Donelan et al, 2012, JGR)
         %% - we want the stress on the ice

      % if chk==1
      %    chk   = 0;
      %    i,j
      %    [theta,cos(theta),sin(theta)]
      %    [S_th,S_th.*cos(theta),S_th.*sin(theta),source]
      %    {tau_x(i,j),tau_y(i,j)}
      %    GEN_pause
      % end

      %%do attenuation
      S(i,j,:) = S(i,j,:)*...
         exp(-atten_dim(i,j)*ag_eff(i,j)*dt);
   end

   %% INTEGRATE SPECTRUM OVER DIRECTION;
   S_freq(i,j) = wt_theta'*squeeze(S(i,j,:));
end
end
