%% adv_atten_timestep_simple.m
%% Author: Timothy Williams
%% Date: 20141018, 18:04:46 CEST

function [S,S_freq,tau_x,tau_y] = ...
   adv_atten_timestep_isotropic(grid_prams,ice_prams,s1,dt)

ndir        = s1.ndir;
wavdir      = s1.wavdir;
S           = s1.Sdir;%% size(S)=[nx,ny,ndir] - do 1 freq at a time;
ag_eff      = s1.ag_eff;
atten_dim   = s1.atten_dim;
damp_dim    = s1.damp_dim;
ICE_MASK    = s1.ICE_MASK;
clear s1;

% ADV_OPT  = 0;%%zeros outside real domain
% ADV_OPT  = 1;%%periodic in x,y
ADV_OPT  = 2;%%periodic in y only

theta = -pi/180*(90+wavdir);%%waves-to, anti-clockwise, radians
% S0 = S;
% {wavdir,theta}

%%advection;
for jth  = 1:ndir
   u           = ag_eff*cos(theta(jth));
   v           = ag_eff*sin(theta(jth));
   S(:,:,jth)  = waveadv_weno(S(:,:,jth),u,v,grid_prams,dt,ADV_OPT);
end
% S0-S
% [min(S0(:)),max(S0(:))]
% [min(S(:)),max(S(:))]
% [min(u(:)),max(u(:))]
% [min(v(:)),max(v(:))]
% GEN_pause

%%attenuation
nx = grid_prams.nx;
ny = grid_prams.ny;
clear grid_prams;

%% weights for integral over directions
%% NB using radians for mwd;
if ndir>1
   dth      = 2*pi/ndir;
   wt_theta = ones(ndir,1)*dth;
else
   wt_theta = 1;
end
S_freq   = zeros(nx,ny);
tau_x    = zeros(nx,ny);
tau_y    = zeros(nx,ny);
%%
oo       = ones(ndir,ndir);
id       = eye(ndir);

for i = 1:nx
for j = 1:ny
   %% atten_dim = ENERGY attenuation coeff [m^{-1}]
   if ICE_MASK(i,j)>0
      %%scattered energy
      q_scat   = atten_dim(i,j);
      %%absorbed energy:
      q_abs    = damp_dim(i,j);
      %%
      %M_bolt   = ( q_scat*(oo/ndir-id)-q_abs*id )/dth;% [m^{-1}]
      M_bolt   = ( q_scat*(oo/ndir-id)-q_abs*id );% [m^{-1}] - dth shouldn't be here
      if 1
         %% get eigenvalues analytically
         %% eigenvalues are -q_scat*[1/dth,1/dth,...,1/dth,0]-q_abs/dth
         %% - should be the same as solving in Fourier space
         %% - TODO: check this
         % dd       = -q_scat./wt_theta;
         % dd(end)  = 0;
         % dd       = dd-q_abs/dth;
         dd       = -q_scat*wt_theta/(2*pi);
         dd(end)  = 0;
         dd       = dd-q_abs;

         %% last eigenvector corresponds to same scattering in all directions
         uu = oo(:,1)/sqrt(ndir);

         %% rest of eigenvectors are in the orthogonal hyper-plane
         %% - can find wih null
         U  = [null(uu'),uu];
      else
         [U,D]    = eig(M_bolt);
         dd       = diag(D);
         %GEN_pause;
      end
      S_th        = squeeze(S(i,j,:));
      source      = ag_eff(i,j)*M_bolt*squeeze(S_th);%% m^{-1}*[m/s]*[m^2s] = m^2
      tau_x(i,j)  = -(cos(theta).*wt_theta)'*source;      %% [m^2]
      tau_y(i,j)  = -(sin(theta).*wt_theta)'*source;      %% [m^2]
         %% tau_x,tau_y need to be multiplied by rho_wtr*g/phase_vel
         %%  and integrated over frequency as well;
         %% units: [m^2]*[kg/m^3]*[m/s^2]*[s/m]*s^{-1}
         %%         = kg/m/s^2 = Pa

         %% NB take '-' because here we have calc'd the stress
         %% ON the waves (cf Donelan et al, 2012, JGR)
         %% - we want the stress on the ice

      cc       = U'*squeeze(S(i,j,:));                %%expand in terms of eigenvectors
      S(i,j,:) = U*diag(exp(dd*ag_eff(i,j)*dt))*cc;   %%exponential attenuation
   end

   %% INTEGRATE SPECTRUM OVER DIRECTION;
   %% TODO: add stress calculation here;
   S_freq(i,j) = wt_theta'*squeeze(S(i,j,:));
end
end

