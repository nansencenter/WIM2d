%% adv_atten_timestep_simple.m
%% Author: Luke Bennetts
%% Date: 20150903

function [S,S_freq,tau_x,tau_y] = ...
   adv_atten_timestep_simple_conserve(grid_prams,ice_prams,s1,dt,adv_options)

nx = grid_prams.nx;
ny = grid_prams.ny;

ndir        = s1.ndir;
wavdir      = s1.wavdir;
S           = s1.Sdir;%% size(S)=[nx,ny,ndir] - do 1 freq at a time;
ag_eff      = s1.ag_eff;
%%atten_dim   = s1.atten_dim + s1.damp_dim;%%treat scattered E and damped E in same way
atten_dim   = s1.atten_dim;%%for scattering
damp_dim    = s1.damp_dim; %%for damping
ICE_MASK    = s1.ICE_MASK;
clear s1;

if mod(ndir,2)
 disp('>>> adv_atten_timestep_simple_conserve requires ndir to be even')
 return
end

% ADV_OPT  = 0;%%zeros outside real domain
% ADV_OPT  = 1;%%periodic in x,y
% ADV_OPT  = 2;%%periodic in y only

theta = -pi/180*(90+wavdir);%%waves-to, anti-clockwise, radians
% S0 = S;
% {wavdir,theta}

%%advection;
if adv_options.ADV_DIM==2
   %%2d advection
   for jth  = 1:ndir
      u           = ag_eff*cos(theta(jth));
      v           = ag_eff*sin(theta(jth));
      S(:,:,jth)  = waveadv_weno(S(:,:,jth),u,v,grid_prams,dt,adv_options);
   end
else
   %%1d advection - 1 row at a time
   for jy=1:ny
      for jth  = 1:ndir
         u           = ag_eff*cos(theta(jth));
         S(:,jy,jth) = waveadv_weno_1d(S(:,jy,jth),u,grid_prams,dt,adv_options);
      end
   end
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
   wt_theta = ones(ndir,1)*(2*pi/ndir);
else
   wt_theta = 1;
end
S_freq   = zeros(nx,ny);
tau_x    = zeros(nx,ny);
tau_y    = zeros(nx,ny);

chk   = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For attenuation
%%% - all the lost energy gets added to reverse dirn
%%% - for ndir=2, the kernel matrix is [0,1;1,0],
%%%      which has eigenvectors and eigen values below
evec0 = [1;1]; % corresponds to eval=0
evec1 = [1;-1]; % corresponds to eval=-2*alpha
Mat0 = [evec0,evec1];
inv_Mat0 = 0.5*Mat0;

%% do everything at once for source calc
M_bolt0  = zeros(ndir,ndir);
for lp=1:ndir/2
   M_bolt0(lp,lp+ndir/2)   = 1;
   M_bolt0(lp+ndir/2,lp)   = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nx
for j = 1:ny
   %% atten_dim = ENERGY attenuation coeff [m^{-1}]
   if ICE_MASK(i,j)>0
      S_th        = squeeze(S(i,j,:));
      cg          = ag_eff(i,j);
      qtot        = cg*(atten_dim(i,j)+damp_dim(i,j));
      source      = -qtot*S_th;
         %% this part of source should be both scat & atten
         %% units: m^{-1}*[m/s]*[m^2s] = m^2

      qscat       = cg*atten_dim(i,j);
      source      = source+qscat*M_bolt0*S_th;%% m^{-1}*[m/s]*[m^2s] = m^2
         %% TODO: this part of source should be just scat
         %% units: m^{-1}*[m/s]*[m^2s] = m^2
         %% NB overall transfer matrix is M_bolt=-qtot*eye(ndir)+qscat*M_bolt0

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
      for lp=1:ndir/2
       dumS = squeeze(S(i,j,[lp,lp+ndir/2]));
       c_vec = inv_Mat0*dumS;%%expand in terms of eigen vectors
       dumS = c_vec(1)*evec0 + ...
        c_vec(2)*exp(-2*atten_dim(i,j)*cg*dt);%%do scattering
       %%S(i,j,[lp,lp+ndir/2]) = dumS;
       S(i,j,[lp,lp+ndir/2]) = dumS*exp(-damp_dim(i,j)*cg*dt);
         %% do damping
         %% -both get damped equally by damp_dim
      end
      clear dumS c_vec
      
   end

   %% INTEGRATE SPECTRUM OVER DIRECTION;
   S_freq(i,j) = wt_theta'*squeeze(S(i,j,:));
end
end

