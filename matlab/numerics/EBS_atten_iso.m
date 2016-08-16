%% EBS_atten_iso.m
%% Author: Timothy Williams
%% Date: 20160812

function [S,S_scattered,S_freq,tau_x,tau_y] = ...
   EBS_atten_iso(grid_prams,ice_prams,s1,dt)

nx = grid_prams.nx;
ny = grid_prams.ny;
clear grid_prams;

ndir        = s1.ndir;
wavdir      = s1.wavdir;
S           = s1.Sdir;%% size(S)=[nx,ny,ndir] - do 1 freq at a time;
S_scattered = s1.Sdir_scattered;%% size(S)=[nx,ny,ndir] - do 1 freq at a time;
ag_eff      = s1.ag_eff;
atten_dim   = s1.atten_dim;
damp_dim    = s1.damp_dim;
ICE_MASK    = s1.ICE_MASK;
clear s1;

theta = -pi/180*(90+wavdir);%%waves-to, anti-clockwise, radians

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
oo    = ones(ndir,ndir);
zz    = zeros(ndir,ndir);
id    = eye(ndir);
id2   = eye(2*ndir);

if 1
   %% choose filter so back-scattered waves
   M_filter = EBS_filter_step(wavdir);
else
   %%for testing - this should give the same answer as SCATMOD==0
   M_filter = zz;
   disp('warning: filter turned off for testing')
end
M_K      = oo/ndir;%% q_scat==1 - can scale eigenvalues & M_bolt inside step_EBS
M11      =  -id+(1-M_filter).*M_K;%% E that stays in "normal" spectrum
M21      = M_filter.*M_K;%% E that is moved to "already-scattered" spectrum
M_bolt   = [M11,zz;M21,zz];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U,D]          = eig(M_bolt);
inputs.e_vals  = diag(D);
inputs.e_vecs  = U;
inputs.M_bolt  = M_bolt;

e_tol = 1e-12;
if ~isempty(find(inputs.e_vals>e_tol))
   disp('Eigen-values:')
   disp(inputs.e_vals)
   error('\nShouldn''t be any positive eigenvalues')
else
   inputs.e_vals   = min(inputs.e_vals,0);
end

ev    = inputs.e_vals(inputs.e_vals<-e_tol);
ev    = max(ev);%%slowest decaying
e_fac = 1;
if 1
   %%scale matrix/evals by e_fac to give correct max e-val of -1
   e_fac = 1./abs(ev);
end
inputs.e_vals  = e_fac*inputs.e_vals;
inputs.M_bolt  = e_fac*inputs.M_bolt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:nx
for j = 1:ny
   %% atten_dim = ENERGY attenuation coeff [m^{-1}]
   if ICE_MASK(i,j)>0

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%scattering cross-section (m^{-1}):
      q_scat   = atten_dim(i,j);

      %%total absorbed energy (m^{-1}):
      q_dis    = damp_dim(i,j);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      cg    = ag_eff(i,j);
      j_E   = (1:ndir)';
      j_S   = ndir+(1:ndir)';
      %%
      S2 = [squeeze(S(i,j,:));squeeze(S_scattered(i,j,:))];

      %% source has units m^{-1}*[m/s]*[m^2s] = m^2
      %% 1st line comes from D_t(E);
      %% 2nd line comes from D_t(E_scattered);
      M_bolt_ij   = q_scat*M_bolt-q_dis*id2;%%proper Boltzmann matrix (scaled & with dissipation) 
      source      = cg*M_bolt_ij(j_E,:)*S2+...
                   +cg*M_bolt_ij(j_S,:)*S2;

      %% ==========================================================
      %% radiation stresses (momentum flux)
      %% NB:
      %% tau_x,tau_y need to be multiplied by rho_wtr*g/phase_vel
      %%  and integrated over frequency as well;
      %% units: [m^2]*[kg/m^3]*[m/s^2]*[s/m]*s^{-1}
      %%         = kg/m/s^2 = Pa

      tau_x(i,j)     = -(cos(theta).*wt_theta)'*source;      %% [m^2]
      tau_y(i,j)     = -(sin(theta).*wt_theta)'*source;      %% [m^2]
         %% NB take '-' because here we have calc'd the stress
         %% ON the waves (cf Donelan et al, 2012, JGR)
         %% - we want the stress on the ice
      %% ==========================================================


      %% ==========================================================
      %%step fwd & apply dissipation afterwards
      Dx = cg*dt;
      S2 = EBS_step(S2,q_scat*Dx,inputs)*exp(-q_dis*Dx);
      %% ==========================================================

      S(i,j,:)             = S2(j_E);
      S_scattered(i,j,:)   = S2(j_S);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   end%%ice present

   %% INTEGRATE SPECTRUM OVER DIRECTION;
   %% (need this even if only water)
   %% - this is to get H_s and other integrals to calculate breaking
   S2_         = S(i,j,:)+S_scattered(i,j,:);%% add the 2 energies together
   S_freq(i,j) = wt_theta'*squeeze(S2_);
end
end
