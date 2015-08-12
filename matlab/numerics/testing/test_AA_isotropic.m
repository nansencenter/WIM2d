%% /Users/timill/GITHUB-REPOSITORIES/WIM2d/matlab/numerics/testing/test_AA_isotropic.m
%% Author: Timothy Williams
%% Date: 20141112, 15:14:20 CET

%%test outputs from adv_atten_isotropic.m
Tfil  = load('test_fou.mat')
%     tau_x: [600x1 double]
%     tau_y: [600x1 double]
%      ndir: 16
%    ag_eff: [600x1 double]
%     theta: [16x1 double]
%        dt: 298.0
%         i: 91
%         j: 1
%    q_scat: 5.206834290781534e-04
%     q_tot: 5.206834290781534e-04
%     K_fou: [16x1 double]
%     S_pre: [16x1 double]
%    S_post: [16x1 double]
%   evals_x: [16x1 double]
%    M_bolt: [16x16 double]
tau_x    = Tfil.tau_x;
tau_y    = Tfil.tau_y;
ndir     = Tfil.ndir;
ag_eff   = Tfil.ag_eff;
theta    = Tfil.theta;
dt       = Tfil.dt;
i        = Tfil.i;
j        = Tfil.j;
q_scat   = Tfil.q_scat;
q_tot    = Tfil.q_tot;
K_fou    = Tfil.K_fou;
S_fou    = Tfil.S_fou;
evals_x  = Tfil.evals_x;
%%
S_pre    = Tfil.S_pre;  %%spec before scattering
S_post   = Tfil.S_post; %%spec after scattering

oo = ones(ndir,ndir);
id = eye(ndir);
%%
dtheta   = 2*pi/ndir;
wt_theta = dtheta*oo(:,1);
cg       = ag_eff(i,j);

%% source term matrix in position space
M_bolt   = ( q_scat*oo/ndir-q_tot*id );% [m^{-1}]
[U,DD]   = eig(M_bolt);
dd       = diag(DD);

if 1%%compare tau_x,tau_y
   %% position space results
   source1  = cg*M_bolt*S_pre;
   tau_x1   = -(cos(theta).*wt_theta)'*source1;%% [m^2]
   tau_y1   = -(sin(theta).*wt_theta)'*source1;%% [m^2]

   %% compare to Fourier space results
   taux_tst_ij = [tau_x1,tau_x(i,j)]
   tauy_tst_ij = [tau_y1,tau_y(i,j)]
end

if 0
   %%compare eigenvalues
   evals_tst   = [evals_x,dd];
else
   cc          = U'*S_pre;                         %% expand in terms of eigenvectors
   tst_atten   = [S_post,U*diag(exp(dd*cg*dt))*cc] %% do exponential attenuation, then go back to position space
end

if 1
   %% is energy conserved?
   Hs_pre      = 4*sqrt(sum(wt_theta.*S_pre ));
   Hs_post     = 4*sqrt(sum(wt_theta.*S_post));
   Econ_tst_Hs = [Hs_pre,Hs_post]
end
