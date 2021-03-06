%% atten_isotropic.m
%% Author: Timothy Williams
%% Date: 20141018, 18:04:46 CEST

function [S,S_freq,tau_x,tau_y] = ...
   atten_isotropic(grid_prams,ice_prams,s1,dt)

nx = grid_prams.nx;
ny = grid_prams.ny;
clear grid_prams;

ndir        = s1.ndir;
wavdir      = s1.wavdir;
S           = s1.Sdir;%% size(S)=[nx,ny,ndir] - do 1 freq at a time;
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
   S_freq0  = dth*sum(S,3);%%init E for testing
else
   wt_theta = 1;
   S_freq0  = S;
end
S_freq   = zeros(nx,ny);
tau_x    = zeros(nx,ny);
tau_y    = zeros(nx,ny);
%%
oo = ones(ndir,ndir);
id = eye(ndir);

%% method of solving \partial_t E=S (2nd part of split step)
%% NB E is variance spectrum [m^2/s], S are sources [m^2]
TIME_SOLN_METHOD  = 1;  %% 1: Fourier method (generalises to non-isotropic case)
                        %% 2: Position space, exact calculation of evals/evecs (specific to isotropic case)
                        %% 3: Position space, numerical calculation of evals/evecs
                        %%    (generalises to non-isotropic case and non-linear scattering)

CHECK_STOP  = 0;%% stop Fourier method run and save results for checking

for i = 1:nx
for j = 1:ny
   %% atten_dim = ENERGY attenuation coeff [m^{-1}]
   if ICE_MASK(i,j)>0
      %%scattering cross-section (m^{-1}):
      q_scat   = atten_dim(i,j);

      %%total absorbed energy (m^{-1}):
      q_tot    = q_scat + damp_dim(i,j);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if TIME_SOLN_METHOD==1
         %%use Fourier method

         %%fourier coeff's of Boltzmann kernel
         K_fou = q_scat*eye(ndir,1);%%K(theta)=q_scat/2/pi: isotropic scattering

         %%fourier coeff's of directional spectrum
         S_th     = squeeze(S(i,j,:));
         %S_fou    = 2*pi*ifft(S_th);%% - for some reason this is different

         nvec  = (0:ndir-1)';
         Ex    = exp(1i*nvec*theta');
         Mft   = Ex*diag(wt_theta);

         %%S_n = \int e^{n i\theta} S(\theta) d\theta
         S_fou = Mft*S_th;

         %% get eigenvalues to do attenuation,
         %% & transform back to get attenuated directional spectrum
         evals_x  = min(0,K_fou-q_tot);%%stop it >0 due to subtraction errors if q_abs=0
         cg       = ag_eff(i,j);
         Mift     = (1/2/pi)*Ex';%%S(\theta_m)=(1/2/pi)*\sum_{n=0} e^{n i\theta_m} S_n
         S(i,j,:) = real( Mift*(S_fou.*exp(evals_x*cg*dt)) );

         %%stresses
         jp1         = 2;
         jm1         = ndir;
         src_p1      = cg*( -q_tot*S_fou(jp1) + K_fou(jp1)*S_fou(jp1) );%%n=+1 coefficient of source term
         src_m1      = cg*( -q_tot*S_fou(jm1) + K_fou(jm1)*S_fou(jm1) );%%n=-1 coefficient of source term
         src_cos     = real(   src_p1+src_m1 )/2;
         src_sin     = real( ( src_p1-src_m1 )/2i );
         tau_x(i,j)  = -src_cos;%%-\int[src*cos(\theta)]d\theta
         tau_y(i,j)  = -src_sin;%%-\int[src*sin(\theta)]d\theta
            %% tau_x,tau_y need to be multiplied by rho_wtr*g/phase_vel
            %%  and integrated over frequency as well;
            %% units: [m^2]*[kg/m^3]*[m/s^2]*[s/m]*s^{-1}
            %%         = kg/m/s^2 = Pa

            %% NB take '-' because here we have calc'd the stress
            %% ON the waves (cf Donelan et al, 2012, JGR)
            %% - we want the stress on the ice

         if CHECK_STOP
            %% Give a test and exit
            %% (also save results to ../numerics/testing/test_fou.mat
            %%  and test with testing/test_AA_isotropic.m)
            if abs(S_fou(1))>1e-4
               % %%"inputs"
               % theta
               % S_th
               % S_fou

               S_pre       = S_th;
               S_post      = squeeze(S(i,j,:));
               %%
               ofil  = '../numerics/testing/test_fou.mat';
               save(ofil,'S_post','dt','evals_x','S_pre','S_fou',...
                    'K_fou','q_scat','q_tot','ndir','ag_eff','i','j',...
                    'tau_x','tau_y','theta');
               error(['Saving test results (to ',ofil,') and stopping']);
            end
         end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      else
         M_bolt   = ( q_scat*oo/ndir-q_tot*id );% [m^{-1}]
         if TIME_SOLN_METHOD==2
            %% get eigenvalues by inspection
            %% eigenvalues are -q_scat*[1/dth,1/dth,...,1/dth,0]-q_abs/dth
            %% - should be the same as solving in Fourier space
            %% - TODO: check this
            %%   (have same eigenvalues, & 1st eigenvector is the same)
            evals_x    = -q_tot*oo(:,1);
            evals_x(1) = evals_x(1)+q_scat;

            %% 1st eigenvector corresponds to same scattering in all directions
            uu = oo(:,1)/sqrt(ndir);

            %% rest of eigenvectors are in the orthogonal hyper-plane
            %% - can find with null
            U  = [uu,null(uu')];
         else
            [U,D]    = eig(M_bolt);
            evals_x  = diag(D);
         end
         S_th        = squeeze(S(i,j,:));
         cg          = ag_eff(i,j);
         source      = cg*M_bolt*S_th;%% m^{-1}*[m/s]*[m^2s] = m^2
         tau_x(i,j)  = -(cos(theta).*wt_theta)'*source;      %% [m^2]
         tau_y(i,j)  = -(sin(theta).*wt_theta)'*source;      %% [m^2]
            %% tau_x,tau_y need to be multiplied by rho_wtr*g/phase_vel
            %%  and integrated over frequency as well;
            %% units: [m^2]*[kg/m^3]*[m/s^2]*[s/m]*s^{-1}
            %%         = kg/m/s^2 = Pa

            %% NB take '-' because here we have calc'd the stress
            %% ON the waves (cf Donelan et al, 2012, JGR)
            %% - we want the stress on the ice

         cc       = U'*S_th;                       %% expand in terms of eigenvectors
         S(i,j,:) = U*diag(exp(evals_x*cg*dt))*cc; %% do exponential attenuation then transform back
      end%%choice of method
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   end%%ice present

   %% INTEGRATE SPECTRUM OVER DIRECTION;
   %% (need this even if only water)
   S_freq(i,j) = wt_theta'*squeeze(S(i,j,:));
end%j
end%i

%E_test   = [sum(S_freq(:)),sum(S_freq0(:))]/nx/ny
