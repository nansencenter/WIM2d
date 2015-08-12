%% /Users/timill/GITHUB-REPOSITORIES/WIM2d/matlab/numerics/testing/test_AA_isotropic.m
%% Author: Timothy Williams
%% Date: 20141112, 15:14:20 CET

%%test outputs from adv_atten_isotropic.m
load test_fou
oo = ones(ndir,ndir);
id = eye(ndir);
%%
dtheta   = 2*pi/ndir;
wt_theta = dtheta*oo(:,1);
Hs0   = 4*sqrt(sum(S_th.*wt_theta))
S_th

%% old way
cg       = ag_eff(i,j);
source1  = cg*M_bolt*S_th;
tau_x1   = -(cos(theta).*wt_theta)'*source1%% [m^2]
tau_y1   = -(sin(theta).*wt_theta)'*source1%% [m^2]

%%compare to results computed the new way
[tau_x1,tau_x(i,j)]
[tau_y1,tau_y(i,j)]

S_th2,dd

Hs2   = 4*sqrt(sum(wt_theta.*S_th2))
return

%%new way
nvec     = (0:ndir-1)';
Ex       = exp(1i*nvec*theta');
Mft      = Ex*diag(wt_theta);%%\int e^{n i\theta} [] d\theta
S_fou2   = Mft*S_th;
Mift     = Ex'/2/pi;%%S(\theta_m)=(1/2/pi)*\sum_{n=0} e^{-n i\theta_m} S_n

if 0
   [S_fou,2*pi*ifft(S_th),S_fou2]
   if 0
      tstv  = oo(:,1)/2/pi;
   else
      tstv  = exp(-1i*theta)/2/pi;
   end
   [Mft*tstv,2*pi*ifft(tstv)]
end

jp1   = 2;
jm1   = ndir;
S_cos = real( S_fou2(jp1)+S_fou2(jm1) )/2;
S_sin = real((S_fou2(jp1)-S_fou2(jm1))/2i);
[S_cos,sum(wt_theta.*cos(theta).*S_th)]
[S_sin,sum(wt_theta.*sin(theta).*S_th)]

src2_fou = cg*( -q_tot*S_fou2+K_fou.*S_fou2 );
%[source2_fou,-cg*q_tot*S_fou2]
src2_cos = real(  src2_fou(jp1)+src2_fou(jm1) )/2;
src2_sin = real( (src2_fou(jp1)-src2_fou(jm1) )/2i);
[tau_x1,-src2_cos]
[tau_y1,-src2_sin]
