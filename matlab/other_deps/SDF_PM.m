%% SDF_PM.m
%% Author: Timothy Williams
%% Date: 20161209
%
% modiefied from SDF_Pierson_Moskowitz.m by LB (Jan 2021)

function S = SDF_PM(omega,sdf_prams,moment)

omega = abs(omega);
Tp    = sdf_prams{1};
% om_m  = 2*pi/Tm;
Hs    = sdf_prams{2};

g        = 9.81;

if 0 % Tim's original code
 
 if ~exist('moment','var')
  moment   = 0;
 end
 
 %% 1st calc parameters
 alp      = 8.1e-3;
 bet      = .74;
 
 om_0        = 2*pi/Tp/.877;
 U_19pt5     = g/om_0;
 U10         = U_19pt5/1.026;
 Hs_0        = 2/g*U_19pt5^2*sqrt(alp/bet);
 Hs_fac      = Hs/Hs_0;
 
 DO_TEST  = 0;
 if ~exist('omega','var')
  DO_TEST  = 1;
  T        = (25:-1:1)';
  omega    = 2*pi./T;
 end
 omega = abs(omega);
 
 f1 = alp*g^2*omega.^(-5);
 f2 = exp(-bet*(om_0./omega).^4);
 S  = Hs_fac*f1.*f2;
 
 S(omega==0) = 0;
 
 % if DO_TEST
 %    params
 %    S_B   = SDF_Bretschneider(omega,{params.Tp,params.Hs,moment});
 %    plot(2*pi./omega,S);
 %    hold on;
 %    plot(2*pi./omega,S_B,'--r');
 %    hold off;
 % end
 
else                        % Consistent with Aalto tests (Luke)
 
 f  = omega/2/pi;
 fp = 1/Tp;
 
 f1 = ((1/2/pi)^5)*g^2*f.^(-5);
 f2 = exp(-(5/4)*(fp./f).^4);
 
 wt_simp            = 2+0*omega;
 wt_simp([1 end])   = 1;
 wt_simp(2:2:end-1) = 4;
 %%NB om_vec needs to be equally spaced;
 dom    = abs(omega(2)-omega(1));
 wt_om  = dom/3*wt_simp;
 
 S  = f1.*f2;
 S(omega==0) = 0;
 
 mom0 = sum(wt_om.*S);
 
 Hs_0 = 4*sqrt(mom0);
 
 alp_J = (Hs/Hs_0)^2; % Phillip's parameter
 
 S = alp_J*S;
 
end

return
