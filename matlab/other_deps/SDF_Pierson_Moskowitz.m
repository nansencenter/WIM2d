%% SDF_Pierson_Moskowitz.m
%% Author: Timothy Williams
%% Date: 20161209

function [params,S]  = SDF_Pierson_Moskowitz(params,omega,moment)

if ~exist('params','var')
   params.U10  = 10;
end
if ~exist('moment','var')
   moment   = 0;
end

%% 1st calc parameters
alp      = 8.1e-3;
bet      = .74;
g        = 9.81;
if isfield(params,'Hs')
   U_19pt5     = (g^2*params.Hs^2*bet/(4*alp))^.25;
   params.U10  = U_19pt5/1.026;
   om_0        = g/U_19pt5;
   params.Tp   = 2*pi/(.877*om_0);%%peak period
elseif isfield(params,'U10')
   U_19pt5     = 1.026*params.U10;
   params.Hs   = 2/g*U_19pt5^2*sqrt(alp/bet);
   om_0        = g/U_19pt5;
   params.Tp   = 2*pi/(.877*om_0);%%peak period
elseif isfield(params,'Tp')
   om_0        = 2*pi/params.Tp/.877;
   U_19pt5     = g/om_0;
   params.U10  = U_19pt5/1.026;
   params.Hs   = 2/g*U_19pt5^2*sqrt(alp/bet);
end

%% only calc spectrum if needed
if nargout==1
   return%comment to do testing
end

DO_TEST  = 0;
if ~exist('omega','var')
   DO_TEST  = 1;
   T        = (25:-1:1)';
   omega    = 2*pi./T;
end
omega = abs(omega);

f1 = alp*g^2*omega.^(-5);
f2 = exp(-bet*(om_0./omega).^4);
S  = f1.*f2;

S(omega==0) = 0;

if DO_TEST
   params
   S_B   = SDF_Bretschneider(omega,{params.Tp,params.Hs,moment});
   plot(2*pi./omega,S);
   hold on;
   plot(2*pi./omega,S_B,'--r');
   hold off;
end
