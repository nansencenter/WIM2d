function s1 = fn_cohesion(s1,INPUT_TYPE)
%% call: ice_prams = fn_cohesion(ice_prams,INPUT_TYPE)
%% Author: Timothy Williams
%% Date: 20161111
%% convert between cohesion,critical strain and critical stress
%% - makes strain_c, stress_c and cohesion consistent with each other
%%
%% INPUTS:
%% if INPUT_TYPE = 'strain':
%%    ice_prams = struct eg:
%%      poisson: 0.300000000000000
%%     friction: 0.700000000000000
%%        young: 5.000000000000000e+09
%%     strain_c: 5.000000000000000e-05
%%
%% if INPUT_TYPE = 'stress':
%%    ice_prams = struct eg:
%%      poisson: 0.300000000000000
%%     friction: 0.700000000000000
%%        young: 5.000000000000000e+09
%%      stress_c: 2.747252747252747e+05
%%
%% if INPUT_TYPE = 'cohesion':
%%    ice_prams = struct eg:
%%      poisson: 0.300000000000000
%%     friction: 0.700000000000000
%%        young: 5.000000000000000e+09
%%      cohesion: 2.458791208791209e+05
%%
%%
%% OUTPUTS:
%% ice_prams = struct eg:
%%       poisson: 0.300000000000000
%%      friction: 0.700000000000000
%%      strain_c: 5.000000000000000e-05
%%         young: 5.000000000000000e+09
%%      stress_c: 2.747252747252747e+05
%%      cohesion: 2.458791208791209e+05
%%  ref_cohesion: 1100000
%%     ref_scale: 1.000000000000000e-03
%%         scale: 0.020014376580007

alpha = (1-s1.poisson)/(1+s1.poisson);

switch INPUT_TYPE
case 'strain'
   s1.stress_c = s1.young/(1-s1.poisson^2)*s1.strain_c; %%convert to plate breaking stress
      %%=E/(1-nu^2)*strain_c = thin plate (plane stress: \sigma_33=0)
   sig_N       = (1-s1.poisson)/2*s1.stress_c;
   s1.cohesion = (alpha+s1.friction)*sig_N;

case 'stress'
   s1.strain_c = (1-s1.poisson^2)*s1.stress_c/s1.young;
   sig_N       = (1-s1.poisson)/2*s1.stress_c;
   s1.cohesion = (alpha+s1.friction)*sig_N;

case 'cohesion'
   sig_N       = s1.cohesion/(alpha+s1.friction);
      %compressive stress  = .5*(sig_11+sig_22) = .5*(1-poisson)*sig11 at breaking point
   s1.stress_c = 2*sig_N/(1-s1.poisson);
      %%=E/(1-nu^2)*strain_c = thin plate (plane stress: \sigma_33=0)
   s1.strain_c = (1-s1.poisson^2)*s1.stress_c/s1.young;

otherwise
   disp('INPUT_TYPE:');
   disp(INPUT_TYPE);
   error('Unknown option for INPUT_TYPE');
end

%%make reference cohesion the lab scale (1.1MPa,1mm)
s1.ref_cohesion   = 1.1e6;
s1.ref_scale      = 1.0e-3;
s1.scale          = s1.ref_scale*(s1.ref_cohesion/s1.cohesion)^2;
