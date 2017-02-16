%% fn_fill_ice_prams.m
%% Author: Timothy Williams
%% Date: 20141016, 12:00:41 CEST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s1 = fn_fill_iceprams(ice_prams)

if exist('ice_prams','var')
   s1 = ice_prams;%%shorten name for convenience
   clear ice_prams;
end

%% Model Parameters
s1.rhowtr   = 1025;  % Ice density      [kg/m^3]
s1.rhoice   = 922.5; % Ice density      [kg/m^3]
s1.g        = 9.81;  % Gravity          [m/s^2]
s1.poisson  = .3;    % Poisson's ratio

if ~isfield(s1,'drag_rp')
 s1.drag_rp  = 13;    % Robinson-Palmer drag coefficient [Pa/(m/s)]
end
if ~isfield(s1,'visc_ws')
 s1.visc_ws  = 0;    % Wang-Shen viscosity coefficient [m^2/s]
end

%%Brine vol fraction -> Young's modulus and flexural strength
if 1%%just set
   s1.vbf   = .1;          %% brine volume fraction 
   s1.vb    = s1.vbf*1e3;  %% brine volume [ppt]
else%%F&G (1967)
   if 1%%Dany's estimate
      salt  = 5;   % Ice salinity    [psu]
      temp  = -10; % Ice temperature [oC]
   else%%vbf=.1
      temp  = -10;                %[oC]
      salt  = 18.346940647647006; %[psu]
   end

   % Brine volume (Frankenstein and Gardner 1967)
   s1.vb    = salt.*(49.185./abs(temp) + 0.532); % [ppt]
   s1.vbf   = s1.vb/1e3;                         % brine volume fraction 
end

if ~isfield(s1,'young');

   if ~isfield(s1,'young_opt');
      s1.young_opt   = 2;
   end

   if s1.young_opt==0%%just set it
      s1.young = 2e9;      %%lower ~ Marchenko
   elseif s1.young_opt==1%%Vernon's est from vbf
      s1.young = 10e9*(1-3.51*s1.vbf)-1e9;          % Young's modulus [Pa]
   elseif s1.young_opt==2%%just set it
      s1.young = 5.49e9;   %%higher ~ Vernon's guess (s1.vbf=.1)
   end

end

% Flexural strength (Timco and O'Brien 1994)
s1.sigma_c  = 1.76e6.*exp(-5.88.*sqrt(s1.vbf)); % [Pa]


%%Breaking criterion option
if ~isfield(s1,'BRK_OPT');
   s1.BRK_OPT = 0;
end

if s1.BRK_OPT==0|s1.BRK_OPT==1%%beam test
   %%not used if BRK_OPT==0
   s1.strain_c = s1.sigma_c/s1.young;
   s1          = fn_cohesion(s1,'strain');
elseif s1.BRK_OPT==2%% Marchenko's stress criterion
    %% - convert to strain criterion
   s1.stress_c = 2.6*s1.sigma_c;
   s1          = fn_cohesion(s1,'stress');
elseif s1.BRK_OPT==3%% Mohr-Coulomb stress criterion
   %% - convert to strain criterion
   s1 = fn_cohesion(s1,'cohesion');
end

%% flex rigidity = s1.flex_rig_coeff*h^3
s1.flex_rig_coeff = s1.young/12/(1-s1.poisson^2);

% Parameters for the floe size distribution
fsd_params     = set_FSD_params();
s1.Dmax_min    = fsd_params.Dmax_min;  % min floe size [m]
s1.xi          = fsd_params.xi;        % [-]
s1.fragility   = fsd_params.fragility; % [-]
s1.Dthresh     = fsd_params.Dthresh;   % change from power law to uniform FSD here [m]
s1.cice_min    = fsd_params.cice_min;  % min conc where atten happens
