function simul_out   = initialise_wim2sim(simul_out,nextsim_timestep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% integer parameters
SCATMOD     = 1;%%scattering model: 0, no directional scattering; 1, directional scattering;
ADV_DIM     = 2;%%2: 2d advection; 1,1d advection (keep as 2 in most cases);
CHECK_FINAL = 1;%%1/0: do/don't dump binary files after call to WIM;
CHECK_PROG  = 1;%%1/0: do/don't dump binary files regularly during call to WIM;
CHECK_INIT  = 1;%%1/0: do/don't dump binary files after entering WIM;
DO_BREAKING = 1;%%1/0: do/don't do breaking (change Dmax);
%%
simul_out.wim.int_prams = [SCATMOD,ADV_DIM,...
                           CHECK_FINAL,CHECK_PROG,CHECK_INIT,...
                           DO_BREAKING];

if CHECK_FINAL+CHECK_PROG+CHECK_INIT>0
   %%make directories for binaries
   dirs  = {'out_io',
            'out_io/log',
            'out_io/binaries',
            'out_io/binaries/prog'};

   for j=1:length(dirs)
      dd = dirs{j};
      if ~exist(dd)
         mkdir(dd);
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% real parameters
young          = 5.49e9;   %% Young's modulus [Pa]
visc_rp        = 13;       %% Robinson-Palmer damping parameter [Pa/(m/s)]

%% coupling frequency between WIM and neXtSIM
if 0
  %%couple every 2nd time-step
  duration  = 2*nextsim_timestep;
else
  duration_hours = 3;                   %% length of each call to WIM [h]
  duration       = duration_hours*60*60; %% convert duration to seconds
end

simul_out.wim.real_prams         = [young,visc_rp,duration];
simul_out.wim.coupling_freq      = duration;          %%useful to have this explicitly
simul_out.wim.nextsim_timestep   = nextsim_timestep;  %% useful to have this explicitly
                                                      %% (kept in simul_in,not simul_out)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% other parameters
simul_out.wim.other_prams  = struct('Dmax_pack',300.,... %Dmax when ice is unbroken
                                    'cice_min',0.2);     %conc criteria for no ice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get WIM's grid (needs to be regular at the moment)
%% NB WIM2d mex function needs to be compiled with this grid
WIM2d_path  = getenv('WIM2D_PATH');
if strcmp(WIM2d_path,'')==0
   gdir  = [WIM2d_path,'/fortran/run/inputs/']%%directory with grid files
else
   error('$WIM2D_PATH environment not set: add to ~/.bash_profile or ~/.bashrc and launch matlab from terminal')
end
simul_out.wim.gridprams = fn_get_grid(gdir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% want to call WIM at first time step,
%% so set last_call to NaN
simul_out.wim.last_call = NaN;
simul_out.wim.INIT_DMAX = 1;     %need to initialise Dmax (no restart available)
