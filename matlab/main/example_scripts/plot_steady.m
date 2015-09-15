%% plot_steady.m
%% Author: Timothy Williams
%% Date: 20150820, 10:12:34 CEST
%% Script to plot W_miz vs T_p

params   = read_infile_matlab();

%% run pure matlab code
params.MEX_OPT          = 0;
params.DO_BREAKING      = 0;
params.duration_hours   = 72;
params.visc_rp          = 0;

%% cancel some saving and plotting inside WIM2d.m
params.SV_BIN     = 1;
params.PLOT_INIT  = 0;
params.PLOT_FINAL = 0;
params.PLOT_PROG  = 0;

%% initial ice and waves
params.Hs_init    = 3;
params.T_init     = 12;
params.dir_init   = -90;
params.h_init     = 2;
params.conc_init  = .7;
params.Dmax_init  = 100;
params.nw         = 1;
params.ndir       = 2^7;

%% set grid here and pass into run_WIM2d.m
grid_prams.x0  = 0;
grid_prams.y0  = 0;
grid_prams.dx  = 1000;
grid_prams.dy  = 1000;
grid_prams.nx  = 300;
grid_prams.ny  = 1;
%%
params.ADV_DIM = 1;%%ny=1 so need 1d advection
params.OPT     = 3;%%no land
grid_prams     = get_grid(grid_prams,params.OPT);
%%        y0: 0
%%        nx: 51
%%        ny: 51
%%        dx: 4000
%%        dy: 4000
%%         X: [51x51 double]
%%         Y: [51x51 double]
%%  LANDMASK: [51x51 double]
%%      scuy: [51x51 double]
%%      scvx: [51x51 double]
%%      scp2: [51x51 double]
%%     scp2i: [51x51 double]
fnames   = fieldnames(grid_prams);
for j=1:length(fnames)
   fname = fnames{j};
   eval([fname,' = grid_prams.',fname,';']);
end

%% set ice
ice_fields.WTR_MASK  = (X<50e3)|(X>200e3);
ice_fields.ICE_MASK  = (1-ice_fields.WTR_MASK).*(1-LANDMASK);%%0 on land & water
%%
ice_fields.cice   = params.conc_init*ice_fields.ICE_MASK;
ice_fields.hice   = params.h_init   *ice_fields.ICE_MASK;
ice_fields.Dmax   = params.Dmax_init*ice_fields.ICE_MASK;
%% ice_fields  = structure eg:
%%      cice: [51x51 double]
%%      hice: [51x51 double]
%%      Dmax: [51x51 double]
%%  WTR_MASK: [51x51 logical]
%%  ICE_MASK: [51x51 double]

%%set waves
wave_fields.WAVE_MASK   = (X<30e3);
wave_fields.Hs          = params.Hs_init *wave_fields.WAVE_MASK;
wave_fields.mwd         = params.dir_init*wave_fields.WAVE_MASK;
wave_fields.Tp          = params.T_init*wave_fields.WAVE_MASK;
%% structure eg:
%%         Hs: [150x10 double]
%%         Tp: [150x10 double]
%%        mwd: [150x10 double]
%%  WAVE_MASK: [150x10 logical]

out_fields,diagnostics  = run_WIM2d(params,grid_prams,ice_fields,wave_fields);
