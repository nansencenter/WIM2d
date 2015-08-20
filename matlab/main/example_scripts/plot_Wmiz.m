%% plot_Wmiz.m
%% Author: Timothy Williams
%% Date: 20150820, 10:12:34 CEST
%% Script to plot W_miz vs T_p

params   = read_infile_matlab();

%% run pure matlab code
params.MEX_OPT          = 0;
params.DO_BREAKING      = 1;
params.duration_hours   = 6;
params.visc_rp          = 13;

%% cancel some saving and plotting inside WIM2d.m
params.SV_BIN     = 0;
params.PLOT_INIT  = 0;
params.PLOT_FINAL = 0;
params.PLOT_PROG  = 0;

%% initial ice and waves
params.Hs_init    = 2;
params.dir_init   = -90;
params.h_init     = 2;
params.conc_init  = .7;
params.Dmax_init  = 300;
params.nw         = 25;
params.ndir       = 16;

%% set grid here and pass into run_WIM2d.m
grid_prams.x0  = 0;
grid_prams.y0  = 0;
grid_prams.dx  = 2000;
grid_prams.dy  = 2000;
grid_prams.nx  = 150;
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
xav                  = mean(X(:));
xm                   = max(X(:)-xav);
ice_fields.WTR_MASK  = ((X-xav)<-.7*xm);
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
wave_fields.WAVE_MASK   = ((X-xav)<-.8*xm);
wave_fields.Hs          = params.Hs_init *wave_fields.WAVE_MASK;
wave_fields.mwd         = params.dir_init*wave_fields.WAVE_MASK;
%% structure eg:
%%         Hs: [150x10 double]
%%         Tp: [150x10 double]
%%        mwd: [150x10 double]
%%  WAVE_MASK: [150x10 logical]
nTp      = 9;
Tp_vec   = linspace(8,12,nTp)';

%% fields to plot
vbls  = {'Dmax','Hs','tau_x'};
lab3  = {'{\itD}_{\rm max}, m','{\itH}_{\rm s}, s','\tau_x, Pa'};
Nv    = length(vbls);
for j=1:Nv
   vbl      = vbls{j};
   V.(vbl)  = zeros(nx,nTp);
end

for n=1:nTp
   params.T_init  = Tp_vec(n)
   wave_fields.Tp = params.T_init*wave_fields.WAVE_MASK;
   out_fields     = run_WIM2d(params,grid_prams,ice_fields,wave_fields);
   for j=1:Nv
      vbl            = vbls{j};
      V.(vbl)(:,n)   = out_fields.(vbl);
   end
end

%%3d plots of variables
for j=1:Nv
   figure(100+j);
   vbl   = vbls{j};
   labs  = {'\itx, \rmkm','\itT_{\rmp}, \rms',lab3{j}};
   fn_pcolor(X/1e3,Tp_vec,V.(vbl),labs);
   fn_fullscreen;
end

%% plot Wmiz
Wmiz  = 0*Tp_vec;
for n=1:nTp
   Dm       = V.Dmax(:,n);
   jMIZ     = find(Dm>0&Dm<300);
   Wmiz(n)  = dx*length(jMIZ)/1e3;
end

figure(201);
fn_fullscreen;
plot(Tp_vec,Wmiz);
GEN_proc_fig('\itT_{\rmp}, \rms','W_{\rmMIZ}, km');
