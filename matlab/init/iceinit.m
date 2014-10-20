%% iceinit.m
%% Author: Timothy Williams
%% Date: 20141016, 12:00:41 CEST
function [ice_fields,ice_prams] = iceinit(ice_prams,grid_prams,OPT)
%% ice_prams is struct eg:
%%   c: 0.750000000000000  (conc)
%%   h: 1                  (thickness)
%% optional fields 'young'  (Young's modulus)
%%             and 'bc_opt' [breaking criteria = 0(beam) or 1(Marchenko)]
%% grid_prams is struct eg:
%%   nx: 15                (size of grid in x dirn)
%%   ny: 15                (size of grid in y dirn)
%% opts is struct eg:
%%   ICE_SHAPE: 1          (0/1/2: different options for initial conditions)

do_test  = 0;
if nargin==0
   do_test     = 1;
   h           = 2;
   c           = 0.75;
   ice_prams   = struct('c'         ,c,...
                        'h'         ,h,...
                        'bc_opt'    ,0,...
                        'young_opt' ,1)
   %%
   nx          = 51;
   dx          = 4e3;
   %%
   OPT         = 1;
   grid_prams  = struct('nx',nx,'ny',nx,...
                        'dx',dx,'dy',dx);
   grid_prams  = get_grid(grid_prams,OPT)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get rest of ice_prams
ice_prams   = fn_fill_iceprams(ice_prams);
%% structure eg:
%%               c: 0.750000000000000
%%               h: 2
%%           young: 2.000000000000000e+09
%%          bc_opt: 0
%%          rhoice: 9.225000000000000e+02
%%               g: 9.810000000000000
%%         poisson: 0.300000000000000
%%             vbf: 0.100000000000000
%%              vb: 100
%%         sigma_c: 2.741429878818372e+05
%%        strain_c: 1.370714939409186e-04
%%  flex_rig_coeff: 1.831501831501831e+08
%%            Dmin: 20
%%              xi: 2
%%       fragility: 0.900000000000000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cice/hice/Dmax in ice_fields
nx       = grid_prams.nx;
ny       = grid_prams.ny;
dx       = grid_prams.dx;
dy       = grid_prams.dy;
LANDMASK = grid_prams.LANDMASK;
X        = grid_prams.X;
Y        = grid_prams.Y;

cice  = zeros(nx,ny);
hice  = zeros(nx,ny);
Dmax  = zeros(nx,ny);
jL    = find(LANDMASK);

if OPT==0
   %%island in corner, with some ice around it
   x0    = mean(X(jL));%%centre
   y0    = mean(Y(jL));%%centre
   R0    = 70e3;       %%radius
   Rsq   = (X-x0).^2+(Y-y0).^2;
   %%
   WTR_MASK = ( Rsq>=R0^2 );%%outside a circle
   jW       = find(WTR_MASK==1);
   ICE_MASK = (1-WTR_MASK).*(1-LANDMASK);%%0 on land & water
   jI       = find( ICE_MASK==1 );

elseif OPT==1
   %% south-east corner land surrounded by ice;
   xm       = max(X(:));
   WTR_MASK = (X<-.4*xm);
   jW       = find(WTR_MASK==1);
   ICE_MASK = (1-WTR_MASK).*(1-LANDMASK);%%0 on land & water
   jI       = find( ICE_MASK==1 );

elseif OPT==2
   %%TEST_1d
   %%land to right,water to left
   WTR_MASK = (X<0)|(Y>0);
   jW       = find(WTR_MASK==1);
   ICE_MASK = (1-WTR_MASK).*(1-LANDMASK);%%0 on land & water
   jI       = find( ICE_MASK==1 );

end

cice(jI) = ice_prams.c;%% ice
cice(jL) = NaN;%% land
%%
Dmax(jI) = 250;
Dmax(jL) = NaN;
%%
hice(jI) = ice_prams.h;
hice(jL) = NaN;

%% outputs:
ice_fields  = struct('cice',cice,...
                     'hice',hice,...
                     'Dmax',Dmax,...
                     'WTR_MASK',WTR_MASK,...
                     'ICE_MASK',ICE_MASK);

if do_test
   fn_plot_ice(grid_prams,ice_fields);
end
