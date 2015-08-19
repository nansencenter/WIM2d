%% iceinit.m
%% Author: Timothy Williams
%% Date: 20141016, 12:00:41 CEST
function ice_fields = iceinit(ice_prams,grid_prams)
%% ice_prams is struct eg:
%%   conc_init : 0.750000000000000  (conc)
%%   h_init    : 1                  (thickness)
%%   Dmax_init : 1                  (max floe size)
%%   OPT       : 1
%% OPT can be 1,2,3 (different options for initial conditions)
%%
%% grid_prams  = structure, eg:
%%        x0: 0
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
%%

do_test  = 0;
if nargin==0
   do_test     = 1;
   %%
   ice_prams.conc_init  = .75;
   ice_prams.h_init     = 2;
   ice_prams.Dmax_init  = 300;
   ice_prams.OPT        = 1
   %%
   nx = 150;
   ny = 40;
   dx = 2e3;
   %%
   grid_prams  = struct('x0',0,'y0',0,...
                        'nx',nx,'ny',ny,...
                        'dx',dx,'dy',dx);
   grid_prams  = get_grid(grid_prams,ice_prams.OPT)
end
OPT   = ice_prams.OPT;

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
   %% column of land on right
   xav      = mean(X(:));
   xm       = max(X(:)-xav);
   WTR_MASK = ((X-xav)<-.7*xm);
   jW       = find(WTR_MASK==1);
   ICE_MASK = (1-WTR_MASK).*(1-LANDMASK);%%0 on land & water
   jI       = find( ICE_MASK==1 );

elseif OPT==2
   %%ice in lower-left corner
   xav      = mean(X(:));
   yav      = mean(Y(:));
   WTR_MASK = ((X-xav)<0)|((Y-yav)>0);
   jW       = find(WTR_MASK==1);
   ICE_MASK = (1-WTR_MASK).*(1-LANDMASK);%%0 on land & water
   jI       = find( ICE_MASK==1 );

elseif OPT==3
   %% ice strip in middle
   xav         = mean(X(:));
   xmin        = min(X(:));
   xmax        = max(X(:));
   xm          = .5*(dx+xmax-xmin);
   xe          = xav-.7*xm;
   strip_width = 100.0e3;
   xe2         = xe+strip_width;
   %%
   ICE_MASK = (X>xe)&(X<xe2);
   WTR_MASK = (1-ICE_MASK).*(1-LANDMASK);%%0 on land & ice
   jW       = find(WTR_MASK==1);
   jI       = find(ICE_MASK==1);
end

cice(jI) = ice_prams.conc_init;%% ice
cice(jL) = NaN;%% land
%%
Dmax(jI) = ice_prams.Dmax_init;
Dmax(jL) = NaN;
%%
hice(jI) = ice_prams.h_init;
hice(jL) = NaN;

%% outputs:
ice_fields  = struct('cice'      ,cice    ,...
                     'hice'      ,hice    ,...
                     'Dmax'      ,Dmax    ,...
                     'WTR_MASK'  ,WTR_MASK,...
                     'ICE_MASK'  ,ICE_MASK);

if do_test
   fn_plot_ice(grid_prams,ice_fields);
   fn_fullscreen;
end
