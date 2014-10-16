%% iceinit.m
%% Author: Timothy Williams
%% Date: 20141016, 12:00:41 CEST
function ice_fields = iceinit(ice_prams,grid_prams,OPT)
%% ice_prams is struct eg:
%%   c: 0.750000000000000  (conc)
%%   h: 1                  (thickness)
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
   ice_prams   = struct('c',c,'h',h)
   %%
   nx          = 51;
   dx          = 4e3;
   %%
   OPT         = 1;
   grid_prams  = struct('nx',nx,'ny',nx,...
                        'dx',dx,'dy',dx);
   grid_prams  = get_grid(grid_prams,OPT)
end

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

   if 0
      pcolor((1-LANDMASK).*Rsq);colorbar;daspect([1 1 1]);pause
   elseif 0
      pcolor(LANDMASK*1.0);colorbar;daspect([1 1 1]);pause
   elseif 0
      pcolor(WTR_MASK+LANDMASK);colorbar;daspect([1 1 1]);pause
   end

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
cice(jL) = -.5;%% land
%%
Dmax(jI) = 500;
Dmax(jL) = -250;
%%
hice(jI) = ice_prams.h;
hice(jL) = -1;

%% outputs:
ice_fields  = struct('cice',cice,...
                     'hice',hice,...
                     'Dmax',Dmax,...
                     'WTR_MASK',WTR_MASK,...
                     'ICE_MASK',ICE_MASK);

if do_test
   subplot(1,3,1);
   %GEN_plot_matrix(X/1e3,Y/1e3,cice,[-1 1]);
   H  = pcolor(X/1e3,Y/1e3,cice);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   colorbar;
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   ttl   = title('Concentration');
   GEN_font(ttl);
   colorbar;
   %%
   subplot(1,3,3);
   %GEN_plot_matrix(X/1e3,Y/1e3,Dmax,[-1 500]);
   H  = pcolor(X/1e3,Y/1e3,Dmax);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   ttl   = title('\itD_{\rmmax}');
   GEN_font(ttl);
   colorbar;
   %%
   subplot(1,3,2);
   %GEN_plot_matrix(X/1e3,Y/1e3,S(:,:,jpp,jmwd));
   H  = pcolor(X/1e3,Y/1e3,hice);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   colorbar;
   ttl   = title('\ith, \rmm');
   GEN_font(ttl);
   %%
end
