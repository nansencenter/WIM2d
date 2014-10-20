function grid_prams = get_grid(grid_prams,OPT)
%% get_grid.m
%% Author: Timothy Williams
%% Date:   20140826, 03:08:57 CEST

do_test  = 0;
if nargin==0
   do_test  = 1;
   OPT      = 1;
   %%
   ii = 51;
   jj = 51;
   dx = 4e3;%m
   dy = 4e3;%m
   %%
   grid_prams  = struct('nx',ii,'ny',jj,...
                        'dx',dx,'dy',dy);
end
ii = grid_prams.nx;
jj = grid_prams.ny;
dx = grid_prams.dx;
dy = grid_prams.dy;

xm = (ii+1)/2*dx;%%x\in[-xm,xm]
ym = (jj+1)/2*dy;%%y\in[-ym,ym]
xx = -xm+dx*(1:ii)';
yy = -ym+dy*(1:jj)';
%
%[Y,X] = meshgrid(yy,xx);
[X,Y] = meshgrid(xx,yy);
R     = sqrt(X.^2+Y.^2);
Theta = atan2(Y,X);
%%
scuy     = 0*X+dy;
scvx     = 0*X+dx;
scp2     = scuy.*scvx;%%cell area
LANDMASK = 0*X;

if OPT==0%%make an island;
   x0             = xm-.25*(2*xm); %%centre
   y0             = -ym+.25*(2*ym); %%centre
   R0             = 30e3;           %%radius
   j0             = find( (X-x0).^2+(Y-y0).^2<R0^2 );
   %%
   LANDMASK(j0)   = 1;
elseif OPT==1
   x0             = .5*xm;
   j0             = find(X>x0);
   LANDMASK(j0)   = 1;
elseif OPT==2
   x0             = .5*xm;
   y0             = -.5*ym;
   j0             = find((X>x0)&(Y<y0));
   LANDMASK(j0)   = 1;
end

grid_prams.X         = X;
grid_prams.Y         = Y;
grid_prams.LANDMASK  = LANDMASK;
grid_prams.scuy      = scuy;
grid_prams.scvx      = scvx;
grid_prams.scp2      = scp2;%%cell area
grid_prams.scp2i     = 1./scp2;   %%1/[cell area]

if do_test
   H  = pcolor(X/1e3,Y/1e3,LANDMASK);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   %GEN_proc_fig('\ity, \rmkm','\itx, \rmkm');
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   ttl   = title('Land mask');
   GEN_font(ttl);
   colorbar;
end
