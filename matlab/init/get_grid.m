function grid_prams = get_grid(grid_prams,OPT)
%% get_grid.m
%% Author: Timothy Williams
%% Date:   20140826, 03:08:57 CEST

do_test  = 0;
if nargin==0
   do_test  = 1;
   OPT      = 1;
   %%
   x0 = 0.0;
   y0 = 0.0;
   ii = 49;
   jj = 51;
   dx = 4e3;%m
   dy = 4e3;%m
   %%
   grid_prams  = struct('x0',x0,'y0',y0,...
                        'nx',ii,'ny',jj,...
                        'dx',dx,'dy',dy);
end
x0 = grid_prams.x0;
y0 = grid_prams.y0;
ii = grid_prams.nx;
jj = grid_prams.ny;
dx = grid_prams.dx;
dy = grid_prams.dy;

%% centres of grid cells
xx = x0-.5*dx+dx*(1:ii)';
yy = y0-.5*dy+dy*(1:jj)';
xmax  = max(xx);
ymax  = max(yy);
xm = .5*(xmax+.5*dx-x0);
ym = .5*(ymax+.5*dy-y0);
%
[Y,X] = meshgrid(yy,xx);%%x~i,y~j
R     = sqrt(X.^2+Y.^2);
Theta = atan2(Y,X);
%%
scuy     = 0*X+dy;
scvx     = 0*X+dx;
scp2     = scuy.*scvx;%%cell area
LANDMASK = 0*X;

if OPT==1
   %% land on right
   x1             = x0+1.6*xm;
   j0             = find(X>x1);
   LANDMASK(j0)   = 1;
elseif OPT==2
   %% land in lower right corner
   x1             = x0+1.5*xm;
   y1             = y0+.5*ym;
   j0             = find((X>x1)&(Y<y1));
   LANDMASK(j0)   = 1;
elseif OPT==3
   %%make an island;
   %%(circle in lower right)
   x1             = x0+1.5*xm;   %%centre
   y1             = y0+.5*ym;    %%centre
   R0             = 30e3;        %%radius
   j0             = find( (X-x1).^2+(Y-y1).^2<R0^2 );
   %%
   LANDMASK(j0)   = 1;
elseif OPT==4
   %%land on upper,lower and rhs
   LANDMASK(:,end)   = 1;% upper
   LANDMASK(:,1)     = 1;% lower
   LANDMASK(end,:)   = 1;% rhs
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
