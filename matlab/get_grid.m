function [X,Y,scuy,scvx,scp2,scp2i,LANDMASK] = get_grid()
%% get_grid.m
%% Author: Timothy Williams
%% Date:   20140826, 03:08:57 CEST

%%testing:
ii = 51;
jj = 51;
dx = 4e3;%m
dy = 4e3;%m
%%
xm = (ii+1)/2*dx;
ym = (jj+1)/2*dy;
xx = -xm+dx*(1:ii)';
yy = -ym+dy*(1:jj)';
%
[Y,X] = meshgrid(yy,xx);
R     = sqrt(X.^2+Y.^2);
Theta = atan2(Y,X);

scuy     = 0*X+dy;
scvx     = 0*X+dx;
scp2     = scuy.*scvx;
scp2i    = 1./scp2;
LANDMASK = 0*X;

%%make an island;
x0             = -xm+.25*(2*xm); %%centre
y0             = -ym+.25*(2*ym); %%centre
R0             = 30e3;           %%radius
j0             = find( (X-x0).^2+(Y-y0).^2<R0^2 );
LANDMASK(j0)   = 1;

if 1
   pcolor(LANDMASK);
end
