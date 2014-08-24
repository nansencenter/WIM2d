function testFBwave_advect(dirn)
%% dirn is degrees clockwise from north

if 0
   adv_fn   = @FBwave_advect2D
elseif 0
   adv_fn   = @FBwave_advect2D_v1
else
   adv_fn   = @FBwave_advect2D_v2
end

dx       = 3000;%% m
dy       = dx;
dt       = 240;%% s
nsteps   = 20;
%%
if 0
   ag = dx/dt;
else
   g  = 9.81;
   T  = 16;
   om = 2*pi/T;
   k  = om^2/g;
   ag = om/k/2;
end

%%define grid
nx    = 30;
ny    = 30;
xx    = dx/2+(0:nx-1)'*dx;
yy    = flipud(dy/2+(0:ny-1)'*dy);
Dx    = nx*dx;
Dy    = ny*dy;
%%
conc  = zeros(ny,nx);%% water
S     = conc;
%%
c     = .75;
conc(end-round(ny/3):end,1:round(.8*nx)) = c;%% ice
conc(end-round(ny/5):end,1:round(.3*nx)) = -1;%% land

%% define waves
jj       = find(conc==0);
if 0
   nn = max(jj);
   j0 = jj(round(nn/2));
   S(j0:j0+10)   = 1;
   E  = sum(sum(S))
else
   S(jj)    = 1;
end
[YY,XX]  = meshgrid(yy/1e3,xx/1e3);
%plot3(XX,YY,S);
subplot(1,2,1), GEN_plot_matrix(xx/1e3,yy/1e3,conc);
daspect([1 1 1]);
%%
subplot(1,2,2), GEN_plot_matrix(xx/1e3,yy/1e3,S,[0 1]);
daspect([1 1 1]);
disp('push any key to start advection'),pause

wave_props  = {ag,dirn};
grid_props  = {dt,dx,dy};


for n = 1:nsteps
   S  = feval(adv_fn,S,wave_props,conc,grid_props);
   E  = sum(sum(S))
   subplot(1,2,2), GEN_plot_matrix(xx/1e3,yy/1e3,S,[0 1]);
   daspect([1 1 1]);
   disp([n nsteps]),pause
end

function GEN_plot_matrix(x,y,M,Mlims)

if nargin==3
   m0 = min(min(M));
   m1 = max(max(M));
else
   m0 = Mlims(1);
   m1 = Mlims(2);
end
C  = [.5 0.5 0];
%%
nx = length(x);
ny = length(y);
%%
dx    = x(2)-x(1);
dy    = y(2)-y(1);
dx0   = dx/2;
dy0   = dy/2;
%%
for i=1:ny
   for j=1:nx
      %% vertices of grid cell;
      X     = [x(j)-dx0,x(j)-dx0,x(j)+dx0,x(j)+dx0];
      Y     = [y(i)-dy0,y(i)+dy0,y(i)+dy0,y(i)-dy0];
      C(3)  = (M(i,j)-m0)/(m1-m0);%% colour;
      patch(X,Y,C), hold on;
   end
end
