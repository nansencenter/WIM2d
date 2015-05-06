%% test_advection_weno_1d.m
%% Author: Timothy Williams
%% Date:   20150505
clear;

%%boundary conditions:
%adv_options.ADV_OPT  = 0;%waves escape domain
adv_options.ADV_OPT  = 1;%waves periodic in i

%%testing:
ii = 49;
jj = 1;
dx = 4e3;%m
dy = 4e3;%m - not real, just used in 1/dx=scp2i*dy
%%
xm = (ii+1)/2*dx;
X  = -xm+dx*(1:ii)';
%
s1.nx       = ii;
s1.ny       = jj;
s1.scuy     = 0*X+dy;
s1.scvx     = 0*X+dx;
s1.scp2     = s1.scuy.*s1.scvx;
s1.scp2i    = 1./s1.scp2;
s1.LANDMASK = 0*X;
grid_prams  = s1;
clear s1;

CFL   = .4;

uc          = -30;%const speed m/s
xc          = 2*xm/3;
u           = 0*X+uc;
h           = 0*X;
h(X>xc)     = 1;
%%
dt          = abs(CFL*dx/uc   );
nt          = abs(2*xm/(uc*dt));
if adv_options.ADV_OPT==1
   nt = 2*nt;
end

for j=1:2%%plot initial h
   subplot(2,1,j);
   ax = plot(X/1e3,h);
   ylim([0 2]);
   GEN_proc_fig('x, km','h(0), m');
   %GEN_pause;
end

for n = 1:nt
   [n,nt]
   h_    = h;
   h     = waveadv_weno_1d(h,u,grid_prams,dt,adv_options);
   hmax  = max(h)
   %%
   subplot(2,1,2);
   plot(X/1e3,h);
   ylim([0 2]);
   GEN_proc_fig('x, km','h(t), m');
   hold on;

   x1 = xc+uc*n*dt;
   x2 = xm+uc*n*dt;
   y3 = 0*X;
   y3((X>=x1)&(X<=x2))   = 1;
   plot(X/1e3,y3,'--k');
   hold off;
   pause(.1);
end
