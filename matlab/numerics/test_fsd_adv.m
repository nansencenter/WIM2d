%% test_advection_weno.m
%% Author: Timothy Williams
%% Date:   20140821, 12:22:17 CEST
clear;

%%boundary conditions:
%bc_opt   = 0; %waves escape domain
bc_opt   = 1; %waves periodic in i,j
%bc_opt   = 2; %waves periodic in j (y) only

%%testing:
ii = 49;
jj = 51;
dx = 4e3;%m
dy = 4e3;%m
%%
xm = (ii+1)/2*dx;
ym = (jj+1)/2*dy;
xx = -xm+dx*(1:ii)';
yy = -ym+dy*(1:jj)';
%
%[X,Y] = meshgrid(xx,yy);
[Y,X] = meshgrid(yy,xx);%%x~i,y~j
R     = sqrt(X.^2+Y.^2);
Theta = atan2(Y,X);

s1.nx       = ii;
s1.ny       = jj;
s1.scuy     = 0*X+dy;
s1.scvx     = 0*X+dx;
s1.scp2     = s1.scuy.*s1.scvx;
s1.scp2i    = 1./s1.scp2;
s1.LANDMASK = 0*X;
grid_prams  = s1;
clear s1;

OPT   = 1;
CFL   = .4;

if OPT==1
   uc          = 10;%% const speed [km/day]
   uc          = uc*1e3/(24*3600);%const speed [m/s]
   xc          = xm/3;
   %theta       = 180;%wave-to direction [deg] - to left
   %theta       = 135;%wave-to direction [deg] - up and to the left
   theta       = 0;%wave-to direction [deg] - to right
   u           = 0*X+uc*cos(pi/180*theta);
   v           = 0*X+uc*sin(pi/180*theta);
   %%
   ICE_MASK             = 0*X;
   ICE_MASK(abs(X)<xc)  = 1;
   fice                 = .8*ICE_MASK;
   Dfloe                = 250*ICE_MASK;
   %%
   dt          = CFL*dx/uc;
   nt          = 2*xm/(uc*dt);
   if bc_opt==1
      nt = 2*nt;
   end
elseif OPT==2
   uc          = 30;%const speed m/s
   xc          = 2*xm/3;
   u           = -uc*X/xm;
   v           = 0*X;
   h           = 0*X;
   h(X>xc) = 1;
   %%
   dt          = CFL*dx/uc;
   nt          = 2*xm/(uc*dt);
elseif OPT==3
   Rm    = xm/3;
   Ym    = ym/12;
   %%
   h        = 0*X;
   jwave    = find((R<Rm)&(X<=0));%%lhs of pacman
   h(jwave) = 1;
   jwave    = find((R<Rm)&(X>0)&(abs(Y)>Ym));%%rhs of pacman
   h(jwave) = 1;
   %%
   angrot   = (1/20)*pi/180;%%radian/s
   u        = 0*X;
   v        = 0*X;
   juv      = find(R<(Rm*1.45));
   u(juv)   = -Y(juv).*angrot;
   v(juv)   =  X(juv).*angrot;
   %%
   %%
   max_speed   = Rm*angrot
   dt          = CFL*dx/max_speed
   nt          = round(2*pi/(angrot*dt))
   dtheta      = dt*angrot;
end

Dfloe2   = Dfloe;
Nfloe    = ICE_MASK;
for i=1:ii
for j=1:jj
   if ICE_MASK(i,j)==1
      Nfloe(i,j)  = fice(i,j)/Dfloe(i,j)^2;
   end
end
end

%%ice edge
col         = 'm';
[C,H]       = contour(X/1e3,Y/1e3,ICE_MASK,[0,0]);
nc          = C(2,1);
x1          = C(1,(1:nc)+1);
y1          = C(2,(1:nc)+1);
[C,H]       = contour(X/1e3,Y/1e3,ICE_MASK,[1,1]);
C(:,1:nc+1) = [];
nc          = C(2,1);
x2          = C(1,(1:nc)+1);
y2          = C(2,(1:nc)+1);
%%
subplot(2,2,1);
ax = pcolor(X/1e3,Y/1e3,fice);
set(ax, 'EdgeColor', 'none');
colorbar;
daspect([1 1 1]);
ttl   = title('Concentration');
GEN_font(ttl);
hold on;
plot(x1,y1,col);
plot(x2,y2,col);
hold off;
GEN_proc_fig('x, km','y, km');
%%
subplot(2,2,2);
ax = pcolor(X/1e3,Y/1e3,Dfloe);
set(ax, 'EdgeColor', 'none');
colorbar;
daspect([1 1 1]);
ttl   = title('Dfloe, m');
GEN_font(ttl);
hold on;
plot(x1,y1,col);
plot(x2,y2,col);
hold off;
caxis([0,300]);
GEN_proc_fig('x, km','y, km');
%%
subplot(2,2,3);
ax = pcolor(X/1e3,Y/1e3,Nfloe);
set(ax, 'EdgeColor', 'none');
colorbar;
daspect([1 1 1]);
ttl   = title('Nfloe, m^{-2}');
GEN_font(ttl);
hold on;
plot(x1,y1,col);
plot(x2,y2,col);
hold off;
GEN_proc_fig('x, km','y, km');
%%
subplot(2,2,4);
ax = pcolor(X/1e3,Y/1e3,Dfloe2);
set(ax, 'EdgeColor', 'none');
colorbar;
daspect([1 1 1]);
ttl   = title('Dfloe2, m');
GEN_font(ttl);
hold on;
plot(x1,y1,col);
plot(x2,y2,col);
hold off;
caxis([0,300]);
GEN_proc_fig('x, km','y, km');
%%
%pause(.1);
GEN_pause;

fmin  = .05;
for n = 1:nt
   [n,nt]
   im1         = 0*fice;
   jI          = find(fice>fmin);
   im1(jI)     = 1;
   Nfloe2      = 0*fice;
   Nfloe2(jI)  = fice(jI)./Dfloe2(jI).^2;
   %%
   fice     = waveadv_weno(fice,u,v,grid_prams,dt,bc_opt);
   Dfloe    = waveadv_weno(Dfloe,u,v,grid_prams,dt,bc_opt);
   Nfloe    = waveadv_weno(Nfloe,u,v,grid_prams,dt,bc_opt);
   Nfloe2   = waveadv_weno(Nfloe2,u,v,grid_prams,dt,bc_opt);
   %%
   Dfloe2         = 0*Nfloe2;
   im2            = 0*Nfloe2;
   im2(fice>fmin) = 1;
   Dfloe2         = im2;
   jI             = find(im2);
   if 1
      Dfloe2(jI)     = sqrt(fice(jI)./Nfloe(jI));
      max(Dfloe2(:))
   else
      Dfloe2(jI)     = sqrt(fice(jI)./Nfloe2(jI));
      max(Dfloe2(:))
   end

   %%ice edge
   [C,H]       = contour(X/1e3,Y/1e3,im2,[0,0]);
   nc          = C(2,1);
   x1          = C(1,(1:nc)+1);
   y1          = C(2,(1:nc)+1);
   [C,H]       = contour(X/1e3,Y/1e3,im2,[1,1]);
   C(:,1:nc+1) = [];
   nc          = C(2,1);
   x2          = C(1,(1:nc)+1);
   y2          = C(2,(1:nc)+1);

   %%plot
   subplot(2,2,1);
   ax = pcolor(X/1e3,Y/1e3,fice);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   daspect([1 1 1]);
   ttl   = title('Concentration');
   GEN_font(ttl);
   hold on;
   plot(x1,y1,col);
   plot(x2,y2,col);
   hold off;
   GEN_proc_fig('x, km','y, km');
   %%
   subplot(2,2,2);
   ax = pcolor(X/1e3,Y/1e3,Dfloe);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   daspect([1 1 1]);
   ttl   = title('Dfloe, m');
   GEN_font(ttl);
   hold on;
   plot(x1,y1,col);
   plot(x2,y2,col);
   hold off;
   caxis([0,300]);
   GEN_proc_fig('x, km','y, km');
   %%
   subplot(2,2,3);
   ax = pcolor(X/1e3,Y/1e3,Nfloe);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   daspect([1 1 1]);
   ttl   = title('Nfloe, m^{-2}');
   GEN_font(ttl);
   hold on;
   plot(x1,y1,col);
   plot(x2,y2,col);
   hold off;
   GEN_proc_fig('x, km','y, km');
   %%
   subplot(2,2,4);
   if 0
      plot(X(:,10)/1e3,Dfloe2(:,10));
      yl = get(gca,'ylim');
      hold on;
      plot(x1(10)+0*yl,yl,col);
      plot(x2(10)+0*yl,yl,col);
      hold off;
      ylim([0,300]);
      GEN_proc_fig('x, km','Dfloe2, m');
   else
      ax = pcolor(X/1e3,Y/1e3,Dfloe2);
      set(ax, 'EdgeColor', 'none');
      colorbar;
      daspect([1 1 1]);
      ttl   = title('Dfloe2, m');
      GEN_font(ttl);
      hold on;
      plot(x1,y1,col);
      plot(x2,y2,col);
      hold off;
      caxis([0,300]);
      GEN_proc_fig('x, km','y, km');
   end
   pause(.1);
end
