%% test_advection_weno.m
%% Author: Timothy Williams
%% Date:   20140821, 12:22:17 CEST
clear;

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
   uc          = 30;%const speed m/s
   xc          = 2*xm/3;
   %theta       = 180;%deg straight across
   theta       = 135;%deg
   u           = 0*X+uc*cos(pi/180*theta);
   v           = 0*X+uc*sin(pi/180*theta);
   h           = 0*X;
   h(X>xc) = 1;
   %%
   dt          = CFL*dx/uc;
   nt          = 2*xm/(uc*dt);
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

if 1%%plot u,v,h
   subplot(2,2,1);
   ax = pcolor(X/1e3,Y/1e3,u);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   daspect([1 1 1]);
   ttl   = title('u, m/s');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   %%
   subplot(2,2,2);
   ax = pcolor(X/1e3,Y/1e3,v);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   daspect([1 1 1]);
   ttl   = title('v, m/s');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   %%
   subplot(2,2,3);
   %[C,H] = contour(X,Y,h,[1,1]);
   [C,H] = contour(X,Y,h,[0,0]);
   nc    = C(2,1);
   x1    = C(1,(1:nc)+1);
   y1    = C(2,(1:nc)+1);
   %%
   [C,H] = contour(X,Y,h,[1,1]);
   nc    = C(2,1);
   x2    = C(1,(1:nc)+1);
   y2    = C(2,(1:nc)+1);
   %%
   cla;
   ax = pcolor(X/1e3,Y/1e3,h);
   set(ax, 'EdgeColor', 'none');
   caxis([0 2]);
   colorbar;
   daspect([1 1 1]);
   ttl   = title('h(0), m');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   %%
   if OPT==3
      hold on;
      plot(x1/1e3,y1/1e3,'r');
      plot(x2/1e3,y2/1e3,'m');
      hold off;
   end
   %%
   subplot(2,2,4);
   ax = pcolor(X/1e3,Y/1e3,h);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   caxis([0 2]);
   daspect([1 1 1]);
   ttl   = title('h(t), m');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   %%
   if OPT==3
      hold on;
      plot(x1/1e3,y1/1e3,'r');
      plot(x2/1e3,y2/1e3,'m');
      hold off;
   end
   GEN_pause;
end


for n = 1:nt
   [n,nt]
   %h     = waveadv_weno(h,u,v,scuy,scvx,scp2i,scp2,dt,LANDMASK);
   h     = waveadv_weno(h,u,v,grid_prams,dt);
   hmax  = max(h(:))
   %%
   subplot(2,2,4);
   ax = pcolor(X/1e3,Y/1e3,h);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   caxis([0 2]);
   daspect([1 1 1]);
   ttl   = title('h(t), m');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');

   if OPT==1
      x1 = xc+uc*cos(pi*theta/180)*n*dt;
      y1 = -ym+uc*sin(pi*theta/180)*n*dt
      x2 = xm+uc*cos(pi*theta/180)*n*dt;
      hold on;
      plot(x1/1e3+0*yy(yy>y1),yy(yy>y1)/1e3,'r');
      plot(x2/1e3+0*yy(yy>y1),yy(yy>y1)/1e3,'r');
      plot([x1,x2]/1e3,y1/1e3*[1 1],'r');
      if 1
         y3 = 0*xx;
         y3((xx>=x1)&(xx<=x2))   = ym/2;
         plot(xx/1e3,y3/1e3,'k');
         plot(xx/1e3,h(:,46)*ym/2/1e3,'--k');
      end
      hold off;

   elseif OPT==3
      hold on;
      x3 = x1*cos(n*dtheta)-y1*sin(n*dtheta);
      y3 = y1*cos(n*dtheta)+x1*sin(n*dtheta);
      plot(x3/1e3,y3/1e3,'r');
      x4 = x2*cos(n*dtheta)-y2*sin(n*dtheta);
      y4 = y2*cos(n*dtheta)+x2*sin(n*dtheta);
      plot(x4/1e3,y4/1e3,'m');
      hold off;
   end
   pause(.1);
end
