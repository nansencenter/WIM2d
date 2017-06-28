%% test_iceadv_weno.m
%% Author: Timothy Williams
%% Date:   20170217
clear;
format long;

MEX      = 1
CFL      = .7;
DO_FIG2  = 1;%2x2 plot: landmask,u,v,h

%%testing:
if 0
   OPT   = 1;
   %OPT   = 2;
   ii = 150;
   jj = 10;
   iland = 80:90;
   %jland = 1:jj;
   jland = 3:6;
else
   OPT   = 3;
   ii    = 70;
   jj    = 70;
   %SHIFT = -12;
   SHIFT = 0;
   iland = (40:50)+SHIFT;
   jland = iland;
end
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
if 1
   s1.LANDMASK(iland,jland)   = 1;
end
s1.dx       = dx;
s1.dy       = dy;
s1.X        = X;
s1.Y        = Y;
grid_prams  = s1;
clear s1;


nbdy     = 4;
masks    = make_masks_puv(1-grid_prams.LANDMASK,nbdy);
masks_i  = make_masks_puv(1-grid_prams.LANDMASK,nbdy,'rows');

if OPT==1
   uc      = 30;%const speed m/s
   xc      = 2*xm/3;
   %theta   = 170;%wave-to direction [deg] - to left & up
   %theta   = 180;%wave-to direction [deg] - to left
   %theta   = 135;%wave-to direction [deg] - up and to the left
   %theta   = 0;%wave-to direction [deg] - to right
   theta   = 45;%wave-to direction [deg] - up and to the right
   %theta    = 90;%up
   u       = 0*X+uc*cos(pi/180*theta);
   v       = 0*X+uc*sin(pi/180*theta);
   h       = 0*X;
   h(X>xc) = 1;
   %%
   dt = CFL*dx/uc;
   nt = 2*xm/(uc*dt);
elseif OPT==2
   uc      = 30;%const speed m/s
   xc      = 2*xm/3;
   u       = -uc*X/xm;
   v       = 0*X;
   h       = 0*X;
   h(X>xc) = 1;
   %%
   dt = CFL*dx/uc;
   nt = 2*xm/(uc*dt);
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
   fac      = -1.5;
   angrot   = fac*(1/20)*pi/180;%%radian/s
   u        = 0*X;
   v        = 0*X;
   juv      = find(R<(Rm*1.45));
   u(juv)   = -Y(juv).*angrot;
   v(juv)   =  X(juv).*angrot;
   %%
   max_speed   = abs(Rm*angrot)
   dt          = CFL*dx/max_speed
   nt          = round(2*pi/(abs(angrot)*dt))
   dtheta      = dt*angrot;
end

h           = h.*(1-grid_prams.LANDMASK);
u           = u.*(1-grid_prams.LANDMASK);
v           = v.*(1-grid_prams.LANDMASK);
Diag.hmax   = max(h(:));
Diag.mass   = sum(h(:))

%extend for plotting with pcolor (ignores last row/col)
X(end+1,:)  = X(end,:)+dx;
Y(end+1,:)  = Y(end,:);
Y(:,end+1)  = Y(:,end)+dy;
X(:,end+1)  = X(:,end);

if DO_FIG2%%plot u,v,h
   figure(2);
   fn_fullscreen;
   subplot(2,2,1);
   Z  = 0*X;
   Z(1:end-1,1:end-1)   = u;
   ax = pcolor(X/1e3,Y/1e3,Z);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   daspect([1 1 1]);
   ix = iland(1):iland(end)+1;
   iy = jland(1):jland(end)+1;
   x  = grid_prams.X(ix,1)/1e3;
   y  = grid_prams.Y(1,iy)/1e3;
   hold on;
   plot(x,y(1)+0*x,'m','linewidth',1.5);
   plot(x,y(end)+0*x,'m','linewidth',1.5);
   plot(x(1)+0*y,y,'m','linewidth',1.5);
   plot(x(end)+0*y,y,'m','linewidth',1.5);
   hold off;
   ttl   = title('u, m/s');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   %%
   subplot(2,2,2);
   Z  = 0*X;
   Z(1:end-1,1:end-1)   = v;
   ax = pcolor(X/1e3,Y/1e3,Z);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   daspect([1 1 1]);
   ttl   = title('v, m/s');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   %%
   subplot(2,2,3);
   Z  = 0*X;
   Z(1:end-1,1:end-1)   = h;
   %%
   cla;
   ax = pcolor(X/1e3,Y/1e3,Z);
   set(ax, 'EdgeColor', 'none');
   caxis([0 2]);
   colorbar;
   daspect([1 1 1]);
   ttl   = title('h(0), m');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   %%
   subplot(2,2,4);
   Z  = 0*X;
   Z(1:end-1,1:end-1)   = h;
   ax = pcolor(X/1e3,Y/1e3,Z);
   set(ax, 'EdgeColor', 'none');
   daspect([1 1 1]);
   %colorbar;
   hold on;
   plot(x,y(1)+0*x,'m','linewidth',1.5);
   plot(x,y(end)+0*x,'m','linewidth',1.5);
   plot(x(1)+0*y,y,'m','linewidth',1.5);
   plot(x(end)+0*y,y,'m','linewidth',1.5);
   hold off;
   caxis([0 2]);
   ttl   = title('h(t), m');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   %%
   GEN_pause;
end

h0 = h;
for n = 1:nt
   [n,nt]
   if MEX==0
      'iceadv_weno'
      h           = iceadv_weno(h,u,v,grid_prams,masks,dt,nbdy);
      Diag.hmax   = max(h(:));
      Diag.mass   = sum(h(:))
   else
      if MEX==1
         'iceadv_weno_mex'
         %[h,test]    =...
         h  =...
               iceadv_weno_mex(grid_prams.nx,grid_prams.ny,dt,nbdy,h, u, v,...
                  grid_prams.scp2, grid_prams.scp2i, grid_prams.scuy, grid_prams.scvx,...
                  masks_i.pmask,masks_i.umask,masks_i.vmask,masks_i.isp,masks_i.isu,masks_i.isv,...
                  masks_i.ifp,masks_i.ifu,masks_i.ifv,masks_i.ilp,masks_i.ilu,masks_i.ilv);
      else
         'iceadv_weno_matmem_mex'
         %[h,test]  =...
         h  =...
               iceadv_weno_matmem_mex(grid_prams.nx,grid_prams.ny,dt,nbdy,h, u, v,...
                  grid_prams.scp2, grid_prams.scp2i, grid_prams.scuy, grid_prams.scvx,...
                  masks.pmask,masks.umask,masks.vmask,masks.isp,masks.isu,masks.isv,...
                  masks.ifp,masks.ifu,masks.ifv,masks.ilp,masks.ilu,masks.ilv);
      end

      Diag.hmax   = max(h(:));
      Diag.mass   = sum(h(:))
      if 0
         test
         return
      elseif ~DO_FIG2
         if 1
            'compare to matlab function'
            [h1,test1]  = iceadv_weno(h0,u,v,grid_prams,masks,dt,nbdy);
         elseif MEX==1
            'compare to original mex function (matmem)'
            h1 =...
                  iceadv_weno_matmem_mex(grid_prams.nx,grid_prams.ny,dt,nbdy,h0, u, v,...
                     grid_prams.scp2, grid_prams.scp2i, grid_prams.scuy, grid_prams.scvx,...
                     masks.pmask,masks.umask,masks.vmask,masks.isp,masks.isu,masks.isv,...
                     masks.ifp,masks.ifu,masks.ifv,masks.ilp,masks.ilu,masks.ilv);
         else
            'compare to new mex function'
            h1 =...
                  iceadv_weno_mex(grid_prams.nx,grid_prams.ny,dt,nbdy,h0, u, v,...
                     grid_prams.scp2, grid_prams.scp2i, grid_prams.scuy, grid_prams.scvx,...
                     masks_i.pmask,masks_i.umask,masks_i.vmask,masks_i.isp,masks_i.isu,masks_i.isv,...
                     masks_i.ifp,masks_i.ifu,masks_i.ifv,masks_i.ilp,masks_i.ilu,masks_i.ilv);
         end
         if 1
            PADDED   = 0;
            tst      = h1-h;
            Z1       = h1;
            Z        = h;
            tst_rng  = [min(tst(:)),max(tst(:))]
         else
            PADDED   = 1;
            tst      = test1-test;
            Z1       = test1;
            Z        = test;
            tst_rng  = [min(tst(:)),max(tst(:))]
         end
         %%
         figure;
         fn_fullscreen;
         if 1
            subplot(2,2,1)
            ax = pcolor(grid_prams.LANDMASK.');
            set(ax, 'EdgeColor', 'none');
            colorbar;
            hold on;
            x  = (iland(1):iland(end)+1);
            y  = (jland(1):jland(end)+1);
            plot(x,jland(1)+0*x,'m','linewidth',1.5);
            plot(x,jland(end)+1+0*x,'m','linewidth',1.5);
            plot(iland(1)+0*y,y,'m','linewidth',1.5);
            plot(iland(end)+1+0*y,y,'m','linewidth',1.5);
            daspect([1 1 1]);
            hold off;
            %%
            subplot(2,2,2)
            ax = pcolor(tst.');
            set(ax, 'EdgeColor', 'none');
            daspect([1 1 1]);
            colorbar;
            hold on;
            x  = PADDED*nbdy+(iland(1):iland(end)+1);
            y  = PADDED*nbdy+(jland(1):jland(end)+1);
            plot(x,PADDED*nbdy+jland(1)+0*x,'m','linewidth',1.5);
            plot(x,PADDED*nbdy+jland(end)+1+0*x,'m','linewidth',1.5);
            plot(PADDED*nbdy+iland(1)+0*y,y,'m','linewidth',1.5);
            plot(PADDED*nbdy+iland(end)+1+0*y,y,'m','linewidth',1.5);
            daspect([1 1 1]);
            hold off;
            %%
            subplot(2,2,3)
            ax = pcolor(Z1.');
            set(ax, 'EdgeColor', 'none');
            daspect([1 1 1]);
            colorbar;
            hold on;
            x  = PADDED*nbdy+(iland(1):iland(end)+1);
            y  = PADDED*nbdy+(jland(1):jland(end)+1);
            plot(x,PADDED*nbdy+jland(1)+0*x,'m','linewidth',1.5);
            plot(x,PADDED*nbdy+jland(end)+1+0*x,'m','linewidth',1.5);
            plot(PADDED*nbdy+iland(1)+0*y,y,'m','linewidth',1.5);
            plot(PADDED*nbdy+iland(end)+1+0*y,y,'m','linewidth',1.5);
            daspect([1 1 1]);
            hold off;
            %%
            subplot(2,2,4)
            ax = pcolor(Z.');
            set(ax, 'EdgeColor', 'none');
            daspect([1 1 1]);
            colorbar;
            hold on;
            x  = PADDED*nbdy+(iland(1):iland(end)+1);
            y  = PADDED*nbdy+(jland(1):jland(end)+1);
            plot(x,PADDED*nbdy+jland(1)+0*x,'m','linewidth',1.5);
            plot(x,PADDED*nbdy+jland(end)+1+0*x,'m','linewidth',1.5);
            plot(PADDED*nbdy+iland(1)+0*y,y,'m','linewidth',1.5);
            plot(PADDED*nbdy+iland(end)+1+0*y,y,'m','linewidth',1.5);
            daspect([1 1 1]);
            hold off;
         end
         return;
      end
   end
   %%
   if 0
      if OPT==1
         figure(3);
         [hmax,imax] = max(h(:,1));
         yp          = Y(imax,:);
         hp          = h(imax,:);
         plot(yp/1.0e3,hp-mean(hp));
         GEN_proc_fig('y, km','h');
         ttl   = title(['max h = ',num2str(max(h(:))),'; x = ',num2str(X(imax,1)/1.0e3),'km']);
         GEN_font(ttl);
      end
   else
      figure(2);
      subplot(2,2,4);
      if 0%n==1
         colorbar;
      end
      %cla;
      Z  = 0*X;
      Z(1:end-1,1:end-1)   = h;
      ax = pcolor(X/1e3,Y/1e3,Z);
      set(ax, 'EdgeColor', 'none');
      daspect([1 1 1]);
      caxis([0 2]);
      hold on;
      plot(x,y(1)+0*x,'m','linewidth',1.5);
      plot(x,y(end)+0*x,'m','linewidth',1.5);
      plot(x(1)+0*y,y,'m','linewidth',1.5);
      plot(x(end)+0*y,y,'m','linewidth',1.5);
      hold off;
      ttl   = title('h(t), m');
      GEN_font(ttl);
      GEN_proc_fig('x, km','y, km');
   end
   %GEN_pause;
   %pause(.1);
   drawnow;
end

if OPT==1
   figure(3);
   [hmax,imax] = max(h(:,1));
   Y_          = Y(1:end-1,1:end-1);
   yp          = Y_(imax,:);
   hp          = h(imax,:);
   plot(yp/1.0e3,hp-mean(hp));
   GEN_proc_fig('y, km','h');
   ttl   = title(['max h = ',num2str(max(h(:))),'; x = ',num2str(X(imax,1)/1.0e3),'km']);
   GEN_font(ttl);
end
