%% test_advection_weno.m
%% Author: Timothy Williams
%% Date:   20140821, 12:22:17 CEST
clear;

%%boundary conditions:
%ADV_OPT  = 0; %waves escape domain
ADV_OPT  = 1; %waves periodic in i,j
%ADV_OPT  = 2; %waves periodic in j (y) only
adv_options.ADV_OPT  = ADV_OPT;

%%testing:
ii = 150;
jj = 20;
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
   uc      = 30;%const speed m/s
   xc      = 2*xm/3;
   %theta   = 180;%wave-to direction [deg] - to left
   theta   = 135;%wave-to direction [deg] - up and to the left
   %theta   = 0;%wave-to direction [deg] - to right
   %theta   = 45;%wave-to direction [deg] - up and to the right
   u       = 0*X+uc*cos(pi/180*theta);
   v       = 0*X+uc*sin(pi/180*theta);
   h       = 0*X;
   h(X>xc) = 1;
   %%
   dt = CFL*dx/uc;
   nt = 2*xm/(uc*dt);
   if ADV_OPT==1
      nt = 2*nt;
   end
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

% location of fortran outputs
w2d   = getenv('WIM2D_PATH')
Fdir  = [w2d,'/fortran/numerics/test_advection/out'];
fmt   = 'float32';%%single precision

figure(2);
fn_fullscreen;

MaxDiff  = zeros(nt+1,1);
for n = 0:nt
   [n,nt]

   if n>0
      %% call matlab advection routine
      h     = waveadv_weno(h,u,v,grid_prams,dt,adv_options);
      hmax  = max(h(:))
   end

   %% load fortran binary files
   nnn      = num2str(n,'%3.3d');
   afile    = [Fdir,'/ADVweno',nnn,'.a']
   if exist(afile)
      aid      = fopen(afile,'rb');
      XX       = reshape( fread(aid,ii*jj,fmt), ii,jj );
      YY       = reshape( fread(aid,ii*jj,fmt), ii,jj );
      hF       = reshape( fread(aid,ii*jj,fmt), ii,jj );
      fclose(aid);
   else
      MaxDiff  = MaxDiff(1:n);
      break;
   end

   %% plot h
   subplot(3,1,1);
   ax = pcolor(X/1e3,Y/1e3,h);
   set(ax, 'EdgeColor', 'none');
   daspect([1 1 1]);
   colorbar;
   caxis([0 2]);
   ttl   = title('h(t), m (.m)');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');

   %% plot hF
   subplot(3,1,2);
   ax = pcolor(X/1e3,Y/1e3,hF);
   set(ax, 'EdgeColor', 'none');
   daspect([1 1 1]);
   colorbar;
   caxis([0 2]);
   ttl   = title('h(t), m (.F)');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');

   %% plot difference
   dh             = hF-h;
   maxdiff        = max(abs(dh(:)));
   MaxDiff(n+1)   = maxdiff;
   %%
   subplot(3,1,3);
   ax = pcolor(X/1e3,Y/1e3,dh);
   set(ax, 'EdgeColor', 'none');
   daspect([1 1 1]);
   colorbar;
   caxis(1e-6*[-1 1]);
   ttl   = title(['{\Delta}h(t), m (.F - .m), |{\Delta}h|<',num2str(maxdiff)]);
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   %%
   drawnow;
end

figure(3);
plot(MaxDiff);
GEN_proc_fig('Time step','max |{\Delta}h|');
