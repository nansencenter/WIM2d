%% test_waveadv_weno_F.m
%% Author: Timothy Williams
%% Date:   20140821, 12:22:17 CEST

clear;
%%testing:

%%set precision of binary outputs
fmt   = 'float32';%%single precision
%fmt   = 'float64';%%single precision

ADV_DIM  = 1;
nbdy     = 3;

ii = 150;
jj = 50;
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

OPT   = 1;
CFL   = .4;

if OPT==1
   uc    = 30;%const speed m/s
   xc    = 2*xm/3;
   %theta = 180;%deg straight across
   theta = 135;%deg
   u     = 0*X+uc*cos(pi/180*theta);
   v     = 0*X+uc*sin(pi/180*theta);
   %%
   dt = CFL*dx/uc;
   nt = 2*xm/(uc*dt);
elseif OPT==2
   uc = 30;%const speed m/s
   xc = 2*xm/3;
   u  = -uc*X/xm;
   v  = 0*X;
   %%
   dt = CFL*dx/uc;
   nt = 2*xm/(uc*dt);
elseif OPT==3
   Rm    = xm/3;
   Ym    = ym/12;
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

if 0%%test outputs from mod_waveadv_weno.F
   %afile       = 'test_out/scp2i.a';
   afile       = 'test_out/h.a';
   aid         = fopen(afile,'rb');
   fad         = ADV_DIM-1;
   test_pad    = reshape( fread(aid,fmt), ii+2*nbdy,jj+2*nbdy*fad );
   fclose(aid);
   %%
   ip = (1-nbdy:ii+nbdy);
   jp = (1-nbdy*fad:jj+nbdy*fad);
   [J,I] = meshgrid(jp,ip);
   ax    = pcolor(I,J,test_pad);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   GEN_proc_fig('I','J');
   return;
end

if 0%%test outputs from mod_waveadv_weno.F

   fad   = ADV_DIM-1;
   Nf    = 4+2*fad;

   if 1
      afile    = 'test_out/all.a';
      aid      = fopen(afile,'rb');
      test_pad = fread(aid,fmt);
      fclose(aid);
      test_pad = reshape( test_pad, ii+2*nbdy,jj+2*nbdy*fad,Nf );
      ip       = (1-nbdy:ii+nbdy);
      jp       = (1-nbdy*fad:jj+nbdy*fad);
   else
      afile    = 'test_out/all0.a';
      aid      = fopen(afile,'rb');
      test_pad = fread(aid,fmt);
      fclose(aid);
      test_pad = reshape( test_pad, ii,jj,Nf );
      ip       = (1:ii);
      jp       = (1:jj);
   end
   %%
   [J,I] = meshgrid(jp,ip);
   if ADV_DIM==2
      ttls  = {'u','v','scp2','scp2i','scuy','scvx'};
      tst   = [NaN,NaN,dx*dy,1/dx/dy,dy,dx];
   else
      ttls  = {'u','scp2','scp2i','scuy'};
      tst   = [NaN,dx*dy,1/dx/dy,dy];
   end
   %%
   for j=1:Nf
      tp    = test_pad(:,:,j);
      disp([ttls{j},' : ',num2str(tst(j))]);
      rng   = [min(tp(:)),max(tp(:))]
      ax    = pcolor(I,J,tp);
      set(ax, 'EdgeColor', 'none');
      GEN_proc_fig('I','J');
      ttl   = title(ttls{j});
      GEN_font(ttl);
      colorbar;
      GEN_pause;
   end
   return;
end

%%read in initial X,Y,h
outfile  = 'out/ADVweno';
n        = 0;
nnn      = num2str(n,'%3.3d');
afile    = [outfile,nnn,'.a'];
aid      = fopen(afile,'rb');
XX       = reshape( fread(aid,ii*jj,fmt), ii,jj );
YY       = reshape( fread(aid,ii*jj,fmt), ii,jj );
hh       = reshape( fread(aid,ii*jj,fmt), ii,jj );
fclose(aid);

if 0%%test initial conditions
   ax = pcolor(XX/1e3,YY/1e3,hh);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   daspect([1 1 1]);
   ttl   = title('Advected quantity (h)');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   return
end

if 1%%plot u,v,h
   subplot(2,2,1);
   ax = pcolor(XX/1e3,YY/1e3,u);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   daspect([1 1 1]);
   ttl   = title('u, m/s');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   %%
   subplot(2,2,2);
   ax = pcolor(XX/1e3,YY/1e3,v);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   daspect([1 1 1]);
   ttl   = title('v, m/s');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');
   %%
   subplot(2,2,3);
   %[C,H] = contour(X,Y,h,[1,1]);
   [C,H] = contour(XX,YY,hh,[0,0]);
   nc    = C(2,1);
   x1    = C(1,(1:nc)+1);
   y1    = C(2,(1:nc)+1);
   %%
   [C,H] = contour(XX,YY,hh,[1,1]);
   nc    = C(2,1);
   x2    = C(1,(1:nc)+1);
   y2    = C(2,(1:nc)+1);
   %%
   cla;
   ax = pcolor(XX/1e3,YY/1e3,hh);
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
   ax = pcolor(XX/1e3,YY/1e3,hh);
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

nt = length(dir('out/*.a'))-1;
if strcmp(fmt,'float32')
   element_size   = 4;
elseif strcmp(fmt,'float64')
   element_size   = 8;
end

for n = 1:20:nt
   %% open output file
   nnn      = num2str(n,'%3.3d');
   afile    = [outfile,nnn,'.a'];
   aid      = fopen(afile,'rb');
   fseek(aid,2*element_size*ii*jj,'bof');%%skip ahead 2 records (2*8B*(ii*jj))
   hh       = reshape( fread(aid,ii*jj,fmt), ii,jj );
   fclose(aid);
   %%
   [n,nt]
   hmax  = max(hh(:))

   %%plot advected thing
   subplot(2,2,4);
   ax = pcolor(XX/1e3,YY/1e3,hh);
   set(ax, 'EdgeColor', 'none');
   colorbar;
   caxis([0 2]);
   daspect([1 1 1]);
   ttl   = title('h(t), m');
   GEN_font(ttl);
   GEN_proc_fig('x, km','y, km');

   %%add test plots
   if OPT==1
      x1 = xc+uc*cos(pi*theta/180)*n*dt;
      y1 = -ym+uc*sin(pi*theta/180)*n*dt;
      x2 = xm+uc*cos(pi*theta/180)*n*dt;
      hold on;
      plot(x1/1e3+0*yy(yy>y1),yy(yy>y1)/1e3,'r');
      plot(x2/1e3+0*yy(yy>y1),yy(yy>y1)/1e3,'r');
      plot([x1,x2]/1e3,y1/1e3*[1 1],'r');
      if 1
         y3 = 0*xx;
         y3((xx>=x1)&(xx<=x2))   = ym/2;
         plot(xx/1e3,y3/1e3,'k');
         plot(xx/1e3,hh(:,46)*ym/2/1e3,'--k');
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
   %GEN_pause
end
