function S = fn_pixel_edges(x,y,z)
%% fn_pcolor.m
%% Author: Timothy Williams
%% Date: 20150820, 13:32:28 CEST

DO_TEST  = 0;
if nargin==0
   %%test inputs
   z        = hadamard(12);
   x        = 1:6;
   y        = 1:12;
   z        = z(x,:);
   z(z==-1) = 0;
   DO_TEST  = 1;
end

if size(x,2)>1
   x  = x';
end
if size(y,2)>1
   y  = y';
end

dx    = mean(diff(x));
dy    = mean(diff(y));
cval  = 1;%% contour to look for

if DO_TEST==1
   %% pcolor ignores last if shading is faceted
   [nx,ny]        = size(z);
   Z              = zeros(nx+1,ny+1);
   Z(1:nx,1:ny)   = z;

   %% want x,y to refer to edges of cells
   xx = cen2edges(x);
   yy = cen2edges(y);

   figure(101);
   fn_fullscreen;
   H  = pcolor(xx,yy,Z');
   set(H,'EdgeColor', 'none');
   hold on;
   [C,H] = contour(x,y,z',cval);
else
   C  = contourc(x,y,z',cval);
end

M  = 0;
while ~isempty(C);
   M  = M+1;
   N  = C(2,1);
   xm = C(1,2:N+1);
   ym = C(2,2:N+1);
   %%
   C(:,1:N+1)        = [];
   [S(M).x,S(M).y]   = conts2pixel_borders(xm,ym,dx,dy);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xx = cen2edges(x);
dfx   = diff(x);
xx    = (x(1:end-1)+x(2:end))/2;
xx    = [x(1)-.5*dfx(1);xx;x(end)+.5*dfx(end)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,yy] = conts2pixel_borders(x,y,dx,dy);

DO_TEST  = 0;
xx       = [];
yy       = [];

for j=1:length(x)
   M  = length(xx);
   xj = x(j)-x(1);
   yj = y(j)-y(1);

   if iseq(round(xj/dx),xj/dx,1e-7*dx )
      %% x is integer: horizontal boundary
      if DO_TEST
         disp('horizontal boundary');
      end
      Xj = x(1)+[xj-.5*dx,xj+.5*dx]';
      Yj = y(1)+[yj,yj]';
   else
      %% x is half-integer: vertical boundary
      if DO_TEST
         disp('vertical boundary');
      end
      Xj = x(1)+[xj,xj]';
      Yj = y(1)+[yj-.5*dy,yj+.5*dy]';
   end

   %% fill up output vectors
   if M==0

      %% check sense of traversing contour
      uj = x(j+1);
      vj = y(j+1);
      s1 = rdist(uj,vj,Xj(1),Yj(1));
      s2 = rdist(uj,vj,Xj(2),Yj(2));
      if s1<s2
         Xj = flipud(Xj);
         Yj = flipud(Yj);
      end

      xx = Xj;
      yy = Yj;

   else
      s1 = rdist(xx(end),yy(end),Xj(1),Yj(1));
      if s1<1e-7*min(dx,dy) 
         xx = [xx;Xj(2)];
         yy = [yy;Yj(2)];
      else
         xx = [xx;Xj(1)];
         yy = [yy;Yj(1)];
      end
   end

   if DO_TEST
      xj/dx,yj/dy,xx/dx,yy/dy
      plot(xx,yy,'c','linewidth',2);
      pause
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=iseq(a,b,tol)

y  = ( abs(b-a)<tol );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=rdist(x,y,u,v)

s  = sqrt((x-u)^2+(y-v)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
