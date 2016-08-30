%% interp1d.m
%% Author: Timothy Williams
%% Date: 20160830
function M_interp = interp1d(nodes_out,nodes_in,method)

DO_TEST  = 0;
if nargin==0
   DO_TEST     = 1;
   nodes_in    = pi*linspace(0,2,500)';
   nodes_out   = pi*linspace(.5,1.5,50)';
   %method      = 'NodeToElement'
   method      = 'ElementToElement'
   %method      = 'ElementToNode' %TODO
   %method      = 'NodeToNode'    %TODO
end

nn1   = length(nodes_in);
nn2   = length(nodes_out);
ne1   = nn1-1;
ne2   = nn2-1;
tol   = 1e-12;

switch method;
case('ElementToElement')
   %% average value over each x interval
   %% - matrix is just fractional overlap

   M_interp = zeros(ne2,ne1);
   for i=1:ne2
      ym = nodes_out(i);
      yp = nodes_out(i+1);
      dy = yp-ym;
      %%
      xL = max(ym,nodes_in(1:end-1));
      xU = min(yp,nodes_in(2:end));
      Hm = heaviside(nodes_in(2:end)-ym).*heaviside(yp-nodes_in(1:end-1));
      Im = Hm.*(xU-xL);
      %%
      M_interp(i,:)  = Im/dy;
   end

case('ElementToNode')
   %% Get best straight-line fit
   error('method  = ''ElementToNode'' not implemented');

   M_interp = zeros(nn2,ne1);
   for i=1:nn2
      xn = nodes_out(i);
      dx = abs(xn-nodes_in);
      jn = find(min(dx)==dx);
      if (xn<nodes_in(1)-tol)|(xn>nodes_in(end)+tol)
         continue;
      end
      if (dx>tol)&(jn>1)
         M_interp(i,jn)    = 1;
      else
         M_interp(i,jn-1)  = .5;
         M_interp(i,jn)    = .5;
      end
   end

case('NodeToNode')
   %% linear interpolation
   error('method = ''NodeToNode'' not implemented');

   M_interp = zeros(nn2,nn1);

case('NodeToElement')
   %% average value over each x interval

   M_interp = zeros(ne2,nn1);
   for i=1:ne2
      ym = nodes_out(i);
      yp = nodes_out(i+1);
      dy = yp-ym;
      %%
      xL = max(ym,nodes_in(1:end-1));
      xU = min(yp,nodes_in(2:end));
      Dx = nodes_in(2:end)-nodes_in(1:end-1);
      Hm = heaviside(nodes_in(2:end)-ym).*heaviside(yp-nodes_in(1:end-1));
      %%
      Im1   = .5*Hm./Dx.*...
               ( (xU-nodes_in(1:end-1)).^2 - (xL-nodes_in(1:end-1)).^2 );
      Im0   = (xU-xL).*Hm-Im1;
      %%
      M_interp(i,1:end-1)  = Im0/dy;
      M_interp(i,2:end)    = M_interp(i,2:end).'+Im1/dy;
   end

end


if DO_TEST
   ss = strsplit(method,'To');
   f  = cos(nodes_in);

   if strcmp(ss{1},'Node')
      %% interp from nodes
      plot(nodes_in,f,'-x');
      hold on;
   else
      f     = .5*(f(2:end)+f(1:end-1));
      X     = zeros(2*ne1,1);
      Y     = zeros(2*ne1,1);
      X(1:2:2*ne1)   = nodes_in(1:end-1);
      X(2:2:2*ne1)   = nodes_in(2:end);
      Y(1:2:2*ne1)   = f;
      Y(2:2:2*ne1)   = f;
      plot(X,Y,'-b');
      hold on;
   end

   fap   = M_interp*f;
   if strcmp(ss{2},'Node')
      plot(nodes_out,fap,'-xr');
   else
      X     = zeros(2*ne2,1);
      Y     = zeros(2*ne2,1);
      X(1:2:2*ne2)   = nodes_out(1:end-1);
      X(2:2:2*ne2)   = nodes_out(2:end);
      Y(1:2:2*ne2)   = fap;
      Y(2:2:2*ne2)   = fap;
      plot(X,Y,'--r');
   end
   hold off;
   fn_fullscreen;

end

function y=heaviside(x)
y  = max(0,sign(x));
