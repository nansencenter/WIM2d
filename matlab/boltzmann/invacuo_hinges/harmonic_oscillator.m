%% harmonic_oscillator.m
%% Author: Timothy Williams
%% Date: 20171010

classdef harmonic_oscillator

%% ====================================
%% properties
properties

nu,om,k;
eta0,deta0,d2eta0;
soln_type;

end
%%end properties
%% ====================================


%% ====================================
methods

function obj = harmonic_oscillator(k,nu,eta0,deta0)

   if ~exist('k','var')
      k  = 1;
   end
   if ~exist('nu','var')
      nu = .5;
   end
   if ~exist('eta0','var')
      eta0  = 1;
   end
   if ~exist('deta0','var')
      deta0 = 0;
   end

   %% class constructor
   obj.k       = k;
   obj.nu      = nu;
   obj.eta0    = eta0;
   obj.deta0   = deta0;
   obj.d2eta0  = -k^2*eta0-2*nu*deta0;

   %% dispersion relation
   if k==nu
      obj.soln_type  = 'double_pole';
      obj.om         = -1i*nu*[1;1];
   elseif k==0
      obj.soln_type  = 'zero_pole'
      obj.om         = [-2i*nu,0];
   elseif nu==0
      obj.soln_type  = 'real_poles'
      obj.om         = k*[1,-1];
   else
      obj.soln_type  = 'simple_poles'
      obj.om         = -1i*nu+sqrt(k^2-nu^2)*[-1,1];
   end

   %%
end

function out = psi(obj,omega)
   out   = obj.deta0 - 1i*omega*obj.eta0;
end

function out = disp_exact(obj,tvec,ndiff)

   if ~exist('ndiff','var')
      ndiff = 0;
   end

   if strcmp(obj.soln_type,'double_pole')

      %% function
      y  = ( obj.eta0*(1+obj.nu*tvec) + obj.deta0 ).*exp(-obj.nu*tvec);
      if ndiff==0
         out   = y;
         return
      end

      %% 1st deriv
      dy = -obj.nu*y + obj.eta0*obj.nu*exp(-obj.nu*tvec);
      if ndiff==1
         out   = dy;
         return
      end
      
      %% 2nd deriv
      out   = -obj.nu*dy - obj.eta0*obj.nu^2*exp(-obj.nu*tvec);
      return

   elseif strcmp(obj.soln_type,'zero_pole')

      b  = -obj.deta0/obj.nu;
      a  = obj.deta0/obj.nu+obj.eta0;

      if ndiff==0
         out   = a+b*exp(-obj.nu*tvec);
      else
         out   = b*(-obj.nu)^ndiff*exp(-obj.nu*tvec);
      end
      return
      
   elseif strcmp(obj.soln_type,'real_poles')
      
      a     = .5*(obj.eta0-1i/obj.k*obj.deta0);%e^{ikt}
      b     = .5*(obj.eta0+1i/obj.k*obj.deta0);%e^{-ikt}
      out   = a*(1i*obj.k)^ndiff*exp(1i*obj.k*tvec)+...
               b*(-1i*obj.k)^ndiff*exp(-1i*obj.k*tvec);
      out   = real(out);
      return
   else
      %% simple_poles
      out   = 0*tvec;
      for r=1:2
         om_r  = obj.om(r);
         a     = .5i*obj.psi(om_r+2i*obj.nu)/(om_r+1i*obj.nu);
         out   = out+a*(-1i*om_r)^ndiff*exp(-1i*om_r*tvec);
      end
      out   = real(out);
      return
   end
end

function out = disp_corrector(obj,tvec,mu)
   a     = obj.eta0;
   b     = obj.deta0+mu*obj.eta0;
   c     = obj.d2eta0+2*mu*obj.deta0+mu^2*obj.eta0;
   out   = (a+b*tvec+c/2*tvec.^2).*exp(-mu*tvec);
end

function out = transform_corrector(obj,omvec,mu)
   a     = obj.eta0;
   b     = obj.deta0+mu*obj.eta0;
   c     = obj.d2eta0+2*mu*obj.deta0+mu^2*obj.eta0;
   %%
   w     = omvec+1i*mu;
   t0    = 1i./w;
   t1    = 1i./w.*t0;
   t2    = 1i./w.*t1;
   out   = a*t0+b*t1+c*t2;
end

function out = transform(obj,omega)
   num   = -obj.psi(omega+2i*obj.nu);
   out   = num./(omega.^2+2i*obj.nu*omega-obj.k^2);
end

function out = disp_approx(obj,tvec,W,N,DO_CORR)
   
   if ~exist('W','var')
      W  = 1e3;
   end
   if ~exist('N','var')
      N  = 2^12
   end
   if ~exist('DO_CORR','var')
      DO_CORR  = 1;
   end

   % get omega values
   dom   = 2*W/N;%resolution in omega
   omega = -W-dom/2+dom*(1:N)';%centres (avoid omega=0)
   Mift  = (dom/2/pi)*exp(-1i*tvec*omega');%matrix to invert the fourier transform

   % fourier transform
   ft    = obj.transform(omega);
   corr  = 0;

   if DO_CORR
      %correct for jump at t=0
      mu    = max(obj.nu,obj.k);
      ft    = ft - obj.transform_corrector(omega,mu);
      corr  = obj.disp_corrector(tvec,mu);
   end

   % calc displacement
   out   = real(Mift*ft + corr);

end

function tvec = plot_exact(obj)
   if strcmp(obj.soln_type,'real_poles')
      T     = 2*pi/obj.k;
      tmax  = 4*T;
   else
      le    = log(1e-4);
      tmax  = -le/obj.nu;
   end

   tvec  = linspace(0,tmax,100)';
   eta   = obj.disp_exact(tvec);
   plot(tvec,eta)
end

function tvec = plot_approx(obj)

   tvec  = obj.plot_exact();
   hold on;

   eta   = obj.disp_approx(tvec);
   plot(tvec,eta,'--r');
   hold off;
end

end
%%end methods
%% ====================================

end
%%end classdef harmonic_oscillator
