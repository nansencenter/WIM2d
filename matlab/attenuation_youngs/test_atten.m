addpath ../init;
addpath ../misc;

nw    = 25;    % no of frequencies
Tmin  = 2.5;   % min period to test
Tmax  = 25.0;  % max period to test
wtest = 1;     % extra printouts for this frequency
     
h           = 2.0;%thickness
gravity     = 9.81;
ice_prams   = fn_fill_iceprams()%other ice parameters

fmin  = 1/Tmax;
fmax  = 1/Tmin;
df    = (fmax-fmin)/(nw-1.0);
for w=1:nw
   freq  = fmin+(w-1)*df;
   om    = 2*pi*freq;
   T     = 1.0/freq;

   %% initial guess for ice wavenumber
   kw_inf   = om^2/gravity;
   wlng     = 2*pi/kw_inf;
   if (w==1)
      guess = kw_inf;
   end

   %% get atten coeff
   [damping,kice,kwtr,int_adm,NDprams,...
      alp_scat,modT,argR,argT] =...
         RT_param_outer(om,h,ice_prams,guess);

   %% update guess for ice wavenumber
   guess = kice;

   %%check attenuation coeff's, wavenumbers
   if (w==wtest)
      disp('***********************************************');
      disp('check outputs from RTparam_outer:');
      disp(sprintf('T,h %f %f',T,iceh));
      disp(sprintf('atten,damping %f %f' ,alp_scat,damping));
      disp(sprintf('ki,kw,2pi/wlng_wtr,guess %f %f %f %f',kice,kwtr,kw_inf,guess));
      disp(sprintf('|T|,argRT,s %f %f %f %f',modT,argR,argT,int_adm));
      disp('***********************************************');
      disp(' ');
   end

   disp(sprintf('T (s), atten (per floe), damping (/m): %f %f %f',T,alp_scat,damping));
end
