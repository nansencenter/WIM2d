function [Dave,FSD] = floe_scaling(Dmax,prams,moment)
%% CALL [Dave,FSD] = floe_scaling(Dmax,prams,moment)
%% This function computes the average floe size within a grid cell as a
%% function of the maximum floe size using a bounded fractal renormalization
%% group method.
%%
%% We suggest using Dmin >= 20. Below that value, there is no scattering
%% by the floes for periods larger than 6 s (it's probably viscous however).
%%
%% Dave  = average floe size
%% Dmax  = max floe size
%% prams = [structure] eg
%%            xi: 2
%%     fragility: 0.900000000000000
%%      Dmax_min: 20
%%       Dthresh: 200 (optional)
%% moment=1,2 -> <D> or <D^2>

if ~exist('moment','var')
   moment   = 1;
end
f        = prams.fragility;
xi       = prams.xi;
Dmin     = prams.Dmax_min;
Dave     = max(Dmax.^moment,Dmin^moment);
if isfield(prams,'Dthresh')
   Dthresh  = prams.Dthresh;
else
   Dthresh  = 200;
end
F        = f*xi^2;

want_fsd = 0;
if nargout==2
   want_fsd = 1;
end

for j=1:length(Dmax(:));
   dmax     = Dmax(j);
   unifdist = (dmax<xi*Dmin)|(dmax>Dthresh);

   if ~unifdist
      r  = dmax/Dmin;
      M  = 0;
      while ( r >= xi ) 
         %% if r<xi,no more breaking
         %% - don't change Dave
         r  = r/xi;
         M  = M +1;
      end
      
      if want_fsd
         %%store FSD
         Nvec  = zeros(M+1,1);
         Dvec  = dmax./(xi.^(0:M))';
         %%
         Nvec(1:M)   = (1-f)*F.^(0:M-1);
         Nvec(M+1)   = F^M;
         %%
         Nsum  = sum(Nvec);
         Ndsum = sum(Nvec.*Dvec.^moment);
         
         %% FSD output
         FSD(j).N             = Nvec;
         FSD(j).probability   = Nvec/Nsum;
         FSD(j).floe_sizes    = Dvec;
         FSD(j).type          = 'RG';
      else
         %% loop, don't store FSD
         %% - how it is in fortran and C++
         Nm1   = 1;
         Dm    = dmax;
         Nsum  = 0;
         Ndsum = 0;

         for m=0:M-1
            Nm    = (1-f)*Nm1;%(1-f)*F^m
            Nsum  = Nsum+Nm;
            Ndsum = Ndsum+Nm*Dm^moment;
            %
            Nm1   = Nm1*F;%F^(m+1)
            Dm    = Dm/xi;%Dmax/(xi^(m+1))
         end

         %m=M
         Nsum  = Nsum+Nm1;
         Ndsum = Ndsum+Nm1*Dm^moment;
      end

      %calc average
      Dave(j)  = Ndsum/Nsum;
   else
      Dm       = max(dmax,Dmin);
      Dave(j)  = Dm^moment;
      if want_fsd
         FSD(j).floe_sizes    = Dm;
         FSD(j).probability   = 1;
         FSD(j).type          = 'uniform';
      end
   end
end
