function [Dave,FSD] = floe_scaling(Dmax,prams,moment)
% This function computes the average floe size within a grid cell as a
% function of the maximum floe size using a bounded fractal renormalization
% groud method.

% We suggest to use Dmin >= 20. Below that value, there is no scattering
% by the floes for periods larger than 6 s (it's probably viscous however).

if ~exist('moment','var')
   moment   = 1;
end
f        = prams.f;
xi       = prams.xi;
Dmin     = prams.Dmin;
Dave     = max(Dmax.^moment,Dmin^moment);
Dthresh  = 200;
F        = f*xi^2

want_fsd = 0;
if nargout==2
   want_fsd = 1;
end

for j=1:length(Dmax(:));
   dmax     = Dmax(j);
   unifdist = (dmax<xi*Dmin)|(dmax>=Dthresh);

   if ~unifdist
      M  = floor(log2(dmax/Dmin));%should be >=1
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
         FSD(j).probability   = Nvec/Nsum;
         FSD(j).floe_sizes    = Dvec;
         FSD(j).type          = 'RG';
      else
         %% loop, don't store FSD
         %% - how it is in fortran and C++
         Nm1   = 1;
         Nsum  = 0;
         Ndsum = 0;
         Dm    = dmax;

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
