%% FSD_wrapper.m
%% Author: Timothy Williams
%% Date: 20161213
function output = FSD_wrapper(inputs,params,output_type)

DO_TEST  = 0;
if nargin==0
   inputs.Dmax = [0;(10:10:300)'];
   inputs.c    = [0;.9+0*inputs.Dmax];
   %%
   params         = set_FSD_params();
   params.FSD_OPT = 1;
   %%
   output_type = 'cohesion_factor';
   DO_TEST     = 1;
end

if ~exist('output_type','var')
   output_type = 'Dmean';
end

if isfield(inputs,'Dmax')
   Dmax  = inputs.Dmax;

elseif isfield(inputs,'Nfloes')&isfield(inputs,'c')
   %Dmax  = Nfloes_to_Dmax(inputs.Nfloes,inputs.c,params);
   Dmax        = 0*inputs.Nfloes;
   jpos        = find(inputs.Nfloes>0);
   Dmax(jpos)  = sqrt(inputs.c(jpos)./inputs.Nfloes(jpos));

elseif isfield(inputs,'cDmax')&isfield(inputs,'c')
   Dmax        = 0*inputs.cDmax;
   jpos        = find(inputs.c>0);
   Dmax(jpos)  = inputs.cDmax(jpos)./inputs.c(jpos);

elseif isfield(inputs,'Dmean')
   inverse  = 1;
   if params.FSD_OPT==1
      Dmax  = floe_scaling_smooth(Dmean,params,1,inverse);
   else
      error('<D> -> D_max mapping not available for FSD_OPT==0 (RG method)')
   end
end


%% check Dmax hasn't got too big on grid:
jbig        = find(Dmax>params.Dmax_pack_thresh);
Dmax(jbig)  = params.Dmax_pack;

%% if conc too low set Dmax=0
jlow        = find(inputs.c<params.cice_min);
Dmax(jlow)  = 0;


if strcmp(output_type,'Dmax')
   output   = Dmax;

   if DO_TEST
      plot(Dmax,output);
   end

   return
end

if strcmp(output_type,'Dmean')
   if params.FSD_OPT==1
      output   = floe_scaling_smooth(Dmax,params,1,0);
   else
      output   = floe_scaling(Dmax,params,1);
   end

   if DO_TEST
      plot(Dmax,output);
   end

   return
end

if strcmp(output_type,'Nfloes')
   output         = 0*inputs.c;
   output(Dmax>0) = inputs.c(Dmax>0)./(Dmax(Dmax>0).^2);

   if DO_TEST
      plot(Dmax,output);
   end

   return
end
if strcmp(output_type,'cDmax')
   output   = inputs.c.*Dmax;

   if DO_TEST
      plot(Dmax,output);
   end

   return
end

if strcmp(output_type,'lat_surface_fraction')
   if params.FSD_OPT==0
      Dmean    = floe_scaling(Dmax,params,1);
      Dsq_mean = floe_scaling(Dmax,params,2);
   else
      Dmean    = floe_scaling_smooth(Dmax,params,1,0);
      Dsq_mean = floe_scaling_smooth(Dmax,params,2,0);
   end
   output      = 0*Dmean;
   jp          = find(Dsq_mean>0);
   output(jp)  = (pi*Dmean(jp).*inputs.c(jp).*inputs.h(jp))./(pi/4*Dsq_mean(jp));

   if DO_TEST
      plot(Dmax,output);
   end

   return
end

if strcmp(output_type,'cohesion_factor')

   if params.FSD_OPT==0
      Dmean       = floe_scaling(Dmax,params,1);
      Dsq_mean    = floe_scaling(Dmax,params,2);
      Dmean_L     = floe_scaling(params.Dthresh,params,1);
      Dsq_mean_L  = floe_scaling(params.Dthresh,params,2);
   else
      Dmean       = floe_scaling_smooth(Dmax,params,1,0);
      Dsq_mean    = floe_scaling_smooth(Dmax,params,2,0);
      Dmean_L     = floe_scaling_smooth(params.Dthresh,params,1,0);
      Dsq_mean_L  = floe_scaling_smooth(params.Dthresh,params,2,0);
   end

   %% P=perimeter per unit area
   %% P=N*(pi*Dmean), N=floes per unit area
   %% N=c*(pi/4*Dsq_mean)
   %% cohesion_factor=sqrt(P_pack/P) -> 1^+ as Dmax->Dmax_pack
   output      = 1+0*Dmean;
   jp          = find(Dmean>0&inputs.c>0);
   P_inv       = (pi/4*Dsq_mean(jp))./(pi*Dmean(jp).*inputs.c(jp));%1/P 
   P_pack      = (pi*Dmean_L*inputs.c(jp))./(pi/4*Dsq_mean_L);
   output(jp)  = min(1,sqrt(P_pack.*P_inv));%<1 if Dmax<Dmax_pack

   if DO_TEST
      plot(Dmax,output);
   end

   return

end
