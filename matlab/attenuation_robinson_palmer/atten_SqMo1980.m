function [T,att,err,floe_density] = atten_SqMo1980(cice)
%%points from Squire & Moore (1980)

T     = [12.2,9.4,7.6,6.4,5.5];
att   = [.272 .438 .855 1.087 1.214]*1e-4;%m^{-1}
err   = [54 36 49 37 192]*1e-7;%m^{-1}

if nargin==0
   cice  = .5;%%concentration - not in paper
end

d_av           = 50;%% [m] approximate from Fig 1
floe_density   = cice/d_av;%%  floes per metre (1d model)
%att_nd   = att/floe_density;%% attenuation per floe
