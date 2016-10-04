function [T,att,err,floe_density] = atten_WaEtal1988(cice)
%% points from Wadhams et al (1988) - also cf Bennetts et al (2010, JGR)
%% thickness = 3.1m

T     = [14.03,11.88,10.31,9.10,8.14];
att   = [.29,    .73, 1.23,2.01,2.66]*1e-4;%m^{-1}
err   = [29 ,     25,   19,  17,  22]*1e-6;%m^{-1}

if nargin==0
   cice  = .3;%% concentration
end

d_av           = 65;%% [m]
floe_density   = cice/d_av;%%  floes per metre (1d model)
%att_nd   = att/floe_density;%% attenuation per floe
