%% set_FSD_prams.m
%% Author: Timothy Williams
%% Date: 20161213

function params = set_FSD_params()

params   = ...
   struct('Dmax_pack',         300.,... % Dmax when ice is unbroken [m]
          'Dmax_pack_thresh',  400.,... % Dmax should be changed if it grows larger than this [m]
          'Dmax_min',          20. ,... % Dmax should be changed if it grows smaller than this [m]
          'cice_min',          0.05,... % min conc
          'fragility',         .9,...   % fragility for FSD
          'xi',                2,...    % no of pieces from each breaking for FSD
          'Dthresh',           200);    % change from power law to uniform FSD above this value [m]
          
