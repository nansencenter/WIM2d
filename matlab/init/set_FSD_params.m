%% set_FSD_prams.m
%% Author: Timothy Williams
%% Date: 20161213

function params = set_FSD_params()

params.Dmax_pack        = 300;   % Dmax when ice is unbroken [m]
params.Dmax_pack_thresh = 400;   % Dmax should be changed if it grows larger than this [m]
params.Dmax_min         = 20;    % Dmax should be changed if it grows smaller than this [m]
params.cice_min         = 0.05;  % min conc
params.fragility        = .9;    % fragility for FSD
params.xi               = 2;     % no of pieces from each breaking for FSD
params.Dthresh          = 200;   % change from power law to uniform FSD above this value [m]
%params.FSD_OPT          = 1;     % model for FSD - 0: RG; 1: smooth RG
