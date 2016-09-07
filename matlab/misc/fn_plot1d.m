function H  = fn_plot1d(x,y,labs,col,STEPPING)
%% fn_pcolor.m
%% Author: Timothy Williams
%% Date: 20150820, 13:32:28 CEST

if ~exist('labs','var')
   labs  = {'\itx, \rmkm',' '};
end

if ~exist('col','var')
   col   = '-';
end

if ~exist('STEPPING')
   STEPPING = 1;
end

if STEPPING
   [x,y]  = step_1d(x,y);
end

H  = plot(x,y,col);

GEN_proc_fig(labs{1},labs{2});
