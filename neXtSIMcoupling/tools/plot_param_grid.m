function figname = plot_param_grid(param,simul_outfile,varargin)
%% CALL: plot_param_grid(param,simul_outfile,domain,colormap_name,manual_axis_range,
%%          box_bound,figure_format,visible,textstring,remove_outer)



% plot_param  plots the variable param of the simul_outfile with the given colormap. 
%   plot_param_v2 was created to avoid issues with global variables and can be called directly from the command line.
%   plot_param_v2 calls plot_tricorner for the actual plotting.
%   The necessary domain file to be loaded is detected by comparing the simul_outfile dyting with a list of known domains. 
%   plot_param_v2 requires the domian file and arctic_coasts_light.mat file to be in the matlab path.
%   Example to plot the concentraiton after 5 steps:
%   >> plot_param_v2('c','path/to/data/simul_out_bigarctic10km_test_step5.mat')
%   >> plot_param_v2('c','path/to/data/simul_out_bigarctic10km_test_step5.mat','bigarctic')
%   >> plot_param_v2('c','path/to/data/simul_out_bigarctic10km_test_step5.mat','bigarctic','rev_gris')
%   >> plot_param_v2('c','path/to/data/simul_out_bigarctic10km_test_step5.mat','bigarctic','rev_gris',[0 1])
%   >> plot_param_v2('c','path/to/data/simul_out_bigarctic10km_test_step5.mat','bigarctic','rev_gris',[0 1],[-2000,-1000,-300,700])
%   >> plot_param_v2('c','path/to/data/simul_out_bigarctic10km_test_step5.mat','bigarctic','rev_gris',[0 1],[-2000,-1000,-300,700],'pdf')
%   >> plot_param_v2('c','path/to/data/simul_out_bigarctic10km_test_step5.mat','bigarctic','rev_gris',[0 1],[-2000,-1000,-300,700],'pdf',1)
%   Created 2014-07-01 by Philipp 

% default values
figname  = '';
%domain   = 'bigarctic'; 
domain   = ''; 
colormap_name='jet';
manual_axis_range = []; 
figure_format=''; % with '' the plot is not saved 
box_bound=[];
date_flag=1; %Set to 1 to plot the date, 0 to disable Phil
plot_coastline  =   1; %0: A coarse Arctic mesh is plotted which is loaded from artic_coasts_light.mat
                       %1: The actual domain boundaries are plotted, black for closed and white for open
visible     = 1;
textstring  = '';
box on;

remove_outer   = 0;
nVarargs = length(varargin);
if nVarargs >= 1, domain = varargin{1}; end
if nVarargs >= 2, colormap_name = varargin{2}; end
if nVarargs >= 3, manual_axis_range = varargin{3}; end
if nVarargs >= 4, box_bound = varargin{4}; end
if nVarargs >= 5, figure_format = varargin{5}; end
if nVarargs >= 6, visible = varargin{6}; end
if nVarargs >= 7, textstring = varargin{7}; end
if nVarargs >= 8, remove_outer = varargin{8}; end
if nVarargs >= 9, error('Too many inputs'), end

% %Detecting domain name from simul_outfile string
% str1     = strsplit(simul_outfile,'_');
% str2     = {'simplesquaresplit2', 'MITgcmsplit2', 'topazrefined', 'topazsplit2', 'topazsplit4', 'topazsplit8', 'arctic50km',     'bigkara20km', 'kara1-5km',  'kara4-5km',  'arctic10km',  'arctic15km', 'bigarctic50km',  'bigkara2km',   'kara15km',  'kara7-5km', 'bigarctic10km',  'kara45km ', 'bk4km', 'bk2km', 'bk1km','bk8km','smallkara1km','smallkara2km','smallkara4km','smallkara8km','smallkara16km','square1km','squarebig1km','squarebig10km','squaresmall10km','squarebig10km','squaresmall2km','squaresmall1km','bigarctic-5km','squaresmall5km','tinykara2km','tinykara1km'};
% meshfile = str1{ismember(str1,str2)};

%Loading data
load arctic_coasts_light.mat;
if isstr(simul_outfile)
   load(simul_outfile);
elseif isstruct(simul_outfile)
   simul_out   = simul_outfile;
else
   error('pass in filename as string or simul_out structure')
end

% Load/read mesh
X  = simul_out.wim.gridprams.X(:,1)/1.e3;%km
Y  = simul_out.wim.gridprams.Y(1,:)/1.e3;%km

% for param == 'speed' we plot the speed in km/day
tstr  = param;
if strcmp(param,'taux_waves')
   Z  = log10(simul_out.wim.waves_for_nodes.tau_x);
   tstr  = 'log_{10}({\it\tau_x}), Pa';
elseif strcmp(param,'tauy_waves')
   Z     = log10(simul_out.wim.waves_for_nodes.tau_y);
   tstr  = 'log_{10}({\it\tau_y}), Pa';
elseif strcmp(param,'Hs')
   Z     = simul_out.wim.wave_fields.Hs;
   tstr  = '{\itH}_{s}, m';
   if nVarargs < 2, colormap_name = 'gray2red'; end
elseif strcmp(param,'Tp')
   Z     = simul_out.wim.wave_fields.Tp;
   tstr  = '{\itT}_{p}, s';
elseif strcmp(param,'mwd')
   Z     = simul_out.wim.wave_fields.mwd;
   tstr  = 'mwd, degrees';
elseif strcmp(param,'cice')|strcmp(param,'c')
   Z     = simul_out.wim.ice_on_grid.cice;
   tstr  = 'concentration';
   if nVarargs < 2, colormap_name = 'rev_gris'; end
elseif strcmp(param,'hice')|strcmp(param,'h')
   Z     = simul_out.wim.ice_on_grid.hice;
   tstr  = 'thickness, m';
elseif strcmp(param,'thick')
   c        = simul_out.wim.ice_on_grid.cice;
   Z        = 0*c;
   jpos     = find(c>0);
   Z(jpos)  = simul_out.wim.ice_on_grid.hice(jpos)./c(jpos);
   tstr     = 'thickness, m';
   clear c;
elseif strcmp(param,'Dmax')
   Z     = Nfloes_to_Dmax( simul_out.wim.ice_for_elements.Nfloes,...
                           simul_out.wim.ice_on_grid.cice,...
                           simul_out.wim.other_prams);
   tstr  = '{\itD}_{max}, m';
   clear Nfloes;
elseif strcmp(param,'Nfloes')
   Z     = 1e6*simul_out.wim.ice_for_elements.Nfloes;
   tstr  = '{\itN}_{floes}, km^{-2}';
end

%%make sure things are 0 on land
Z  = Z.*(1-simul_out.wim.gridprams.LANDMASK);

if 1
   labs  = {'\itx, \rmkm','\ity, \rmkm',tstr};
else
   labs  = {'\itx, \rmkm','\ity, \rmkm',[]};
end
P  = fn_pcolor(X,Y,Z,labs,manual_axis_range,remove_outer);
colormap(colormap_name);



if nVarargs < 6
end
if date_flag == 1 & isfield(simul_out,'current_time') & strcmp(textstring,'')==0
   textstring = datestr(simul_out.current_time,'yyyy/mm/dd   HH:MM');
end

if strcmp(figure_format,'eps')
   ss = strsplit(simul_outfile,'.mat');
   figname  = [ss{1},'_',param,'.',figure_format];
   disp(['saving to ',figname]);
   print(gcf,figname,'-painters','-r300','-depsc');
elseif strcmp(figure_format,'tiff')
   ss = strsplit(simul_outfile,'.mat');
   figname  = [ss{1},'_',param,'.',figure_format];
   disp(['saving to ',figname]);
   print(gcf,figname,'-dtiff','-r300');
   %saveas(gcf,figname);
elseif ~strcmp(figure_format,'')
   ss = strsplit(simul_outfile,'.mat');
   figname  = [ss{1},'_',param,'.',figure_format];
   disp(['saving to ',figname]);
   saveas(gcf,figname);
end
