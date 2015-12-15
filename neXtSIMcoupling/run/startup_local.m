%%define location of /Data/sim:
if exist('data_sim','var');
   clear data_sim
end

SHOW_WARNING   = 0;
try_extHD      = 1;
try_local      = 1;

data_locs   = {'/Volumes/sim',...
               '/Volumes/Tim_Ext_HD2/WORK/neXtSIM',...
               '../local_data'};

for loop_i=1:length(data_locs)
   dloc     = data_locs{loop_i+1};
   %%
   disp(['Checking for ',dloc]);
   if exist(dloc)
      dd = dir(dloc);
      if length(dd)>2
         %% sometimes dir can show up even if it's not loaded
         %% more than only ".", ".." in dir
         data_sim = dloc;
         disp(['Using SIM data from',data_sim]);
         break;
      end
   end
end

if ~exist('data_sim','var')
   disp({'No path to SIM data';
         '- may need to:';
         ' > connect to johansen (/Data/sim) with cmd+k';
         ' > attach external HD';
         ' > mkdir ../local_data and add/link basic data:';
         '     mesh, m_map'});
   return;
end

rmpaths;

%% if we need to make a new grid
addpath('../grid_setup');
addpath('../tools');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define neXtSIM paths;
gitdir   = getenv('GIT_REPOS');
if strcmp(gitdir,'')
   error('$GIT_REPOS not set - launch matlab from terminal')
end
xsimdir        = [gitdir,'/neXtSIM'];
nextsim_path   = [xsimdir,'/neXtSIM-trunk-SourceTree'];
bamg_path      = [xsimdir,'/ISSM-trunk-jpl-svn/lib'];

%% add local paths
local_dirs     = {nextsim_path,...
                  [nextsim_path,'/code'],...
                  [nextsim_path,'/WIM/code'],...
                  [nextsim_path,'/tools'],...
                  bamg_path};

for loop_i=1:length(local_dirs)
   addpath(local_dirs{loop_i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% need to compile SuiteSparse3 locally (neXtSIM solvers)
ss3   = [getenv('HOME'),'/MATHS/programs/matlab/SuiteSparse3']
addpath([ss3,'/CHOLMOD/MATLAB/']);
addpath([ss3,'/AMD/MATLAB/']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add wim2d paths
%%
wimdir1        = [gitdir,'/WIM2d/fortran'];
wimdirs{1}     = [wimdir1,'/bin'];%%mex funs
wimdirs{end+1} = [wimdir1,'/run'];%%interface to mex funs
wimdirs{end+1} = [wimdir1,'/matlab_funs'];
wimdirs{end+1} = [wimdir1,'/../matlab/other_deps'];
wimdirs{end+1} = [gitdir,'/matlab/Semi-Infinite-Elastic-Plate/GEN_progs'];
wimdirs{end+1} = [gitdir,'/WIM2d/matlab/misc'];
for loop_i=1:length(wimdirs)
   addpath(wimdirs{loop_i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rest of SIM data
data_paths = [data_sim,'/data'];%%+all subdirs
topaz_path     = [data_paths,'/TOPAZ4/198910_201112'];%topaz data
%topaz_path     = [data_paths,'/TOPAZ4/198910_201312'];%topaz data
amsre_path     = [data_paths,'/AMSRE_ice_conc/2008'];%ice conc
etopo_path     = [data_paths,'/BATHYMETRY/etopo1_ice_c_i2'];%%bathymetry

%% add paths
joh_dirs = {topaz_path,...
            amsre_path,...
            etopo_path};
for loop_i=1:length(joh_dirs)
   disp(joh_dirs{loop_i})
   addpath(joh_dirs{loop_i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% other "MATLAB" paths
mat_path = [data_sim,'/MATLAB'];
addpath([mat_path]);
addpath([mat_path,'/age/']);
addpath([mat_path,'/defo_rgps/']);
addpath([mat_path,'/from An/']);
addpath([mat_path,'/m_map/']);


%% get all the tools,data:
%% - need all the sub-directories in these folders
lookin_dirs{2} = data_paths;
lookin_dirs{3} = [mat_path,'/m_map/'];

lookin_dirs{1} = [nextsim_path,'/tools'];

for j=1:length(lookin_dirs)

   lookin   = lookin_dirs{j};

   subdirs     = dir(lookin);
   isub        = [subdirs(:).isdir]; %# returns logical vector
   nameFolds   = {subdirs(isub).name}';

   nameFolds(ismember(nameFolds,{'.','..'})) = [];

   for loop_i=1:length(nameFolds)
      addpath([lookin '/' nameFolds{loop_i}]);
   end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showpaths;

ll = '*****************************************************************';
disp(' ');
disp(ll);
disp('Remember to open "Parallel>Manage Cluster Profiles..." menu,');
disp('then run "matlabpool".');
disp(ll);
disp(' ');

clear;
