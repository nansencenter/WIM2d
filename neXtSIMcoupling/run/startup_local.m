%%define location of /Data/sim:
SHOW_WARNING   = 0;
if exist('/Volumes/sim')
   %% johansen
   data_sim = '/Volumes/sim'
elseif exist('/Volumes/Tim_Ext_HD2/WORK/neXtSIM')
   %% external hard drive
   data_sim = '/Volumes/Tim_Ext_HD2/WORK/neXtSIM'
else
   SHOW_WARNING   = 1;
   %%
   warnings = {'No paths to neXtSIM data present';
               '- eg may need to connect to johansen (/Data/sim) with cmd+k';
               '  or attach external HD.'};
end
rmpaths;

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
%% add wim2d paths
%%
wimdir1     = [gitdir,'/WIM2d/fortran'];
wimdirs{1}  = [wimdir1,'/bin'];%%mex funs
wimdirs{2}  = [wimdir1,'/run'];%%interface to mex funs
wimdirs{3}  = [wimdir1,'/matlab_funs'];
wimdirs{4}  = [gitdir,'/matlab/Semi-Infinite-Elastic-Plate/GEN_progs'];
for loop_i=1:length(wimdirs)
   addpath(wimdirs{loop_i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('data_sim','var')

   johansen_paths = [data_sim,'/data'];%%+all subdirs
   topaz_path     = [johansen_paths,'/TOPAZ4/200709_201102'];%topaz data
   amsre_path     = [johansen_paths,'/AMSRE_ice_conc/2008/mar'];%ice conc
   etopo_path     = [johansen_paths,'/BATHYMETRY/etopo1_ice_c_i2'];%%bathymetry

   %% add paths
   joh_dirs = {topaz_path,...
               amsre_path,...
               etopo_path};
   for loop_i=1:length(joh_dirs)
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
   addpath([mat_path,'/SuiteSparse/CHOLMOD/MATLAB/']);
   addpath([mat_path,'/m_map/']);

   %% get all the tools,data:
   %% - need all the sub-directories in these folders
   lookin_dirs{2} = johansen_paths;
   lookin_dirs{3} = [mat_path,'/m_map/'];
end

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

if SHOW_WARNING==1
   disp(' ');
   ll = '***********************************************************************';
   disp(ll);
   for j=1:length(warnings)
      disp(warnings{j});
   end
   disp(ll);
   disp(' ');
end

clear;