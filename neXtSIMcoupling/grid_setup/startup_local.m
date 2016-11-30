%%define location of /Data/sim:
if exist('data_sim','var');
   clear data_sim
end

SHOW_WARNING   = 0;
try_extHD      = 1;
try_local      = 1;

gitdir      = getenv('GIT_REPOS');
w2d_path    = [gitdir,'/WIM2d/'];
data_locs   = {'/Volumes/sim',...
               '/Volumes/Tim_Ext_HD2/WORK/neXtSIM',...
               [w2d_path,'/neXtSIMcoupling/local_data']};

for loop_i=1:length(data_locs)
   dloc     = data_locs{loop_i};
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
   txt   = {'No path to SIM data';...
         '- may need to:';...
         ' > connect to johansen (/Data/sim) with cmd+k';...
         ' > attach external HD';...
         [' > mkdir ',data_locs{3},' and add/link basic data:'];...
         '     mesh, m_map'}
   error(txt);
end

rmpaths;
ALLDIRS  = {};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define neXtSIM paths;

if strcmp(gitdir,'')
   error('$GIT_REPOS not set - launch matlab from terminal')
end
xsimdir        = [gitdir,'/neXtSIM'];
nextsim_path   = [xsimdir,'/neXtSIM-trunk-SourceTree'];
bamg_path      = [xsimdir,'/ISSM-trunk-jpl-svn/lib'];

%% add local paths
ALLDIRS(1:6)   = {nextsim_path,...
                  [nextsim_path,'/code'],...
                  [nextsim_path,'/bin'],...
                  [nextsim_path,'/code/WIM'],...
                  [nextsim_path,'/tools'],...
                  bamg_path};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% need to compile SuiteSparse3 locally (neXtSIM solvers)
ss3   = [getenv('HOME'),'/MATHS/programs/matlab/SuiteSparse3'];
ALLDIRS{end+1} = [ss3,'/CHOLMOD/MATLAB/'];
ALLDIRS{end+1} = [ss3,'/AMD/MATLAB/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add wim2d paths
ALLDIRS{end+1} = 'bin';%mex functions
%ALLDIRS{end+1} = [w2d_path,'/neXtSIMcoupling/run'];%scripts to launch function
ALLDIRS{end+1} = [w2d_path,'/neXtSIMcoupling/grid_setup'];%scripts to launch function
ALLDIRS{end+1} = [w2d_path,'/neXtSIMcoupling/tools'];%scripts to launch function
ALLDIRS{end+1} = [w2d_path,'/matlab/main'];%%interface to mex funs
%%
wimdir1        = [gitdir,'/WIM2d/fortran'];
ALLDIRS{end+1} = [wimdir1,'/matlab_funs'];
ALLDIRS{end+1} = [wimdir1,'/../matlab/other_deps'];
ALLDIRS{end+1} = [gitdir,'/matlab/Semi-Infinite-Elastic-Plate/GEN_progs'];
ALLDIRS{end+1} = [gitdir,'/WIM2d/matlab/misc'];

if 1
   %%if want to run pure matlab code (not mex)
   ALLDIRS{end+1} = [w2d_path,'/matlab/misc'];
   ALLDIRS{end+1} = [w2d_path,'/matlab/main'];
   ALLDIRS{end+1} = [w2d_path,'/matlab/init'];
   ALLDIRS{end+1} = [w2d_path,'/matlab/numerics'];
   ALLDIRS{end+1} = [w2d_path,'/matlab/boltzmann'];
   ALLDIRS{end+1} = [w2d_path,'/matlab/attenuation_youngs'];
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
ALLDIRS(end+1:end+3) = {topaz_path,...
                        amsre_path,...
                        etopo_path};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% other "MATLAB" paths
mat_path = [data_sim,'/MATLAB'];
ALLDIRS{end+1} = [mat_path];
ALLDIRS{end+1} = [mat_path,'/age/'];
ALLDIRS{end+1} = [mat_path,'/defo_rgps/'];
ALLDIRS{end+1} = [mat_path,'/from An/'];
ALLDIRS{end+1} = [mat_path,'/m_map/'];


%% get all the tools,data:
%% - need all the sub-directories in these folders
lookin_dirs{2} = data_paths;
lookin_dirs{3} = [mat_path,'/m_map/'];
lookin_dirs{1} = [nextsim_path,'/tools'];

for j=1:length(lookin_dirs)

   lookin      = lookin_dirs{j};
   subdirs     = dir(lookin);
   isub        = [subdirs(:).isdir]; %# returns logical vector
   nameFolds   = {subdirs(isub).name}';

   nameFolds(ismember(nameFolds,{'.','..'})) = [];
   for loop_i=1:length(nameFolds)
      subdir   = nameFolds{loop_i};
      if ~strcmp(subdir,'private')
         ALLDIRS{end+1} = [lookin '/' subdir];
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bad   = {};
for loop_i=1:length(ALLDIRS)
   dd = ALLDIRS{loop_i};
   if exist(dd,'dir')
      addpath(dd);
   else
      bad{end+1,1}  = dd;
   end
end

ll = '*****************************************************************';
disp(' ');
disp(ll);
showpaths;
disp(ll);
disp(' ');
disp(ll);
disp('Directories not present:');
for j=1:length(bad)
   disp(bad{j});
end
disp(ll);

disp(' ');
disp(ll);
disp('Remember to open "Parallel>Manage Cluster Profiles..." menu,');
disp('then run "matlabpool".');
disp(ll);
disp(' ');

clear;
