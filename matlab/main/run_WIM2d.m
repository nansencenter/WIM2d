function out_fields = run_WIM2d(params_in,grid_prams,...
                                ice_fields,wave_fields,wave_stuff)
%% CALL: out_fields = run_WIM2d(params_in,grid_prams,...
%%                              ice_fields,wave_fields,wave_stuff)
%%
%% ALL INPUTS ARE OPTIONAL
%% *params_in is a structure that can be used to overwrite particular parameters
%% from infile_matlab.txt (eg in a loop over 1 parameter)
%% *ice_fields  = structure eg:
%%      cice: [51x51 double]
%%      hice: [51x51 double]
%%      Dmax: [51x51 double]
%%  WTR_MASK: [51x51 logical]
%%  ICE_MASK: [51x51 double]
%% *wave_fields  = structure eg:
%%         Hs: [150x10 double]
%%         Tp: [150x10 double]
%%        mwd: [150x10 double]
%%  WAVE_MASK: [150x10 logical]
%% *wave_stuff = structure eg
%%        nfreq: 25
%%         ndir: 16
%%         freq: [25x1 double]
%%         dirs: [16x1 double]
%%     dir_spec: [150x20x16x2 double]
%%
%% OUTPUTS:
%% *out_fields = structure eg
%%  tau_x: [150x20 double]
%%  tau_y: [150x20 double]
%%   Dmax: [150x20 double]
%%     Hs: [150x20 double]
%%     Tp: [150x20 double]
%% *wave_stuff - like the input, but modified by WIM2d.m

infile         = 'infile_matlab.txt';
infile_version = 6;%%latest infile version

if ~exist(infile)
   %% now need infile to run code
   error([infile,' not present - get example from "matlab/main/infiles" directory'])
else
   disp('********************************************************')
   disp('reading options from infile:')
   disp(infile)
   disp('********************************************************')
   disp(' ')
   fid   = fopen(infile);

   %%check infile version:
   infile_version_   = read_next(fid);
   if infile_version_~=infile_version
      error(['Infile version number is: ',num2str(infile_version_),' - should be: ',num2str(infile_version)]);
   end

   %%read in rest of variables:
   while ~feof(fid)
      [x,name] = read_next(fid);
      if ~isempty(x)
         cmd   = ['params.',name,' = ',num2str(x),';'];
         disp(cmd);
         eval(cmd);
      end
   end
   fclose(fid);
end

%%Override parameters in infile with these
if exist('params_in','var')
   if ~isempty(params_in)
      fnames   = fieldnames(params_in);
      for loop_i  = 1:length(fnames)
         vbl   = fnames{loop_i};
         eval(['params.',vbl,' = ','params_in.',vbl]);
      end
   end
   clear params_in;
end

TEST_INC_SPEC  = 0;

HAVE_GRID   = 0;
HAVE_ICE    = 0;
HAVE_WAVES  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.MEX_OPT>0
   %% mex functions need to have pre-saved grid,
   %% already compiled into fortran code
   if exist('grid_prams')
      if ~isempty(grid_prams)
         error('Using mex function (MEX_OPT>0) but trying to pass in grid')
      end
   end

   infile_dirs = 'infile_dirs.txt';
   if ~exist(infile_dirs)
      error([infile_dirs,' not present (needed by mex functions)'])
   else
      disp(' ');
      disp('*******************************************************');
      disp(['Reading ',infile_dirs,'...']);
      fid      = fopen(infile_dirs);
      indir    = strtrim(fgets(fid));
      outdir   = strtrim(fgets(fid));
      fclose(fid);
   end

   %% get grid
   disp(' ');
   disp(['Getting grid from ',indir,'...']);
   grid_prams     = fn_get_grid(indir);
   grid_prams.x0  = min(grid_prams.X(:));
   grid_prams.y0  = min(grid_prams.Y(:));

   %% get nw,ndir,Tmin,Tmax from wave_info.h
   fdir  = '../../fortran';
   hfil  = [fdir,'/header_files/wave_info.h'];
   disp(' ');
   disp(['Getting wave grid from ',hfil,'...']);
   fid         = fopen(hfil);
   lin         = strsplit(strtrim(fgets(fid)));
   params.ndir = find_num(lin{5});
   lin         = strsplit(strtrim(fgets(fid)));
   params.nw   = find_num(lin{5});
   lin         = strsplit(strtrim(fgets(fid)));
   lin         = strsplit(strtrim(fgets(fid)));
   params.Tmin = find_num(lin{5});
   lin         = strsplit(strtrim(fgets(fid)));
   params.Tmax = find_num(lin{5});
   fclose(fid);

   disp(' ');
   disp(['Will save results in ',outdir]);
   disp('*******************************************************');
   disp(' ');

   %% create dirs if necessary;
   prep_mex_dirs(outdir);

   if 0
      %% get initial conditions from fortran run
      [grid_prams,ice_fields,wave_fields] = fn_check_init(outdir);
      ice_prams.c          = 'given';
      ice_prams.h          = 'given';
      ice_prams.Dmax       = 'given';
      ice_prams.break_opt  = 0;
      if ~isnan(young)
         ice_prams.young_opt  = NaN;
      else
         ice_prams.young_opt  = 1;
      end
      if ~isnan(visc_rp)
         ice_prams.visc_rp = visc_rp;
      end
   end
elseif exist('grid_prams','var');
   if ~isempty(grid_prams)
      HAVE_GRID   = 1;
      fnames      = {'x0','y0','nx','ny','dx','dy'};
      for loop_i=1:length(fnames)
         fname = fnames{loop_i};
         eval(['params.',fname,' = grid_prams.',fname,';']);
      end
   end
end

if HAVE_GRID==0
   fnames      = {'x0','y0','nx','ny','dx','dy'};
   for loop_i=1:length(fnames)
      fname = fnames{loop_i};
      eval(['grid_prams.',fname,' = params.',fname,';']);
   end
   grid_prams  = get_grid(grid_prams,params.OPT);
   %% grid_prams  = structure, eg:
   %%        x0: 0
   %%        y0: 0
   %%        nx: 51
   %%        ny: 51
   %%        dx: 4000
   %%        dy: 4000
   %%         X: [51x51 double]
   %%         Y: [51x51 double]
   %%  LANDMASK: [51x51 double]
   %%      scuy: [51x51 double]
   %%      scvx: [51x51 double]
   %%      scp2: [51x51 double]
   %%     scp2i: [51x51 double]
end

%% ice
if exist('ice_fields','var')
   if ~isempty(ice_fields)
      HAVE_ICE    = 1;

      %% check size of arrays
      fnames   = {'cice','hice','Dmax','WTR_MASK','ICE_MASK'};
      for loop_i=1:length(fnames)
         fname = fnames{loop_i};
         eval(['tst_arr = params.',fname,';']);
         check_size(tst_arr,[params.nx,params.ny]);
      end
   end
end

if HAVE_ICE==0
   %%
   ice_fields  = iceinit(params,grid_prams);
   %% ice_fields  = structure eg:
   %%      cice: [51x51 double]
   %%      hice: [51x51 double]
   %%      Dmax: [51x51 double]
   %%  WTR_MASK: [51x51 logical]
   %%  ICE_MASK: [51x51 double]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%waves
if exist('wave_fields','var');
   if ~isempty(wave_fields)
      HAVE_WAVES  = 1;

      %% check size of arrays
      fnames   = {'Hs','Tp','mwd'};
      for loop_i=1:length(fnames)
         fname = fnames{loop_i};
         eval(['tst_arr = params.',fname,';']);
         check_size(tst_arr,[params.nx,params.ny]);
      end
   end
end

if HAVE_WAVES==0
   wave_prams  = struct('Hs' ,params.Hs_init,...
                        'Tp' ,params.T_init,...
                        'mwd',params.dir_init,...
                        'OPT',params.OPT);

   wave_fields = waves_init(grid_prams,wave_prams,ice_fields);
   %% structure eg:
   %%         Hs: [150x10 double]
   %%         Tp: [150x10 double]
   %%        mwd: [150x10 double]
   %%  WAVE_MASK: [150x10 logical]
end

if ~exist('wave_stuff','var')
   wave_stuff  = set_incident_waves(grid_prams,wave_fields,params);
   %% wave_stuff = structure eg
   %%        nfreq: 1
   %%         ndir: 16
   %%         freq: 0.0833
   %%         dirs: [16x1 double]
   %%     dir_spec: [150x20x16 double]
end

if TEST_INC_SPEC==1
   nw       = wave_stuff.nfreq;     %% number of frequencies
   om_vec   = 2*pi*wave_stuff.freq; %% radial freq
   ndir     = wave_stuff.ndir;      %% number of directions
   wavdir   = wave_stuff.dirs;      %% wave from, degrees, clockwise
   Sdir     = wave_stuff.dir_spec;  %% initial directional spectrum

   disp(' ');
   disp('Testing initial spectrum...');
   [wf.Hs,wf.Tp,wf.mwd] =...
      fn_spectral_integrals(om_vec,wavdir,Sdir);

   vbls  = fieldnames(wf);
   for n=1:length(vbls)
      vbl   = vbls{n};
      disp(' ');
      disp(['comparing field: ',vbl]);
      v1    = wf.(vbl);
      v2    = wave_fields.(vbl);
      diff  = abs(v2-v1);
      disp(['max diff: ',num2str(max(diff(:)))]);
      disp(' ');
   end

   return
end

[out_fields,wave_stuff] = WIM2d(params,grid_prams,ice_fields,wave_fields,wave_stuff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,name]  = read_next(fid)
%% read next line in text file

lin   = strtrim(fgets(fid));  %% trim leading white space
lin2  = strsplit(lin);        %%split using spaces
x     = lin2{1};              %%get 1st thing in line

if strcmp(x,'')
   % blank line
   disp(' ');
   x     = [];
   name  = [];
elseif strcmp(x,'#')
   % comment
   disp(lin);
   x     = [];
   name  = [];
else
   % proper variable
   x     = str2num(x);
   name  = lin2{3};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=find_num(txt)

x0 = strsplit(txt,'!');
x  = str2num(x0{1});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prep_mex_dirs(outdir)
%% mex function saves results as binary files:
%% needs some directories to exist - otherwise it crashes

if ~exist(outdir,'dir')
   eval(['!mkdir ',outdir]);
end

odirs = {[outdir,'/binaries'],...
         [outdir,'/binaries/prog'],...
         [outdir,'/log'],...
         [outdir,'/figs'],...
         [outdir,'/figs/prog']};
for j=1:length(odirs)
   outdir   = odirs{j};
   if ~exist(outdir,'dir')
      eval(['!mkdir ',outdir]);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_size(tst_arr,sz0);

sz = size(tst_arr);
ss = ['Array is wrong size (should be',...
      num2str(sz0(1)),'x',num2str(sz0(2)),')'];
if length(sz)~=length(sz0)
   error(ss);
end
for j=1:length(sz)
   if sz(j)~=sz0(j)
      error(ss);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
