function [ice_fields,wave_fields] = WIM2d()

%DO_SAVE     = 0;
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
         cmd   = [name,' = ',num2str(x),';'];
         disp(cmd);
         eval(cmd);
      end
   end
   fclose(fid);
end

if (ADV_DIM==2)&(ny<4)
   error('incompatible values of ADV_DIM and ny: increase ny or use ADV_DIM=1');
end

%% TURN ON/OFF PLOTTING:
PLOT_OPT    = 2;%%plot option (only used if doing plotting)
adv_options = struct('ADV_DIM',ADV_DIM,...
                     'ADV_OPT',ADV_OPT);
TEST_INC_SPEC     = 0;
TEST_FINAL_SPEC   = 0;

%% make a log file similar to fortran file
log_dir  = 'log';
if ~exist('log','dir')
   eval(['!mkdir ',log_dir]);
end
log_file    = [log_dir,'/WIM2d_matlab.log'];
this_subr   = mfilename();

%%open log file for writing (clear contents)
logid       = fopen(log_file,'w');
fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s\n','Outer subroutine:');
fprintf(logid,'%s%s%s\n','>> ',this_subr,'.m');
fprintf(logid,'%s\n','***********************************************');

fprintf(logid,'%s\n',' ');
fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s\n','Main parameters:');
fprintf(logid,'%s%2.2d\n','SCATMOD:                          ',SCATMOD);
fprintf(logid,'%s%2.2d\n','ADV_DIM:                          ',ADV_DIM);
fprintf(logid,'%s%2.2d\n','ADV_OPT:                          ',ADV_OPT);
fprintf(logid,'%s%2.2d\n','DO_BREAKING:                      ',DO_BREAKING);
fprintf(logid,'%s%2.2d\n','STEADY:                           ',STEADY);
fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s\n',' ');
fclose(logid);

%if DO_SAVE
%   filename=['wim2d_out',num2str(Tm),...
%            's_',num2str(mwd_dim),...
%            'deg_spread','.mat'];
%else
%   disp('not saving');
%   %disp('not saving - push any key');
%   %pause;
%end

%%important settings
itest = 24;
jtest = 1;
format long

disp('Initialization')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%check what inputs we are given;
if MEX_OPT>0
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
   fid   = fopen(hfil);
   lin   = strsplit(strtrim(fgets(fid)));
   ndir  = find_num(lin{5});
   lin   = strsplit(strtrim(fgets(fid)));
   nw    = find_num(lin{5});
   lin   = strsplit(strtrim(fgets(fid)));
   lin   = strsplit(strtrim(fgets(fid)));
   Tmin  = find_num(lin{5});
   lin   = strsplit(strtrim(fgets(fid)));
   Tmax  = find_num(lin{5});
   fclose(fid);

   disp(' ');
   disp(['Will save results in ',outdir]);
   disp('*******************************************************');
   disp(' ');

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
end
HAVE_GRID   = exist('grid_prams','var');

HAVE_ICE    = exist('ice_fields','var');
if HAVE_ICE
   if ~isfield(ice_fields,'cice')
      ice_prams   = ice_fields;
      HAVE_ICE    = 0;
   else
      %% have ice_fields (arrays) already
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
      ice_prams   = fn_fill_iceprams(ice_prams);
   end
end

HAVE_WAVES  = exist('wave_prams','var');
if HAVE_WAVES
   %%check if full info is given or just partial info;
   if isfield(wave_prams,'dir_spec')
      wave_fields = wave_prams;
   else
      HAVE_WAVES  = 0;
   end
end

HAVE2 = HAVE_ICE+HAVE_WAVES;
if (0<HAVE2)&(HAVE2<3)
   disp(' ');
   disp('*********************************************************************');
   disp('Please specify both of ''ice_fields'' and ''wave_fields''');
   disp('(or specify none, and let configuration be chosen through local variable ''OPT'')');
   disp(['HAVE_ICE    = ',num2str(HAVE_ICE  )]);
   disp(['HAVE_WAVES  = ',num2str(HAVE_WAVES)]);
   disp('*********************************************************************');
   error('WIM2d initialisation');
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Grid
if HAVE_GRID==0

   %%can pass in all/none of grid_prams,ice_fields,wave_fields
   %%but not some
   if exist('ice_fields','var')|exist('wave_fields','var')
      disp('WIM2d.m');
      disp('Please specify all 3 of ''grid_prams'',''ice_fields'',''wave_fields''');
      disp('(at least ''grid_prams'' not given)');
      return;
   end

   grid_prams  = struct('x0',x0,'y0',y0,...
                        'nx',nx,'ny',ny,...
                        'dx',dx,'dy',dy);
   grid_prams  = get_grid(grid_prams,OPT);
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
nx       = grid_prams.nx;
ny       = grid_prams.ny;
dy       = grid_prams.dy;
dy       = grid_prams.dy;
X        = grid_prams.X;
Y        = grid_prams.Y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ice/water;
if HAVE_ICE==0
   ice_prams.c    = conc_init;
   ice_prams.h    = h_init;
   ice_prams.Dmax = Dmax_init;

   break_opt   = 0;%% breaking condition (0=beam;1=Marchenko)
   ice_prams.break_opt  = break_opt;
   if ~isnan(young)
      ice_prams.young      = young;%%Young's modulus [Pa]
      ice_prams.young_opt  = NaN;%%not used
   else
      %%Option for selecting Young's modulus
      %% - only used if Young's modulus not given
      %% 0: 2e9 [Marchenko estimate]
      %% 1: Vernon's estimate from brine volume fraction (vbf)
      %% 2: same as above but with vbf=.1 -> young = 5.49e9 Pa
      ice_prams.young_opt  = 0;
   end
   if ~isnan(visc_rp)
      ice_prams.visc_rp = visc_rp;%%Robinson-Palmer damping [Pa/(m/s)]
   end
   %%
   [ice_fields,ice_prams]  = iceinit(ice_prams,grid_prams,OPT);
   %% ice_fields  = structure:
   %%      cice: [51x51 double]
   %%      hice: [51x51 double]
   %%      Dmax: [51x51 double]
   %%  WTR_MASK: [51x51 logical]
   %%  ICE_MASK: [51x51 double]

   %% ice_prams = structure eg:
   %%               c: 0.750000000000000
   %%               h: 2
   %%            Dmax: 300
   %%           young: 2.000000000000000e+09
   %%          bc_opt: 0
   %%         visc_rp: 13
   %%          rhowtr: 1.025000000000000e+03
   %%          rhoice: 9.225000000000000e+02
   %%               g: 9.810000000000000
   %%         poisson: 0.300000000000000
   %%             vbf: 0.100000000000000
   %%              vb: 100
   %%         sigma_c: 2.741429878818372e+05
   %%        strain_c: 1.370714939409186e-04
   %%  flex_rig_coeff: 1.831501831501831e+08
   %%            Dmin: 20
   %%              xi: 2
   %%       fragility: 0.900000000000000
end

cice     = ice_fields.cice;
hice     = ice_fields.hice;
Dmax     = ice_fields.Dmax;
WTR_MASK = ice_fields.WTR_MASK;
ICE_MASK = ice_fields.ICE_MASK;

%% add wave stress computation
ice_fields.tau_x  = 0*cice;
ice_fields.tau_y  = 0*cice;

%% Model Parameters
rho_wtr  = ice_prams.rhowtr;   % Water density   [m/s^2]
g        = ice_prams.g;        % Gravity         [m/s^2]
strain_c = ice_prams.strain_c; % breaking strain [-] 
young    = ice_prams.young;    % Young's modulus [Pa]
visc_rp  = ice_prams.visc_rp;  % Robinson-Palmer coeff [Pa/(m/s)]

% Parameters for the floe size distribution
Dmin        = ice_prams.Dmin;      % [m]
xi          = ice_prams.xi;        % [-]
fragility   = ice_prams.fragility; % [-]

if 0
   fn_plot_ice(grid_prams,ice_fields);
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%append to log file
logid = fopen(log_file,'a');
fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s\n','WIM parameters:');
fprintf(logid,'%s%4.2f\n','Brine volume fraction:       ' ,ice_prams.vbf);
fprintf(logid,'%s%10.3e\n','Youngs modulus (Pa):        ' ,ice_prams.young);
fprintf(logid,'%s%10.3e\n','Flexural strength (Pa):     ' ,ice_prams.sigma_c);
fprintf(logid,'%s%10.3f\n','Breaking strain:            ' ,ice_prams.strain_c);
fprintf(logid,'%s%5.2f\n','Damping (Pa.s/m):           '  ,ice_prams.visc_rp);
fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s\n','');
fclose(logid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%waves
if HAVE_WAVES==0
   if ~exist('wave_prams','var');
      wave_prams  = struct('Hs',Hs_init,...
                           'Tp',T_init,...
                           'mwd',dir_init);
   end
   wave_fields = waves_init(grid_prams,wave_prams,ice_fields,OPT);
   %% structure eg:
   %%         Hs: [150x10 double]
   %%         Tp: [150x10 double]
   %%        mwd: [150x10 double]
   %%  WAVE_MASK: [150x10 logical]
   WAVE_MASK   = wave_fields.WAVE_MASK;
   %%
   nw,ndir,Tmin,Tmax
   inc_options = struct('nw',nw,...
                        'ndir',ndir,...
                        'Tmin',Tmin,...
                        'Tmax',Tmax,...
                        'DIRSPEC_INC_OPT',DIRSPEC_INC_OPT);
   wave_stuff  = set_incident_waves(grid_prams,wave_fields,inc_options);
end

nw       = wave_stuff.nfreq;     %% number of frequencies
om_vec   = 2*pi*wave_stuff.freq; %% radial freq
ndir     = wave_stuff.ndir;      %% number of directions
wavdir   = wave_stuff.dirs;      %% wave from, degrees, clockwise
Sdir     = wave_stuff.dir_spec;  %% initial directional spectrum

if TEST_INC_SPEC==1
   disp(' ');
   disp('Testing initial spectrum...');
   [wf.Hs,wf.Tp,wf.mwd] = fn_spectral_integrals(om_vec,wavdir,Sdir);
   vbls  = {'Hs','Tp','mwd'};
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

if STEADY==1
   S_inc    = Sdir;
   theta    = -pi/180*(90+wavdir);
   j_fwd    = find(cos(theta)>=0);
   WAVE_MASK2        = 0*WAVE_MASK;
   WAVE_MASK2(1:3,:) = 1;
   %[i1,j1]  = find(WAVE_MASK==1,1,'first');
   %J_STEADY = find(S_inc(i1:i1+2,j1:j1+,j_fwd,:)>0);%%only get "forward" directions
   % squeeze(S_inc(i1,j1,:))
   % GEN_pause;
end
%%
T  = 2*pi./om_vec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Water wavelength and wave speed
%% is a function only of wave period
wlng     = g.*T.^2./(2.*pi);
ap       = sqrt(g.*wlng./(2.*pi)); % Phase speed
ag       = ap./2;                  % Group speed

ag_eff      = zeros(nx,ny,nw);
ap_eff      = zeros(nx,ny,nw);
wlng_ice    = zeros(nx,ny,nw);
disp_ratio  = ones (nx,ny,nw);
atten_nond  = zeros(nx,ny,nw);
damping     = zeros(nx,ny,nw);

% Display some parameters here (since initialisation can be slow)
Info  = { '------------------------------------';
         ['h          = ' num2str(ice_prams.h) ' m const'];
         ['Tp         = ' num2str(wave_prams.Tp) ' s'];
         ['Hs         = ' num2str(wave_prams.Hs) ' m'];
         ['CFL        = ' num2str(CFL)];
         ['nfreq      = ' num2str(nw)];
         ['ndir       = ' num2str(ndir)];
         ['SCATMOD    = ' num2str(SCATMOD)];
         '------------------------------------';
         ' '};
disp(strvcat(Info));

for i = 1:nx
for j = 1:ny

%  if j==1
%     %%progress report - can be slow
%     disp([' - initialised ',num2str(i),' rows out of ',num2str(nx)])
%  end

   if ICE_MASK(i,j)==1
      %% if ice is present:
      %% get ice wavelengths, group velocities,
      %% displacement ratios (convert wtr disp to ice),
      %% attenuation;
      if DO_ATTEN==1
         [damping_rp,kice,kwtr,int_adm,NDprams,...
            alp_scat,modT,argR,argT] =...
               RT_param_outer(hice(i,j),om_vec,young,visc_rp);
         %%
         if CHK_ATTEN==1
            %%check with old version
         end
         atten_nond(i,j,:) = alp_scat;
         damping(i,j,:)    = damping_rp;
      else
         [damping_rp,kice,kwtr,int_adm,NDprams] =...
            RT_param_outer(hice(i,j),om_vec,young,visc_rp);
         modT  = 1;
      end

      if USE_ICE_VEL==0
         %%use wtr group vel;
         ag_eff(i,j,:)  = ag;
         ap_eff(i,j,:)  = ap;
      else
         %%TODO check if this is correct
         %%weighted avg of ice and wtr group vel;
         ag_ice         = GEN_get_ice_groupvel(hice(i,j),T,Inf,young);
         ag_eff(i,j,:)  = cice(i,j)*ag_ice+...
                           +(1-cice(i,j))*ag;

         %%weighted avg of ice and wtr phase vel;
         ap_ice         = om_vec./k_ice;
         ap_eff(i,j,:)  = cice(i,j)*ap_ice+...
                           +(1-cice(i,j))*ap;
      end

      wlng_ice(i,j,:)   = 2*pi./kice;
      disp_ratio(i,j,:) = (kice./kwtr).*modT;
      %%
%     if (i==itest)&(j==jtest)
%        disp('om_vec,T,h')
%        disp([om_vec(1),T(1),hice(i,j)])
%        disp('atten')
%        ss = [num2str(atten_nond(i,j,1),'%7.7e'),'   ',...
%              num2str(damping(i,j,1),'%7.7e')];
%        disp(ss);
%        disp('ki,kw,2pi/wlng_wtr')
%        ss = [num2str(kice,'%7.7e'),'   ',...
%              num2str(kwtr,'%7.7e'),'   '  ,...
%              num2str(2*pi./wlng,'%7.7e')];
%        disp(ss);
%        disp('lam,|T|,disp_rat')
%        ss = [num2str(wlng_ice(i,j,1),'%4.4f'),'   ',...
%              num2str(modT,'%7.7f'),'   ',...
%              num2str(disp_ratio(i,j,1),'%7.7f')];
%        disp(ss);
%        disp('argRT,s')
%        ss = [num2str(argR,'%7.7e'),'   ',...
%              num2str(argT,'%7.7e'),'   ',...
%              num2str(int_adm,'%7.7e')];
%        disp(ss);
%        %return
%     end
   else
      ag_eff(i,j,:)     = ag;
      ap_eff(i,j,:)     = ap;
      wlng_ice(i,j,:)   = wlng;
   end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amax     = max(ag_eff(:))
dt       = CFL*dx/max(ag_eff(:)); 

if 0
   nt    = 50;
   L     = nt*dt*amax
elseif 0
   L     = max(X(:))-min(X(:));
   amin  = min(ag_eff(:));
   uc    = amin+.7*(amax-amin);
   nt    = round( L/uc/dt );
else
   L     = max(X(:))-min(X(:));
   amin  = min(ag_eff(:));
   uc    = amin+.7*(amax-amin);
   %%
   nt = floor(duration_hours*3600/dt);
end
duration = nt*dt%%duration in seconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display parameters
Info  = { '------------------------------------';
         ['Young''s modulus    = ' num2str(ice_prams.young,'%5.5e') ' Pa'];
         ['sigma_c            = '  num2str(ice_prams.sigma_c,'%5.5e') ' Pa'];
         ['fragility          = '  num2str(fragility)];
         ['strain_c           = '  num2str(strain_c,'%5.5e')];
         ['h                  = '  num2str(ice_prams.h) ' m const'];
         ['c                  = '  num2str(ice_prams.c) ' const'];
         ['Damping            = '  num2str(visc_rp) ' Pa.s/m'];
         [' '];
         ['Tp                 = '  num2str(wave_prams.Tp) ' s'];
         ['Hs                 = '  num2str(wave_prams.Hs) ' m'];
         ['nfreq              = '  num2str(nw)];
         ['ndir               = '  num2str(ndir)];
         ['SCATMOD            = '  num2str(SCATMOD)];
         [' '];
         ['CFL                = '  num2str(CFL)];
         ['dt                 = '  num2str(dt,'%1.1f')];
         ['nt                 = '  num2str(nt)];
         ['Time interval      = '  num2str(duration/3600,'%1.1f') 'h'];
         [' '];
         ['nx                 = '  num2str(nx)];
         ['ny                 = '  num2str(ny)];
         ['dx                 = '  num2str(dx/1e3) ' km'];
         ['dy                 = '  num2str(dy/1e3) ' km'];
         ['x extent           = '  num2str(nx*dx/1e3) ' km'];
         ['y extent           = '  num2str(ny*dy/1e3) ' km'];
         '------------------------------------';
         ' '};

disp(strvcat(Info));

%% Integration
t0       = now;   %%days
t0_fac   = 24*60; %%days to minutes
brkcrt  = zeros(nx,ny,nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%append to log file
logid = fopen(log_file,'a');
fprintf(logid,'%s\n',' ');
fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s\n','Other Parameters:');
fprintf(logid,'%s%6.1f\n','Time step (s):                    ' ,dt);
fprintf(logid,'%s%4.3f\n','CFL number:                       ' ,CFL);
fprintf(logid,'%s%5.2f\n','Maximum wave group velocity (m/s):' ,amax);
fprintf(logid,'%s%6.1f\n','Time step (s):                    ' ,dt);
fprintf(logid,'%s%4.4d\n','Number of time steps:             ',nt);
fprintf(logid,'%s%5.2f\n','Time interval (h):                ',nt*dt/3600 );
fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s\n',' ');

fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s%4.4d%s%4.4d\n','Grid dimensions:                  ' ,...
   nx,' ',ny);
fprintf(logid,'%s%4.1f%s%4.1f\n','Spatial resolution (km):          ' ,...
   dx/1.0e3,' ',dy/1.0e3);
fprintf(logid,'%s%4.1f%s%4.1f\n','Extent of domain   (km):          ' ,...
   nx*dx/1.0e3,' ',nx*dy/1.0e3);

fprintf(logid,'%s\n',' ');
fprintf(logid,'%s%5.2f\n','Minimum period (s):               ',1/max(wave_stuff.freq) );
fprintf(logid,'%s%5.2f\n','Maximum period (s):               ',1/max(wave_stuff.freq) );
fprintf(logid,'%s%4.4d\n','Number of wave frequencies:       ',nw);
fprintf(logid,'%s%4.4d\n','Number of wave directions:        ',ndir);
fprintf(logid,'%s%5.2f\n','Directional resolution (degrees): ',360.0/(1.0*ndir) );
fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s\n',' ');
fclose(logid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define weights for numerical quadrature;
if nw>1%% weights for integral over frequency
   %% (Simpson's rule);
   wt_simp            = 2+0*om_vec;
   wt_simp([1 end])   = 1;
   wt_simp(2:2:end-1) = 4;

   %%NB om_vec needs to be equally spaced;
   dom    = abs(om_vec(2)-om_vec(1));
   wt_om  = dom/3*wt_simp;

   %if 0
   %   Hs2      = 4*sqrt(wt_om*S00');
   %   testHs   = [Hs2,Hs0]
   %   pause;
   %end
else
   wt_om = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DIAG1d==1
   cols     = {'-k','-c','-m','-r','-g','-b'};
   labs1d_1 = {'\itx, \rmkm','{\itH}_{\rm s}, m'};
   labs1d_2 = {'\itx, \rmkm','{\itH}_{\rm s}^+, m'};
   labs1d_3 = {'\itx, \rmkm','{\itH}_{\rm s}^-, m'};
   labs1d_4 = {'\itx, \rmkm','c'};
end

if PLOT_INIT
   %%
   figure(1),clf;
   fn_fullscreen;
   fn_plot_ice(grid_prams,ice_fields);
   pause(0.1);
   %%
   figure(2),clf;
   fn_fullscreen;
   Tc    = 12;%check this period
   jchq  = find(abs(T-Tc)==min(abs(T-Tc)));
   jdir  = round(ndir/2);
   s1 = struct('dir',wave_stuff.dirs(jdir),...
               'period',Tc,...
               'Sdir',Sdir(:,:,jdir,jchq));
   %%
   if PLOT_OPT==1
      fn_plot_spec(X,Y,wave_fields.Hs,wave_fields.Tp,Dmax,s1);
   else
      fn_plot_spec_2(X,Y,wave_fields.Hs,ice_fields.tau_x,...
         Dmax,ice_fields.tau_y);
   end
   %if OPT==1
   %   subplot(2,2,1);
   %   hold on;
   %   x0 = min(X(:));
   %   x1 = X(find(WAVE_MASK(:,1)==0,1,'first'),1);
   %   % {x0/1e3,x1/1e3}
   %   yc = .3*max(Y(:))/1e3;
   %   x_ = [min(X(:)),x0,x0,x1,x1,max(X(:))]/1e3;
   %   y_ = [0,0,yc*[1,1],0,0];
   %   plot(x_,y_,'k');
   %   plot(X(:,1)/1e3,wave_fields.Hs(:,1)*yc,'--k');
   %   xlabel('$x$, km','interpreter','latex','fontsize',20); 
   %   ylabel('$\hat{H}_{s}$, m','interpreter','latex','fontsize',20)
   %   hold off;
   %end
   %%
   clear s1;
   %%
   if DIAG1d==1
      %% initial
      figure(4),clf;
      fn_fullscreen;
      loop_col = 1;
      %%
      subplot(4,1,4);
      fn_plot1d(X(:,1)/1e3,ice_fields.cice(:,1),labs1d_4,cols{loop_col});
      %%
      subplot(4,1,1);
      %fn_plot1d(X(:,1)/1e3,wave_fields.Hs(:,1),labs1d_1,cols{loop_col});
      fn_plot1d(X(:,1)/1e3,mean(wave_fields.Hs,2),labs1d_1,cols{loop_col});%%average over y (columns)
      hold on;
      %%
      [Ep,Em,Et1,Et2]   = fn_split_energy(om_vec,wavdir,Sdir);
      %Hp                = 4*sqrt(Ep(:,1));
      %Hm                = 4*sqrt(Em(:,1));
      Hp                = 4*sqrt(mean(Ep,2));
      Hm                = 4*sqrt(mean(Em,2));
      if DIAG1d_OPT==0
         %Hs2   = 4*sqrt(Ep(:,1)+Em(:,1));%%add Ep + Em
         Hs2   = 4*sqrt(mean(Ep,2)+mean(Em,2));
      elseif DIAG1d_OPT==1
         %Hs2   = 4*sqrt(Et1(:,1));%%check const panel integration
         Hs2   = 4*sqrt(mean(Et1,2));
      elseif DIAG1d_OPT==2
         %Hs2   = 4*sqrt(Et2(:,1));%%check Simpson's rule integration
         Hs2   = 4*sqrt(mean(Et2,2));
      end
      fn_plot1d(X(:,1)/1e3,Hs2,labs1d_1,['-',cols{loop_col}]);
      hold on;
      %%
      subplot(4,1,2);
      fn_plot1d(X(:,1)/1e3,Hp,labs1d_2,cols{loop_col});
      hold on;
      fn_plot1d(X(:,1)/1e3,Hs2,labs1d_2,['-',cols{loop_col}]);
      hold on;
      %%
      subplot(4,1,3);
      fn_plot1d(X(:,1)/1e3,Hm,labs1d_3,cols{loop_col});
      hold on;
      %%
      loop_col = loop_col+1;
   end
   pause(0.1);
   %pause;
end

CSUM  = DO_CHECK_INIT+DO_CHECK_PROG+DO_CHECK_FINAL;
if (SV_BIN==1) & (MEX_OPT==0) & (CSUM>0)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% save some fields as binaries to m_out
   !mkdir -p  m_out
   !mkdir -p  m_out/binaries
   !rm    -rf m_out/binaries/prog
   Bdims = [nx,ny,nw,ndir];

   reps_ab  = 10;%%save every 10 time-steps if DO_CHECK_PROG==1
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% save grid files
   Fdir  = 'm_out/binaries';
   Froot = [Fdir,'/wim_grid'];
   %%
   pairs = {};
   pairs{end+1}   = {'X'         ,grid_prams.X};
   pairs{end+1}   = {'Y'         ,grid_prams.Y};
   pairs{end+1}   = {'scuy'      ,grid_prams.scuy};
   pairs{end+1}   = {'scvx'      ,grid_prams.scvx};
   pairs{end+1}   = {'scp2'      ,grid_prams.scp2};
   pairs{end+1}   = {'scp2i'     ,grid_prams.scp2i};
   pairs{end+1}   = {'LANDMASK'  ,grid_prams.LANDMASK};
   %%
   fn_save_binary(Froot,Bdims,[],pairs);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   SV_BIN   = 0;
end


if (SV_BIN==1) & (DO_CHECK_INIT==1)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% save init files
   Fdir  = 'm_out/binaries';
   Froot = [Fdir,'/wim_init'];
   %%
   pairs = {};
   pairs{end+1}   = {'cice',ice_fields.cice};
   pairs{end+1}   = {'hice',ice_fields.hice};
   pairs{end+1}   = {'Dmax',ice_fields.Dmax};
   pairs{end+1}   = {'Hs',wave_fields.Hs};
   pairs{end+1}   = {'Tp',wave_fields.Tp};
   pairs{end+1}   = {'mwd',wave_fields.mwd};
   %%
   fn_save_binary(Froot,Bdims,[],pairs);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if (SV_BIN==1) & (DO_CHECK_PROG==1)

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% 1st prog file
   !mkdir -p  m_out/binaries/prog
   Fdir  = 'm_out/binaries/prog';
   cn    = num2str(nt);
   cn(:) = '0';
   Froot = [Fdir,'/wim_prog',cn];
   %%
   pairs = {};
   pairs{end+1}   = {'Dmax',ice_fields.Dmax};
   pairs{end+1}   = {'tau_x',ice_fields.tau_x};
   pairs{end+1}   = {'tau_y',ice_fields.tau_y};
   pairs{end+1}   = {'Hs',wave_fields.Hs};
   pairs{end+1}   = {'Tp',wave_fields.Tp};
   %%
   fn_save_binary(Froot,Bdims,0,pairs);
   %eval(['!cat ',Froot,'.b'])
   %pause;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

disp('beginning main integration...');
if MEX_OPT==1

   disp(' ');
   disp('*****************************************************************');
   disp('Running fortran code with mex function: run_WIM2d_io_mex_v2');
   disp('*****************************************************************');
   disp(' ');

   % real parameters
   real_prams  = [ice_prams.young,ice_prams.visc_rp,duration,CFL];

   % integer parameters
   int_prams   = [SCATMOD,ADV_DIM,ADV_OPT,...
                  DO_CHECK_FINAL,DO_CHECK_PROG,DO_CHECK_INIT,...
                  DO_BREAKING,STEADY];

   in_arrays   = zeros(nx,ny,6);
   in_arrays(:,:,1)  = ice_fields.cice;
   in_arrays(:,:,2)  = ice_fields.hice;
   in_arrays(:,:,3)  = ice_fields.Dmax;
   in_arrays(:,:,4)  = wave_fields.Hs;
   in_arrays(:,:,5)  = wave_fields.Tp;
   in_arrays(:,:,6)  = wave_fields.mwd;
   %[min(ice_fields.cice(:)),max(ice_fields.cice(:))]
   %[min(ice_fields.hice(:)),max(ice_fields.hice(:))]
   %[min(ice_fields.Dmax(:)),max(ice_fields.Dmax(:))]
   %[min(wave_fields.Hs(:)) ,max(wave_fields.Hs(:)) ]
   %[min(wave_fields.Tp(:)) ,max(wave_fields.Tp(:)) ]
   %[min(wave_fields.mwd(:)),max(wave_fields.mwd(:))]
   %pause

   %% make the call!
   prep_mex_dirs(outdir);
   tic;
   out_arrays  = WIM2d_run_io_mex_v2(in_arrays(:),int_prams,real_prams);
   toc;

   %% extract outputs
   fldnames    = {'Dmax','tau_x','tau_y','Hs','Tp'};
   Nout        = length(fldnames);
   out_arrays  = reshape(out_arrays,[nx,ny,Nout]);
   for j=1:3
      ice_fields.(fldnames{j})   = out_arrays(:,:,j);
   end
   for j=4:5
      wave_fields.(fldnames{j})  = out_arrays(:,:,j);
   end

   % delete annoying file
   !rm -f fort.6

elseif MEX_OPT==2

   disp(' ');
   disp('*****************************************************************');
   disp('Running fortran code with mex function: run_WIM2d_io_mex_vSdir');
   disp('*****************************************************************');
   disp(' ');

   % real parameters
   real_prams  = [ice_prams.young,ice_prams.visc_rp,duration,CFL];

   % integer parameters
   int_prams   = [SCATMOD,ADV_DIM,ADV_OPT,...
                  DO_CHECK_FINAL,DO_CHECK_PROG,DO_CHECK_INIT,...
                  DO_BREAKING,STEADY];

   in_arrays         = zeros(nx,ny,3);
   in_arrays(:,:,1)  = ice_fields.cice;
   in_arrays(:,:,2)  = ice_fields.hice;
   in_arrays(:,:,3)  = ice_fields.Dmax;

   %% make the call!
   prep_mex_dirs(outdir);
   tic;
   [Sdir,out_arrays] = WIM2d_run_io_mex_vSdir(...
      Sdir(:),in_arrays(:),int_prams,real_prams,T_init,dir_init);
   Sdir  = reshape(Sdir,nx,ny,ndir,nw);
   toc;

   %% extract outputs
   fldnames    = {'Dmax','tau_x','tau_y','Hs','Tp'};
   Nout        = length(fldnames);
   out_arrays  = reshape(out_arrays,[nx,ny,Nout]);
   for j=1:3
      ice_fields.(fldnames{j})   = out_arrays(:,:,j);
   end
   for j=4:5
      wave_fields.(fldnames{j})  = out_arrays(:,:,j);
   end
  
   % delete annoying file
   !rm -f fort.6

else
   disp('Running pure matlab code');
   COMP_F   = 0;
   if COMP_F==1
      %% load prog binaries and compare saved wim_prog*.[ab] files
      %% to matlab results
      
      Fdir     = 'out_2/binaries/prog/';
      FF       = dir([Fdir,'wim_prog*.a'])
      if length(FF)==0
         error('COMP_F==1, but no fortran files to compare to');
      end
      %%
      F0 = FF(1).name;
      cn = F0(9:end-2);
      Ln = length(cn);
      fmt_n = sprintf('%%%d.%dd',Ln,Ln);
      afile    = [Fdir,F0];
      %%
      fmt            = 'float32';
      aid            = fopen(afile,'rb');
      F_fields.Dmax  = reshape( fread(aid,nx*ny,fmt), nx,ny );
      F_fields.tau_x = reshape( fread(aid,nx*ny,fmt), nx,ny );
      F_fields.tau_y = reshape( fread(aid,nx*ny,fmt), nx,ny );
      F_fields.Hs    = reshape( fread(aid,nx*ny,fmt), nx,ny );
      F_fields.Tp    = reshape( fread(aid,nx*ny,fmt), nx,ny );
      fclose(aid);
      %%
      if 1
         if 0
            %% plot Hs
            vc = {'Hs','H_s, m'};
            v1 = wave_fields.(vc{1});
         else
            %% plot Dmax
            vc = {'Dmax','D_{max}, m'};
            v1 = ice_fields.(vc{1});
         end
         v2    = F_fields.(vc{1});
         subplot(2,1,1);
         fn_plot_vbl(X,Y,v1,vc{2});
         subplot(2,1,2);
         fn_plot_vbl(X,Y,v2,vc{2});
         maxes = {max(v1(:)),max(v1(:))}
      elseif 1
         %% plot relative diff's
         vlist = {'Hs','Tp','tau_x'};
         lbl   = {'{\Delta}H_s/H_s','{\Delta}T_p/T_p','\Delta{\tau}_x/{\tau}_x'};
         for j=1:2
            subplot(3,1,j);
            v1    = wave_fields.(vlist{j});
            v2    = F_fields.(vlist{j});
            Z     = 0*v2;
            jn    = find(abs(v2)>0);
            Z(jn) = 1-v1(jn)./v2(jn);%%relative difference
            fn_plot_vbl(X,Y,Z,lbl{j})
         end
         for j=3:3
            subplot(3,1,j);
            v1    = ice_fields.(vlist{j});
            v2    = F_fields.(vlist{j});
            Z     = 0*v2;
            jn    = find(abs(v2)>0);
            Z(jn) = 1-v1(jn)./v2(jn);%%relative difference
            fn_plot_vbl(X,Y,Z,lbl{j})
         end
      end
   end

   %% also give progress report every 'reps' time steps;
   %reps  = nt+1;%%go straight through without reporting back or plotting
   reps     = 50;
   GET_OUT  = 1;
   if GET_OUT
      Dmax_all         = zeros(nx,ny,1+floor(nt/reps));
      Dmax_all(:,:,1)  = Dmax;
   end

   %nt = 13%%stop straight away for testing
   for n = 1:nt
      disp([n nt])

      %% spectral moments;
      mom0  = zeros(nx,ny);
      mom2  = zeros(nx,ny);
      mom0w = zeros(nx,ny);
      mom2w = zeros(nx,ny);

      %% wave stresses;
      tau_x = zeros(nx,ny);
      tau_y = zeros(nx,ny);

      %% variances of stress and strain;
      var_stress     = zeros(nx,ny);
      var_strain     = zeros(nx,ny);

      % %% test integrals;
      % var_boundary   = cell(1,length(Jy_boundary));
      % for r=1:length(Jy_boundary)
      %    var_boundary{r}   = 0;
      % end
      % var_boundary0  = var_boundary;  
      %%
      if STEADY==1
         for i = 1:nx
         for j = 1:ny
            %%top-up waves in wave mask if STEADY==1
            %%(steady-state solution);
            if WAVE_MASK2(i,j)>0 & STEADY==1
               Sdir(i,j,j_fwd,:)  = S_inc(i,j,j_fwd,:);
            end
         end
         end
      end
         
      for jw   = 1:nw

         %% CALC DIMENSIONAL ATTEN COEFF;
         atten_dim   = 0*X;
         damp_dim    = 0*X;
         for i = 1:nx
         for j = 1:ny

            
            if ICE_MASK(i,j)>0 & DO_ATTEN==1

               %% get expected no of floes met per unit
               %%  distance if travelling in a line;
               if Dmax(i,j) < 200
                  %%power law distribution
                  Dave  = floe_scaling(fragility,xi,...
                           Dmin,Dmax(i,j));
               else
                  %% uniform lengths
                  Dave  = Dmax(i,j);
               end
               c1d = cice(i,j)/Dave;%% floes per unit length;

               %% ENERGY attenuation coeff;
               atten_dim(i,j) = atten_nond(i,j,jw)*c1d;%%scattering
               damp_dim(i,j)  = 2*damping(i,j,jw)*cice(i,j);%%damping

   %           if (i==itest)&(j==jtest)
   %              disp(['Hs (pre)   = ',num2str(wave_fields.Hs(i,j)),' m']);
   %              disp(['Tp (pre)   = ',num2str(wave_fields.Tp(i,j)),' s']);
   %              disp(['Dmax       = ',num2str(Dmax(i,j)),' m']);
   %              disp(['Dave       = ',num2str(Dave),' m']);
   %              disp(['c1d        = ',num2str(c1d)]);
   %              disp(['q_scat     = ',num2str(atten_dim(i,j),'%7.7e'),' /m']);
   %              disp(['q_abs      = ',num2str(damp_dim(i,j),'%7.7e'),' /m']);
   %           end

            end
         end% j
         end% i
         % max(ag_eff(:))
         % max(atten_dim(:))
         % max(damp_dim(:))
         %1e3*max(ag_eff(:))*max(atten_dim(:))
         % GEN_pause

         s1.ndir        = ndir;
         s1.wavdir      = wavdir;
         %s1.Sdir        = reshape( Sdir(:,:,jw,:), nx,ny,ndir);
         s1.Sdir        = Sdir(:,:,:,jw);
         s1.ag_eff      = ag_eff(:,:,jw);
         s1.atten_dim   = atten_dim;
         s1.damp_dim    = damp_dim;
         s1.ICE_MASK    = ice_fields.ICE_MASK;

         if ndir==1
            if SCATMOD~=0
               disp('warning: changing SCATMOD option as not enough directions');
               disp(['(ndir = ',num2str(ndir)]);
            end
            SCATMOD  = 0;
         end

         if SCATMOD==0
            %% Simple attenuation scheme - doesn't conserve scattered energy
            [Sdir(:,:,:,jw),S_freq,tau_x_om,tau_y_om] = ...
               adv_atten_simple(grid_prams,ice_prams,s1,dt,adv_options);
            clear s1 S_out;
         elseif SCATMOD==1
            %% same as SCATMOD==0, but scattered energy
            %% is distributed isotropically
            [Sdir(:,:,:,jw),S_freq,tau_x_om,tau_y_om] = ...
               adv_atten_isotropic(grid_prams,ice_prams,s1,dt,adv_options);
            clear s1 S_out;
         end

         %% integrate stress densities over frequency
         %% TODO: check if this is correct for ice-covered water
         tmp   = rho_wtr*g*tau_x_om./ap_eff(:,:,jw);  %%[Pa*s]
         tau_x = tau_x+wt_om(jw)*tmp;                 %%[Pa]
         tmp   = rho_wtr*g*tau_y_om./ap_eff(:,:,jw);  %%[Pa*s]
         tau_y = tau_y+wt_om(jw)*tmp;                 %%[Pa]
         clear tmp;
         %GEN_pause

         %% INTEGRALS FOR BREAKING PROB:
         for i = 1:nx
         for j = 1:ny

            %% INTEGRATE SPECTRUM OVER DIRECTION;
            %S_freq   = wt_theta'*squeeze(Sdir(i,j,jw,:));

            %% convert from water amp's to ice amp's;
            F     = disp_ratio(i,j,jw);%%|T| also
            k_ice = 2*pi/wlng_ice(i,j,jw);

            %% SPECTRAL MOMENTS;
            %%take abs as small errors can make S_freq negative
            mom0(i,j)   = mom0(i,j)+abs( wt_om(jw)*S_freq(i,j)*F^2 );
            mom2(i,j)   = mom2(i,j)+abs( wt_om(jw)*S_freq(i,j)*F^2*om_vec(jw)^2 );

            %% reference values in water (F=1)
            mom0w(i,j)  = mom0w(i,j)+abs( wt_om(jw)*S_freq(i,j) );
            mom2w(i,j)  = mom2w(i,j)+abs( wt_om(jw)*S_freq(i,j)*om_vec(jw)^2 );

            if ICE_MASK(i,j)==1
               %% VARIANCE OF STRAIN;
               strain_density    = abs( S_freq(i,j)*F^2*...
                                     (k_ice^2*hice(i,j)/2)^2 );
               var_strain(i,j)   = var_strain(i,j)+...
                                    + wt_om(jw)*strain_density;
            end
         end%% end spatial loop x;
         end%% end spatial loop y;

      end%% end spectral loop;
      %mom0,mom2,wlng_ice,return

      %%calc Hs, Tw (into wave_fields.Tp)
      if REF_Hs_ICE==1
         %% diagnostic variables Hs/Tp related to ice displacment in ice-covered areas
         wave_fields.Hs       = 4*sqrt(mom0);%%diagnostic variable - for waves in water
         wave_fields.Tp       = 0*X;
         jnz                  = find(mom2>0);
         wave_fields.Tp(jnz)  = 2*pi*sqrt(mom0(jnz)./mom2(jnz));
      else
         %% diagnostic variables Hs/Tp related to water displacment in ice-covered areas
         wave_fields.Hs       = 4*sqrt(mom0w);%%diagnostic variable - for waves in water
         wave_fields.Tp       = 0*X;
         jnz                  = find(mom2w>0);
         wave_fields.Tp(jnz)  = 2*pi*sqrt(mom0w(jnz)./mom2w(jnz));
      end

      %% wave stresses
      ice_fields.tau_x  = tau_x;
      ice_fields.tau_y  = tau_y;

      %% FINALLY DO FLOE BREAKING;
      for i=1:nx
      for j=1:ny
         if ICE_MASK(i,j)==1 & mom0(i,j)>0
            %% only try breaking if ice is present
            %%  & some waves have arrived;

            %% significant strain amp
            sig_strain  = 2*sqrt(var_strain(i,j));

            %%  probability of critical strain
            %%  being exceeded from Rayleigh distribution;
            Pstrain  = exp( -strain_c^2/(2*var_strain(i,j)) );
            P_crit   = (1-DO_BREAKING)+exp(-1);%%this is critical prob if monochromatic wave

            %% FLOE BREAKING:
            BREAK_CRIT     = ( Pstrain>=P_crit );%%breaks if larger than this
            brkcrt(i,j,n)  = BREAK_CRIT;

            if BREAK_CRIT
               %% use crest period to work out wavelength
               %% - half this is max poss floe length;
               T_crit   = wave_fields.Tp(i,j);
               if 0
                  %%get wavelength directly
                  wlng_crest  = ...
                     GEN_get_ice_wavelength(hice(i,j),T_crit,Inf,young);
               else
                  %% interpolate (to check fortran code)
                  %% NB slight difference due to RT_param_outer
                  %% using finite depth wavelength approx to
                  %% inifinit depth one
                  om       = 2*pi/T_crit;
                  om_min   = om_vec(1);
                  om_max   = om_vec(end);
                  if (om<=om_min)
                     wlng_crest  = wlng_ice(i,j,1);
                     %tst_crest   = [wlng_crest,...
                     %   GEN_get_ice_wavelength(hice(i,j),T_crit,Inf,young)]
                  elseif (om>=om_max)
                     wlng_crest  = wlng_ice(i,j,nw);
                  else
                     jcrest      = floor((om-om_min+dom)/dom);
                     om1         = 2*PI*freq_vec(jcrest);
                     lam1        = wlng_ice(i,j,jcrest);
                     lam2        = wlng_ice(i,j,jcrest+1);
                     wlng_crest  = lam1+(om-om1)*(lam2-lam1)/dom;
                  end
               end

               Dc          = max(Dmin,wlng_crest/2);
               Dmax(i,j)   = min(Dc,Dmax(i,j));
            end%% end breaking action;

            if 0%i==11 & j==1
               BREAK_CRIT
               Dmax(i,j)
            end

         elseif WTR_MASK(i,j)==1%% only water present
            Dmax(i,j)   = 0;
         end

      end%% end spatial loop j in y;
      end%% end spatial loop i in x;

      ice_fields.Dmax   = Dmax;
      jmiz              = find((Dmax>0)&(Dmax<250));
      % {jmiz}

      %% progress report;
      if round(n/reps)==(n/reps)

         if GET_OUT
            Dmax_all(:,:,n/reps) = Dmax;
         end

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         disp('##############################################################');
         t1 = now;
         disp([num2str(n),' time steps done, out of ',num2str(nt)]);
         disp(['Time taken (mins)      : ' ,num2str(t0_fac*(t1-t0))]);
         disp(['Model time passed (h)  : ' ,num2str(n*dt/3600.)]);
         disp('##############################################################');
         disp(' ');
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         disp(strvcat(Info));

         if PLOT_PROG
            if DIAG1d==0
               figure(2),clf;
               fn_fullscreen;
               %%
               if PLOT_OPT==1
                  s1 = struct('dir',wave_stuff.dirs(jdir),...
                              'period',Tc,...
                              'Sdir',Sdir(:,:,jdir,jchq));
                  fn_plot_spec(X,Y,wave_fields.Hs,wave_fields.Tp,Dmax,s1);
               else
                  fn_plot_spec_2(X,Y,wave_fields.Hs,ice_fields.tau_x,...
                     Dmax,ice_fields.tau_y);
               end

               %if OPT==1
               %   subplot(2,2,1);
               %   hold on;
               %   x0 = min(X(:))+uc*n*dt;
               %   x1 = X(find(WAVE_MASK(:,1)==0,1,'first'),1)+uc*n*dt;
               %   % {uc,x0/1e3,x1/1e3}
               %   yc = .3*max(Y(:))/1e3;
               %   x_ = [min(X(:)),x0,x0,x1,x1,max(X(:))]/1e3;
               %   y_ = [0,0,yc*[1,1],0,0];
               %   plot(x_,y_,'k');
               %   plot(X(:,1)/1e3,wave_fields.Hs(:,1)*yc,'--k');
               %   hold off;
               %   %%
               %   subplot(2,2,2);
               %   hold on;
               %   plot(x_,y_,'k');
               %   plot(X(:,1)/1e3,wave_fields.Hs(:,1)*yc,'--k');
               %   hold off;
               %end
            else
               %% DIAG1d==1
               %% during run
               if 1
                  figure(4);
                  %% check symmetry
                  [hmax,imax] = max(wave_fields.Hs(:,1));
                  hp          = wave_fields.Hs(imax,:);
                  yp          = grid_prams.Y(imax,:);
                  plot(yp/1e3,hp);
                  ttl   = title(['max h = ',num2str(max(wave_fields.Hs(:))),'; x = ',num2str(X(imax,1)/1.0e3),'km']);
                  GEN_font(ttl);
                  GEN_proc_fig('y, km','H_s, m')
               else
                  %%check partition of fwd and back energy
                  figure(4);
                  subplot(4,1,1);
                  hold on;
                  %%
                  %fn_plot1d(X(:,1)/1e3,wave_fields.Hs(:,1),labs1d_1,cols{loop_col});
                  fn_plot1d(X(:,1)/1e3,mean(wave_fields.Hs,2),labs1d_1,cols{loop_col});%%average over y (columns)
                  hold on;
                  %%
                  [Ep,Em,Et1,Et2]   = fn_split_energy(om_vec,wavdir,Sdir);
                  %Hp                = 4*sqrt(Ep(:,1));
                  %Hm                = 4*sqrt(Em(:,1));
                  Hp                = 4*sqrt(mean(Ep,2));
                  Hm                = 4*sqrt(mean(Em,2));
                  if DIAG1d_OPT==0
                     %Hs2   = 4*sqrt(Ep(:,1)+Em(:,1));%%add Ep + Em
                     Hs2   = 4*sqrt(mean(Ep,2)+mean(Em,2));
                  elseif DIAG1d_OPT==1
                     %Hs2   = 4*sqrt(Et1(:,1));%%check const panel integration
                     Hs2   = 4*sqrt(mean(Et1,2));
                  elseif DIAG1d_OPT==2
                     %Hs2   = 4*sqrt(Et2(:,1));%%check Simpson's rule integration
                     Hs2   = 4*sqrt(mean(Et2,2));
                  end
                  fn_plot1d(X(:,1)/1e3,Hs2,labs1d_1,['-',cols{loop_col}]);
                  hold on;
                  %%
                  subplot(4,1,2);
                  fn_plot1d(X(:,1)/1e3,Hp,labs1d_2,cols{loop_col});
                  hold on;
                  fn_plot1d(X(:,1)/1e3,Hs2,labs1d_2,['-',cols{loop_col}]);
                  hold on;
                  %%
                  subplot(4,1,3);
                  fn_plot1d(X(:,1)/1e3,Hm,labs1d_3,cols{loop_col});
                  hold on;
                  %%
                  loop_col = loop_col+1;
                  if loop_col>length(cols)
                     loop_col = 1;
                  end
               end
            end

            if COMP_F==1
               %% load prog binaries and compare saved wim_prog*.[ab] files
               %% to matlab results
               Fdir           = 'out_2/binaries/prog';
               nnn            = num2str(n,'%3.3d');
               afile          = [Fdir,nnn,'.a'];
               aid            = fopen(afile,'rb');
               F_fields.Dmax  = reshape( fread(aid,nx*ny,fmt), nx,ny );
               F_fields.tau_x = reshape( fread(aid,nx*ny,fmt), nx,ny );
               F_fields.tau_y = reshape( fread(aid,nx*ny,fmt), nx,ny );
               F_fields.Hs    = reshape( fread(aid,nx*ny,fmt), nx,ny );
               F_fields.Tp    = reshape( fread(aid,nx*ny,fmt), nx,ny );
               fclose(aid);
            end

            clear s1;
            pause(0.1);
         end
         %%
         if 0
            Jlook = nx-16+(1:16);
            Dmax(edge:edge+9,Jlook)
         end
      end


      for jw=[]%11
         [n,T(jw)],%[ag_ice(jw,1:11,1)]
         testSatt=S(1:11,1,jw,jmwd)
         pause
      end

      if (SV_BIN==1)&(mod(n,reps_ab)==0)&(DO_CHECK_PROG==1)
         %% save matlab files as binaries
         %% to matlab results
         Fdir              = 'm_out/binaries/prog';
         cnt               = num2str(nt);
         cnt(:)            = '0';
         cn                = num2str(n);
         lc                = length(cn);
         cnt(end+1-lc:end) = cn;
         Froot = [Fdir,'/wim_prog',cnt];
         %%
         pairs = {};
         pairs{end+1}   = {'Dmax' ,ice_fields.Dmax};
         pairs{end+1}   = {'tau_x',ice_fields.tau_x};
         pairs{end+1}   = {'tau_y',ice_fields.tau_y};
         pairs{end+1}   = {'Hs'   ,wave_fields.Hs};
         pairs{end+1}   = {'Tp'   ,wave_fields.Tp};
         %%
         fn_save_binary(Froot,Bdims,n*dt,pairs);
      end

   end%% end time loop
end%%MEX_OPT==0 option

if (SV_BIN==1)&(DO_CHECK_FINAL==1)
   %% save matlab files as binaries
   %% to matlab results
   Fdir  = 'm_out/binaries';
   Froot = [Fdir,'/wim_out'];
   %%
   pairs = {};
   pairs{end+1}   = {'Dmax' ,ice_fields.Dmax};
   pairs{end+1}   = {'tau_x',ice_fields.tau_x};
   pairs{end+1}   = {'tau_y',ice_fields.tau_y};
   pairs{end+1}   = {'Hs'   ,wave_fields.Hs};
   pairs{end+1}   = {'Tp'   ,wave_fields.Tp};
   %%
   fn_save_binary(Froot,Bdims,duration,pairs);
end

t1 = now;

if (OPT==1)|(OPT==3)
   Dmax_j   = ice_fields.Dmax(:,1);
   jmiz     = find((Dmax_j>0)&(Dmax_j<250));
   Wmiz     = dx/1e3*length(jmiz);
   %%
   Dmax_min = min(Dmax_j);
   Dmax_max = max(Dmax_j);
   %%
   disp(' ');
   disp(['MIZ width = ',num2str(Wmiz),' km']);
end

taux_min = min(ice_fields.tau_x(:))
taux_max = max(ice_fields.tau_x(:))
tauy_min = min(ice_fields.tau_x(:))
tauy_max = max(ice_fields.tau_y(:))
disp(['max tau_x = ',num2str(taux_max),' Pa']);
disp(['max tau_y = ',num2str(tauy_max),' Pa']);
disp(' ');

%%append to log file
logid = fopen(log_file,'a');
fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s\n','Diagnostics:');
if (OPT==1)|(OPT==3)
   fprintf(logid,'%s%6.1f\n','MIZ width (km): ',Wmiz);
   fprintf(logid,'%s%6.1f%s%6.1f\n','Dmax range in MIZ (m): ',...
      Dmax_min,' ',Dmax_max);
end
fprintf(logid,'%s%10.3e%s%10.3e\n','tau_x range (Pa): ',...
   taux_min,' ',taux_max);
fprintf(logid,'%s%10.3e%s%10.3e\n','tau_y range (Pa): ',...
   tauy_min,' ',tauy_max);
fprintf(logid,'%s\n','***********************************************');

fprintf(logid,'%s,\n',' ');
fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s,%7.1f\n','Elapsed time (min):',t0_fac*(t1-t0));
fprintf(logid,'%s\n','***********************************************');
fprintf(logid,'%s\n',' ');
fclose(logid);

if SV_SPEC
   %% save final directional spectrum
   !mkdir -p m_out
   freq_vec = om_vec/2/pi;
   S_inc    = wave_stuff.dir_spec;
   Hs       = wave_fields.Hs;
   save('m_out/Sdir.mat','Sdir','wavdir','freq_vec','grid_prams','S_inc','cice','Hs');
   clear freq_vec S_inc;
end

if TEST_FINAL_SPEC==1
   disp(' ');
   disp('Testing final spectrum...');

   if 0
      disp('(Check integrals of output spectrum vs output Hs,Tp,etc)');
      %% check consistency of Sdir with Hs,Tp
      [wf.Hs,wf.Tp,wf.mwd] = fn_spectral_integrals(om_vec,wavdir,Sdir);

      vbls  = {'Hs','Tp'};%,'mwd'};%mwd currently not updated
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
   elseif MEX_OPT>0
      disp('(Check outputs vs values in binary files)');
      of2         = fn_check_final(outdir);%%set in infile_dirs.txt
      of1.Hs      = wave_fields.Hs;
      of1.Tp      = wave_fields.Tp;
      of1.tau_x   = ice_fields.tau_x;
      of1.tau_y   = ice_fields.tau_y;
      of1.Dmax    = ice_fields.Dmax;
      %%
      vbls  = {'Hs','Tp','tau_x','tau_y','Dmax'};%,'mwd'};%mwd currently not updated
      for n=1:length(vbls)
         vbl   = vbls{n};
         disp(' ');
         disp(['comparing field: ',vbl]);
         v1    = of1.(vbl);
         v2    = of2.(vbl);
         diff  = abs(v2-v1);
         disp(['max diff: ',num2str(max(diff(:)))]);
         disp(' ');
      end
   end

   return
end

if PLOT_FINAL%%check exponential attenuation
   figure(2),clf;
   fn_fullscreen;
   if PLOT_OPT==1
      s1 = struct('dir',wave_stuff.dirs(jdir),...
                  'period',Tc,...
                  'Sdir',Sdir(:,:,jdir,jchq));
      fn_plot_spec(X,Y,wave_fields.Hs,wave_fields.Tp,Dmax,s1);
   else
      fn_plot_spec_2(X,Y,wave_fields.Hs,ice_fields.tau_x,...
         Dmax,ice_fields.tau_y);
   end
   %%
   if 0
      %% figure testing how 1d results are (only appropriate for 1d geometries)
      figure(3),clf;
      fn_fullscreen;
      xx = X(:,1);
      if 0
         vbl   = 'Hs';
         Vbl   = wave_fields.(vbl);
      else
         vbl   = 'tau_x';
         Vbl   = ice_fields.(vbl);
      end

      if 1
         subplot(2,1,2)
         yy    = Y(1,:);
         xp    = 110e3;
         ix    = find(abs(xx-xp)==min(abs(xx-xp)));
         ix    = ix(1);
         xp    = xx(ix);
         Vy    = Vbl(ix,:);
         plot(yy/1e3,Vy);
         if strcmp(vbl,'Hs')
            GEN_proc_fig('{\ity}, km','{\itH}_s, m');
         elseif strcmp(vbl,'tau_x')
            GEN_proc_fig('{\ity}, km','{\tau}_x, Pa');
         end
         ttl   = title(['Profile at x = ',num2str(xp/1e3,'%7.2f'),'km']);
         GEN_font(ttl);
         if Vy>0
            ylim(sort([0,1.1*max(Vy)]));
         end
         %%
         subplot(2,1,1)
      end
      plot(xx/1e3,mean(Vbl,2),'-k');
      hold on;
      plot(xx/1e3,min(Vbl,[],2),'--r');
      plot(xx/1e3,max(Vbl,[],2),'--c');
      %set(gca,'yscale','log');
      if strcmp(vbl,'Hs')
         GEN_proc_fig('{\itx}, km','{\itH}_s, m');
      elseif strcmp(vbl,'tau_x')
         GEN_proc_fig('{\itx}, km','{\tau}_x, Pa');
      end
      legend('Mean','Min','Max');
   end
   %%
   if DIAG1d==1
      %% final
      figure(5);
      fn_fullscreen;
      clf;

      COMP_STEADY = 0;%%compare to fortran results
      dfiles      = {};
      leg_text    = {};
      fortcols    = {};
      if COMP_STEADY
         frun        = '../../fortran/run/';
         fdir        = [frun,'/fig_scripts/figs/TC2S/'];%%location of text files with fortran results
         mdir        = '../boltzmann/out/';

         % time-dep results (from fortran code)
         dfiles  {end+1} = [fdir,'/test_steady1.dat']; % file name
         leg_text{end+1} = 'F77 (time-dep)';          % legend text
         fortcols{end+1} = 'b';                       % colour

         % steady-state results (from python code)
         dfiles  {end+1} = [fdir,'/test_steady2.dat']; % file name
         leg_text{end+1} = 'python (steady)';         % legend text
         fortcols{end+1} = '--g';                     % colour

         % steady-state results (from python code v2)
         dfiles  {end+1} = [fdir,'/test_steady2_FT.dat']; % file name
         leg_text{end+1} = 'python (steady, v2)';        % legend text
         fortcols{end+1} = '--c';                        % colour

         % steady-state results (from matlab code)
         dfiles  {end+1} = [mdir,'/test_steady_mat.dat']; % file name
         leg_text{end+1} = 'matlab (steady)';            % legend text
         fortcols{end+1} = '--m';                        % colour

         leg_text_used  = {};
         fortcols_used  = {};
         for k=1:length(dfiles)
            dfil  = dfiles{k};
            if exist(dfil)
               leg_text_used{end+1} = leg_text{k};
               fortcols_used{end+1} = fortcols{k};
               %%
               disp(['opening ',dfil,'']);
               fid   = fopen(dfil,'r');

               %% search for hash lines
               found_hash  = 0;
               while ~found_hash
                  lin   = fgets(fid);
                  if length(lin>=5)
                     found_hash  = strcmp(lin(1:5),'#####');
                  end
               end

               %%skip one more line, then get data;
               fgets(fid);
               columns  = textscan(fid,'%f %f');
               xx_f     = columns{1};
               Hs_f     = columns{2};
               fclose(fid);

               plot(xx_f/1e3,Hs_f,fortcols{k},'linewidth',2);
               hold on;

               %% check y-dependance
               if k==1
                  figure(3);
                  ix = find(abs(xx_f-xp)==min(abs(xx-xp)));
                  subplot(2,1,2)
                  hold on;
                  for loop_ix=1:length(ix)
                     plot(yy/1e3,Hs_f(ix(loop_ix))+0*yy,'--m');
                  end
                  figure(5);
               end
            else
               disp([dfil,' not present']);
               disp('To create, run ../../fortran/run/fig_scripts/fig_test_convergence2steady.py');
               disp('or ../boltmann/fig_Boltzmann_Steady.m');
            end
         end
         leg_text = leg_text_used;
         fortcols = fortcols_used;
      end

      fcols    = cols;
      fcols{1} = '-r';
      %fn_plot1d(X(:,1)/1e3,wave_fields.Hs(:,1),labs1d_1,cols{1});
      fn_plot1d(X(:,1)/1e3,mean(wave_fields.Hs,2),labs1d_1,fcols{1});
      leg_text{end+1}   = 'Total';
      hold on;
      %%
      if 0
         [Ep,Em,Et1,Et2]   = fn_split_energy(om_vec,wavdir,Sdir);
         %Hp                = 4*sqrt(Ep(:,1));
         %Hm                = 4*sqrt(Em(:,1));
         Hp                = 4*sqrt(mean(Ep,2));
         Hm                = 4*sqrt(mean(Em,2));
         if DIAG1d_OPT==0
            %Hs2               = 4*sqrt(Ep(:,1)+Em(:,1));%%add Ep + Em
            Hs2               = 4*sqrt(mean(Ep,2)+mean(Em,2));
            leg_text{end+1}   = 'Total (test 1)';
         elseif DIAG1d_OPT==1
            %Hs2               = 4*sqrt(Et1(:,1));%%check const panel integration
            Hs2               = 4*sqrt(mean(Et1,2));
            leg_text{end+1}   = 'Total (test 2)';
         elseif DIAG1d_OPT==2
            %Hs2               = 4*sqrt(Et2(:,1));%%check Simpson's rule integration
            Hs2               = 4*sqrt(mean(Et2,2));
            leg_text{end+1}   = 'Total (Simpson''s)';
         end

         fn_plot1d(X(:,1)/1e3,Hs2,labs1d_1,['-',fcols{1}]);
         hold on;

         fn_plot1d(X(:,1)/1e3,Hp,labs1d_1,fcols{2});
         leg_text{end+1}   = 'Fwd';
         hold on;

         fn_plot1d(X(:,1)/1e3,Hm,labs1d_1,fcols{3});
         leg_text{end+1}   = 'Back';
      end
      %%

      %% make legend
      cmd   = 'legend(';
      for k=1:length(leg_text)
         cmd   = [cmd,'''',leg_text{k},''','];
      end
      cmd(end:end+1) = ');';
      eval(cmd);
   end

   if SV_FIG%%save figures

      if nw==1
         if SCATMOD==1
            fig_dir  = 'out/isotropic_1freq';  %%use this for monochromatic wave
         elseif SCATMOD==0
            fig_dir  = 'out/simple_1freq';  %%use this for monochromatic wave
         end
      else
         if SCATMOD==1
            fig_dir  = 'out/isotropic_spec';  %%use this for spectrum
         elseif SCATMOD==0
            fig_dir  = 'out/simple_spec';  %%use this for spectrum
         end
      end

      if ~exist(fig_dir)
         mkdir(fig_dir)
      end

      if ~exist([fig_dir,'/att_fig'])
         mkdir([fig_dir,'/att_fig']);
         mkdir([fig_dir,'/att_png']);
         mkdir([fig_dir,'/fig']);
         mkdir([fig_dir,'/png']);
      end

      figure(2);
      saveas(gcf,[fig_dir,'/fig/B',num2str(ndir,'%3.3d'),'.fig']);
      saveas(gcf,[fig_dir,'/png/B',num2str(ndir,'%3.3d'),'.png']);

      if 0
         figure(3);
         if 0
            %%fix position so that comparison is easier between computers
            pos   = [0.13   0.121428571428571   0.775   0.803571428571429];
            set(gca,'position',pos);
         end
         saveas(gcf,[fig_dir,'/att_fig/B',num2str(ndir,'%3.3d'),'_atten.fig']);
         saveas(gcf,[fig_dir,'/att_png/B',num2str(ndir,'%3.3d'),'_atten.png']);
      end
   end
end

%%display info again
disp(strvcat(Info));

%% save time-stepped Dmax;
%if DO_SAVE
%   save(filename,'out');
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn_plot_spec(X,Y,Hs,Tw,Dmax,s1)
%% plot Hs, Dmax, Tw &
%% S for 1 particular freq and dir

[nx,ny]  = size(X);

subplot(2,2,1);
if ny==1
   [xx,yy]  = step_1d(X/1e3,Hs);
   plot(xx,yy);
   GEN_proc_fig('\itx, \rmkm','{\itH}_{\rm s}, m');
else
   H  = pcolor(X/1e3,Y/1e3,Hs);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   colorbar;
   ttl   = title('{\itH}_{\rm s}, m');
   GEN_font(ttl);
end

subplot(2,2,2);
if ny==1
   [xx,yy]  = step_1d(X/1e3,Dmax);
   plot(xx,yy);
   GEN_proc_fig('\itx, \rmkm','{\itD}_{\rm max}, m');
else
   H  = pcolor(X/1e3,Y/1e3,Dmax);
   caxis([0 250]);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   colorbar;
   ttl   = title('{\itD}_{\rm max}, m');
   GEN_font(ttl);
end

subplot(2,2,3);
if ny==1
   [xx,yy]  = step_1d(X/1e3,Tw);
   plot(xx,yy);
   GEN_proc_fig('\itx, \rmkm','{\itT}_{\rm w}, s');
else
   H  = pcolor(X/1e3,Y/1e3,Tw);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   colorbar;
   ttl   = title('{\itT}_{\rm w}, s');
   GEN_font(ttl);
end

subplot(2,2,4);
detls = ['S: ',num2str(s1.period),'s, ',num2str(s1.dir),'^o'];
if ny==1
   [xx,yy]  = step_1d(X/1e3,s1.Sdir);
   plot(xx,yy);
   GEN_proc_fig('\itx, \rmkm',detls);
else
   H  = pcolor(X/1e3,Y/1e3,s1.Sdir);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   colorbar;
   ttl   = dtls;
   ttl   = title({'\itS, \rmm^2s',ttl});
   GEN_font(ttl);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn_plot_spec_2(X,Y,Hs,tau_x,Dmax,tau_y)
%% plot Hs, Dmax, Tw &
%% S for 1 particular freq and dir

[nx,ny]  = size(X);

%%fix positions so figures can be compared more easily between computers
pos1  = [0.130000000000000   0.583837209302326   0.334659090909091   0.341162790697674];
pos2  = [0.570340909090909   0.583837209302326   0.334659090909091   0.341162790697674];
pos3  = [0.130000000000000   0.110000000000000   0.334659090909091   0.341162790697674];
pos4  = [0.570340909090909   0.110000000000000   0.334659090909091   0.341162790697674];

%subplot(2,2,1);
subplot('position',pos1);
if ny==1
   [xx,yy]  = step_1d(X/1e3,Hs);
   plot(xx,yy);
   GEN_proc_fig('\itx, \rmkm','{\itH}_{\rm s}, m');
else
   H  = pcolor(X/1e3,Y/1e3,Hs);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   colorbar;
   GEN_font(gca);
   ttl   = title('{\itH}_{\rm s}, m');
   GEN_font(ttl);
end

%subplot(2,2,2);
subplot('position',pos2);
if ny==1
   [xx,yy]  = step_1d(X/1e3,Dmax);
   plot(xx,yy);
   GEN_proc_fig('\itx, \rmkm','{\itD}_{\rm max}, m');
else
   H  = pcolor(X/1e3,Y/1e3,Dmax);
   caxis([0 250]);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   colorbar;
   GEN_font(gca);
   ttl   = title('{\itD}_{\rm max}, m');
   GEN_font(ttl);
end

%subplot(2,2,3);
subplot('position',pos3);
if ny==1
   [xx,yy]  = step_1d(X/1e3,tau_x);
   plot(xx,yy);
   GEN_proc_fig('\itx, \rmkm','{\tau}_{x}, Pa');
else
   H  = pcolor(X/1e3,Y/1e3,tau_x);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   colorbar;
   GEN_font(gca);
   ttl   = title('{\tau}_{x}, Pa');
   GEN_font(ttl);
end

%subplot(2,2,4);
subplot('position',pos4);
if ny==1
   [xx,yy]  = step_1d(X/1e3,tau_y);
   plot(xx,yy);
   GEN_proc_fig('\itx, \rmkm','{\tau}_{y}, Pa');
else
   H  = pcolor(X/1e3,Y/1e3,tau_y);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   colorbar;
   GEN_font(gca);
   ttl   = title('{\tau}_{y}, Pa');
   GEN_font(ttl);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = fn_plot1d(x,y,labs,col)
if ~exist('col','var')
   col   = '-k';
end
H  = plot(x,y,col);
GEN_proc_fig(labs{1},labs{2});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
function fn_plot_vbl(X,Y,Z,lbl)

H  = pcolor(X/1e3,Y/1e3,Z);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title(lbl);
GEN_font(ttl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
