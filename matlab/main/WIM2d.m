function [out_fields,wave_stuff,diagnostics,mesh_e] =...
   WIM2d(params_in,gridprams,ice_fields,wave_fields,wave_stuff,mesh_e)

%% ============================================================
%% Inputs:
%%
%% params_in = structure eg:
%%           SCATMOD: 1
%%           ADV_DIM: 2
%%           ADV_OPT: 2
%%     DO_CHECK_INIT: 1
%%     DO_CHECK_PROG: 1
%%    DO_CHECK_FINAL: 1
%%       BRK_OPT: 0
%%            STEADY: 1
%%          DO_ATTEN: 1
%%             young: 5.4900e+09
%%           drag_rp: 0
%%               CFL: 0.7000
%%    duration_hours: 24
%%           DO_DISP: 0
%%         PLOT_INIT: 1
%%         PLOT_PROG: 1
%%        PLOT_FINAL: 1
%%         CHK_ATTEN: 0
%%       USE_ICE_VEL: 0
%%            DIAG1d: 0
%%        DIAG1d_OPT: 1
%%
%%
%% ice_fields  = structure eg:
%%      cice: [51x51 double]
%%      hice: [51x51 double]
%%      Dmax: [51x51 double]
%%
%%
%% wave_fields = structure eg
%%           Hs: [150x10 double]
%%           Tp: [150x10 double]
%%          mwd: [150x10 double]
%%  STEADY_MASK: [150x10 logical]
%%
%%
%% wave_stuff = structure eg
%%        nfreq: 2
%%         ndir: 16
%%         freq: [23x1 double]
%%         dirs: [16x1 double]
%%     dir_spec: [150x20x16x23 double]
%%
%%
%% optional:
%% - give coordinates of FEM mesh centres
%% - do breaking on mesh as well to try to reduce numerical diffusion
%%   caused by interpolation between grid and mesh
%% mesh_e = structure eg
%%         xe: [760x1 double]
%%         ye: [760x1 double]
%%          c: [760x1 double]
%%          h: [760x1 double]
%%     Nfloes: [760x1 double]
%% DAMAGE_OPT: 1
%% ============================================================

format long;

%%check params_in has the needed fields
check_params_in(params_in);

%% check if we want to do breaking on the mesh also
INTERP_MESH = 0;
if exist('mesh_e','var') & params_in.MEX_OPT==0
   INTERP_MESH    = 1;
   mesh_e.broken  = 0*mesh_e.c;%0/1 if ice was broken by waves this time
else
   mesh_e  = [];
end

USE_EBS  = 0;
if params_in.SCATMOD==3
   USE_EBS  = 1;
end

%% get date & time;
date_vector    = [params_in.start_year,...
                  params_in.start_month,...
                  params_in.start_day,...
                  params_in.start_hour,...
                  params_in.start_minute,...
                  params_in.start_second];
year_info      = datevec2year_info(date_vector)
model_day      = year_info.model_day;
model_seconds  = year_info.model_seconds;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (params_in.ADV_DIM==2)&(gridprams.ny<4)
   error({'incompatible values of params_in.ADV_DIM and gridprams.ny:';
          'increase gridprams.ny or use params_in.ADV_DIM=1'});
else
   adv_options = struct('ADV_DIM',params_in.ADV_DIM,...
                        'ADV_OPT',params_in.ADV_OPT);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% TURN ON/OFF PLOTTING:
PLOT_OPT          = 2;%%plot option (only used if doing plotting)
                      %%(if params_in.PLOT_INIT==1|params_in.PLOT_PROG==1|params_in.PLOT_FINAL==1)
TEST_INC_SPEC     = 0;
TEST_FINAL_SPEC   = 0;

TEST_IJ  = (params_in.itest>0)&(params_in.jtest>0);

CSUM     = params_in.DO_CHECK_INIT+params_in.DO_CHECK_PROG+params_in.DO_CHECK_FINAL;
SV_BIN   = (CSUM>0);

FSUM     = params_in.PLOT_INIT+params_in.PLOT_PROG+params_in.PLOT_FINAL;
SV_FIGS  = (FSUM>0);

COMP_F   = 0;
compFdir = 'out_2/binaries/prog/';

%% make a log file similar to fortran file
if TEST_IJ | params_in.SV_LOG | SV_BIN | SV_FIGS
   if params_in.outdir==0
   end
   log_dir  = params_in.outdir;
   eval(['!mkdir -p ',log_dir]);
end

if TEST_IJ | params_in.SV_LOG
   log_dir  = [params_in.outdir,'/diagnostics'];
   eval(['!mkdir -p ',log_dir]);
end

if SV_BIN
   log_dir  = [params_in.outdir,'/binaries'];
   eval(['!mkdir -p ',log_dir]);
   log_dir  = [params_in.outdir,'/binaries/prog'];
   eval(['!mkdir -p ',log_dir]);
end

if SV_FIGS
   log_dir  = [params_in.outdir,'/figs'];
   eval(['!mkdir -p ',log_dir]);
end

if TEST_IJ
   log_dir2 = [params_in.outdir,'/diagnostics/local'];
   eval(['!mkdir -p ',log_dir]);
end

if params_in.SV_LOG
   log_dir  = [params_in.outdir,'/diagnostics/global'];
   eval(['!mkdir -p ',log_dir]);

   log_file    = [log_dir,'/WIM2d_diagnostics',year_info.date_string,'.txt'];
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
   fprintf(logid,'%s%2.2d\n','SCATMOD:                          ',params_in.SCATMOD);
   fprintf(logid,'%s%2.2d\n','ADV_DIM:                          ',params_in.ADV_DIM);
   fprintf(logid,'%s%2.2d\n','ADV_OPT:                          ',params_in.ADV_OPT);
   fprintf(logid,'%s%2.2d\n','BRK_OPT:                          ',params_in.BRK_OPT);
   fprintf(logid,'%s%2.2d\n','STEADY:                           ',params_in.STEADY);
   fprintf(logid,'%s%2.2d\n','DO_ATTEN:                         ',params_in.DO_ATTEN);
   fprintf(logid,'%s\n','***********************************************');
   fprintf(logid,'%s\n',' ');
   fclose(logid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if params_in.DO_DISP; disp('Initialization'); end

cmin     = 0;
ICE_MASK = zeros(size(ice_fields.cice));
ICE_MASK(ice_fields.cice>cmin)   = 1;
WTR_MASK = (1-ICE_MASK).*(1-gridprams.LANDMASK);

%% add wave stress computation
out_fields.tau_x  = 0*ice_fields.cice;
out_fields.tau_y  = 0*ice_fields.cice;
out_fields.Dmax   = ice_fields.Dmax;

if params_in.DO_DISP;
   Dice  = out_fields.Dmax(ice_fields.cice>0);
   Drng  = [min(Dice),max(Dice(Dice<300)),max(out_fields.Dmax(:))]
   clear Dice
end
%% add max distance at which waves could have broken the ice
if ~params_in.BRK_OPT
 diagnostics.break_max = nan;%%TW move scalar outputs to new structure "diagnostics"
end

% %% WAVES
om_vec   = 2*pi*wave_stuff.freq; %% radial freq

%% set mask
WAVE_MASK   = 0*wave_fields.Tp;
WAVE_MASK(wave_fields.Tp>0)   = 1;

%% get rest of ice_prams
if ~isnan(params_in.young)
   ice_prams.young      = params_in.young;
   ice_prams.young_opt  = NaN;
else
   ice_prams.young_opt  = 1;
end
if ~isnan(params_in.drag_rp)
   ice_prams.drag_rp = params_in.drag_rp;
end
if ~isnan(params_in.visc_ws)
   ice_prams.visc_ws  = params_in.visc_ws;
end
ice_prams.BRK_OPT = params_in.BRK_OPT;
ice_prams         = fn_fill_iceprams(ice_prams);
%% ice_prams = structure eg:
%%               c: 0.750000000000000
%%               h: 2
%%            Dmax: 300
%%           young: 2.000000000000000e+09       % Young's modulus                [Pa]
%%          bc_opt: 0
%%         drag_rp: 13                          % Robinson-Palmer drag coeff     [Pa/(m/s)]
%% visc_ws: 0                           % Wang-Shen viscoelastic coeff   [Pa/(m/s)]
%%          rhowtr: 1.025000000000000e+03       % Water density                  [kg/m^3]
%%          rhoice: 9.225000000000000e+02       % Ice density                    [kg/m^3]
%%               g: 9.810000000000000           % Gravitational acceleration     [m/s^2]
%%         poisson: 0.300000000000000           % Poisson's ratio
%%             vbf: 0.100000000000000           % brine volume fraction
%%              vb: 100                         % ppt (1e3*vbf)
%%         sigma_c: 2.741429878818372e+05       % breaking stress                [Pa] 
%%        strain_c: 1.370714939409186e-04       % breaking strain                [-] 
%%  flex_rig_coeff: 1.831501831501831e+08
%%            Dmin: 20                          % minimum floe size              [m]
%%              xi: 2                           % no of pieces floes break into
%%       fragility: 0.900000000000000           % probability that floes break

if 0
   figure,fn_fullscreen;
   fn_plot_ice(gridprams,ice_fields);
   figure,fn_fullscreen;
   fn_plot_waves(gridprams,wave_fields);
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%append to log file
if params_in.SV_LOG
   logid = fopen(log_file,'a');
   fprintf(logid,'%s\n','***********************************************');
   fprintf(logid,'%s\n','WIM parameters:');
   fprintf(logid,'%s%4.2f\n','Brine volume fraction:       ' ,ice_prams.vbf);
   fprintf(logid,'%s%10.3e\n','Youngs modulus (Pa):        ' ,ice_prams.young);
   fprintf(logid,'%s%10.3e\n','Flexural strength (Pa):     ' ,ice_prams.sigma_c);
   fprintf(logid,'%s%10.3e\n','Breaking stress (Pa):       ' ,ice_prams.stress_c);
   fprintf(logid,'%s%10.3f\n','Breaking strain:            ' ,ice_prams.strain_c);
   fprintf(logid,'%s%5.2f\n','Drag RP (Pa.s/m):            ' ,ice_prams.drag_rp);
   fprintf(logid,'%s%5.2f\n','Viscoelastic WS (m^2/s):     ' ,ice_prams.visc_ws);
   fprintf(logid,'%s\n','***********************************************');
   fprintf(logid,'%s\n','');
   fclose(logid);
end

if TEST_INC_SPEC==1
   if params_in.DO_DISP; disp(' ');
   disp('Testing initial spectrum...'); end
   [wf.Hs,wf.Tp,wf.mwd] = fn_spectral_integrals(om_vec,wave_stuff.dirs,Sdir);
   vbls  = {'Hs','Tp','mwd'};
   for n=1:length(vbls)
      vbl   = vbls{n};
      if params_in.DO_DISP; disp(' ');
      disp(['comparing field: ',vbl]); end
      v1    = wf.(vbl);
      v2    = wave_fields.(vbl);
      diff  = abs(v2-v1);
      if params_in.DO_DISP; disp(['max diff: ',num2str(max(diff(:)))]);
      disp(' '); end
   end

   return
end

if params_in.STEADY==1
   S_inc = wave_stuff.dir_spec;
   theta = -pi/180*(90+wave_stuff.dirs);
   j_fwd = find(cos(theta)>=0);

   if isfield(wave_fields,'STEADY_MASK');
      WAVE_MASK2  = wave_fields.STEADY_MASK;
   else
      WAVE_MASK2        = 0*WAVE_MASK;
      WAVE_MASK2(1:3,:) = 1;
   end
end
%%
T  = 2*pi./om_vec;
if USE_EBS==1
   %% enhanced back-scatter method
   %% - backwards-scattered waves are left alone
   wave_stuff.dir_spec_scattered = 0*wave_stuff.dir_spec;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Water wavelength and wave speed
%% is a function only of wave period
wlng     = ice_prams.g.*T.^2./(2.*pi);
ap       = sqrt(ice_prams.g.*wlng./(2.*pi)); % Phase speed
ag       = ap./2;                  % Group speed

diagnostics.phase_speed = ap;
diagnostics.group_speed = ag;

ag_eff      = zeros(gridprams.nx,gridprams.ny,wave_stuff.nfreq);
ap_eff      = zeros(gridprams.nx,gridprams.ny,wave_stuff.nfreq);
wlng_ice    = zeros(gridprams.nx,gridprams.ny,wave_stuff.nfreq);
disp_ratio  = ones (gridprams.nx,gridprams.ny,wave_stuff.nfreq);
atten_nond  = zeros(gridprams.nx,gridprams.ny,wave_stuff.nfreq);
damping     = zeros(gridprams.nx,gridprams.ny,wave_stuff.nfreq);

% Display some parameters here (since initialisation can be slow)
Nice  = sum(ICE_MASK(:));
h_av  = sum(ice_fields.hice(find(ICE_MASK)))/Nice;
c_av  = sum(ice_fields.cice(find(ICE_MASK)))/Nice;
%%
Nwav  = sum(WAVE_MASK(:));
Hs_av = sum(wave_fields.Hs(:))/Nwav;
Tp_av = sum(wave_fields.Tp(:))/Nwav;

Info  = { '------------------------------------';
         ['c          = ' num2str(c_av)  ' const'];
         ['h          = ' num2str(h_av)  ' m const'];
         ['Hs         = ' num2str(Hs_av) ' m'];
         ['Tp         = ' num2str(Tp_av) ' s'];
         ['CFL        = ' num2str(params_in.CFL)];
         ['nfreq      = ' num2str(wave_stuff.nfreq)];
         ['ndir       = ' num2str(wave_stuff.ndir)];
         ['SCATMOD    = ' num2str(params_in.SCATMOD)];
         '------------------------------------';
         ' '};
if params_in.DO_DISP; disp(strvcat(Info)); end


for i = 1:gridprams.nx
for j = 1:gridprams.ny

%  if j==1
%     %%progress report - can be slow
%     disp([' - initialised ',num2str(i),' rows out of ',num2str(gridprams.nx)])
%  end

   if ICE_MASK(i,j)==1
      %% if ice is present:
      %% get ice wavelengths, group velocities,
      %% displacement ratios (convert wtr disp to ice),
      %% attenuation;
      if params_in.DO_ATTEN==1
         [damping_rp,kice,kwtr,int_adm,NDprams,...
            alp_scat,modT,argR,argT] =...
               RT_param_outer(om_vec,ice_fields.hice(i,j),ice_prams);
         %%
         if params_in.CHK_ATTEN==1
            %%check with old version
            %%TODO remove from infile
         end
         atten_nond(i,j,:) = alp_scat;
         damping(i,j,:)    = damping_rp;
      else
         [damping_rp,kice,kwtr,int_adm,NDprams] =...
            RT_param_outer(om_vec,ice_fields.hice(i,j),ice_prams);
         modT  = 1;
      end

      if params_in.USE_ICE_VEL==0
         %%use wtr group vel;
         ag_eff(i,j,:)  = ag;
         ap_eff(i,j,:)  = ap;
      else
         %%TODO check if this is correct
         %%weighted avg of ice and wtr group vel;
         ag_ice         = GEN_get_ice_groupvel(ice_fields.hice(i,j),T,Inf,ice_prams.young);
         ag_eff(i,j,:)  = ice_fields.cice(i,j)*ag_ice+...
                           +(1-ice_fields.cice(i,j))*ag;

         %%weighted avg of ice and wtr phase vel;
         ap_ice         = om_vec./k_ice;
         ap_eff(i,j,:)  = ice_fields.cice(i,j)*ap_ice+...
                           +(1-ice_fields.cice(i,j))*ap;
      end

      wlng_ice(i,j,:)   = 2*pi./kice;
      disp_ratio(i,j,:) = (kice./kwtr).*modT;
      %%
%     if (i==itest)&(j==jtest)
%        disp('om_vec,T,h')
%        disp([om_vec(1),T(1),ice_fields.hice(i,j)])
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

amax     = max(ag_eff(:));
dt       = params_in.CFL*gridprams.dx/max(ag_eff(:)); 
duration = params_in.duration_hours*3600;%%duration in seconds;
nt       = ceil(duration/dt);
dt       = duration/nt;%%reduce dt so it divides duration perfectly

% L     = max(gridprams.X(:))-min(gridprams.X(:));
% amin  = min(ag_eff(:));
% uc    = amin+.7*(amax-amin);

diagnostics.wave_travel_dist = duration*ag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display parameters
Info  = { '------------------------------------';
         ['Young''s modulus    = ' num2str(ice_prams.young,'%5.5e') ' Pa'];
         ['sigma_c            = '  num2str(ice_prams.sigma_c,'%5.5e') ' Pa'];
         ['fragility          = '  num2str(ice_prams.fragility)];
         ['strain_c           = '  num2str(ice_prams.strain_c,'%5.5e')];
         ['h                  = '  num2str(h_av) ' m const'];
         ['c                  = '  num2str(c_av) ' const'];
         ['Damping            = '  num2str(ice_prams.drag_rp) ' Pa.s/m'];
         [' '];
         ['Tp                 = '  num2str(Tp_av) ' s'];
         ['Hs                 = '  num2str(Hs_av) ' m'];
         ['nfreq              = '  num2str(wave_stuff.nfreq)];
         ['ndir               = '  num2str(wave_stuff.ndir)];
         ['SCATMOD            = '  num2str(params_in.SCATMOD)];
         [' '];
         ['FSD_OPT            = '  num2str(params_in.FSD_OPT)];
         ['BRK_OPT            = '  num2str(params_in.BRK_OPT)];
         [' '];
         ['CFL                = '  num2str(params_in.CFL)];
         ['dt                 = '  num2str(dt,'%1.1f') ' s'];
         ['nt                 = '  num2str(nt)];
         ['Time interval      = '  num2str(duration/3600,'%1.1f') ' h'];
         [' '];
         ['nx                 = '  num2str(gridprams.nx)];
         ['ny                 = '  num2str(gridprams.ny)];
         ['dx                 = '  num2str(gridprams.dx/1e3)    ' km'];
         ['dy                 = '  num2str(gridprams.dy/1e3)    ' km'];
         ['x extent           = '  num2str(gridprams.nx*gridprams.dx/1e3) ' km'];
         ['y extent           = '  num2str(gridprams.ny*gridprams.dy/1e3) ' km'];
         '------------------------------------';
         ' '};

if params_in.DO_DISP; disp(strvcat(Info)); end

%% Integration
t0       = now;   %%days
t0_fac   = 24*60; %%days to minutes
brkcrt  = zeros(gridprams.nx,gridprams.ny,nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%append to log file
if params_in.SV_LOG
   logid = fopen(log_file,'a');
   fprintf(logid,'%s\n',' ');
   fprintf(logid,'%s\n','***********************************************');
   fprintf(logid,'%s\n','Other Parameters:');
   fprintf(logid,'%s%6.1f\n','Time step (s):                    ',dt);
   fprintf(logid,'%s%4.3f\n','CFL number:                       ',params_in.CFL);
   fprintf(logid,'%s%5.2f\n','Maximum wave group velocity (m/s):',amax);
   fprintf(logid,'%s%4.4d\n','Number of time steps:             ',nt);
   fprintf(logid,'%s%5.2f\n','Time interval (h):                ',nt*dt/3600 );
   fprintf(logid,'%s\n','***********************************************');
   fprintf(logid,'%s\n',' ');

   fprintf(logid,'%s\n','***********************************************');
   fprintf(logid,'%s%4.4d%s%4.4d\n','Grid dimensions:                  ' ,...
      gridprams.nx,' ',gridprams.ny);
   fprintf(logid,'%s%4.1f%s%4.1f\n','Spatial resolution (km):          ' ,...
      gridprams.dx/1.0e3,' ',gridprams.dy/1.0e3);
   fprintf(logid,'%s%4.1f%s%4.1f\n','Extent of domain   (km):          ' ,...
      gridprams.nx*gridprams.dx/1.0e3,' ',gridprams.ny*gridprams.dy/1.0e3);

   fprintf(logid,'%s\n',' ');
   fprintf(logid,'%s%5.2f\n','Minimum period (s):               ',1/max(wave_stuff.freq) );
   fprintf(logid,'%s%5.2f\n','Maximum period (s):               ',1/min(wave_stuff.freq) );
   fprintf(logid,'%s%4.4d\n','Number of wave frequencies:       ',wave_stuff.nfreq);
   fprintf(logid,'%s%4.4d\n','Number of wave directions:        ',wave_stuff.ndir);
   fprintf(logid,'%s%5.2f\n','Directional resolution (degrees): ',360.0/(1.0*wave_stuff.ndir) );
   fprintf(logid,'%s\n','***********************************************');
   fprintf(logid,'%s\n',' ');
   fclose(logid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define weights for numerical quadrature;
if wave_stuff.nfreq>1%% weights for integral over frequency
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

if params_in.DIAG1d==1
   cols     = {'-k','-c','-m','-r','-g','-b'};
   labs1d_1 = {'\itx, \rmkm','{\itH}_{\rm s}, m'};
   labs1d_2 = {'\itx, \rmkm','{\itH}_{\rm s}^+, m'};
   labs1d_3 = {'\itx, \rmkm','{\itH}_{\rm s}^-, m'};
   labs1d_4 = {'\itx, \rmkm','c'};
end

if params_in.PLOT_INIT
   %%
   figure(1),clf;
   fn_fullscreen;
   fn_plot_ice(gridprams,ice_fields);
   if params_in.SV_FIG==1
      fig_dir  = [params_in.outdir,'/figs/init'];
      eval(['!mkdir -p ',fig_dir]);
      figname  = [fig_dir,'/ice_init.png'];
      saveas(gcf,figname);
   end
   pause(0.1);
   %%
   figure(2),clf;
   fn_fullscreen;
   fn_plot_waves(gridprams,wave_fields);
   if params_in.SV_FIG==1
      figname  = [fig_dir,'/waves_init.png'];
      saveas(gcf,figname);
   end
   pause(0.1);
   %%
   figure(3),clf;
   fn_fullscreen;
   Tc    = 12;%check this period
   jchq  = find(abs(T-Tc)==min(abs(T-Tc)));
   jdir  = round(wave_stuff.ndir/2);
   s1 = struct('dir',wave_stuff.dirs(jdir),...
               'period',Tc,...
               'Sdir',wave_stuff.dir_spec(:,:,jdir,jchq));
   %%
   if PLOT_OPT==1
      fn_plot_spec(gridprams.X,gridprams.Y,wave_fields.Hs,wave_fields.Tp,out_fields.Dmax,s1);
   else
      fn_plot_spec_2(gridprams.X,gridprams.Y,wave_fields.Hs,out_fields.tau_x,...
         out_fields.Dmax,out_fields.tau_y);
   end

   if params_in.SV_FIG==1
      figname  = [fig_dir,'/wim_init2d.png'];
      saveas(gcf,figname);
   end
   %if params_in.OPT==1
   %   subplot(2,2,1);
   %   hold on;
   %   x0 = min(gridprams.X(:));
   %   x1 = gridprams.X(find(WAVE_MASK(:,1)==0,1,'first'),1);
   %   % {x0/1e3,x1/1e3}
   %   yc = .3*max(gridprams.Y(:))/1e3;
   %   x_ = [min(gridprams.X(:)),x0,x0,x1,x1,max(gridprams.X(:))]/1e3;
   %   y_ = [0,0,yc*[1,1],0,0];
   %   plot(x_,y_,'k');
   %   plot(gridprams.X(:,1)/1e3,wave_fields.Hs(:,1)*yc,'--k');
   %   xlabel('$x$, km','interpreter','latex','fontsize',20); 
   %   ylabel('$\hat{H}_{s}$, m','interpreter','latex','fontsize',20)
   %   hold off;
   %end
   %%
   clear s1;
   %%
   if params_in.DIAG1d==1
      %% initial
      figure(4),clf;
      fn_fullscreen;
      loop_col = 1;
      %%
      subplot(4,1,4);
      fn_plot1d(gridprams.X(:,1)/1e3,ice_fields.cice(:,1),labs1d_4,cols{loop_col});
      %%
      subplot(4,1,1);
      %fn_plot1d(gridprams.X(:,1)/1e3,wave_fields.Hs(:,1),labs1d_1,cols{loop_col});
      fn_plot1d(gridprams.X(:,1)/1e3,mean(wave_fields.Hs,2),labs1d_1,cols{loop_col});%%average over y (columns)
      hold on;
      %%
      [Ep,Em,Et1,Et2]   = fn_split_energy(om_vec,wave_stuff.dirs,wave_stuff.dir_spec);
      %Hp                = 4*sqrt(Ep(:,1));
      %Hm                = 4*sqrt(Em(:,1));
      Hp                = 4*sqrt(mean(Ep,2));
      Hm                = 4*sqrt(mean(Em,2));
      if params_in.DIAG1d_OPT==0
         %Hs2   = 4*sqrt(Ep(:,1)+Em(:,1));%%add Ep + Em
         Hs2   = 4*sqrt(mean(Ep,2)+mean(Em,2));
      elseif params_in.DIAG1d_OPT==1
         %Hs2   = 4*sqrt(Et1(:,1));%%check const panel integration
         Hs2   = 4*sqrt(mean(Et1,2));
      elseif params_in.DIAG1d_OPT==2
         %Hs2   = 4*sqrt(Et2(:,1));%%check Simpson's rule integration
         Hs2   = 4*sqrt(mean(Et2,2));
      end
      fn_plot1d(gridprams.X(:,1)/1e3,Hs2,labs1d_1,['-',cols{loop_col}]);
      hold on;
      %%
      subplot(4,1,2);
      fn_plot1d(gridprams.X(:,1)/1e3,Hp,labs1d_2,cols{loop_col});
      hold on;
      fn_plot1d(gridprams.X(:,1)/1e3,Hs2,labs1d_2,['-',cols{loop_col}]);
      hold on;
      %%
      subplot(4,1,3);
      fn_plot1d(gridprams.X(:,1)/1e3,Hm,labs1d_3,cols{loop_col});
      hold on;
      %%
      loop_col = loop_col+1;
      if params_in.SV_FIG==1
         fig_dir  = [params_in.outdir,'/figs/init'];
         figname  = [fig_dir,'/wim_diag1d.png'];
      end
   end
   pause(0.1);
   %pause;
end

if SV_BIN & (params_in.MEX_OPT==0)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% save some fields as binaries to [params_in.outdir]/binaries
   Bdims = [gridprams.nx,gridprams.ny,wave_stuff.nfreq,wave_stuff.ndir];

   reps_ab  = params_in.dumpfreq;%%save every 10 time-steps if params_in.DO_CHECK_PROG==1
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% save grid files
   Fdir  = [params_in.outdir,'/binaries'];
   Froot = [Fdir,'/wim_grid'];
   %%
   pairs = {};
   pairs{end+1}   = {'X'         ,gridprams.X};
   pairs{end+1}   = {'Y'         ,gridprams.Y};
   pairs{end+1}   = {'scuy'      ,gridprams.scuy};
   pairs{end+1}   = {'scvx'      ,gridprams.scvx};
   pairs{end+1}   = {'scp2'      ,gridprams.scp2};
   pairs{end+1}   = {'scp2i'     ,gridprams.scp2i};
   pairs{end+1}   = {'LANDMASK'  ,gridprams.LANDMASK};
   %%
   fn_save_binary(Froot,Bdims,[],pairs);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   %% don't save binaries if MEX_OPT >0 :
   %% - saved by fortran instead
   SV_BIN   = 0;
end


if (SV_BIN==1) & (params_in.DO_CHECK_INIT==1)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% save init files
   Fdir  = [params_in.outdir,'/binaries'];
   Froot = [Fdir,'/wim_init'];
   %%
   pairs = {};
   pairs{end+1}   = {'cice',ice_fields.cice};
   pairs{end+1}   = {'hice',ice_fields.hice};
   pairs{end+1}   = {'Dmax',out_fields.Dmax};
   pairs{end+1}   = {'Hs'  ,wave_fields.Hs};
   pairs{end+1}   = {'Tp'  ,wave_fields.Tp};
   pairs{end+1}   = {'mwd' ,wave_fields.mwd};
   %%
   fn_save_binary(Froot,Bdims,year_info,pairs);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if (SV_BIN==1) & (params_in.DO_CHECK_PROG==1)

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% 1st prog file
   Fdir  = [params_in.outdir,'/binaries/prog'];
   Froot = [Fdir,'/wim_prog'];
   %%
   pairs = {};
   pairs{end+1}   = {'Dmax',out_fields.Dmax};
   pairs{end+1}   = {'tau_x',out_fields.tau_x};
   pairs{end+1}   = {'tau_y',out_fields.tau_y};
   pairs{end+1}   = {'Hs',wave_fields.Hs};
   pairs{end+1}   = {'Tp',wave_fields.Tp};
   pairs{end+1}   = {'mwd',wave_fields.mwd};
   %%
   fn_save_binary(Froot,Bdims,year_info,pairs);
   %eval(['!cat ',Froot,'.b'])
   %pause;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if params_in.DO_DISP; disp('BEGINNING MAIN INTEGRATION...'); end

%% =========================================================================
%% Different versions of mex functions
if params_in.MEX_OPT==1

   params_mex  = get_params_mex(params_in,duration,ice_prams,year_info);
   out_fields  = WIM2d_mex(params_mex,gridprams,ice_fields,wave_fields);

elseif params_in.MEX_OPT==2

   params_mex  = get_params_mex(params_in,duration,ice_prams,year_info);
   %%
   [out_fields,wave_stuff ]   =...
      WIM2d_mex(params_mex,gridprams,ice_fields,wave_fields,...
                  wave_stuff);

elseif params_in.MEX_OPT==3

   params_mex  = get_params_mex(params_in,duration,ice_prams,year_info);
   %%
   [out_fields,wave_stuff,mesh_e]   =...
      WIM2d_mex(params_mex,gridprams,ice_fields,wave_fields,...
                  wave_stuff,mesh_e);

%% end of options for mex functions
%% =========================================================================

else

%% =========================================================================
%% do time-stepping in matlab

   if params_in.DO_DISP; disp('Running pure matlab code'); end
   if COMP_F==1
      %% load prog binaries and compare saved wim_prog*.[ab] files
      %% to matlab results

      FF       = dir([compFdir,'wim_prog*.a'])
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
      F_fields.Dmax  = reshape( fread(aid,gridprams.nx*gridprams.ny,fmt), gridprams.nx,gridprams.ny );
      F_fields.tau_x = reshape( fread(aid,gridprams.nx*gridprams.ny,fmt), gridprams.nx,gridprams.ny );
      F_fields.tau_y = reshape( fread(aid,gridprams.nx*gridprams.ny,fmt), gridprams.nx,gridprams.ny );
      F_fields.Hs    = reshape( fread(aid,gridprams.nx*gridprams.ny,fmt), gridprams.nx,gridprams.ny );
      F_fields.Tp    = reshape( fread(aid,gridprams.nx*gridprams.ny,fmt), gridprams.nx,gridprams.ny );
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
         fn_pcolor(gridprams.X,gridprams.Y,v1,{'\itx, \rmkm','\ity, \rmkm',vc{1}});
         subplot(2,1,2);
         fn_pcolor(gridprams.X,gridprams.Y,v2,{'\itx, \rmkm','\ity, \rmkm',vc{2}});
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
            fn_pcolor(gridprams.X,gridprams.Y,Z,{'\itx, \rmkm','\ity, \rmkm',lbl{j}});
         end
         for j=3:3
            subplot(3,1,j);
            v1    = ice_fields.(vlist{j});
            v2    = F_fields.(vlist{j});
            Z     = 0*v2;
            jn    = find(abs(v2)>0);
            Z(jn) = 1-v1(jn)./v2(jn);%%relative difference
            fn_pcolor(gridprams.X,gridprams.Y,Z,{'\itx, \rmkm','\ity, \rmkm',lbl{j}});
         end
      end
   end

   %% also give progress report every 'reps' time steps;
   reps     = params_in.dumpfreq;
   GET_OUT  = 1;
   if GET_OUT
      Dmax_all         = zeros(gridprams.nx,gridprams.ny,1+floor(nt/reps));
      Dmax_all(:,:,1)  = out_fields.Dmax;
   end

   %nt = 13%%stop straight away for testing
   for n = 1:nt

      %%determine if we need to dump local diagnostics
      DUMP_DIAG   = (mod(n-1,reps)==0)&TEST_IJ;
      if DUMP_DIAG
         logfile2 = [log_dir2,'/WIMdiagnostics_local',year_info.date_string,'.txt'];
         logid2   = fopen(logfile2,'w');
         fprintf(logid2,'%s\n',[year_info.cdate,' # date']);
         fprintf(logid2,'%s\n',[year_info.ctime,' # time']);
         fprintf(logid2,'%d%s\n',model_day,' # model day');
         fprintf(logid2,'%13.5f%s\n',model_seconds,' # model seconds');
         fprintf(logid2,'%d%s\n',params_in.itest,' # itest');
         fprintf(logid2,'%d%s\n',params_in.jtest,' # jtest');
         fprintf(logid2,'%d%s\n',ICE_MASK(params_in.itest,params_in.jtest),' # ICE_MASK');
         fprintf(logid2,'%s\n',' ');
      end

      if params_in.DO_DISP
         disp([n nt]);
         disp(' ');
         Dice  = out_fields.Dmax(ice_fields.cice>0);
         Drng  = [min(Dice),max(Dice(Dice<300)),max(out_fields.Dmax(:))]
         clear Dice
      end

      %% spectral moments;
      mom0  = zeros(gridprams.nx,gridprams.ny);
      mom2  = zeros(gridprams.nx,gridprams.ny);
      mom0w = zeros(gridprams.nx,gridprams.ny);
      mom2w = zeros(gridprams.nx,gridprams.ny);

      %% wave stresses;
      tau_x = zeros(gridprams.nx,gridprams.ny);
      tau_y = zeros(gridprams.nx,gridprams.ny);

      %% mwd
      mwd   = zeros(gridprams.nx,gridprams.ny);

      %% variances of stress and strain;
      var_stress  = zeros(gridprams.nx,gridprams.ny);
      var_strain  = zeros(gridprams.nx,gridprams.ny);

      % %% test integrals;
      % var_boundary   = cell(1,length(Jy_boundary));
      % for r=1:length(Jy_boundary)
      %    var_boundary{r}   = 0;
      % end
      % var_boundary0  = var_boundary;  
      %%
      if params_in.STEADY==1
         for i = 1:gridprams.nx
         for j = 1:gridprams.ny
            %%top-up waves in wave mask if params_in.STEADY==1
            %%(steady-state solution);
            if WAVE_MASK2(i,j)>0 & params_in.STEADY==1
               wave_stuff.dir_spec(i,j,j_fwd,:)  = S_inc(i,j,j_fwd,:);
            end
         end
         end
      end

      % Dmean
      Dave  = 0*ICE_MASK;
      for i = 1:gridprams.nx
      for j = 1:gridprams.ny
         if out_fields.Dmax(i,j) < 200
            if params_in.FSD_OPT==0
               %%renormalisation group method
               Dave(i,j)  = floe_scaling(out_fields.Dmax(i,j),ice_prams,1);
            else
               %%power law distribution
               Dave(i,j)  = floe_scaling_smooth(out_fields.Dmax(i,j),ice_prams,1);
            end
         else
            %% uniform lengths
            Dave(i,j)  = out_fields.Dmax(i,j);
         end

         test_ij  = (i==params_in.itest)&(j==params_in.jtest);
         if DUMP_DIAG&test_ij&(ICE_MASK(i,j)>0)
            fprintf(logid2,'%s\n','Ice info: pre-breaking');
            fprintf(logid2,'%8.4f%s\n',ice_fields.cice(i,j),' # conc');
            fprintf(logid2,'%8.4f%s\n',ice_fields.hice(i,j),' # h, m');
            fprintf(logid2,'%8.4f%s\n',Dave(i,j),' # D_av, m');
            fprintf(logid2,'%8.4f%s\n',out_fields.Dmax(i,j),' # D_max, m');
            fprintf(logid2,'%s\n',' ');
            fprintf(logid2,'%s\n','# period, s | atten_dim, m^{-1} | damp_dim, m^{-1}');
         end
      end
      end
         
      for jw   = 1:wave_stuff.nfreq

         %% CALC DIMENSIONAL ATTEN COEFF;
         atten_dim   = 0*gridprams.X;
         damp_dim    = 0*gridprams.X;
         for i = 1:gridprams.nx
         for j = 1:gridprams.ny

            
            if ICE_MASK(i,j)>0 & params_in.DO_ATTEN==1

               %% get expected no of floes met per unit
               %%  distance if travelling in a line;
               c1d = ice_fields.cice(i,j)/Dave(i,j);%% floes per unit length;

               %% ENERGY attenuation coeff;
               atten_dim(i,j) = atten_nond(i,j,jw)*c1d;%%scattering
               damp_dim(i,j)  = 2*damping(i,j,jw)*ice_fields.cice(i,j);%%damping

   %           if (i==itest)&(j==jtest)
   %              disp(['Hs (pre)   = ',num2str(wave_fields.Hs(i,j)),' m']);
   %              disp(['Tp (pre)   = ',num2str(wave_fields.Tp(i,j)),' s']);
   %              disp(['Dmax       = ',num2str(out_fields.Dmax(i,j)),' m']);
   %              disp(['Dave       = ',num2str(Dave),' m']);
   %              disp(['c1d        = ',num2str(c1d)]);
   %              disp(['q_scat     = ',num2str(atten_dim(i,j),'%7.7e'),' /m']);
   %              disp(['q_abs      = ',num2str(damp_dim(i,j),'%7.7e'),' /m']);
   %           end
               test_ij  = (i==params_in.itest)&(j==params_in.jtest);
               if DUMP_DIAG&test_ij
                  fprintf(logid2,'%8.4f%s%13.6e%s%13.6e\n',...
                     1/wave_stuff.freq(jw),' | ',...
                     atten_dim(i,j),' | ',...
                     damp_dim(i,j));
               end

            end
         end% j
         end% i
         % max(ag_eff(:))
         % max(atten_dim(:))
         % max(damp_dim(:))
         %1e3*max(ag_eff(:))*max(atten_dim(:))
         % GEN_pause

         s1.ndir        = wave_stuff.ndir;
         s1.wavdir      = wave_stuff.dirs;
         %s1.Sdir        = reshape( wave_stuff.dir_spec(:,:,jw,:), gridprams.nx,gridprams.ny,ndir);
         s1.Sdir        = wave_stuff.dir_spec(:,:,:,jw);
         if USE_EBS == 1;
            %% "already-scattered" energy - doesn't get scattered anymore
            s1.Sdir_scattered = wave_stuff.dir_spec_scattered(:,:,:,jw);
         end
         s1.ag_eff      = ag_eff(:,:,jw);
         s1.atten_dim   = atten_dim;
         s1.damp_dim    = damp_dim;
         s1.ICE_MASK    = ICE_MASK;

         % ==============================================================
         %%advection;
         % ADV_OPT  = 0;%%zeros outside real domain
         % ADV_OPT  = 1;%%periodic in x,y
         % ADV_OPT  = 2;%%periodic in y only
         theta = -pi/180*(90+wave_stuff.dirs);
         if adv_options.ADV_DIM==2
            %%2d advection
            for jth  = 1:s1.ndir
               %% set the velocities
               u  = s1.ag_eff*cos(theta(jth));
               v  = s1.ag_eff*sin(theta(jth));

               %% call advection routine
               s1.Sdir(:,:,jth)  = waveadv_weno(...
                  s1.Sdir(:,:,jth),u,v,gridprams,dt,adv_options);
               if USE_EBS == 1
                  s1.Sdir_scattered(:,:,jth)  = waveadv_weno(...
                     s1.Sdir_scattered(:,:,jth),u,v,gridprams,dt,adv_options);
               end
            end
         else
            %%1d advection - 1 row at a time
            for jy=1:gridprams.ny
               for jth  = 1:s1.ndir
                  %% set the velocity
                  u  = s1.ag_eff(:,jy)*cos(theta(jth));

                  %% call advection routine
                  s1.Sdir(:,jy,jth) = waveadv_weno_1d(...
                     s1.Sdir(:,jy,jth),u,gridprams,dt,adv_options);
                  if USE_EBS == 1
                     s1.Sdir_scattered(:,jy,jth) = waveadv_weno_1d(...
                        s1.Sdir_scattered(:,jy,jth),u,gridprams,dt,adv_options);
                  end
               end
            end
         end
         % ==============================================================

         if wave_stuff.ndir==1
            if params_in.SCATMOD~=0
               if params_in.DO_DISP; disp('warning: changing params_in.SCATMOD option as not enough directions');
               disp(['(ndir = ',num2str(wave_stuff.ndir)]); end
            end
            params_in.SCATMOD  = 0;
         end

         if params_in.SCATMOD==0
            %% Simple attenuation scheme - doesn't conserve scattered energy
            [wave_stuff.dir_spec(:,:,:,jw),S_freq,tau_x_om,tau_y_om] = ...
               atten_simple(gridprams,ice_prams,s1,dt);
            clear s1 S_out;
         elseif params_in.SCATMOD==1
            %% same as params_in.SCATMOD==0, but scattered energy
            %% is distributed isotropically
            [wave_stuff.dir_spec(:,:,:,jw),S_freq,tau_x_om,tau_y_om] = ...
               atten_isotropic(gridprams,ice_prams,s1,dt);
            clear s1 S_out;
         elseif floor(params_in.SCATMOD)==2
            %% same as params_in.SCATMOD==1, but scattered energy
            %% is distributed non-isotropically
            %% pass in SHPOPT as decimal place eg 2.1 or 2.2 or 2.3
            SHPOPT   = round(10*(params_in.SCATMOD-2));
            %SHPOPT   = 1;%% cos^2
            %SHPOPT   = 2;%% cos^4
            %SHPOPT   = 3;%% cos^8
            [wave_stuff.dir_spec(:,:,:,jw),S_freq,tau_x_om,tau_y_om] = ...
               atten_noniso(gridprams,ice_prams,s1,dt,SHPOPT);
            clear s1 S_out;    
         elseif params_in.SCATMOD==-1
            %% Simple attenuation scheme - does conserve scattered energy
            [wave_stuff.dir_spec(:,:,:,jw),S_freq,tau_x_om,tau_y_om] = ...
               atten_simple_conserve(gridprams,ice_prams,s1,dt);
            clear s1 S_out;
         elseif params_in.SCATMOD==3
            %% isotropic scattering, enhanced backscatter (EBS)
            [wave_stuff.dir_spec(:,:,:,jw),...
             wave_stuff.dir_spec_scattered(:,:,:,jw),...
             S_freq,...
             tau_x_om,...
             tau_y_om] = ...
                EBS_atten_iso(gridprams,ice_prams,s1,dt);
            clear s1 S_out;
         end

         %% integrate stress densities over frequency
         %% TODO: check if this is correct for ice-covered water
         tmp   = ice_prams.rhowtr*ice_prams.g*tau_x_om./ap_eff(:,:,jw);  %%[Pa*s]
         tau_x = tau_x+wt_om(jw)*tmp;                 %%[Pa]
         tmp   = ice_prams.rhowtr*ice_prams.g*tau_y_om./ap_eff(:,:,jw);  %%[Pa*s]
         tau_y = tau_y+wt_om(jw)*tmp;                 %%[Pa]

         %% mwd integrals
         tmp3  = wave_stuff.dir_spec(:,:,:,jw);
         if USE_EBS
            tmp3  = tmp3 + wave_stuff.dir_spec_scattered(:,:,:,jw);
         end
         mwd_om   = calc_mwd_1freq(wave_stuff.dirs,tmp3);
         mwd      = mwd+wt_om(jw)*mwd_om;
         %GEN_pause

         %% INTEGRALS FOR BREAKING PROB:
         for i = 1:gridprams.nx
         for j = 1:gridprams.ny

            %% INTEGRATE SPECTRUM OVER DIRECTION;
            %S_freq   = wt_theta'*squeeze(wave_stuff.dir_spec(i,j,jw,:));

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
                                     (k_ice^2*ice_fields.hice(i,j)/2)^2 );
               var_strain(i,j)   = var_strain(i,j)+...
                                    + wt_om(jw)*strain_density;
            end
         end%% end spatial loop x;
         end%% end spatial loop y;

      end%% end spectral loop;
      %mom0,mom2,wlng_ice,return

      %%calc Hs, Tw (into out_fields.Tp)
      if params_in.REF_Hs_ICE==1
         %% diagnostic variables Hs/Tp related to ice displacment in ice-covered areas
         out_fields.Hs       = 4*sqrt(mom0);%%diagnostic variable - for waves in water
         out_fields.Tp       = 0*gridprams.X;
         jnz                  = find(mom2>0);
         out_fields.Tp(jnz)  = 2*pi*sqrt(mom0(jnz)./mom2(jnz));
      else
         %% diagnostic variables Hs/Tp related to water displacment in ice-covered areas
         out_fields.Hs       = 4*sqrt(mom0w);%%diagnostic variable - for waves in water
         out_fields.Tp       = 0*gridprams.X;
         jnz                  = find(mom2w>0);
         out_fields.Tp(jnz)  = 2*pi*sqrt(mom0w(jnz)./mom2w(jnz));
      end

      %% wave stresses
      out_fields.tau_x  = tau_x;
      out_fields.tau_y  = tau_y;

      %% mwd
      out_fields.mwd    = mwd;

      if DUMP_DIAG
         fprintf(logid2,'%s\n',' ');
         fprintf(logid2,'%13.6e,%s\n',mom0w(params_in.itest,params_in.jtest),' # mom0w, m^2');
         fprintf(logid2,'%13.6e,%s\n',mom2w(params_in.itest,params_in.jtest),' # mom2w, m^2/s^2');
         fprintf(logid2,'%13.6e,%s\n',mom0(params_in.itest,params_in.jtest),' # mom0, m^2');
         fprintf(logid2,'%13.6e,%s\n',mom2(params_in.itest,params_in.jtest),' # mom2, m^2/s^2');
         fprintf(logid2,'%8.4f%s\n',out_fields.Hs(params_in.itest,params_in.jtest),' # Hs, m')
         fprintf(logid2,'%10.4f%s\n',out_fields.Tp(params_in.itest,params_in.jtest),' # Tp, s')
         fprintf(logid2,'%10.4f%s\n',mwd(params_in.itest,params_in.jtest),' # mwd, deg');
         fprintf(logid2,'%13.6e%s\n',tau_x(params_in.itest,params_in.jtest),' # tau_x, Pa');
         fprintf(logid2,'%13.6e%s\n',tau_y(params_in.itest,params_in.jtest),' # tau_y, Pa');
         fprintf(logid2,'%s\n',' ');
      end


      %% FINALLY DO FLOE BREAKING;
      %% - ON GRID
      P_crit0  = 0;
      if params_in.BRK_OPT==0
         P_crit0  = 1;
      end

      for i=1:gridprams.nx
      for j=1:gridprams.ny
         if ICE_MASK(i,j)==1 & mom0(i,j)>0
            %% only try breaking if ice is present
            %%  & some waves have arrived;

            %% significant strain amp
            sig_strain  = 2*sqrt(var_strain(i,j));

            %%  probability of critical strain
            %%  being exceeded from Rayleigh distribution;
            Pstrain  = exp( -ice_prams.strain_c^2/(2*var_strain(i,j)) );
            P_crit   = P_crit0+exp(-1);%%this is critical prob if monochromatic wave

            %% FLOE BREAKING:
            BREAK_CRIT     = ( Pstrain>=P_crit );%%breaks if larger than this
            brkcrt(i,j,n)  = BREAK_CRIT;

            if BREAK_CRIT
               %% use crest period to work out wavelength
               %% - half this is max poss floe length;
               T_crit   = out_fields.Tp(i,j);
               if 0
                  %%get wavelength directly
                  wlng_crest  = ...
                     GEN_get_ice_wavelength(ice_fields.hice(i,j),T_crit,Inf,ice_prams.young);
               else
                  %% interpolate (to check fortran code)
                  %% NB slight difference due to RT_param_outer
                  %% using finite depth wavelength approx to
                  %% inifinite depth one
                  om       = 2*pi/T_crit;
                  om_min   = om_vec(1);
                  om_max   = om_vec(end);
                  if (om<=om_min)
                     wlng_crest  = wlng_ice(i,j,1);
                     %tst_crest   = [wlng_crest,...
                     %   GEN_get_ice_wavelength(ice_fields.hice(i,j),T_crit,Inf,ice_prams.young)]
                  elseif (om>=om_max)
                     wlng_crest  = wlng_ice(i,j,wave_stuff.nfreq);
                  else
                     jcrest      = floor((om-om_min+dom)/dom);
                     om1         = 2*pi*wave_stuff.freq(jcrest);
                     lam1        = wlng_ice(i,j,jcrest);
                     lam2        = wlng_ice(i,j,jcrest+1);
                     wlng_crest  = lam1+(om-om1)*(lam2-lam1)/dom;
                  end
               end

               Dc                   = max(ice_prams.Dmin,wlng_crest/2);
               out_fields.Dmax(i,j) = min(Dc,out_fields.Dmax(i,j));
               
            end%% end breaking action;

            % ==================================================================================

            if 0%i==11 & j==1
               BREAK_CRIT
               out_fields.Dmax(i,j)
            end
            
            if params_in.BRK_OPT==0
             P_crit_        = exp(-1);
             BREAK_CRIT     = ( Pstrain>=P_crit_ );
             brkcrt(i,j,n)  = BREAK_CRIT;
             if BREAK_CRIT
              if isnan(diagnostics.break_max)
               diagnostics.break_max = gridprams.X(i);
              else
               if gridprams.X(i)>diagnostics.break_max
                diagnostics.break_max = gridprams.X(i);
               end
              end
             end
            end

            test_ij  = (i==params_in.itest)&(j==params_in.jtest);
            if DUMP_DIAG&test_ij
               fprintf(logid2,'%s\n',' ');
               fprintf(logid2,'%s\n','Ice info: post-breaking');
               fprintf(logid2,'13.6%e%s\n',Pstrain,' # P_strain');
               fprintf(logid2,'13.6%e%s\n',P_crit,' # P_crit');
               fprintf(logid2,'%10.5f%s\n',wlng_crest,' # peak wavelength, m');
               fprintf(logid2,'%9.5f%s\n',out_fields.Dmax(i,j),' # D_max, m');
            end
            
         elseif WTR_MASK(i,j)==1%% only water present
            out_fields.Dmax(i,j)   = 0;
         end

      end%% end spatial loop j in y;
      end%% end spatial loop i in x;

      if DUMP_DIAG
         fclose(logid2);
      end

%% end of time stepping in matlab
%% =========================================================================

%% ==================================================================================
%% extra task when coupling to neXtSIM:
      if INTERP_MESH==1
         %% FOR neXtSIM COUPLING
         %% - DO FLOE BREAKING ON MESH
         X  = gridprams.X.'/1e3;%take transpose, change to km
         Y  = gridprams.Y.'/1e3;%take transpose, change to km

         %% choose interpolation order
         %meth  = 'nearest';
         %meth  = 'linear';
         %meth  = 'spline';
         meth  = 'cubic';

         %% get ice elements
         jice     = find(mesh_e.c>0);
         thick_e  = mesh_e.h(jice)./mesh_e.c(jice);%absolute thickness

         %% Interp mom0,mom2  -> Tp
         mom0_e      = interp2(X,Y,mom0.',mesh_e.xe(jice),mesh_e.ye(jice),meth);
         mom2_e      = interp2(X,Y,mom2.',mesh_e.xe(jice),mesh_e.ye(jice),meth);
         jwav        = find(mom2_e>0);%where waves are
         Tp_e        = 0*mom0_e;
         Tp_e(jwav)  = 2*pi*sqrt(mom0_e(jwav)./mom2_e(jwav));
         %%

         %% Interp var_strain -> (P_strain>P_crit)
         var_strain_e   = interp2(X,Y,var_strain.',mesh_e.xe(jice),mesh_e.ye(jice),meth);
         %%
         for loop_j=1:length(jice)
            vse      = var_strain_e(loop_j);
            Pstrain  = exp( -ice_prams.strain_c^2/(2*vse) );
            if Pstrain>P_crit&vse>0
               % {var_strain_e,Pstrain,thick_e(loop_j),Tp_e(loop_j)}
               jl          = jice(loop_j);
               wlng_crest  = ...
                  GEN_get_ice_wavelength(thick_e(loop_j),Tp_e(loop_j),Inf,ice_prams.young);
               %%
               Dmax              = sqrt(mesh_e.c(jl)/mesh_e.Nfloes(jl));
               Dmax              = max(ice_prams.Dmin,min(Dmax,wlng_crest/2));
               mesh_e.Nfloes(jl) = mesh_e.c(jl)/Dmax^2;
               if mesh_e.DAMAGE_OPT==1
                  mesh_e.broken(jl) = 1;
               end
               %disp('Breaking on mesh')
            end
         end%%loop over ice elements

         if params_in.DO_DISP
            Dmesh = sqrt(mesh_e.c(jice)./mesh_e.Nfloes(jice));
            DrngM = [min(Dmesh),max(Dmesh(Dmesh<300)),max(Dmesh)]
         end
      end%% INTERP_MESH==1

%% end of work for coupling to neXtSIM
%% ==================================================================================


%% ==================================================================================
%% print info, plotting...

      jmiz  = find((out_fields.Dmax>0)&(out_fields.Dmax<250));
      % {jmiz}

      %% progress report;
      if round(n/reps)==(n/reps)

         if GET_OUT
            Dmax_all(:,:,n/reps) = out_fields.Dmax;
         end

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if params_in.DO_DISP; disp('##############################################################'); end
         t1 = now;
         if params_in.DO_DISP; disp([num2str(n),' time steps done, out of ',num2str(nt)]);
         disp(['Time taken (mins)      : ' ,num2str(t0_fac*(t1-t0))]);
         disp(['Model time passed (h)  : ' ,num2str(n*dt/3600.)]);
         disp('##############################################################');
         disp(' '); end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         if params_in.DO_DISP; disp(strvcat(Info)); end

         if params_in.PLOT_PROG
            if params_in.DIAG1d==0
               figure(3),clf;
               fn_fullscreen;
               %%
               if PLOT_OPT==1
                  s1 = struct('dir',wave_stuff.dirs(jdir),...
                              'period',Tc,...
                              'Sdir',wave_stuff.dir_spec(:,:,jdir,jchq));
                  fn_plot_spec(gridprams.X,gridprams.Y,out_fields.Hs,out_fields.Tp,out_fields.Dmax,s1);
               else
                  fn_plot_spec_2(gridprams.X,gridprams.Y,out_fields.Hs,out_fields.tau_x,...
                     out_fields.Dmax,out_fields.tau_y);
               end

               if params_in.SV_FIG==1

                  cnt               = num2str(nt);
                  cnt(:)            = '0';
                  cn                = num2str(n);
                  lc                = length(cn);
                  cnt(end+1-lc:end) = cn;

                  fig_dir  = [params_in.outdir,'/figs/prog'];
                  figname  = [fig_dir,'/wim_prog2d_',cnt,'.png'];
                  saveas(gcf,figname);
               end

               %if params_in.OPT==1
               %   subplot(2,2,1);
               %   hold on;
               %   x0 = min(gridprams.X(:))+uc*n*dt;
               %   x1 = gridprams.X(find(WAVE_MASK(:,1)==0,1,'first'),1)+uc*n*dt;
               %   % {uc,x0/1e3,x1/1e3}
               %   yc = .3*max(Y(:))/1e3;
               %   x_ = [min(gridprams.X(:)),x0,x0,x1,x1,max(gridprams.X(:))]/1e3;
               %   y_ = [0,0,yc*[1,1],0,0];
               %   plot(x_,y_,'k');
               %   plot(gridprams.X(:,1)/1e3,out_fields.Hs(:,1)*yc,'--k');
               %   hold off;
               %   %%
               %   subplot(2,2,2);
               %   hold on;
               %   plot(x_,y_,'k');
               %   plot(gridprams.X(:,1)/1e3,out_fields.Hs(:,1)*yc,'--k');
               %   hold off;
               %end
            else
               %% params_in.DIAG1d==1
               %% during run
               if 1
                  figure(4);
                  %% check symmetry
                  [hmax,imax] = max(out_fields.Hs(:,1));
                  hp          = out_fields.Hs(imax,:);
                  yp          = gridprams.Y(imax,:);
                  plot(yp/1e3,hp);
                  ttl   = title(['max h = ',num2str(max(out_fields.Hs(:))),...
                                 '; x = ',num2str(gridprams.X(imax,1)/1.0e3),'km']);
                  GEN_font(ttl);
                  GEN_proc_fig('y, km','H_s, m')

                  if params_in.SV_FIG==1

                     cnt               = num2str(nt);
                     cnt(:)            = '0';
                     cn                = num2str(n);
                     lc                = length(cn);
                     cnt(end+1-lc:end) = cn;

                     fig_dir  = [params_in.outdir,'/figs/prog'];
                     figname  = [fig_dir,'/wim_prog1d_symm',cn,'.png'];
                     saveas(gcf,figname);
                  end
               else
                  %%check partition of fwd and back energy
                  figure(4);
                  subplot(4,1,1);
                  hold on;
                  %%
                  %fn_plot1d(gridprams.X(:,1)/1e3,out_fields.Hs(:,1),labs1d_1,cols{loop_col});
                  fn_plot1d(gridprams.X(:,1)/1e3,mean(out_fields.Hs,2),labs1d_1,cols{loop_col});%%average over y (columns)
                  hold on;
                  %%
                  [Ep,Em,Et1,Et2]   = fn_split_energy(om_vec,wave_stuff.dirs,wave_stuff.dir_spec);
                  %Hp                = 4*sqrt(Ep(:,1));
                  %Hm                = 4*sqrt(Em(:,1));
                  Hp                = 4*sqrt(mean(Ep,2));
                  Hm                = 4*sqrt(mean(Em,2));
                  if params_in.DIAG1d_OPT==0
                     %Hs2   = 4*sqrt(Ep(:,1)+Em(:,1));%%add Ep + Em
                     Hs2   = 4*sqrt(mean(Ep,2)+mean(Em,2));
                  elseif params_in.DIAG1d_OPT==1
                     %Hs2   = 4*sqrt(Et1(:,1));%%check const panel integration
                     Hs2   = 4*sqrt(mean(Et1,2));
                  elseif params_in.DIAG1d_OPT==2
                     %Hs2   = 4*sqrt(Et2(:,1));%%check Simpson's rule integration
                     Hs2   = 4*sqrt(mean(Et2,2));
                  end
                  fn_plot1d(gridprams.X(:,1)/1e3,Hs2,labs1d_1,['-',cols{loop_col}]);
                  hold on;
                  %%
                  subplot(4,1,2);
                  fn_plot1d(gridprams.X(:,1)/1e3,Hp,labs1d_2,cols{loop_col});
                  hold on;
                  fn_plot1d(gridprams.X(:,1)/1e3,Hs2,labs1d_2,['-',cols{loop_col}]);
                  hold on;
                  %%
                  subplot(4,1,3);
                  fn_plot1d(gridprams.X(:,1)/1e3,Hm,labs1d_3,cols{loop_col});
                  hold on;
                  %%
                  loop_col = loop_col+1;
                  if loop_col>length(cols)
                     loop_col = 1;
                  end

                  if params_in.SV_FIG==1

                     cnt               = num2str(nt);
                     cnt(:)            = '0';
                     cn                = num2str(n);
                     lc                = length(cn);
                     cnt(end+1-lc:end) = cn;

                     fig_dir  = [params_in.outdir,'/figs/prog'];
                     figname  = [fig_dir,'/wim_prog1d_Epart',cn,'.png']
                     saveas(gcf,figname);
                  end
               end
            end

            if COMP_F==1
               %% load prog binaries and compare saved wim_prog*.[ab] files
               %% to matlab results
               nnn            = num2str(n,'%3.3d');
               afile          = [compFdir,nnn,'.a'];
               aid            = fopen(afile,'rb');
               F_fields.Dmax  = reshape( fread(aid,gridprams.nx*gridprams.ny,fmt), gridprams.nx,gridprams.ny );
               F_fields.tau_x = reshape( fread(aid,gridprams.nx*gridprams.ny,fmt), gridprams.nx,gridprams.ny );
               F_fields.tau_y = reshape( fread(aid,gridprams.nx*gridprams.ny,fmt), gridprams.nx,gridprams.ny );
               F_fields.Hs    = reshape( fread(aid,gridprams.nx*gridprams.ny,fmt), gridprams.nx,gridprams.ny );
               F_fields.Tp    = reshape( fread(aid,gridprams.nx*gridprams.ny,fmt), gridprams.nx,gridprams.ny );
               fclose(aid);
            end

            clear s1;
            pause(0.1);
         end
         %%
         if 0
            Jlook = nx-16+(1:16);
            out_fields.Dmax(edge:edge+9,Jlook)
         end
      end


      for jw=[]%11
         [n,T(jw)],%[ag_ice(jw,1:11,1)]
         testSatt=S(1:11,1,jw,jmwd)
         pause
      end

      if (SV_BIN==1)&(mod(n,reps_ab)==0)&(params_in.DO_CHECK_PROG==1)
         %% save matlab files as binaries
         %% to matlab results
         Fdir  = [params_in.outdir,'/binaries/prog'];
         Froot = [Fdir,'/wim_prog'];
         %%
         pairs = {};
         pairs{end+1}   = {'Dmax' ,out_fields.Dmax};
         pairs{end+1}   = {'tau_x',out_fields.tau_x};
         pairs{end+1}   = {'tau_y',out_fields.tau_y};
         pairs{end+1}   = {'Hs'   ,out_fields.Hs};
         pairs{end+1}   = {'Tp'   ,out_fields.Tp};
         pairs{end+1}   = {'mwd'  ,out_fields.mwd};
         %%
         fn_save_binary(Froot,Bdims,year_info,pairs);
      end

%% end of time step
%% ==================================================================================

      %% update time
      model_seconds  = model_seconds+dt;
      ndays_jump     = floor(model_seconds/24/3600);
      model_seconds  = model_seconds-ndays_jump*24*3600;
      model_day      = model_day+ndays_jump;
      year_info      = model_time_to_year_info(model_day,model_seconds);

   end%% end time loop
end%%params_in.MEX_OPT==0 option


if (SV_BIN==1)&(params_in.DO_CHECK_FINAL==1)
   %% save matlab files as binaries
   %% to matlab results
   Fdir  = [params_in.outdir,'/binaries'];
   Froot = [Fdir,'/wim_out'];
   %%
   pairs = {};
   pairs{end+1}   = {'Dmax' ,out_fields.Dmax};
   pairs{end+1}   = {'tau_x',out_fields.tau_x};
   pairs{end+1}   = {'tau_y',out_fields.tau_y};
   pairs{end+1}   = {'Hs'   ,out_fields.Hs};
   pairs{end+1}   = {'Tp'   ,out_fields.Tp};
   pairs{end+1}   = {'mwd'  ,out_fields.mwd};
   %%
   fn_save_binary(Froot,Bdims,year_info,pairs);
end

t1 = now;

if (params_in.OPT==1)|(params_in.OPT==3)
   %%TODO make an explicit parameter to get diagnostics like Wmiz in this way
   Nchk     = ceil(gridprams.ny/2);%%=1 if ny=1
   Dmax_j   = out_fields.Dmax(:,Nchk);
   jmiz     = find((Dmax_j>0)&(Dmax_j<250));
   Wmiz     = gridprams.dx*length(jmiz);
   %%
   diagnostics.MIZ_width   = Wmiz;
   %%
   Dmax_min = min(Dmax_j);
   Dmax_max = max(Dmax_j);
   %%
   if params_in.DO_DISP; disp(' ');
   disp(['MIZ width = ',num2str(Wmiz/1e3),' km']); end
end

taux_min = min(out_fields.tau_x(:));
taux_max = max(out_fields.tau_x(:));
tauy_min = min(out_fields.tau_x(:));
tauy_max = max(out_fields.tau_y(:));
if params_in.DO_DISP; disp(['max tau_x = ',num2str(taux_max,'%0.10e'),' Pa']);
disp(['max tau_y = ',num2str(tauy_max,'%0.10e'),' Pa']);
disp(' '); end
diagnostics.taux_min = taux_min;
diagnostics.taux_max = taux_max;
diagnostics.tauy_min = tauy_min;
diagnostics.tauy_max = tauy_max;

%%append to log file
if params_in.SV_LOG
   logid = fopen(log_file,'a');
   fprintf(logid,'%s\n','***********************************************');
   fprintf(logid,'%s\n','Diagnostics:');
   if (params_in.OPT==1)|(params_in.OPT==3)
      fprintf(logid,'%s%9.4f\n','MIZ width (km): ',Wmiz);
      fprintf(logid,'%s%9.4f%s%9.4f\n','Dmax range in MIZ (m): ',...
         Dmax_min,' ',Dmax_max);
   end
   fprintf(logid,'%s%13.6e%s%13.6e\n','tau_x range (Pa): ',...
      taux_min,' ',taux_max);
   fprintf(logid,'%s%13.6e%s%13.6e\n','tau_y range (Pa): ',...
      tauy_min,' ',tauy_max);
   fprintf(logid,'%s\n','***********************************************');

   fprintf(logid,'%s\n',' ');
   fprintf(logid,'%s\n','***********************************************');
   fprintf(logid,'%s%7.1f\n','Elapsed time (min): ',t0_fac*(t1-t0));
   fprintf(logid,'%s\n','***********************************************');
   fprintf(logid,'%s\n',' ');
   fclose(logid);
end

if params_in.SV_SPEC
   %% save final directional spectrum
   wavdir   = wave_stuff.dirs;
   freq_vec = om_vec/2/pi;
   Sdir     = wave_stuff.dir_spec;
   Hs       = out_fields.Hs;
   cice     = ice_fields.cice;
   hice     = ice_fields.hice;
   save([params_in.outdir,'/Sdir.mat'],...
         'gridprams','S_inc',...
         'wavdir','freq_vec','Sdir',...
         'Hs','cice','hice');
   clear wavdir freq_vec Sdir Hs cice hice;
end

if TEST_FINAL_SPEC==1
   if params_in.DO_DISP; disp(' ');
   disp('Testing final spectrum...'); end

   if 0
      if params_in.DO_DISP; disp('(Check integrals of output spectrum vs output Hs,Tp,etc)');end
      %% check consistency of Sdir with Hs,Tp
      [wf.Hs,wf.Tp,wf.mwd] = fn_spectral_integrals(om_vec,wave_stuff.dirs,wave_stuff.dir_spec);

      vbls  = {'Hs','Tp'};%,'mwd'};%mwd currently not updated
      for n=1:length(vbls)
         vbl   = vbls{n};
         if params_in.DO_DISP; disp(' ');
         disp(['comparing field: ',vbl]); end
         v1    = wf.(vbl);
         v2    = out_fields.(vbl);
         diff  = abs(v2-v1);
         if params_in.DO_DISP; disp(['max diff: ',num2str(max(diff(:)))]);
         disp(' '); end
      end
   elseif params_in.MEX_OPT>0
      if params_in.DO_DISP; disp('(Check outputs vs values in binary files)'); end
      of2         = fn_check_final(outdir);%%set in infile_dirs.txt
      of1.Hs      = out_fields.Hs;
      of1.Tp      = out_fields.Tp;
      of1.tau_x   = out_fields.tau_x;
      of1.tau_y   = out_fields.tau_y;
      of1.Dmax    = out_fields.Dmax;
      %%
      vbls  = {'Hs','Tp','tau_x','tau_y','Dmax'};%,'mwd'};%mwd currently not updated
      for n=1:length(vbls)
         vbl   = vbls{n};
         if params_in.DO_DISP; disp(' ');
         disp(['comparing field: ',vbl]); end
         v1    = of1.(vbl);
         v2    = of2.(vbl);
         diff  = abs(v2-v1);
         if params_in.DO_DISP; disp(['max diff: ',num2str(max(diff(:)))]);
         disp(' '); end
      end
   end

   return
end

if params_in.PLOT_FINAL%%check exponential attenuation
   figure(3),clf;
   fn_fullscreen;
   if PLOT_OPT==1
      s1 = struct('dir',wave_stuff.dirs(jdir),...
                  'period',Tc,...
                  'Sdir',wave_stuff.dir_spec(:,:,jdir,jchq));
      fn_plot_spec(gridprams.X,gridprams.Y,out_fields.Hs,out_fields.Tp,out_fields.Dmax,s1);
   else
      fn_plot_spec_2(gridprams.X,gridprams.Y,out_fields.Hs,out_fields.tau_x,...
         out_fields.Dmax,out_fields.tau_y);
   end
   %%
   if 0
      %% figure testing how 1d results are (only appropriate for 1d geometries)
      figure(4),clf;
      fn_fullscreen;
      xx = gridprams.X(:,1);
      if 0
         vbl   = 'Hs';
         Vbl   = out_fields.(vbl);
      else
         vbl   = 'tau_x';
         Vbl   = ice_fields.(vbl);
      end

      if 1
         subplot(2,1,2)
         yy    = gridprams.Y(1,:);
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
   if params_in.DIAG1d==1
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
               if params_in.DO_DISP; disp(['opening ',dfil,'']); end
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
                  figure(4);
                  ix = find(abs(xx_f-xp)==min(abs(xx-xp)));
                  subplot(2,1,2)
                  hold on;
                  for loop_ix=1:length(ix)
                     plot(yy/1e3,Hs_f(ix(loop_ix))+0*yy,'--m');
                  end
                  figure(5);
               end
            else
               if params_in.DO_DISP; disp([dfil,' not present']);
               disp('To create, run ../../fortran/run/fig_scripts/fig_test_convergence2steady.py');
               disp('or ../boltmann/fig_Boltzmann_Steady.m'); end
            end
         end
         leg_text = leg_text_used;
         fortcols = fortcols_used;
      end

      fcols    = cols;
      fcols{1} = '-r';
      %fn_plot1d(gridprams.X(:,1)/1e3,out_fields.Hs(:,1),labs1d_1,cols{1});
      fn_plot1d(gridprams.X(:,1)/1e3,mean(out_fields.Hs,2),labs1d_1,fcols{1});
      leg_text{end+1}   = 'Total';
      hold on;
      %%
      if 0
         [Ep,Em,Et1,Et2]   = fn_split_energy(om_vec,wave_stuff.dirs,wave_stuff.dir_spec);
         %Hp                = 4*sqrt(Ep(:,1));
         %Hm                = 4*sqrt(Em(:,1));
         Hp                = 4*sqrt(mean(Ep,2));
         Hm                = 4*sqrt(mean(Em,2));
         if params_in.DIAG1d_OPT==0
            %Hs2               = 4*sqrt(Ep(:,1)+Em(:,1));%%add Ep + Em
            Hs2               = 4*sqrt(mean(Ep,2)+mean(Em,2));
            leg_text{end+1}   = 'Total (test 1)';
         elseif params_in.DIAG1d_OPT==1
            %Hs2               = 4*sqrt(Et1(:,1));%%check const panel integration
            Hs2               = 4*sqrt(mean(Et1,2));
            leg_text{end+1}   = 'Total (test 2)';
         elseif params_in.DIAG1d_OPT==2
            %Hs2               = 4*sqrt(Et2(:,1));%%check Simpson's rule integration
            Hs2               = 4*sqrt(mean(Et2,2));
            leg_text{end+1}   = 'Total (Simpson''s)';
         end

         fn_plot1d(gridprams.X(:,1)/1e3,Hs2,labs1d_1,['-',fcols{1}]);
         hold on;

         fn_plot1d(gridprams.X(:,1)/1e3,Hp,labs1d_1,fcols{2});
         leg_text{end+1}   = 'Fwd';
         hold on;

         fn_plot1d(gridprams.X(:,1)/1e3,Hm,labs1d_1,fcols{3});
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

   if params_in.SV_FIG%%save figures

      %if wave_stuff.nfreq==1
      %   if params_in.SCATMOD==1
      %      fig_dir  = 'out/isotropic_1freq';  %%use this for monochromatic wave
      %   elseif params_in.SCATMOD==0
      %      fig_dir  = 'out/simple_1freq';  %%use this for monochromatic wave
      %   end
      %else
      %   if params_in.SCATMOD==1
      %      fig_dir  = 'out/isotropic_spec';  %%use this for spectrum
      %   elseif params_in.SCATMOD==0
      %      fig_dir  = 'out/simple_spec';  %%use this for spectrum
      %   end
      %end

      fig_dir  = [params_in.outdir,'/figs/final'];
      eval(['!mkdir -p ',fig_dir])

      figure(3);
      eval(['!mkdir -p ',fig_dir,'/fig']);
      saveas(gcf,[fig_dir,'/fig/B',num2str(wave_stuff.ndir,'%3.3d'),'.fig']);

      eval(['!mkdir -p ',fig_dir,'/png']);
      saveas(gcf,[fig_dir,'/png/B',num2str(wave_stuff.ndir,'%3.3d'),'.png']);

      if 0
         figure(3);
         if 0
            %%fix position so that comparison is easier between computers
            pos   = [0.13   0.121428571428571   0.775   0.803571428571429];
            set(gca,'position',pos);
         end
         eval(['!mkdir -p ',fig_dir,'/att_fig']);
         eval(['!mkdir -p ',fig_dir,'/att_png']);
         saveas(gcf,[fig_dir,'/att_fig/B',num2str(ndir,'%3.3d'),'_atten.fig']);
         saveas(gcf,[fig_dir,'/att_png/B',num2str(ndir,'%3.3d'),'_atten.png']);
      end
   end
end

%%display info again
if params_in.DO_DISP; disp('##############################################################'); 
t1 = now;
disp(['Time taken (mins)      : ' ,num2str(t0_fac*(t1-t0))]);
disp(['Model time passed (h)  : ' ,num2str(nt*dt/3600.)]);
disp('##############################################################');
disp(' ');
disp(strvcat(Info)); end

%% save time-stepped Dmax;
%if DO_SAVE
%   save(filename,'out');
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn_plot_spec(X,Y,Hs,Tw,Dmax,s1)
%% plot Hs, Dmax, Tw &
%% S for 1 particular freq and dir

[nx,ny]  = size(X);
vbls     = {'Hs','Dmax','Tw','s1.Sdir'};
lab3     = {'{\itH}_{\rm s}, m','{\itD}_{\rm max}, m','{\itT}_{\rm w}, s'};
detls    = ['S: ',num2str(s1.period),'s, ',num2str(s1.dir),'^o'];
if ny==1
   lab3{4}  = detls;
else
   lab3{4}  = {'\itS, \rmm^2s',detls};
end

for j=1:4
   subplot(2,2,j);
   eval(['Z = ',vbls{j},';']);
   if ny==1
      labs  = {'\itx, \rmkm',lab3{j}};
      fn_plot1d(X/1e3,Z,labs);
   else
      labs  = {'\itx, \rmkm','\ity, \rmkm',lab3{j}};
      fn_pcolor(X(:,1)/1e3,Y(1,:)/1e3,Z,labs);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn_plot_spec_2(X,Y,Hs,tau_x,Dmax,tau_y)
%% plot Hs, Dmax, Tw &
%% S for 1 particular freq and dir

[nx,ny]  = size(X);
vbls     = {'Hs','Dmax','tau_x','tau_y'};
lab3     = {'{\itH}_{\rm s}, m','{\itD}_{\rm max}, m',...
            '{\tau}_{x}, \rmPa','{\tau}_{y}, \rmPa'};

%%fix positions so figures can be compared more easily between computers
pos{1}   = [0.130000000000000   0.583837209302326   0.334659090909091   0.341162790697674];
pos{2}   = [0.570340909090909   0.583837209302326   0.334659090909091   0.341162790697674];
pos{3}   = [0.130000000000000   0.110000000000000   0.334659090909091   0.341162790697674];
pos{4}   = [0.570340909090909   0.110000000000000   0.334659090909091   0.341162790697674];

for j=1:4
   subplot('position',pos{j});
   eval(['Z = ',vbls{j},';']);
   if ny==1
      labs  = {'\itx, \rmkm',lab3{j}};
      fn_plot1d(X/1e3,Z,labs);
   else
      labs  = {'\itx, \rmkm','\ity, \rmkm',lab3{j}};
      fn_pcolor(X(:,1)/1e3,Y(1,:)/1e3,Z,labs);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params_mex  = get_params_mex(params,duration,ice_prams,year_info)
%%           MEX_OPT: 1
%%            DODISP: 2
%%        params_vec: [28x1 double]
%%          - order of parameters same as in fortran/infiles/infile_nonIO.txt
%%          - also see read_params_vec subroutine in fortran/src/main/mod_WIM2d_run.F

i  = 0;

%% ================================================
%% set int_prams
fields   = {...
            'SCATMOD',...
            'ADV_DIM',...
            'ADV_OPT',...
            'DO_CHECK_INIT',...
            'DO_CHECK_PROG',...
            'DO_CHECK_FINAL',...
            'STEADY',...
            'BRK_OPT',...
            'DO_ATTEN',...
            };
            
Ni = length(fields);
for j=1:Ni
   i     = i+1;
   fld   = fields{j};
   %params_mex.int_prams(j) = params.(fld);
   params_vec(i)  = params.(fld);
end
%% ================================================

%% ================================================
%% set real_prams
i  = i+1;
% params_mex.real_prams   = [ice_prams.young,...
params_vec(i:i+4) = [ice_prams.young,...
                     ice_prams.drag_rp,...
                     ice_prams.visc_ws,...
                     duration,...
                     params.CFL];
i  = i+4;
%% ================================================
   

%% ===============================================
%% not in params_vec
%% - can be set in automatically
fields   = {...
            'MEX_OPT',...
            'DO_DISP',...
            };

for j=1:length(fields)
   fld               = fields{j};
   params_mex.(fld)  = params.(fld);
end
%% ===============================================


%% ===============================================
%% other integer parameters
fields   = {'FSD_OPT',...
            'REF_Hs_ICE',...
            'USE_ICE_VEL'};
for j=1:length(fields)
   i              = i+1;
   params_vec(i)  = params.(fields{j});
end
%% ===============================================


%% ===============================================
%% initial conditions
fields   = {'Hs_init',...
            'T_init',...
            'dir_init',...
            'conc_init',...
            'h_init',...
            'Dmax_init'};
for j=1:length(fields)
   i              = i+1;
   params_vec(i)  = params.(fields{j});
end
%% ===============================================


%% ===============================================
%% start time
fields   = {'model_day',...
            'model_seconds'};
for j=1:length(fields)
   i              = i+1;
   params_vec(i)  = year_info.(fields{j});
end
%% ===============================================


%% ===============================================
%% diagnostics
fields   = {'itest',...
            'jtest',...
            'dumpfreq'};
for j=1:length(fields)
   i              = i+1;
   params_vec(i)  = params.(fields{j});
end
%% ===============================================

params_mex.params_vec   = params_vec';
if params.DO_DISP
   ice_prams
   params_mex
   params_mex.params_vec
   %pause;
end
return;
