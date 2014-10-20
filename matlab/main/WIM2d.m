function [ice_fields,wave_fields,ice_prams,grid_prams,Dmax_all,brkcrt] =...
   WIM_2D(grid_prams,wave_prams,ice_prams)

DO_SAVE     = 0;
DO_PLOT     = 1;  %% change this to 0
                  %% if graphics aren't supported;
USE_ICE_VEL = 0   %% if 0, approx ice wlng by water wlng;  
DO_ATTEN    = 1   %% if 0, just advect waves
                  %%  without attenuation;

OPT   = 1;%%ice-water-land configuration;
if ~exist('SHARP_DIST')
   SHARP_DIST     = 1
end

%if DO_SAVE
%   filename=['wim2d_out',num2str(Tm),...
%            's_',num2str(mwd_dim),...
%            'deg_spread',num2str(SHARP_DIST),'.mat'];
%else
%   disp('not saving');
%   %disp('not saving - push any key');
%   %pause;
%end

%% set attenuation model;
%% also give progress report every 'reps' time
%%  steps;
reps  = 10;
format long

disp('Initialization')

%% Environment
% homedir = '/home/nersc/dany';
% wavedir = [homedir,'/waves'];
% mizdir  = [homedir,'/validation/miz'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%check what inputs we are given;
HAVE_GRID   = exist('grid_prams','var');
HAVE_ICE    = exist('ice_fields','var');
if HAVE_ICE
   if ~isfield(ice_fields,'cice')
      ice_prams   = ice_fields;
      HAVE_ICE    = 0;
   else
      ice_prams   = struct('c'        ,'given',...
                           'h'        ,'given',...
                           'bc_opt'   ,0,...
                           'young_opt',1);
      ice_prams   = fn_fill_iceprams(ice_prams);
   end
end
HAVE_WAVES  = exist('wave_fields','var');
if HAVE_WAVES
   %%check if full info is given or just partial info;
   if ~isfield(wave_fields,'dir_spec')
      wave_prams  = wave_fields;
      HAVE_ICE    = 0;
   end
end
HAVE3       = HAVE_GRID+HAVE_ICE+HAVE_WAVES;
if (HAVE3>0)&(HAVE3<3)
   disp(' ');
   disp('*********************************************************************');
   disp('Please specify all 3 of ''grid_prams'',''ice_fields'',''wave_fields''');
   disp('(or specify none, and let configuration be chosen through local variable ''OPT'')');
   disp(['HAVE_GRID   = ',num2str(HAVE_GRID )]);
   disp(['HAVE_ICE    = ',num2str(HAVE_ICE  )]);
   disp(['HAVE_WAVES  = ',num2str(HAVE_WAVES)]);
   disp('*********************************************************************');
   disp(' ');
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

   nx = 51;
   ny = 51;
   dx = 4000; % m
   dy = 4000; % m
   %%
   grid_prams  = struct('nx',nx,'ny',ny,...
                        'dx',dx,'dy',dy);
   grid_prams  = get_grid(grid_prams,OPT);
   %% grid_prams  = structure, eg:
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
   h           = 2;
   c           = 0.75;
   bc_opt      = 0;%%breaking condition (0=beam;1=Marchenko)
   young_opt   = 0;%%young's modulus option
   ice_prams   = struct('c'         ,c,...
                        'h'         ,h,...
                        'bc_opt'    ,bc_opt,...
                        'young_opt' ,young_opt);
   %%
   [ice_fields,ice_prams]  = iceinit(ice_prams,grid_prams,OPT);
   %% ice_fields  = structure:
   %%      cice: [51x51 double]
   %%      hice: [51x51 double]
   %%      Dmax: [51x51 double]
   %%  WTR_MASK: [51x51 logical]
   %%  ICE_MASK: [51x51 double]
   cice     = ice_fields.cice;
   hice     = ice_fields.hice;
   Dmax     = ice_fields.Dmax;
   WTR_MASK = ice_fields.WTR_MASK;
   ICE_MASK = ice_fields.ICE_MASK;

   %% ice_prams = structure eg:
   %%               c: 0.750000000000000
   %%               h: 2
   %%           young: 2.000000000000000e+09
   %%          bc_opt: 0
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

%% Model Parameters
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%waves
if HAVE_WAVES==0
   Hs0         = 2;
   Tp0         = 12;
   wave_prams  = struct('Hs',Hs0,...
                        'Tp',Tp0);
   wave_fields = waves_init(grid_prams,wave_prams,ice_fields,OPT);
   wave_stuff  = set_incident_waves(grid_prams,wave_fields);
end

nw       = wave_stuff.nfreq;     %% number of frequencies
om       = 2*pi*wave_stuff.freq; %% radial freq
ndir     = wave_stuff.ndir;      %% number of directions
wavdir   = wave_stuff.dirs;      %% wave from, degrees, clockwise
Sdir     = wave_stuff.dir_spec;  %% initial directional spectrum
%%
T  = 2*pi./om;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Water wavelength and wave speed
%% is a function only of wave period
wlng     = g.*T.^2./(2.*pi);
ap       = sqrt(g.*wlng./(2.*pi)); % Phase speed
ag       = ap./2;                  % Group speed

ag_eff      = zeros(ny,nx,nw);
wlng_ice    = zeros(ny,nx,nw);
disp_ratio  = ones (ny,nx,nw);
atten_nond  = zeros(ny,nx,nw);
damping     = zeros(ny,nx,nw);

for i = 1:nx
for j = 1:ny
   if cice(i,j)>0
      %% if ice is present:
      %% get ice wavelengths, group velocities,
      %% displacement ratios (convert wtr disp to ice),
      %% attenuation;
      [damping_rp,kice,kwtr,int_adm,NDprams,...
            alp_scat,modT,argR,argT] =...
               RT_param_outer(hice(i,j),om,young,visc_rp);

      if USE_ICE_VEL==0
         %%use wtr group vel;
         ag_eff(i,j,:)  = ag;
      else
         %%weighted avg of ice and wtr group vel;
         ag_ice         = GEN_get_ice_groupvel(hice(i,j),T,Inf,young);
         ag_eff(i,j,:)  = cice(i,j)*ag_ice+...
                           +(1-cice(i,j))*ag;
      end

      wlng_ice(i,j,:)   = 2*pi./kice;
      disp_ratio(i,j,:) = wlng.*(kice/2/pi)*modT;
      atten_nond(i,j,:) = alp_scat;
      damping(i,j,:)    = damping_rp;
   else
      wlng_ice(i,j,:)   = wlng;
   end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CFL   = .7;
dt    = CFL*dx/max(ag_eff(:)); 

if 0
   nt    = 50;
   amax  = max(ag_eff(:))
   L     = nt*dt*amax
else
   L     = max(X(:))-min(X(:));
   amax  = max(ag_eff(:));
   nt    = round( L/amax/dt );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display parameters
disp('------------------------------------')
disp(['sigma_c    = ' num2str(ice_prams.sigma_c) ' Pa']);
disp(['fragility  = ' num2str(fragility)]);
disp(['strain_c   = ' num2str(strain_c)]);
disp(['h          = ' num2str(ice_prams.h) ' m const']);
disp(['Tp         = ' num2str(wave_prams.Tp) ' s']);
disp(['Hs         = ' num2str(wave_prams.Hs) ' m']);
disp(['CFL        = ' num2str(CFL)]);
disp(['dt         = ' num2str(dt)]);
disp(['nt         = ' num2str(nt)]);
disp('------------------------------------')
disp(' ');

%% Integration
t0       = now;   %%days
t0_fac   = 24*60; %%days to minutes
brkcrt  = zeros(nx,ny,nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define weights for numerical quadrature;
if nw>1%% weights for integral over frequency
   %% (Simpson's rule);
   wt_simp            = 2+0*om;
   wt_simp([1 end])   = 1;
   wt_simp(2:2:end-1) = 4;

   %%NB om needs to be equally spaced;
   dom    = abs(om(2)-om(1));
   wt_om  = dom/3*wt_simp;

   %if 0
   %   Hs2      = 4*sqrt(wt_om*S00');
   %   testHs   = [Hs2,Hs0]
   %   pause;
   %end
else
   wt_om = 1;
end

%% weights for integral over directions
%% NB using radians for mwd;
if ndir>1
   wt_theta = ones(ndir,1)*(2*pi/ndir);
else
   wt_theta = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GET_OUT  = 1;
if GET_OUT
   Dmax_all         = zeros(ny,nx,floor(nt/reps));
   Dmax_all(:,:,1)  = Dmax;
end

if DO_PLOT
   %%
   figure(1);
   fn_plot_ice(grid_prams,ice_fields);
   %%
   figure(2);
   Tc    = 12;%check this period
   jchq  = find(abs(T-Tc)==min(abs(T-Tc)));
   jdir  = round(ndir/2);
   s1 = struct('dir',wave_stuff.dirs(jdir),...
               'period',Tc,...
               'Sdir',Sdir(:,:,jchq,jdir));
   %%
   fn_plot_spec(X,Y,wave_fields.Hs,wave_fields.Tp,Dmax,s1);
   %%
   clear s1;
   GEN_pause
end

disp('beginning main integration...');

for n = 2:nt
   disp([n nt])

   %% spectral moments;
   mom0  = zeros(ny,nx);
   mom2  = zeros(ny,nx);

   %% variances of stress and strain;
   var_stress     = zeros(ny,nx);
   var_strain     = zeros(ny,nx);

   % %% test integrals;
   % var_boundary   = cell(1,length(Jy_boundary));
   % for r=1:length(Jy_boundary)
   %    var_boundary{r}   = 0;
   % end
   % var_boundary0  = var_boundary;  
   %%
   for jw   = 1:nw


      %% CALC DIMENSIONAL ATTEN COEFF;
      atten_dim   = 0*X;
      for i = 1:nx
      for j = 1:ny
         if cice(i,j)>0
            Dave  = floe_scaling(fragility,xi,...
                     Dmin,Dmax(i,j));

            %% get expected no of floes met per unit
            %%  distance if travelling in a line;
            if Dmax(i,j) < 200
               c1d = cice(i,j)/Dave;
               %% floes per unit length;
            else
               c1d = cice(i)/Dmax(i,j);
               %% uniform lengths
            end

            %% ENERGY attenuation coeff;
            atten_dim(i,j)   = atten_nond(i,j,jw)*c1d+...
                                +2*damping(i,j,jw)*cice(i,j);
         end
      end% j
      end% i

      s1.ndir        = ndir;
      s1.wavdir      = wavdir;
      s1.Sdir        = reshape( Sdir(:,:,jw,:), nx,ny,ndir);
      s1.ag_eff      = ag_eff(:,:,jw);
      s1.atten_dim   = atten_dim;
      s1.ICE_MASK    = ice_fields.ICE_MASK;

      %% Simple attenuation scheme - doesn't conserve scattered energy
      S_out          = adv_atten_timestep_simple(grid_prams,ice_prams,s1,dt);
      Sdir(:,:,jw,:) = reshape( S_out, nx,ny,1,ndir);
      clear s1 S_out;

      %% DO BREAKING:
      for i = 1:nx
      for j = 1:ny

         %% INTEGRATE SPECTRUM OVER DIRECTION;
         Sint  = wt_theta'*squeeze(Sdir(i,j,jw,:));

         %% convert from water amp's to ice amp's;
         F     = disp_ratio(i,j,jw);%%|T| also
         k_ice = 2*pi/wlng_ice(i,j,jw);

         %% SPECTRAL MOMENTS;
         mom0(i,j)   = mom0(i,j)+wt_om(jw)*Sint*F^2;
         mom2(i,j)   = mom2(i,j)+wt_om(jw)*om(jw)^2*Sint*F^2;

         if ICE_MASK(i,j)>0
            %% VARIANCE OF STRAIN;
            strain_density    = Sint*F^2*...
                                  (k_ice^2*hice(i,j)/2)^2;
            var_strain(i,j)   = var_strain(i,j)+...
                                 + wt_om(jw)*strain_density;
         end
      end%% end spatial loop x;
      end%% end spatial loop y;

   end%% end spectral loop;

   % if 1%%plot sig wave height in each cell;
   %    for r = []%1:length(Jy_boundary)
   %       Hs
   %       [4*sqrt(var_boundary0{r});4*sqrt(var_boundary{r})]
   %    end
   %    %%
   %    Hs_all   = sqrt(mom0)*4;
   %    sig_eps  = sqrt(var_strain)*2;
   %    Tc       = 2*pi*sqrt(mom0./mom2);
   %    %%
   %    chq_Hs   = Hs_all(1:10,1:8)
   %    chq_Tc   = Tc(1:10,1:8)
   %    chq_eps  = sig_eps(1:10,1:8)
   %    chq_Dmax = Dmax(1:10,1:8)
   %    %%
   %    pcolor(X,Y,Dmax);
   %    pause;
   % end

   %%calc Hs, Tw (into wave_fields.Tp)
   wave_fields.Hs       = 4*sqrt(mom0);
   wave_fields.Tp       = 0*X;
   jnz                  = find(mom2>0);
   wave_fields.Tp(jnz)  = 2*pi*sqrt(mom0(jnz)./mom2(jnz));

   %% FINALLY DO FLOE BREAKING;
   for i=1:ny
   for j=1:nx
      if ICE_MASK(i,j)>0 & mom0(i,j)>0
         %% only try breaking if ice is present
         %%  & some waves have arrived;

         %% significant strain amp
         sig_strain  = 2*sqrt(var_strain(i,j));

         %%  probability of critical strain
         %%  being exceeded from Rayleigh distribution;
         Pstrain  = exp( -strain_c^2/(2*var_strain(i,j)) );

         %% FLOE BREAKING:
         BREAK_CRIT     = ( Pstrain>=P_crit );%%breaks if larger than this
         brkcrt(i,j,n)  = BREAK_CRIT;

         if BREAK_CRIT
            disp('breaking')
            %% use crest period to work out wavelength
            %% - half this is max poss floe length;
            if ~USE_ICE_WLNG
               wlng_crest = g.*T_crit.^2./(2.*pi);
            else
               wlng_crest =...
                  GEN_get_ice_wavelength(hice(i,j),T_crit,Inf,young); 
            end

            Dc = wlng_crest/2;
            if Dc >= Dmin & Dmax(i,j)>Dc
               Dmax(i,j)   = Dc;
            end
         end%% end breaking action;

         if 0%i==11 & j==1
            BREAK_CRIT
            Dmax(i,j)
         end

      elseif cice(i,j)<=0%% only water present
         Dmax(i,j)   = 0;
      end
      
   end%% end spatial loop j in y;
   end%% end spatial loop i in x;

   ice_fields.Dmax   = Dmax;
   jmiz              = find((Dmax>0)&(Dmax<250));
   {jmiz}

   %% progress report;
   if round(n/reps)==(n/reps)
      if GET_OUT
         Dmax_all(:,:,n/reps) = Dmax;
      end
      disp([num2str(n),...
            ' time steps done, out of ',...
               num2str(nt)]);
      %%
      t1    = now;
      disp(['time taken (mins): ',...
            num2str(t0_fac*(t1-t0))]);

      if DO_PLOT
         figure(2),clf;
         s1 = struct('dir',wave_stuff.dirs(jdir),...
                     'period',Tc,...
                     'Sdir',Sdir(:,:,jchq,jdir));
         %%
         fn_plot_spec(X,Y,wave_fields.Hs,wave_fields.Tp,Dmax,s1);
         clear s1;
         GEN_pause;
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

end%% end time loop

%% save time-stepped Dmax;
if DO_SAVE
   save(filename,'out');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn_plot_spec(X,Y,Hs,Tw,Dmax,s1)
%% plot Hs, Dmax, Tw &
%% S for 1 particular freq and dir

subplot(2,2,1);
H  = pcolor(X/1e3,Y/1e3,Hs);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = title('{\itH}_{\rm s}, m');
GEN_font(ttl);

subplot(2,2,2);
H  = pcolor(X/1e3,Y/1e3,Dmax);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = title('{\itD}_{\rm max}, m');
GEN_font(ttl);

subplot(2,2,3);
H  = pcolor(X/1e3,Y/1e3,Tw);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = title('{\itT}_{\rm w}, m');
GEN_font(ttl);

subplot(2,2,4);
H  = pcolor(X/1e3,Y/1e3,s1.Sdir);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
ttl   = [num2str(s1.period),'s, ',num2str(s1.dir),'^o'];
ttl   = title({'\itS, \rmm^2s',ttl});
GEN_font(ttl);
