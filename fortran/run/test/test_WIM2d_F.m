%% test_WIM2d_F.m
%% Author: Timothy Williams
%% Date:   20141110, 08:37:21 CET

function test_WIM2d_F()
%clear;

SV_FIG   = 1;
IO_OPT   = 0;

if IO_OPT==0
   %% outputs from 'dumb' code with no in/out
   outdir   = '../out';
   matdir   = 'out';%%where to put outputs from this program
else
   %% outputs from code with in/out
   %% - called from python
   outdir   = '../out_io';
   matdir   = 'out_io';%%where to put outputs from this program
end
if ~exist(matdir)
   mkdir(matdir);
end

%%check initialisation
[grid_prams,ice_fields,wave_fields] = check_init(outdir);
if 0
   nn = 1:26;
   tst_masks   = [nn',wave_fields.WAVE_MASK(nn,1)+ice_fields.ICE_MASK(nn,1)]
end

if 1%%plot and save initial conditions
   fig_dir  = [matdir,'/init_cons/'];
   if ~exist(fig_dir)
      mkdir(fig_dir);
   end
   figure(1),clf;
   fn_fullscreen;
   fn_plot_ice(grid_prams,ice_fields);
   saveas(gcf,[fig_dir,'ice.png']);
   %%
   figure(2),clf;
   fn_fullscreen;
   fn_plot_waves(grid_prams,wave_fields);
   saveas(gcf,[fig_dir,'waves.png']);
   %%
   cmd   = ['!cp ',outdir,'/binaries/wim_grid.* ',fig_dir]
   eval(cmd);
   cmd   = ['!cp ',outdir,'/binaries/wim_init.* ',fig_dir]
   eval(cmd);
end

progdir        = [outdir,'/binaries/prog/'];
D              = dir([progdir,'wim_prog*.a']);
nm             = D(1).name;
n0             = str2num(nm(9:11));
nm             = D(2).name;
nstep          = str2num(nm(9:11))-n0;
nm             = D(end).name;
nt             = str2num(nm(9:11))
%%
binary_final   = [outdir,'/binaries/wim_out'];

nvec  = (n0:nstep:nt);

%%further reduce freq of plotting
nstep_   = 50;
nf       = floor(nstep_/nstep);
nvec     = nvec(1:nf:end);

for r = 1:length(nvec)
   n  = nvec(r);
   disp([n,nt]);
   %%
   figure(3),clf;
   fn_fullscreen;
   out_fields  = plot_prog(grid_prams,n,outdir);
   drawnow;
   %GEN_pause;
end

out_fields  = plot_final(grid_prams,outdir);

if 1
   figure(4),clf;
   plot(grid_prams.X(:,1)/1e3,out_fields.Hs(:,1),'-k');
   set(gca,'yscale','log');
   GEN_proc_fig('{\itx}, km','{\itH}_s, m');
   %%
   xmin  = min(grid_prams.X(:))/1e3;%%km
   xmax  = max(grid_prams.X(:))/1e3;%%km
   ymax  = max(grid_prams.Y(:))/1e3;%%km
   axis([xmin,xmax,1e-3,1e1]);
end

%%these parameters determine where to save figure
SOLVER   = out_fields.SOLVER;
nw       = out_fields.n_wave_freq;
ndir     = out_fields.n_wavdir;
disp(out_fields);

if 1
   %% NB this definition won't
   %% work for all configurations
   dx       = grid_prams.dx;
   D_j      = out_fields.Dmax(:,1);
   MIZ_MASK = ((D_j>0)&(D_j<250));
   Wmiz     = sum(MIZ_MASK)*dx/1e3;
   %%
   disp(' ');
   disp(['MIZ width = ',num2str(Wmiz),' km']);
end

taux_max = max(out_fields.tau_x(:));
tauy_max = max(out_fields.tau_y(:));
disp(['max tau_x = ',num2str(taux_max),' Pa']);
disp(['max tau_y = ',num2str(tauy_max),' Pa']);
disp(' ');

if SV_FIG

   %%determine where to save files from parameters
   if nw==1
      if SOLVER==1
         fig_dir  = [matdir,'/isotropic_1freq'];  %%use this for monochromatic wave
      elseif SOLVER==0
         fig_dir  = [matdir,'out/simple_1freq'];  %%use this for monochromatic wave
      end
   else
      if SOLVER==1
         fig_dir  = [matdir,'out/isotropic_spec'];  %%use this for spectrum
      elseif SOLVER==0
         fig_dir  = [matdir,'out/simple_spec'];  %%use this for spectrum
      end
   end
   if ~exist(fig_dir)
      mkdir(fig_dir)
   end

   %%make subdirectories to separate file types
   Dirs  = {[fig_dir,'/binary/'],
            [fig_dir,'/log/'],
            [fig_dir,'/fig/'],
            [fig_dir,'/png/'],
            [fig_dir,'/att_fig/'],
            [fig_dir,'/att_png/']};
   for j=1:length(Dirs)
      if ~exist(Dirs{j})
         mkdir(Dirs{j});
      end
   end

   %%save binary file
   nd3   = num2str(ndir,'%3.3d');
   fn3   = [Dirs{1},'wim_out',nd3];
   cmd   = ['!cp ',binary_final,'.a ',fn3,'.a']
   eval(cmd);
   cmd   = ['!cp ',binary_final,'.b ',fn3,'.b']
   eval(cmd);

   %%save log file
   fil   = [Dirs{2},'wim2d_',nd3,'.log'];
   cmd   = ['!cp ',outdir,'/log/wim2d.log ',fil];
   eval(cmd);

   %%save main figures
   figure(3);
   saveas(gcf,[Dirs{3},'wim_final',nd3,'.fig']);
   saveas(gcf,[Dirs{4},'wim_final',nd3,'.png']);

   %%save plots of Hs
   figure(4);
   pos   = [0.13   0.121428571428571   0.775   0.803571428571429];
   set(gca,'position',pos);
   saveas(gcf,[Dirs{5},'wim_final',nd3,'_atten.fig']);
   saveas(gcf,[Dirs{6},'wim_final',nd3,'_atten.png']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [grid_prams,ice_fields,wave_fields] = check_init(outdir)

afile    = [outdir,'/binaries/wim_grid.a'];
bfile    = [outdir,'/binaries/wim_grid.b'];
afile2   = [outdir,'/binaries/wim_init.a'];

grid_prams  = struct('nx'        ,[],...
                     'ny'        ,[],...
                     'dx'        ,[],...
                     'dy'        ,[],...
                     'X'         ,[],...
                     'Y'         ,[],...
                     'scuy'      ,[],...
                     'scvx'      ,[],...
                     'scp2'      ,[],...
                     'scp2i'     ,[],...
                     'LANDMASK'  ,[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get basic info from bfile
bid   = fopen(bfile);
C     = textscan(bid,'%2.2d %s %s %s',1);
nrec  = C{1};
%%
C              = textscan(bid,'%3.3d %s %s %s %s %s %s',1);
grid_prams.nx  = C{1};
%%
C              = textscan(bid,'%3.3d %s %s %s %s %s %s',1);
grid_prams.ny  = C{1};
fclose(bid);

nx = grid_prams.nx;
ny = grid_prams.ny;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%read afile
fmt   = 'float32';
aid   = fopen(afile);
%%
grid_prams.X         = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
grid_prams.Y         = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
grid_prams.scuy      = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
grid_prams.scvx      = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
grid_prams.scp2      = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
grid_prams.scp2i     = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
grid_prams.LANDMASK  = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
%%
fclose(aid);

grid_prams.dx  = mean(grid_prams.scvx(:));
grid_prams.dy  = mean(grid_prams.scuy(:))

ice_fields  = struct('cice'      ,[],...
                     'hice'      ,[],...
                     'Dmax'      ,[],...
                     'ICE_MASK'  ,[]);

wave_fields = struct('Hs'        ,[],...
                     'Tp'        ,[],...
                     'mwd'       ,[],...
                     'WAVE_MASK' ,[]);

aid         = fopen(afile2);
%%
ice_fields.cice      = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
ice_fields.hice      = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
ice_fields.Dmax      = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
ice_fields.ICE_MASK  = 1.0*(ice_fields.cice>0)
%%
wave_fields.Hs          = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
wave_fields.Tp          = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
wave_fields.mwd         = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
wave_fields.WAVE_MASK   = 1.0*(wave_fields.Hs>0)
%%
fclose(aid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s1 = plot_prog(grid_prams,n,outdir)

cts   = num2str(n,'%3.3d');
afile = [outdir,'/binaries/prog/wim_prog',cts,'.a'];
bfile = [outdir,'/binaries/prog/wim_prog',cts,'.b'];

%% get basic info from bfile
%% - eg:
%% > 04       Number of records
%% > 150      Record length in x direction (elements)
%% > 050      Record length in y direction (elements)
%% > 01       Option number for solver
%% > 01       Number of wave frequencies
%% > 016      Number of wave directions

bid   = fopen(bfile);
C     = textscan(bid,'%2.2d %s %s %s',1);
nrec  = C{1};
%%
C  = textscan(bid,'%3.3d %s %s %s %s %s %s',1);
nx = C{1};
%%
C  = textscan(bid,'%3.3d %s %s %s %s %s %s',1);
ny = C{1};
%%
C        = textscan(bid,'%2.2d %s %s %s %s',1);
SOLVER   = C{1};
%%
C        = textscan(bid,'%2.2d %s %s %s %s',1);
nw       = C{1};
%%
C     = textscan(bid,'%3.3d %s %s %s %s',1);
ndir  = C{1};
fclose(bid);

X  = grid_prams.X;
Y  = grid_prams.Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%read afile
s1    = struct('Dmax'   ,[],...
               'Hs'     ,[],...
               'tau_x'  ,[],...
               'tau_y'  ,[]);
fmt   = 'float32';
aid   = fopen(afile);
%%
s1.Dmax  = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
s1.tau_x = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
s1.tau_y = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
s1.Hs    = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
%%
s1.SOLVER      = SOLVER;
s1.n_wave_freq = nw;
s1.n_wavdir    = ndir;
%%
fclose(aid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Hs = s1.Hs,GEN_pause
%Dm = s1.Dmax,GEN_pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%do plots

%%fix positions so figures can be compared more easily between computers
pos1  = [0.130000000000000   0.583837209302326   0.334659090909091   0.341162790697674];
pos2  = [0.570340909090909   0.583837209302326   0.334659090909091   0.341162790697674];
pos3  = [0.130000000000000   0.110000000000000   0.334659090909091   0.341162790697674];
pos4  = [0.570340909090909   0.110000000000000   0.334659090909091   0.341162790697674];

%subplot(2,2,1);
subplot('position',pos1);
H  = pcolor(X/1e3,Y/1e3,s1.Hs);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\itH}_{\rm s}, m');
GEN_font(ttl);

%subplot(2,2,2);
subplot('position',pos2);
H  = pcolor(X/1e3,Y/1e3,s1.Dmax);
%caxis([0 250]);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\itD}_{\rm max}, m');
GEN_font(ttl);

%subplot(2,2,3);
subplot('position',pos3);
H  = pcolor(X/1e3,Y/1e3,s1.tau_x);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\tau}_{x}, Pa');
GEN_font(ttl);

%subplot(2,2,4);
subplot('position',pos4);
H  = pcolor(X/1e3,Y/1e3,s1.tau_y);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\tau}_{y}, Pa');
GEN_font(ttl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s1 = plot_final(grid_prams,outdir)

afile = [outdir,'/binaries/wim_out.a'];
bfile = [outdir,'/binaries/wim_out.b'];

%% get basic info from bfile
%% - eg:
%% > 04       Number of records
%% > 150      Record length in x direction (elements)
%% > 050      Record length in y direction (elements)
%% > 01       Option number for solver
%% > 01       Number of wave frequencies
%% > 016      Number of wave directions
bid   = fopen(bfile);
C     = textscan(bid,'%2.2d %s %s %s',1);
nrec  = C{1};
%%
C  = textscan(bid,'%3.3d %s %s %s %s %s %s',1);
nx = C{1};
%%
C  = textscan(bid,'%3.3d %s %s %s %s %s %s',1);
ny = C{1};
%%
C        = textscan(bid,'%2.2d %s %s %s %s',1);
SOLVER   = C{1};
%%
C  = textscan(bid,'%2.2d %s %s %s %s',1);
nw = C{1};
%%
C     = textscan(bid,'%3.3d %s %s %s %s',1);
ndir  = C{1};
fclose(bid);

X  = grid_prams.X;
Y  = grid_prams.Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%read afile
s1    = struct('Dmax'   ,[],...
               'tau_x'  ,[],...
               'tau_y'  ,[],...
               'Hs'     ,[]);
fmt   = 'float32';
aid   = fopen(afile);
%%
s1.Dmax  = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
s1.tau_x = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
s1.tau_y = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
s1.Hs    = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
%%
s1.SOLVER      = SOLVER;
s1.n_wave_freq = nw;
s1.n_wavdir    = ndir;
%%
fclose(aid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Hs = s1.Hs,GEN_pause
%Dm = s1.Dmax,GEN_pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%do plots

%%fix positions so figures can be compared more easily between computers
pos1  = [0.130000000000000   0.583837209302326   0.334659090909091   0.341162790697674];
pos2  = [0.570340909090909   0.583837209302326   0.334659090909091   0.341162790697674];
pos3  = [0.130000000000000   0.110000000000000   0.334659090909091   0.341162790697674];
pos4  = [0.570340909090909   0.110000000000000   0.334659090909091   0.341162790697674];

%subplot(2,2,1);
subplot('position',pos1);
H  = pcolor(X/1e3,Y/1e3,s1.Hs);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\itH}_{\rm s}, m');
GEN_font(ttl);

%subplot(2,2,2);
subplot('position',pos2);
H  = pcolor(X/1e3,Y/1e3,s1.Dmax);
%caxis([0 250]);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\itD}_{\rm max}, m');
GEN_font(ttl);

%subplot(2,2,3);
subplot('position',pos3);
H  = pcolor(X/1e3,Y/1e3,s1.tau_x);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\tau}_{x}, Pa');
GEN_font(ttl);

%subplot(2,2,4);
subplot('position',pos4);
H  = pcolor(X/1e3,Y/1e3,s1.tau_y);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\tau}_{y}, Pa');
GEN_font(ttl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
