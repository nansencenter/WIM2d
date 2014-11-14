%% test_WIM2d_F.m
%% Author: Timothy Williams
%% Date:   20141110, 08:37:21 CET

function test_WIM2d_F()
%clear;

SV_FIG   = 1;

%%check initialisation
[grid_prams,ice_fields,wave_fields] = check_init();

if 1
   figure(1),clf;
   fn_fullscreen;
   fn_plot_ice(grid_prams,ice_fields);
   %%
   figure(2),clf;
   fn_fullscreen;
   fn_plot_waves(grid_prams,wave_fields,ice_fields);
end

D  = dir('../out/wim_prog*.a');
nm = D(end).name;
nt = str2num(nm(9:11))

nvec  = (2:80:nt);
if (max(nvec)<nt)
   nvec  = [nvec,nt];
end

for r = 1:length(nvec)
   n  = nvec(r)
   [n,nt]
   %%
   figure(3),clf;
   fn_fullscreen;
   out_fields  = plot_prog(grid_prams,n);
   drawnow;
   %GEN_pause;
end

if 1
   figure(4),clf;
   plot(grid_prams.X(:,1)/1e3,out_fields.Hs(:,1),'-k');
   set(gca,'yscale','log');
   GEN_proc_fig('{\itx}, km','{\itH}_s, m');
end

if SV_FIG
   SOLVER   = 1;
   nw       = 1;
   ndir     = 16;

   if nw==1
      if SOLVER==1
         fig_dir  = 'out/isotropic_1freq';  %%use this for monochromatic wave
      elseif SOLVER==0
         fig_dir  = 'out/simple_1freq';  %%use this for monochromatic wave
      end
   else
      if SOLVER==1
         fig_dir  = 'out/isotropic_spec';  %%use this for spectrum
      elseif SOLVER==0
         fig_dir  = 'out/simple_spec';  %%use this for spectrum
      end
   end

   if ~exist(fig_dir)
      mkdir(fig_dir)
   end

   if ~exist([fig_dir,'/att_png'])
      mkdir([fig_dir,'/fig']);
      mkdir([fig_dir,'/png']);
      mkdir([fig_dir,'/att_fig']);
      mkdir([fig_dir,'/att_png']);
   end

   figure(3);
   saveas(gcf,[fig_dir,'/fig/B',num2str(ndir,'%3.3d'),'.fig']);
   saveas(gcf,[fig_dir,'/png/B',num2str(ndir,'%3.3d'),'.png']);

   figure(4);
   pos   = [0.13   0.121428571428571   0.775   0.803571428571429];
   set(gca,'position',pos);
   saveas(gcf,[fig_dir,'/att_fig/B',num2str(ndir,'%3.3d'),'_atten.fig']);
   saveas(gcf,[fig_dir,'/att_png/B',num2str(ndir,'%3.3d'),'_atten.png']);
end

return;

function [grid_prams,ice_fields,wave_fields] = check_init()

afile = '../out/wim_grid.a';
bfile = '../out/wim_grid.b';

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

afile       = '../out/wim_init.a';
aid         = fopen(afile);
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
function s1 = plot_prog(grid_prams,n)

cts   = num2str(n,'%3.3d');
afile = ['../out/wim_prog',cts,'.a'];
bfile = ['../out/wim_prog',cts,'.b'];

nx = grid_prams.nx;
ny = grid_prams.ny;
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
