%% test_WIM2d_F.m
%% Author: Timothy Williams
%% Date:   20141110, 08:37:21 CET

function test_WIM2d_F()
%clear;

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

nvec  = (2:40:nt);
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

saveas(gcf,'out/wim2d_F.png');

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

afile       = 'out/wim_init.a';
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
subplot(2,2,1);
H  = pcolor(X/1e3,Y/1e3,s1.Hs);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\itH}_{\rm s}, m');
GEN_font(ttl);

subplot(2,2,2);
H  = pcolor(X/1e3,Y/1e3,s1.Dmax);
%caxis([0 250]);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\itD}_{\rm max}, m');
GEN_font(ttl);

subplot(2,2,3);
H  = pcolor(X/1e3,Y/1e3,s1.tau_x);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\tau}_{x}, Pa');
GEN_font(ttl);

subplot(2,2,4);
H  = pcolor(X/1e3,Y/1e3,s1.tau_y);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\tau}_{y}, Pa');
GEN_font(ttl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
