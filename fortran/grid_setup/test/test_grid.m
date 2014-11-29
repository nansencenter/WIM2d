%% test_grid.m
%% Author: Timothy Williams
%% Date:   20141128, 18:18:43 CET

function grid_prams  = test_grid()
%clear;

SV_FIG   = 1;

outdir   = '../../run/inputs';

%%check grid
grid_prams  = check_grid(outdir);
plot_grid(grid_prams);
fn_fullscreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grid_prams = check_grid(outdir)

afile    = [outdir,'/wim_grid.a']
bfile    = [outdir,'/wim_grid.b']

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
C  = textscan(bid,'%3.3d %s %s %s %s %s %s',1);
nx = C{1};
C  = textscan(bid,'%3.3d %s %s %s %s %s %s',1);
ny = C{1};
fclose(bid);
%%
grid_prams.nx  = nx;
grid_prams.ny  = ny;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_grid(grid_prams)

X           = grid_prams.X;
Y           = grid_prams.Y;
LANDMASK    = grid_prams.LANDMASK;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%do plots

H  = pcolor(X/1e3,Y/1e3,LANDMASK);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('Land mask');
GEN_font(ttl);
