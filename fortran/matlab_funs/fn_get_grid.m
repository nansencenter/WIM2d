function grid_prams  = fn_get_grid(outdir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check if gridfiles exist
bdir     = [outdir];
afile    = [bdir,'/wim_grid.a'];
bfile    = [bdir,'/wim_grid.b'];
afile2   = [bdir,'/wim_init.a'];

if ~exist(afile)
   bdir     = [outdir,'/binaries'];
   afile    = [bdir,'/wim_grid.a'];
   bfile    = [bdir,'/wim_grid.b'];
   afile2   = [bdir,'/wim_init.a'];
end

if ~exist(afile)
   error([afile,' not present.'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid_prams.dx  = mean(grid_prams.scvx(:));
grid_prams.dy  = mean(grid_prams.scuy(:))

return;
