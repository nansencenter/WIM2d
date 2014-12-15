function [grid_prams,ice_fields,wave_fields] = fn_check_init(outdir)

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
