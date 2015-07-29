function grid_prams  = fn_get_grid(outdir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check if gridfiles exist
bdir           = [outdir];
afile          = [bdir,'/wim_grid.a'];

if ~exist(afile)
   bdir     = [outdir,'/binaries'];
   afile    = [bdir,'/wim_grid.a'];
end
if ~exist(afile)
   error([afile,' not present.'])
end

%% read file
grid_prams  = fn_read_general_binary(afile);

%% extra info
[grid_prams.nx,grid_prams.ny] = size(grid_prams.X);
grid_prams.dx                 = mean(grid_prams.scvx(:));
grid_prams.dy                 = mean(grid_prams.scuy(:));
