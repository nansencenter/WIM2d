function s1 = check_mex(n_plot)
%% check binary file outputs of mex functions
%% plot fields at a certain time step (n_plot: integer or 'final')

if ~exist('n_plot','var')
   n_plot   = 'final';%%time step to look at
end

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

if strcmp(n_plot,'final')
   s1 = fn_plot_final(grid_prams,outdir);
   fn_fullscreen;
else
   s1 = fn_plot_prog(grid_prams,n_plot,outdir);
   fn_fullscreen;
end
