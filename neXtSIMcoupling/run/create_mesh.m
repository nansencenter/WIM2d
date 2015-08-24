w2d   = getenv('WIM2D_PATH');
gdir  = [w2d,'/fortran/run/Inputs'];
gfil  = [gdir,'/wim_grid.a'];

meshfile = Myconvert_arctic_mesh(gfil);
%% file is saved to gdir

if 1
   %%save mesh to current directory
   disp(['Moving mesh file to current directory']);
   eval(['!mv ',meshfile,' .']);
else
   %% save mesh to /Data/sim/data/mesh
   %% (or to mirror dir on external HD)

   %%define location of /Data/sim:
   try_extHD      = 1;

   if exist('/Volumes/sim')
      %% johansen
      data_sim    = '/Volumes/sim'
      try_extHD   = 0;

      %% sometimes this dir shows up even if it's not loaded
      dd       = dir(data_sim);
      if length(dd)==2
         %%only ".", ".." in dir
         try_extHD   = 1;
      end
   end
      
   if try_extHD
      %% can't find johansen
      if exist('/Volumes/Tim_Ext_HD2/WORK/neXtSIM')
         %% external hard drive
         data_sim = '/Volumes/Tim_Ext_HD2/WORK/neXtSIM'
      else
         error('Can''t move mesh file to data/mesh directory - load johansen or external HD')
      end
   end

   final_dir   = [data_sim,'/data/mesh/'];
   disp(['Moving mesh file to ',final_dir]);
   eval(['!mv ',meshfile,' ',final_dir]);
end

gp = fn_get_grid(gdir);
disp(gp);
x  = gp.X(:,1)/1e3;
y  = gp.Y(1,:)/1e3;
LM = gp.LANDMASK;

%%plot original grid
figure;
fn_fullscreen;
fn_pcolor(x,y,LM);
