function meshfile = create_mesh(OPT,output_type)

if nargin~=2
   error('create_mesh(OPT,output_type) eg OPT=2, output_type=''msh'' or ''mat''')
end

w2d   = getenv('WIM2D_PATH');
wgd   = getenv('WIMGRIDPATH');
if OPT==1
   gdir  = [w2d,'/fortran/run/Inputs'];
   gfil  = [gdir,'/wim_grid.a']
elseif OPT==2
   gdir  = wgd;
   gfil  = [gdir,'/wim_grid_full_ONR_Oct2015_2km_big.a']
elseif OPT==3
   gdir  = wgd;
   gfil  = [gdir,'/wim_grid_full_ONR_Oct2015_2km_small.a']
elseif OPT==4
   gdir  = wgd;
   gfil  = [gdir,'/wim_grid_full_ONR_Oct2015_4km_big.a']
elseif OPT==5
   gdir  = wgd;
   gfil  = [gdir,'/wim_grid_full_FS_Dec2015_4km_big.a']
end

% test reading of file (plot land-mask)
if OPT==1
   gp = fn_get_grid(gdir);
   disp(gp);
   x  = gp.X(:,1)/1e3;
   y  = gp.Y(1,:)/1e3;
   LM = gp.LANDMASK;

   %%plot original grid
   figure;
   fn_fullscreen;
   fn_pcolor(x,y,LM);
else
   [fields,info]  = fn_read_nextwim_binary(gfil)
   [X,Y]          = mapll(fields.qlat,fields.qlon,60,-45,'N');%km
   x              = X(:,1);
   y              = Y(1,:);
 
   %%plot original grid
   figure;
   fn_fullscreen;
   fn_pcolor(x,y,fields.LANDMASK);
   title('Land Mask');
end

meshfile = Myconvert_arctic_mesh(gfil,output_type);
%% mesh file is saved to gdir

if 0
   %%save mesh to current directory
   disp(['Moving mesh file to current directory']);
   eval(['!mv ',meshfile,' .']);
elseif 0
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

