function meshfile = create_mesh(gfil,output_type)
%% CALL: meshfile = create_mesh(gridfile,output_type)
%% output_type = 'msh' or 'mat'

if ~exist('gfil','var')
   error('create_mesh(gridfile,output_type); output_type=''msh'' or ''mat''');
end
if ~exist('output_type','var')
   output_type  = 'msh';
end

ii    = strfind(gfil,'/');
if ~isempty(ii)
   gdir  = gfil(1:ii(end));
   if ii(1)~=1
      gdir  = [pwd,'/',gdir];
   end
   gfil0 = gfil(ii(end)+1:end);
else
   gdir  = pwd;
   gfil0 = gfil;
end
gfil  = [gdir,'/',gfil0];

% test reading of file (plot land-mask)
if strcmp(gfil0,'wim_grid.a')
   gp = fn_get_grid(gdir);
   disp(gp);
   x  = gp.X(:,1)/1e3;
   y  = gp.Y(1,:)/1e3;
   LM = gp.LANDMASK;

   disp('range of land mask');
   disp([min(LM(:)),max(LM(:))]);
   disp(' ');


   %%plot original grid
   figure;
   fn_fullscreen;
   fn_pcolor(x,y,LM);
   title('WIM grid land mask');
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
