function simul_out   = wim2sim(simul_out,mesh,element)
%% function to call from perform_simul.m (called at end of EB_model_v2.m)
%% do it round "Time interpolation of the forcings"
%%    - when Vair and Voce are calculated
%% - wave stresses should go into system_assemble_mex.c at "Step1"
%%    - similar to Voce, Vair?
%%               coef_Vair=Vair_factor?
%% - Dmax should maybe go into Step1 too (influence damage?)
%%    - also thermodynamic effect (lat melt - thermo_ow_mex.c)

if ~exist('simul_out','var')
   testdir     = 'test_inputs';
   %testfile    = [testdir,'/simul_out_squaresmall1km_test2_step0.mat'];
   %testfile    = [testdir,'/simul_out_squaresmall1km_test2_step10.mat'];
   testfile    = [testdir,'/simul_out_squaresmall1km_test15_step0.mat'];
   tf          = load(testfile);
   simul_out   = tf.simul_out;
   clear tf;
end

if ~exist('element','var')
   %[mesh, element, simul_in.ind_node_fix_bnd, simul_in.ind_node_free_bnd, simul_in.ind_element_free_bnd] =... 
   [mesh,element] =...
      importbamg(simul_out.bamg.mesh, simul_out.bamg.geom);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% determine whether to run WIM:
if isnan(simul_out.wim.last_call)
   %% 1st time step
   %% - run and update last_call
   simul_out.wim.last_call = simul_out.current_time;
else
   time_passed = simul_out.current_time-simul_out.wim.last_call;
   if time_passed>=simul_out.wim.coupling_freq
      %% time since last call is >= coupling frequency
      %% - run and update last_call
      simul_out.wim.last_call = simul_out.current_time;
   else
      %% time since last call is < coupling frequency
      %% - go back to perform_simul
      disp('not running WIM yet');
      return;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate needed stuff from mesh to WIM grid.
[xvert,yvert,xcent,ycent] = get_centres(simul_out,mesh,element);
% [x,y]vert: vertices of FEM mesh [km]
% [x,y]cent: centres of FEM mesh  [km]
Nn    = length(xvert);%%number of nodes
Ne    = length(xcent);%%number of elements
index = element.num_node(:,[1 3 2]);%% row number is element number,
                                    %% columns are indices (in eg xvert,yvert) of its 3 nodes

%%data on FEM mesh (at centres)
data        = zeros(length(xcent),3);%%1 col for each field
data(:,1)   = simul_out.c;           %% conc
data(:,2)   = simul_out.h;           %% thickness
if isfield(simul_out,'Dmax')
   INIT_DMAX   = 0;
   data(:,3)   = simul_out.Dmax;
else
   INIT_DMAX   = 1;
   data(:,3)   = [];
end

gridprams   = simul_out.wim.gridprams;
if 1
   %%interp with ISSM
   [griddata,xWIM,yWIM] = interp_SIM2WIM_ISSM(gridprams,index,xvert,yvert,data);
else
   %%just set it for now
   cice     = 0*gridprams.X;
   hice     = 0*gridprams.X;
   JJ       = find((gridprams.Y>-gridprams.X)&(gridprams.LANDMASK==0));
   cice(JJ) = .7;
   hice(JJ) = 1;
   %%
   griddata(:,:,1)   = cice;
   griddata(:,:,2)   = hice;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0%%test interp from vertices
   data     = simul_out.UM(1:2:end);
   size(data)
   [griddata,xWIM,yWIM] = interp_SIM2WIM_ISSM(gridprams,index,mesh.node.x',mesh.node.y',data);
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set conc, thickness, Dmax on the WIM grid.
ICEMASK  = (griddata(:,:,1)>simul_out.wim.other_prams.cice_min);
if INIT_DMAX
   griddata(:,:,3)   = simul_out.wim.other_prams.Dmax_pack*ICEMASK;
end
%%
ice_fields  = struct('cice'      ,griddata(:,:,1),...
                     'hice'      ,griddata(:,:,2),...
                     'Dmax'      ,griddata(:,:,3),...
                     'ICE_MASK'  ,ICEMASK);
if 1
   figure(101);
   P  = pcolor(gridprams.X.',gridprams.Y.',ice_fields.cice.');
   colorbar;
   title('conc');
   set(P, 'EdgeColor', 'none');
   fn_fullscreen;
   daspect([1 1 1]);
   GEN_proc_fig('x, km', 'y, km');
   drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define wave forcing on WIM grid
%% - could also interpolate from wave model to WIM grid
wave_fields.WAVE_MASK      = 0*gridprams.X;
JJ                         = find((gridprams.Y<(-10.0e3-gridprams.X))&(gridprams.LANDMASK==0));
wave_fields.WAVE_MASK(JJ)  = 1.0;
wave_fields.Hs             =   3*wave_fields.WAVE_MASK;
wave_fields.Tp             =  12*wave_fields.WAVE_MASK;
wave_fields.mwd            = 225*wave_fields.WAVE_MASK;%%clockwise from north
if 1
   figure(102);
   P  = pcolor(gridprams.X.',gridprams.Y.',wave_fields.Hs.');
   colorbar;
   title('H_s');
   set(P, 'EdgeColor', 'none');
   fn_fullscreen;
   daspect([1 1 1]);
   GEN_proc_fig('x, km', 'y, km');
   drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call WIM2d
simul_out.wim.int_prams,simul_out.wim.real_prams
out_fields  = run_WIM2d_io_mex(ice_fields,wave_fields,...
                simul_out.wim.int_prams,simul_out.wim.real_prams)
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate
%% taux,tauy
%% onto VERTICES of FEM grid
data        = zeros(length(xWIM),2);
data(:,1)   = out_fields.taux(:);
data(:,2)   = out_fields.tauy(:);
%%
missing_g2m = 0;%missing value - just set to 0 (one of FEM mesh points is out of WIM grid)
data_mesh   = InterpFromGridToMesh(xWIM,yWIM,data,xvert,yvert,missing_g2m);
%%
simul_out.taux_waves = data_mesh(:,1);
simul_out.tauy_waves = data_mesh(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate
%% Dmax
%% onto CENTRES of FEM grid
data(:,1)   = out_fields.Dmax(:);
data(:,2)   = [];
%%
missing_g2m = 0;%missing value - just set to 0 (one of FEM mesh points is out of WIM grid)
data_mesh   = InterpFromGridToMesh(xWIM,yWIM,data,xvert,yvert,missing_g2m);
%%
simul_out.Dmax = data_mesh(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [node_x,node_y,xc,yc] = get_centres(simul_out,mesh,element);

% Adding displacement
node_x = mesh.node.x' + simul_out.UM(1:2:end)*1e-3;% km: size(Nn,1) Nn=no of nodes
node_y = mesh.node.y' + simul_out.UM(2:2:end)*1e-3;% km: (Nn,1) Nn=no of nodes

% Compute the position of the 3 nodes
xy_tricorner(:,:,1) = node_x(element.num_node)*1000;%m: (Ne,3,1) Ne=no of elements
xy_tricorner(:,:,2) = node_y(element.num_node)*1000;%m: (Ne,3,1) Ne=no of elements

% Compute position of the center
xc = mean(xy_tricorner(:,:,1),2)/1000;%km: (Ne,1)
yc = mean(xy_tricorner(:,:,2),2)/1000;%km: (Ne,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out,x_m,y_m] = interp_SIM2WIM_ISSM(gridprams,index,xvert,yvert,data)
%% (xvert,yvert): coords of vertices
%% data (in columns):
%%  - can be defined at centres or nodes 
%%  - ISSM fxn InterpFromMeshToGrid can tell where it is defined by the length of the data
%% index: row corresponds to element number; columns are indices (in xvert) of 3 nodes

DO_TEST  = 0;

Ne = size(index,1);
Nn = length(xvert);
Nd = size(data,1);

%WIM grid info
xmin     = 1e-3*min(gridprams.X(:)); % WIM grid xmin (stere proj, km)
ymax     = 1e-3*max(gridprams.Y(:)); % WIM grid ymax (stere proj, km)
xposting = 1e-3*gridprams.dy;        % res in y dirn (km)
yposting = 1e-3*gridprams.dx;        % res in x dirn (km)
nlines   = gridprams.ny;             % no of cols in WIM grid
ncols    = gridprams.nx;             % no of rows in WIM grid


%% inputs to interp fxn
inputs.x             = xvert;
inputs.y             = yvert;
inputs.data          = data(:,1);
inputs.index         = index;
inputs.xmin          = xmin;
inputs.ymax          = ymax;
inputs.xposting      = xposting;
inputs.yposting      = yposting;
inputs.nlines        = double(nlines);%%inputs need to be doubles
inputs.ncols         = double(ncols); %%inputs need to be doubles
inputs.default_value = NaN;

if DO_TEST
   if Nd==Ne
      disp('Interpolating from centres');
      ttl   = 'conc';
      mfil  = 'test_outputs/test_InterpFromMeshToGrid_conc.mat';
      ffil  = 'test_outputs/test_InterpFromMeshToGrid_conc.png';
      figure(101);
   else
      disp('Interpolating from nodes');
      ttl   = 'node disp (x)';
      mfil  = 'test_outputs/test_InterpFromMeshToGrid_UMx.mat';
      ffil  = 'test_outputs/test_InterpFromMeshToGrid_UMx.png';
      figure(102);
   end
   %%
   disp(inputs);
   save(mfil,'inputs');
end

Nfields  = size(data,2);
out      = zeros(gridprams.nx,gridprams.ny,Nfields);
for j=1:Nfields
   inputs.data          = data(:,j);
   [x_m,y_m,griddata]   =...
      InterpFromMeshToGrid(inputs.index,inputs.x,inputs.y,inputs.data,...
                           inputs.xmin,inputs.ymax,inputs.xposting,inputs.yposting,inputs.nlines,inputs.ncols,...
                           inputs.default_value);

   %x_m is length ncols
   %y_m is length nlines
   %griddata is size [nlines,ncols]
   out(:,:,j)  = griddata.';
end

if DO_TEST
   %x_m is length ncols
   %y_m is length nlines
   %griddata is size [nlines,ncols]
   size(griddata),size(x_m),size(y_m)
   
   %% plot conc
   [Xm,Ym]  = meshgrid(x_m,y_m);
   P        = pcolor(Xm,Ym,griddata);
   colorbar;
   title(ttl);
   set(P, 'EdgeColor', 'none');
   GEN_proc_fig('x, km', 'y, km')
   saveas(gcf,ffil);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
