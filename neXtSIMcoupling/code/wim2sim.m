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

RUN_WIM              = 0;
TEST_INTERP2GRID     = 1;
TEST_INTERP2NODES    = 1;
TEST_INTERP2CENTRES  = 1;

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
gridprams                  = simul_out.wim.gridprams;
[xnode,ynode,xcent,ycent]  = wim2sim_get_meshpoints(simul_out,mesh,element);
% [x,y]vert: vertices of FEM mesh [km]
% [x,y]cent: centres of FEM mesh  [km]

Nn    = mesh.Nn;                       %% number of nodes
Ne    = mesh.Ne;                       %% number of elements
Ng    = gridprams.nx*gridprams.ny;     %% total number of grid points (WIM grid)
index = element.num_node(:,[1 3 2]);   %% row number is element number,
                                       %% columns are indices (in eg xnode,ynode) of its 3 nodes

%%data on FEM mesh (at centres)
data        = zeros(length(xcent),3);  %%1 col for each field
data(:,1)   = simul_out.c;             %% conc
data(:,2)   = simul_out.h;             %% thickness
if simul_out.wim.INIT_DMAX==0
   data(:,3)   = simul_out.wim.wim2centres.Dmax;
else
   data(:,3)   = [];
end


if 1
   %%interp with ISSM
   grid_data  = interp_SIM2WIM_ISSM(gridprams,index,xnode,ynode,data);
else
   %%just set it for now
   cice     = 0*gridprams.X;
   hice     = 0*gridprams.X;
   JJ       = find((gridprams.Y>-gridprams.X)&(gridprams.LANDMASK==0));
   cice(JJ) = .7;
   hice(JJ) = 1;
   %%
   grid_data(:,:,1)   = cice;
   grid_data(:,:,2)   = hice;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set conc, thickness, Dmax on the WIM grid.
ICEMASK  = (grid_data(:,:,1)>simul_out.wim.other_prams.cice_min);
if simul_out.wim.INIT_DMAX==1
   grid_data(:,:,3)        = simul_out.wim.other_prams.Dmax_pack*ICEMASK;
   simul_out.wim.INIT_DMAX = 0;%%Dmax initialised now
end
%%
ice_fields  = struct('cice'      ,grid_data(:,:,1),...
                     'hice'      ,grid_data(:,:,2),...
                     'Dmax'      ,grid_data(:,:,3),...
                     'ICE_MASK'  ,ICEMASK);

if TEST_INTERP2GRID==1
   figure(91);
   P  = pcolor(gridprams.X.',gridprams.Y.',ice_fields.cice.');
   colorbar;
   title('conc');
   set(P, 'EdgeColor', 'none');
   fn_fullscreen;
   daspect([1 1 1]);
   GEN_proc_fig('x, km', 'y, km');
   drawnow;
   %%
   plot_param(simul_out.c,[],[],'small_square','jet');
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

if TEST_INTERP2GRID==1
   figure(92);
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

tfil_wim2sim   = 'test_outputs/wim2sim.mat';
if RUN_WIM==1
   %% Call WIM2d
   simul_out.wim.int_prams,simul_out.wim.real_prams
   out_fields  = run_WIM2d_io_mex(ice_fields,wave_fields,...
                   simul_out.wim.int_prams,simul_out.wim.real_prams)
   if 1
      save(tfil_wim2sim,'out_fields');
   else
      num_node = element.num_node;
      save(tfil_wim2sim,'out_fields','gridprams','xcent','ycent','xnode','ynode','num_node');
   end
else
   load(tfil_wim2sim);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simul_out.wim.out.tau_x = out_fields.tau_x;
simul_out.wim.out.tau_y = out_fields.tau_y;
simul_out.wim.out.Dmax  = out_fields.Dmax;
%%
simul_out   = wim2sim_update_mesh(simul_out,mesh,element);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if TEST_INTERP2NODES==1
   %%plot original value
   fld_names   = {'tau_x','tau_y'};
   for j=1:2
      figure(100+j);
      fld_name = fld_names{j};
      cmd      = ['P  = pcolor(gridprams.X.''/1e3,gridprams.Y.''/1e3,out_fields.',fld_name,'.'');'];
      eval(cmd);
      colorbar;
      title(fld_name);
      set(P, 'EdgeColor', 'none');
      fn_fullscreen;
      daspect([1 1 1]);
      GEN_proc_fig('x, km', 'y, km');
      drawnow;

      %%plot interpolated value (at NODES)
      cmd      = ['plot_param(simul_out.wim.wim2nodes.',fld_name,'(element.num_node(:,1)),[],[],''small_square'',''jet'')'];
      eval(cmd);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if TEST_INTERP2CENTRES==1
   %%plot original value
   figure(103);
   fld_name = 'Dmax';
   cmd      = ['P  = pcolor(gridprams.X.''/1e3,gridprams.Y.''/1e3,out_fields.',fld_name,'.'');']
   eval(cmd);
   colorbar;
   title(fld_name);
   set(P, 'EdgeColor', 'none');
   fn_fullscreen;
   daspect([1 1 1]);
   GEN_proc_fig('x, km', 'y, km');
   drawnow;

   %%plot interpolated value (at CENTRES)
   cmd      = ['plot_param(simul_out.wim.wim2elements.',fld_name,',[],[],''small_square'',''jet'');']
   eval(cmd);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out,x_m,y_m] = interp_SIM2WIM_ISSM(gridprams,index,xnode,ynode,data)
%% (xnode,ynode): coords of vertices
%% data (in columns):
%%  - can be defined at centres or nodes 
%%  - ISSM fxn InterpFromMeshToGrid can tell where it is defined by the length of the data
%% index: row corresponds to element number; columns are indices (in xnode) of 3 nodes

DO_TEST  = 0;

Ne = size(index,1);
Nn = length(xnode);
Nd = size(data,1);

%WIM grid info
xmin     = 1e-3*min(gridprams.X(:)); % WIM grid xmin (stere proj, km)
ymax     = 1e-3*max(gridprams.Y(:)); % WIM grid ymax (stere proj, km)
xposting = 1e-3*gridprams.dy;        % res in y dirn (km)
yposting = 1e-3*gridprams.dx;        % res in x dirn (km)
nlines   = gridprams.ny;             % no of cols in WIM grid
ncols    = gridprams.nx;             % no of rows in WIM grid


%% inputs to interp fxn
inputs.x             = xnode;
inputs.y             = ynode;
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
   [x_m,y_m,grid_data]   =...
      InterpFromMeshToGrid(inputs.index,inputs.x,inputs.y,inputs.data,...
                           inputs.xmin,inputs.ymax,inputs.xposting,inputs.yposting,inputs.nlines,inputs.ncols,...
                           inputs.default_value);

   %x_m is length ncols
   %y_m is length nlines
   %grid_data is size [nlines,ncols]
   out(:,:,j)  = grid_data.';
end

if DO_TEST
   %x_m is length ncols
   %y_m is length nlines
   %grid_data is size [nlines,ncols]
   size(grid_data),size(x_m),size(y_m)
   
   %% plot conc
   [Xm,Ym]  = meshgrid(x_m,y_m);
   P        = pcolor(Xm,Ym,grid_data);
   colorbar;
   title(ttl);
   set(P, 'EdgeColor', 'none');
   GEN_proc_fig('x, km', 'y, km')
   saveas(gcf,ffil);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
