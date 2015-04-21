function simul_out   = wim2sim_update_mesh(simul_out,mesh,element,INTERP_ICE_ELEMENTS)
%% Interpolate
%% taux,tauy
%% onto VERTICES of FEM grid;
%%
%% Interpolate
%% Dmax
%% onto CENTRES of FEM grid;

G2M_METHOD  = 0;%0: use interp2; 1: ISSM (InterpFromGridToMesh)
if ~exist('INTERP_ICE_ELEMENTS')
   INTERP_ICE_ELEMENTS  = 1;
end
%%
missing_g2m = NaN;%missing value - just set to 0 (one of FEM mesh points is out of WIM grid)
gridprams   = simul_out.wim.gridprams;
Ng          = length(gridprams.X(:));
Nn          = mesh.Nn;
Ne          = mesh.Ne;

%%get mesh nodes/centres
[xnode,ynode,xcent,ycent]  = wim2sim_get_meshpoints(simul_out,mesh,element);

if G2M_METHOD==1
   %Use InterpFromGridToMesh:
   %>/Users/timill/GITHUB-REPOSITORIES/neXtSIM/ISSM-trunk-jpl-svn/src/wrappers/InterpFromGridToMesh/InterpFromGridToMesh.cpp
   %>/Users/timill/GITHUB-REPOSITORIES/neXtSIM/ISSM-trunk-jpl-svn/src/c/modules/InterpFromGridToMeshx/InterpFromGridToMeshx.cpp
   xyWIM       = zeros(Ng,2);
   xyWIM(:,1)  = gridprams.X(:);
   xyWIM(:,2)  = gridprams.Y(:);
   %%
   data        = zeros(Ng,2);
   data(:,1)   = simul_out.wim.out.tau_x(:);
   data(:,2)   = simul_out.wim.out.tau_y(:);
   %%
   save('test_outputs/for_InterpFromGridToMesh','xyWIM','data','xnode','ynode','missing_g2m');
   data_mesh   = InterpFromGridToMesh(xyWIM(:,1),xyWIM(:,2),data,xnode,ynode,missing_g2m);
else
   %%use interp2
   %data_names  = {'tau_x','tau_y'};
   sname_in    = 'simul_out.wim.waves_for_nodes';  %% data to be interpolated is here
   sname_out   = 'simul_out.wim.waves_on_nodes';   %% interpolated data goes here
   %%
   eval(['data_names  = fieldnames(',sname_in,');']);
   Nd          = length(data_names);
   %data_mesh   = zeros(Nn,Nd);
   X           = gridprams.X.'/1e3;%take transpose, change to km
   Y           = gridprams.Y.'/1e3;%take transpose, change to km
   ifxn        = 'interp2';
   %ifxn        = 'griddata';
   for j=1:Nd
      dname = data_names{j};
      cmd1  = ['Z    = ',sname_in,'.',dname,'.'';'];%transpose data also
      cmd2  = [sname_out,'.',dname,' = ',ifxn,'(X,Y,Z,xnode,ynode);'];
      eval(cmd1);
      eval(cmd2);
      %eval(['tmp  = ',ifxn,'(X,Y,Z,xnode,ynode);']);
      %data_mesh(:,j) = tmp;
   end
end
%%
%simul_out.wim.waves_on_nodes.tau_x = data_mesh(:,1);
%simul_out.wim.waves_on_nodes.tau_y = data_mesh(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if INTERP_ICE_ELEMENTS==0
   %% don't interpolate Dmax etc (simul_out.wim.ice_on_elements)
   %% - these are advected with ice,
   %%   and may also be affected by thermodynamics,
   %%   so interpolating from wave grid at the wrong time could delete these changes
   return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate
%% Dmax
%% onto CENTRES of FEM grid
if G2M_METHOD==1
   %Use InterpFromGridToMesh:
   %>/Users/timill/GITHUB-REPOSITORIES/neXtSIM/ISSM-trunk-jpl-svn/src/wrappers/InterpFromGridToMesh/InterpFromGridToMesh.cpp
   %>/Users/timill/GITHUB-REPOSITORIES/neXtSIM/ISSM-trunk-jpl-svn/src/c/modules/InterpFromGridToMeshx/InterpFromGridToMeshx.cpp
   data(:,1)   = simul_out.wim.out.Dmax(:);
   data(:,2)   = [];
   data_mesh   = InterpFromGridToMesh(xyWIM(:,1),xyWIM(:,2),data,xcent,ycent,missing_g2m);
else
   %%use interp2
   %data_names  = {'Dmax'};
   sname_in    = 'simul_out.wim.ice_for_elements';  %% data to be interpolated is here
   sname_out   = 'simul_out.wim.ice_on_elements';   %% interpolated data goes here
   %%
   eval(['data_names  = fieldnames(',sname_in,');']);
   Nd          = length(data_names);
   data_mesh   = zeros(Ne,Nd);
   for j=1:Nd
      dname = data_names{j};
      cmd1  = ['Z    = ',sname_in,'.',dname,'.'';'];%transpose data also
      cmd2  = [sname_out,'.',dname,' = ',ifxn,'(X,Y,Z,xcent,ycent);'];
      eval(cmd1);
      eval(cmd2);
      %eval(['Z    = simul_out.wim.out.',data_names{j},'.'';']);%transpose data also
      %eval(['tmp  = ',ifxn,'(X,Y,Z,xcent,ycent);']);
      %data_mesh(:,j) = tmp;
   end
end
%%
%simul_out.wim.wim2elements.Dmax = data_mesh(:,1);
