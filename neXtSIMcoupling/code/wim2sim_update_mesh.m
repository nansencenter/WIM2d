function simul_out   = wim2sim_update_mesh(simul_out,mesh,element)
%% Interpolate
%% taux,tauy
%% onto VERTICES of FEM grid;
%%
%% Interpolate
%% Dmax
%% onto CENTRES of FEM grid;

G2M_METHOD  = 0;%0: use interp2; 1: ISSM (InterpFromGridToMesh)
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
   data_names  = {'tau_x','tau_y'};
   Nd          = length(data_names);
   data_mesh   = zeros(Nn,Nd);
   X           = gridprams.X.'/1e3;%take transpose, change to km
   Y           = gridprams.Y.'/1e3;%take transpose, change to km
   ifxn        = 'interp2';
   %ifxn        = 'griddata';
   for j=1:Nd
      eval(['Z    = simul_out.wim.out.',data_names{j},'.'';']);%transpose data also
      eval(['tmp  = ',ifxn,'(X,Y,Z,xnode,ynode);']);
      data_mesh(:,j) = tmp;
   end
end
%%
simul_out.wim.wim2nodes.tau_x = data_mesh(:,1);
simul_out.wim.wim2nodes.tau_y = data_mesh(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
   data_names  = {'Dmax'};
   Nd          = length(data_names);
   data_mesh   = zeros(Ne,Nd);
   for j=1:Nd
      eval(['Z    = simul_out.wim.out.',data_names{j},'.'';']);%transpose data also
      eval(['tmp  = ',ifxn,'(X,Y,Z,xcent,ycent);']);
      data_mesh(:,j) = tmp;
   end
end
%%
simul_out.wim.wim2elements.Dmax = data_mesh(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xnode,ynode,xcent,ycent] = get_meshpoints(simul_out,mesh,element);
% 
% % Adding displacement
% xnode = mesh.node.x' + simul_out.UM(1:2:end)*1e-3;% km: size(Nn,1) Nn=no of nodes
% ynode = mesh.node.y' + simul_out.UM(2:2:end)*1e-3;% km: size(Nn,1) Nn=no of nodes
% 
% % Compute the position of the 3 nodes
% xy_tricorner(:,:,1) = xnode(element.num_node)*1000;%m: size(Ne,3,1) Ne=no of elements
% xy_tricorner(:,:,2) = ynode(element.num_node)*1000;%m: size(Ne,3,1) Ne=no of elements
% 
% % Compute position of the center
% xcent = mean(xy_tricorner(:,:,1),2)/1000;%km: (Ne,1)
% ycent = mean(xy_tricorner(:,:,2),2)/1000;%km: (Ne,1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
