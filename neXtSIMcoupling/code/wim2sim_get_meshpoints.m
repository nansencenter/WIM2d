function [xnode,ynode,xcent,ycent] = wim2sim_get_meshpoints(simul_out,mesh,element,GET_EXACT);
%%function to get nodes and centres of neXtSIM mesh

if ~exist('GET_EXACT','var')
   GET_EXACT   = 1;%%default is to get exact position of nodes and centres
end

if GET_EXACT==0
   %% don't worry about taking account of the moving mesh

   xnode = mesh.node.x';% km: size(Nn,1) Nn=no of nodes
   ynode = mesh.node.y';% km: size(Nn,1) Nn=no of nodes
   xcent = element.x;
   ycent = element.y;
else
   %% do worry about taking account of the moving mesh

   % Adding displacement
   xnode = mesh.node.x' + simul_out.UM(1:2:end)*1e-3;% km: size(Nn,1) Nn=no of nodes
   ynode = mesh.node.y' + simul_out.UM(2:2:end)*1e-3;% km: size(Nn,1) Nn=no of nodes

   % Compute the position of the 3 nodes
   xy_tricorner(:,:,1) = xnode(element.num_node)*1000;%m: size(Ne,3,1) Ne=no of elements
   xy_tricorner(:,:,2) = ynode(element.num_node)*1000;%m: size(Ne,3,1) Ne=no of elements

   % Compute position of the center
   xcent = mean(xy_tricorner(:,:,1),2)/1000;%km: (Ne,1)
   ycent = mean(xy_tricorner(:,:,2),2)/1000;%km: (Ne,1)
end
