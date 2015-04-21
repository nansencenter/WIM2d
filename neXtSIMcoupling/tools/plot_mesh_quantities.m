function h=plot_mesh_quantities(param,mesh,element,domain,colormap_name)
%% adapted from plot_param

load arctic_coasts_light.mat;

%global mesh element
node_x     = mesh.node.x;
node_y     = mesh.node.y;
num_node   = element.num_node;

defo_tmp.data.xy_tricorner(:,:,1)=node_x(num_node)*1000;
defo_tmp.data.xy_tricorner(:,:,2)=node_y(num_node)*1000;

defo_tmp.data.variable=param;
caxis_range=[min(param),max(param)];
title_figure='';
masked=0;
visible=1;

%Selecting the mesh boundaries
boundary = mesh.boundary.from_msh;
%Selecting only closed boundaries
b1 = find(1==boundary(:,3));
%Selecting only closed boundaries
b2 = find(0==boundary(:,3));
 
defo_tmp.boundaryX  = node_x(boundary(b1,1:2,1));
defo_tmp.boundaryY  = node_y(boundary(b1,1:2,1));
defo_tmp.boundaryXo = node_x(boundary(b2,1:2,1));
defo_tmp.boundaryYo = node_y(boundary(b2,1:2,1));

plot_tricorner(defo_tmp,'variable',colormap_name,caxis_range,title_figure,domain,[],masked,visible,'','');
