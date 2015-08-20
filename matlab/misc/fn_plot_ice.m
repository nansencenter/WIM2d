function fn_plot_ice(grid_prams,ice_fields);

X        = grid_prams.X;
Y        = grid_prams.Y;
jL       = find(grid_prams.LANDMASK==1);
[nx,ny]  = size(X);
%%
vbls  = {'cice','hice','Dmax'};
if ny==1
   lab3  = {'\itc','\ith, \rmm','\itD_{\rm max}, \rmm'};
else
   lab3  = {'Concentration','Thickness, m','\itD_{\rm max}, \rmm'};
end

for j=1:3
   subplot(3,1,j);
   eval(['Z = ice_fields.',vbls{j},';']);
   Z(jL) = NaN;
   if ny==1
      labs  = {'\itx, \rmkm',lab3{j}};
      fn_plot1d(X/1e3,Z);
   else
      labs  = {'\itx, \rmkm','\ity, \rmkm',lab3{j}};
      fn_pcolor(X(:,1)/1e3,Y(1,:)/1e3,Z,labs);
      daspect([1,1,1]);
   end
end
