function fn_plot_waves(grid_prams,wave_fields);

X        = grid_prams.X;
Y        = grid_prams.Y;
jL       = find(grid_prams.LANDMASK==1);
[nx,ny]  = size(X);
%%
vbls  = {'Hs','Tp','mwd'};
lab3  = {'{\itH}_{\rms}, m','{\itT}_{\rmp}, s','<\theta>, degrees'};

for j=1:3
   subplot(3,1,j);
   eval(['Z = wave_fields.',vbls{j},';']);
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
