if  0
   [X,Y,V] = peaks(50);
   subplot(1,2,1);
   P  = pcolor(X,Y,V);
   %%X should increase along rows,Y should increase down columns
   set(P, 'EdgeColor', 'none');
   daspect([1 1 1]);
   %%
   [Xq,Yq]  = meshgrid(-3:.01:3,-3:.01:3);
   Vq       = 0*Xq;
   Vq(:)    = interp2(X,Y,V,Xq(:),Yq(:));
   subplot(1,2,2);
   P  = pcolor(Xq,Yq,Vq);
   set(P, 'EdgeColor', 'none');
   daspect([1 1 1]);
else
   tf = 'test_outputs/wim2sim.mat';
   load(tf);
   X  = gridprams.X.'/1e3;%transpose (X should increase along rows),   convert to km
   Y  = gridprams.Y.'/1e3;%transpose (Y should increase down columns), convert to km
   V  = out_fields.Dmax.';%transpose also, to follow X,Y
   %%
   P  = pcolor(X,Y,V);
   %%X should increase along rows,Y should increase down columns
   set(P, 'EdgeColor', 'none');
   fn_fullscreen;
   daspect([1 1 1]);
   %%
   if 1
      Vq = interp2(X,Y,V,xcent,ycent);
      plot_param(Vq,[],[],'small_square','jet')
   else
      Vq = interp2(X,Y,V,xnode,ynode);
      plot_param(Vq(num_node(:,1)),[],[],'small_square','jet')
   end
end
