%%testWIM2D_ideal.m

paths;

dmwd     = 180/8;
mwd_dim  = 180:dmwd:270-1; %mean wave direction
na       = length(mwd_dim);

Tm = 5:1:9;%%peak period
nT = length(Tm);

SHP_DST  = [0 1];%%sharp dist==1 => single direction
ns       = length(SHP_DST);

if 1

   for j=1:nT
      for r=1:na
         for s=1%:ns
            Tm0   = Tm(j);
            mwd0  = mwd_dim(r);
            sd    = SHP_DST(s);
            WIM_2D_ideal(Tm0,mwd0,sd);
         end
      end
   end

else

   DO_PLOT  = 1;

   for j=nT-1:nT
      for r=1:na

         for s=1:ns
            Tm0   = Tm(j);
            mwd0  = mwd_dim(r);
            sd    = SHP_DST(s);
            %%
            filename =...
               ['wim2d_out',num2str(Tm0),...
                's_',num2str(mwd0),...
                'deg_spread',num2str(sd),'.mat'];
            if exist(filename)
               load(filename);
               %%
               jj    = find(out(end,25,:));
               Dmax  = out(:,:,jj(end));
               sd
               %%
               ny = size(Dmax,1);
               nx = size(Dmax,2);
               %%
               %jx_test  = 1:10;
               jx_test  = [1:5,nx-4:nx];
               jy_test  = find(Dmax(:,round(nx/2)));
               edge     = jy_test(1);
               jy_test  = jy_test(1:16);
               tstDmax  = Dmax(jy_test,jx_test)
               %%
               dx    = 5;%% km
               x     = (1:nx)'*dx-dx/2;
               dy    = dx;
               y     = flipud((1:ny)')*dy-dy/2;

               disp([Tm0,mwd0,sd]);
               if DO_PLOT
                  if 0
                     subplot(1,2,s);
                     JX = nx-14:nx;
                     GEN_plot_matrix(x(JX)+dx/2,y+dy/2,Dmax(:,JX));
                     daspect([1 1 1]);
                  else
                     subplot(1,2,s);
                     [X,Y] = meshgrid(x,y);
                     %contourf(X,Y,Dmax);
                     surf(X,Y,Dmax);
                     set(gca,'CameraPosition',[mean(x),mean(y),510]);
                     daspect([1 1 1]);
                     %%
                     hold on;
                     plot(X(edge,:),Y(edge,:)+dy/2,'g');
                     hold off;
                  end
               end
            else
               disp('calc not done');
               disp([Tm0,mwd0,sd]);
            end
         end
         pause;

      end
   end
end
