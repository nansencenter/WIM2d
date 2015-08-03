function fn_plot_ice(grid_prams,ice_fields);

X        = grid_prams.X;
Y        = grid_prams.Y;
[nx,ny]  = size(X);
%%
cice  = ice_fields.cice;
hice  = ice_fields.hice;
Dmax  = ice_fields.Dmax;
%%
jL       = find(grid_prams.LANDMASK==1);
cice(jL) = NaN;
hice(jL) = NaN;
Dmax(jL) = NaN;

subplot(3,1,1);
if ny==1
   plot(X/1e3,cice);
   GEN_proc_fig('\itx, \rmkm','\itc');
else
   H  = pcolor(X/1e3,Y/1e3,cice);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   colorbar;
   GEN_font(gca);
   %%
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   ttl   = title('Concentration');
   GEN_font(ttl);
end
%%
subplot(3,1,3);
if ny==1
   plot(X/1e3,Dmax);
   GEN_proc_fig('\itx, \rmkm','\itD_{\rm max}, \rmm');
else
   H  = pcolor(X/1e3,Y/1e3,Dmax);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   ttl   = title('\itD_{\rmmax}, \rmm');
   GEN_font(ttl);
   colorbar;
   GEN_font(gca);
end
%%
subplot(3,1,2);
if ny==1
   plot(X/1e3,hice);
   GEN_proc_fig('\itx, \rmkm','\ith, \rmm');
else
   H  = pcolor(X/1e3,Y/1e3,hice);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   ttl   = title('\ith, \rmm');
   GEN_font(ttl);
   colorbar;
   GEN_font(gca);
end
