function fn_plot_ice(grid_prams,ice_fields);

X  = grid_prams.X;
Y  = grid_prams.Y;
%%
cice  = ice_fields.cice;
hice  = ice_fields.hice;
Dmax  = ice_fields.Dmax;
%%
jL       = find(grid_prams.LANDMASK==1);
cice(jL) = NaN;
hice(jL) = NaN;
Dmax(jL) = NaN;

subplot(1,3,1);
%GEN_plot_matrix(X/1e3,Y/1e3,cice,[-1 1]);
H  = pcolor(X/1e3,Y/1e3,cice);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
colorbar;
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
ttl   = title('Concentration');
GEN_font(ttl);
colorbar;
GEN_font(gca);
%%
subplot(1,3,3);
%GEN_plot_matrix(X/1e3,Y/1e3,Dmax,[-1 500]);
H  = pcolor(X/1e3,Y/1e3,Dmax);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
ttl   = title('\itD_{\rmmax}');
GEN_font(ttl);
colorbar;
GEN_font(gca);
%%
subplot(1,3,2);
%GEN_plot_matrix(X/1e3,Y/1e3,S(:,:,jpp,jmwd));
H  = pcolor(X/1e3,Y/1e3,hice);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
ttl   = title('\ith, \rmm');
GEN_font(ttl);
colorbar;
GEN_font(gca);
