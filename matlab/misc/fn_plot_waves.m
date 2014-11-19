%% waves_init.m
%% Author: Timothy Williams
%% Date: 20141016, 18:13:51 CEST

function fn_plot_waves(grid_prams,wave_fields);

X        = grid_prams.X;
Y        = grid_prams.Y;
Hs_mat   = wave_fields.Hs;
Tp_mat   = wave_fields.Tp;
mwd_mat  = wave_fields.mwd;
%%
jL          = find(grid_prams.LANDMASK==1);
Hs_mat(jL)  = NaN;
Tp_mat(jL)  = NaN;
mwd_mat(jL) = NaN;

subplot(3,1,1);
%GEN_plot_matrix(X/1e3,Y/1e3,cice,[-1 1]);
H  = pcolor(X/1e3,Y/1e3,Hs_mat);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
colorbar;
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
ttl   = title('{\itH}_{\rms}, m');
GEN_font(ttl);
colorbar;
GEN_font(gca);
%%
subplot(3,1,2);
%GEN_plot_matrix(X/1e3,Y/1e3,Dmax,[-1 500]);
H  = pcolor(X/1e3,Y/1e3,Tp_mat);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
ttl   = title('\itT_{\rmp}');
GEN_font(ttl);
colorbar;
GEN_font(gca);
%%
subplot(3,1,3);
%GEN_plot_matrix(X/1e3,Y/1e3,S(:,:,jpp,jmwd));
H  = pcolor(X/1e3,Y/1e3,mwd_mat);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('<theta>, \rmm');
GEN_font(ttl);
