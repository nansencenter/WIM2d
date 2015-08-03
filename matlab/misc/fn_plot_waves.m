%% waves_init.m
%% Author: Timothy Williams
%% Date: 20141016, 18:13:51 CEST

function fn_plot_waves(grid_prams,wave_fields);

X        = grid_prams.X;
Y        = grid_prams.Y;
[nx,ny]  = size(X);
%%
Hs_mat   = wave_fields.Hs;
Tp_mat   = wave_fields.Tp;
mwd_mat  = wave_fields.mwd;
%%
jL          = find(grid_prams.LANDMASK==1);
Hs_mat(jL)  = NaN;
Tp_mat(jL)  = NaN;
mwd_mat(jL) = NaN;

subplot(3,1,1);
if ny==1
   plot(X/1e3,Hs_mat);
   GEN_proc_fig('\itx, \rmkm','{\itH}_{\rms}, m');
else
   H  = pcolor(X/1e3,Y/1e3,Hs_mat);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   ttl   = title('{\itH}_{\rms}, m');
   GEN_font(ttl);
   colorbar;
   GEN_font(gca);
end
%%
subplot(3,1,2);
if ny==1
   plot(X/1e3,Tp_mat);
   GEN_proc_fig('\itx, \rmkm','{\itT}_{\rmp}, s');
else
   H  = pcolor(X/1e3,Y/1e3,Tp_mat);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   ttl   = title('\itT_{\rmp}');
   GEN_font(ttl);
   colorbar;
   GEN_font(gca);
end
%%
subplot(3,1,3);
if ny==1
   plot(X/1e3,mwd_mat);
   GEN_proc_fig('\itx, \rmkm','<\theta>, degrees');
else
   H  = pcolor(X/1e3,Y/1e3,mwd_mat);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   colorbar;
   GEN_font(gca);
   ttl   = title('<\theta>, degrees');
   GEN_font(ttl);
end
