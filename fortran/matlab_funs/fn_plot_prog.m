function out = fn_plot_prog(grid_prams,n,outdir)

out   = fn_check_prog(outdir,n);
fn_fullscreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%do plots

%%fix positions so figures can be compared more easily between computers
pos1  = [0.130000000000000   0.583837209302326   0.334659090909091   0.341162790697674];
pos2  = [0.570340909090909   0.583837209302326   0.334659090909091   0.341162790697674];
pos3  = [0.130000000000000   0.110000000000000   0.334659090909091   0.341162790697674];
pos4  = [0.570340909090909   0.110000000000000   0.334659090909091   0.341162790697674];

%subplot(2,2,1);
subplot('position',pos1);
H  = pcolor(grid_prams.X/1e3,grid_prams.Y/1e3,out.Hs);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\itH}_{\rm s}, m');
GEN_font(ttl);

%subplot(2,2,2);
subplot('position',pos2);
H  = pcolor(grid_prams.X/1e3,grid_prams.Y/1e3,out.Dmax);
%caxis([0 250]);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\itD}_{\rm max}, m');
GEN_font(ttl);

%subplot(2,2,3);
subplot('position',pos3);
H  = pcolor(grid_prams.X/1e3,grid_prams.Y/1e3,out.tau_x);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\tau}_{x}, Pa');
GEN_font(ttl);

%subplot(2,2,4);
subplot('position',pos4);
H  = pcolor(grid_prams.X/1e3,grid_prams.Y/1e3,out.tau_y);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\tau}_{y}, Pa');
GEN_font(ttl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
