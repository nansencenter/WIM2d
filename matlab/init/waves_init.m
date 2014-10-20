%% waves_init.m
%% Author: Timothy Williams
%% Date: 20141016, 18:13:51 CEST

function wave_fields = waves_init(grid_prams,wave_prams,ice_fields,OPT);

do_test  = 0;
if nargin==0
   do_test     = 1;
   %%
   nx          = 51;
   dx          = 4e3;
   %%
   OPT         = 0;
   grid_prams  = struct('nx',nx,'ny',nx,...
                        'dx',dx,'dy',dx);
   grid_prams  = get_grid(grid_prams,OPT)
   %%
   h           = 2;
   c           = 0.75;
   ice_prams   = struct('c',c,'h',h)
   ice_fields  = iceinit(ice_prams,grid_prams,OPT);
   %%
   Hs          = 3;
   Tp          = 11;
   wave_prams  = struct('Hs',Hs,'Tp',Tp);
end

if OPT==0
   mwd   = 135;%%waves-from direction
elseif OPT==1
   mwd   = -90;%%waves-from direction
elseif OPT==2
   mwd   = 135;%%waves-from direction
end

X  = grid_prams.X;
Y  = grid_prams.Y;
nx = grid_prams.nx;
ny = grid_prams.ny;
xm = max(X(:));
ym = max(Y(:));

if OPT==0
   WAVE_MASK   = (X<-.4*xm)|(Y>.4*ym);
   jWV         = find(WAVE_MASK==1);
elseif OPT==1
   WAVE_MASK   = (X<-.6*xm)*1.0;
   jWV         = find(WAVE_MASK==1);
elseif OPT==2
   WAVE_MASK   = (X<-.25*xm)|(Y>.25*ym);
   jWV         = find(WAVE_MASK==1);
end

Hs = wave_prams.Hs;
Tp = wave_prams.Tp;

Hs_mat   = zeros(nx,ny);
Tp_mat   = zeros(nx,ny);
mwd_mat  = zeros(nx,ny);
%%
Hs_mat(jWV)    = Hs;
Tp_mat(jWV)    = Tp;
mwd_mat(jWV)   = mwd;
wave_fields = struct('Hs',Hs_mat,...
                     'Tp',Tp_mat,...
                     'mwd',mwd_mat,...
                     'WAVE_MASK',WAVE_MASK);

if do_test
   jI          = find(ice_fields.ICE_MASK==1);
   Hs_mat(jI)  = NaN;
   Tp_mat(jI)  = NaN;
   mwd_mat(jI) = NaN;

   subplot(1,3,1);
   %GEN_plot_matrix(X/1e3,Y/1e3,cice,[-1 1]);
   H  = pcolor(X/1e3,Y/1e3,Hs_mat);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   colorbar;
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   ttl   = title('{\itH}_{\rms}, m');
   GEN_font(ttl);
   colorbar;
   %%
   subplot(1,3,2);
   %GEN_plot_matrix(X/1e3,Y/1e3,Dmax,[-1 500]);
   H  = pcolor(X/1e3,Y/1e3,Tp_mat);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   ttl   = title('\itT_{\rmp}');
   GEN_font(ttl);
   colorbar;
   %%
   subplot(1,3,3);
   %GEN_plot_matrix(X/1e3,Y/1e3,S(:,:,jpp,jmwd));
   H  = pcolor(X/1e3,Y/1e3,mwd_mat);
   set(H,'EdgeColor', 'none');
   daspect([1 1 1]);
   GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
   colorbar;
   ttl   = title('<theta>, \rmm');
   GEN_font(ttl);
   %%
end
