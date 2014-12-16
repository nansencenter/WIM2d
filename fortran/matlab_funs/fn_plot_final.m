function s1 = fn_plot_final(grid_prams,outdir)

if isstruct(outdir)
   %% already have the data required to plot;
   s1 = outdir;
   X  = grid_prams.X;
   Y  = grid_prams.Y;
else
   %% load from binaries;
   afile = [outdir,'/binaries/wim_out.a'];
   bfile = [outdir,'/binaries/wim_out.b'];

   %% get basic info from bfile
   %% - eg:
   %% > 04       Number of records
   %% > 150      Record length in x direction (elements)
   %% > 050      Record length in y direction (elements)
   %% > 01       Option number for solver
   %% > 01       Number of wave frequencies
   %% > 016      Number of wave directions
   bid   = fopen(bfile);
   C     = textscan(bid,'%2.2d %s %s %s',1);
   nrec  = C{1};
   %%
   C  = textscan(bid,'%3.3d %s %s %s %s %s %s',1);
   nx = C{1};
   %%
   C  = textscan(bid,'%3.3d %s %s %s %s %s %s',1);
   ny = C{1};
   %%
   C        = textscan(bid,'%2.2d %s %s %s %s',1);
   SOLVER   = C{1};
   %%
   C  = textscan(bid,'%2.2d %s %s %s %s',1);
   nw = C{1};
   %%
   C     = textscan(bid,'%3.3d %s %s %s %s',1);
   ndir  = C{1};
   fclose(bid);

   X  = grid_prams.X;
   Y  = grid_prams.Y;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%read afile
   s1    = struct('Dmax'   ,[],...
                  'tau_x'  ,[],...
                  'tau_y'  ,[],...
                  'Hs'     ,[]);
   fmt   = 'float32';
   aid   = fopen(afile);
   %%
   s1.Dmax  = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
   s1.tau_x = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
   s1.tau_y = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
   s1.Hs    = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
   %%
   s1.SOLVER      = SOLVER;
   s1.n_wave_freq = nw;
   s1.n_wavdir    = ndir;
   %%
   fclose(aid);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%Hs = s1.Hs,GEN_pause
%Dm = s1.Dmax,GEN_pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%do plots

%%fix positions so figures can be compared more easily between computers
pos1  = [0.130000000000000   0.583837209302326   0.334659090909091   0.341162790697674];
pos2  = [0.570340909090909   0.583837209302326   0.334659090909091   0.341162790697674];
pos3  = [0.130000000000000   0.110000000000000   0.334659090909091   0.341162790697674];
pos4  = [0.570340909090909   0.110000000000000   0.334659090909091   0.341162790697674];

%subplot(2,2,1);
subplot('position',pos1);
H  = pcolor(X/1e3,Y/1e3,s1.Hs);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\itH}_{\rm s}, m');
GEN_font(ttl);

%subplot(2,2,2);
subplot('position',pos2);
H  = pcolor(X/1e3,Y/1e3,s1.Dmax);
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
H  = pcolor(X/1e3,Y/1e3,s1.tau_x);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\tau}_{x}, Pa');
GEN_font(ttl);

%subplot(2,2,4);
subplot('position',pos4);
H  = pcolor(X/1e3,Y/1e3,s1.tau_y);
set(H,'EdgeColor', 'none');
daspect([1 1 1]);
GEN_proc_fig('\itx, \rmkm','\ity, \rmkm');
colorbar;
GEN_font(gca);
ttl   = title('{\tau}_{y}, Pa');
GEN_font(ttl);
