function s1 = fn_check_final(outdir)

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
SCATMOD  = C{1};
%%
C  = textscan(bid,'%2.2d %s %s %s %s',1);
nw = C{1};
%%
C     = textscan(bid,'%3.3d %s %s %s %s',1);
ndir  = C{1};
fclose(bid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%read afile
fmt   = 'float32';
aid   = fopen(afile);
%%
s1.Dmax  = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
s1.tau_x = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
s1.tau_y = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
s1.Hs    = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
s1.Tp    = reshape( fread(aid,nx*ny,fmt) ,nx,ny );
%%
s1.SCATMOD     = SCATMOD;
s1.n_wave_freq = nw;
s1.n_wavdir    = ndir;
%%
fclose(aid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
