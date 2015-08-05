function [out_fields,info] = fn_read_general_binary(afile)

%% load from binaries;
bfile = [afile(1:end-2),'.b'];
if ~exist(afile)
   error([afile,' not present.'])
end
if ~exist(bfile)
   error([bfile,' not present.'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get basic info from bfile
%% - eg:
%% > 04       Nrecs  # Number of records
%% > 01       Norder # Storage order
%% > 150      nx     # Record length in x direction (elements)
%% > 050      ny     # Record length in y direction (elements)
% > 0.00     t_out  # Model time at output (s)
[info,vlist]   = fn_bfile_info(bfile);
nrec           = info.Nrecs;
Nord           = info.Norder;
nx             = info.nx;
ny             = info.ny;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%read afile
fmt         = 'float32';
aid         = fopen(afile);
%%
for n=1:nrec
   vec   = fread(aid,nx*ny,fmt);
   if Nord==1
      %% mat/fortran order
      mat = reshape(vec,nx,ny);
   else
      %% python/C order (need to transpose)
      mat = reshape(vec,ny,nx).';
   end
   out_fields.(vlist{n})   = mat;
end
%%
fclose(aid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
