function [out_fields,info] = fn_read_nextwim_binary(afile)

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
if info.nbytes==4
   fmt   = 'float32';
else
   fmt   = 'float64';
end
aid   = fopen(afile);
%%
for n=1:nrec
   vbl   = vlist{n}
   if ~isempty(strfind(vbl,'sc'))
      Nx = nx;
      Ny = ny;
   elseif ~isempty(strfind(vbl,'q'))
      Nx = nx+1;
      Ny = ny+1;
   elseif ~isempty(strfind(vbl,'p'))
      Nx = nx;
      Ny = ny;
   elseif ~isempty(strfind(vbl,'u'))
      Nx = nx+1;
      Ny = ny;
   elseif ~isempty(strfind(vbl,'v'))
      Nx = nx;
      Ny = ny+1;
   elseif ~isempty(strfind(vbl,'LAND'))
      Nx = nx;
      Ny = ny;
   end

   vec   = fread(aid,Nx*Ny,fmt);
   if Nord==1
      %% mat/fortran order
      mat = reshape(vec,Nx,Ny);
   else
      %% python/C order (need to transpose)
      mat = reshape(vec,Ny,Nx).';
   end
   out_fields.(vlist{n})   = mat;
end
%%
fclose(aid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
