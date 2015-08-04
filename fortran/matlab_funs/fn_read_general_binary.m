function [out_fields,info] = fn_read_general_binary(afile,WANT_TOUT)

%% load from binaries;
bfile = [afile(1:end-2),'.b'];
if ~exist(afile)
   error([afile,' not present.'])
end
if ~exist(bfile)
   error([bfile,' not present.'])
end

if ~exist('WANT_TOUT')
   WANT_TOUT   = 0;
end

%% get basic info from bfile
%% - eg:
%% > 04       Nrecs  # Number of records
%% > 01       Norder # Storage order
%% > 150      nx     # Record length in x direction (elements)
%% > 050      ny     # Record length in y direction (elements)
%% > 0.00     t_out  # Model time at output (s)

bid   = fopen(bfile);
C     = textscan(bid,'%2.2d',1);
nrec  = C{1};
C     = textscan(bid,'%s',5);
%%
C     = textscan(bid,'%2.2d',1);
Nord  = C{1};
C     = textscan(bid,'%s',12);
%%
C  = textscan(bid,'%4.4d',1);
nx = C{1};
C  = textscan(bid,'%s',8);
%%
C  = textscan(bid,'%4.4d',1);
ny = C{1};
C  = textscan(bid,'%s',8);
%%
if WANT_TOUT==1
   C           = textscan(bid,'%9.1f',1);
   t_out       = C{1};
   C           = textscan(bid,'%s',7);
   info.t_out  = t_out;
end

%% output the above info
info.Nrecs  = nrec;
info.Norder = Nord;
info.nx     = nx;
info.ny     = ny;

%% Find "Record number and name"
C  = textscan(bid,'%s',1);
while ~strcmp(C{1},'Record');
   C  = textscan(bid,'%s',1);
end
C  = textscan(bid,'%s',3);

%% get variable names
vlist = cell(nrec,1);
for n=1:nrec
   C        = textscan(bid,'%2.2d %s',1);
   vlist{n} = C{2}{1};
end
fclose(bid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%read afile
fmt         = 'float32';
aid         = fopen(afile);
out_arrays  = zeros(nx,ny,nrec);
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
