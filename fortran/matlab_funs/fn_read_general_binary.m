function out_fields = fn_read_general_binary(afile)

%% load from binaries;
bfile = [afile(1:end-2),'.b'];
if ~exist(afile)
   error([afile,' not present.'])
end
if ~exist(bfile)
   error([bfile,' not present.'])
end

%% get basic info from bfile
%% - eg:
%% > 04       Number of records
%% > 150      Record length in x direction (elements)
%% > 050      Record length in y direction (elements)
%% > 01       Option number for solver
%% > 01       Number of wave frequencies
%% > 016      Number of wave directions
bid   = fopen(bfile);
C     = textscan(bid,'%2.2d',1);
nrec  = C{1};
C     = textscan(bid,'%s',3);
%%
C     = textscan(bid,'%2.2d',1);
Nord  = C{1};
C     = textscan(bid,'%s',10);
%%
C  = textscan(bid,'%3.3d',1);
nx = C{1};
C  = textscan(bid,'%s',6);
%%
C  = textscan(bid,'%3.3d',1);
ny = C{1};
C  = textscan(bid,'%s',6);

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
