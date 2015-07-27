function fn_save_binary(Froot,dims,varargin);

if nargin<2
   error('not enough inputs');
elseif nargin==2
   disp('nothing to do');
   return
else
   pairs = varargin{:};
   nx    = dims(1);
   ny    = dims(2);
   nw    = dims(3);
   ndir  = dims(4);
end

Nrecs = length(pairs);
vlist = {};

%% save data
afile = [Froot,'.a'];
fmt   = 'float32';
aid   = fopen(afile,'wb');
for loop_i=1:Nrecs
   %pairs{loop_i}
   vlist{end+1}   = pairs{loop_i}{1};%name of vbl
   dat            = pairs{loop_i}{2};%array
   fwrite(aid,dat,fmt,'l');%'l'=little-endian (default)
end
fclose(aid);
disp(['Saved ',afile]);

%% save descriptor file
bfile = [Froot,'.b'];
%%
bid   = fopen(bfile,'w');
fprintf(bid,'%2.2d       Number of records\n',Nrecs);
fprintf(bid,'%1.1d        Storage order [column-major (F/matlab) = 1, row-major (C) = 0]\n\n',1);
fprintf(bid,'%3.3d      Record length in x direction (elements)\n',nx);
fprintf(bid,'%3.3d      Record length in y direction (elements)\n',ny);
fprintf(bid,'%2.2d       Number of wave frequencies\n',nw);
fprintf(bid,'%3.3d      Number of wave directions\n',ndir);
%%
fprintf(bid,'%s','Record number and name:');
for loop_i=1:length(vlist)
   fprintf(bid,'\n%2.2d       %s',loop_i,vlist{loop_i});
end
fclose(bid);
disp(['Saved ',bfile]);
disp(' ');
