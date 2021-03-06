function fn_save_binary(Froot,dims,year_info,varargin);

if nargin<3
   error('not enough inputs');
elseif nargin<=3
   error('nothing to do');
else
   pairs = varargin{:};
   nx    = dims(1);
   ny    = dims(2);
   nw    = dims(3);
   ndir  = dims(4);
end

Nrecs = length(pairs);
vlist = {};

if ~isempty(year_info)
   Froot = [Froot,year_info.date_string];
   stime = sprintf('%s   t_out  # Model time of output',year_info.date_string);
end

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
fprintf(bid,'%2.2d         Nrecs  # Number of records\n',Nrecs);
fprintf(bid,'%1.1d          Norder # Storage order [column-major (F/matlab) = 1, row-major (C) = 0]\n',1);
fprintf(bid,'%3.3d        nx     # Record length in x direction (elements)\n',nx);
if isempty(year_info)
   fprintf(bid,'%3.3d        ny     # Record length in y direction (elements)\n\n',ny);
else
   fprintf(bid,'%3.3d        ny     # Record length in y direction (elements)\n',ny);
   fprintf(bid,'%s\n\n',stime);
end
%%
fprintf(bid,'%s','Record number and name:');
for loop_i=1:length(vlist)
   fprintf(bid,'\n%2.2d       %s',loop_i,vlist{loop_i});
end
fclose(bid);
disp(['Saved ',bfile]);
disp(' ');
