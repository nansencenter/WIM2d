addpath('misc_numeric');
addpath('misc');
addpath('roots');
addpath('invacuo_hinges');
addpath('../misc');
addpath('../attenuation_youngs');
hmdir    = getenv('HOME');
gitdir   = getenv('GIT_REPOS');

%%rm ISSM paths
issmdir  = getenv('ISSM_DIR');
spath    = path;
jissm    = strfind(path,issmdir);
for j=1:length(jissm)

   path2 = spath(jissm(j):end);
   jcol  = strfind(path2,':');
   if length(jcol)>0
      path2 = path2(1:jcol(1)-1);
   end

   disp(' ');
   disp(['rmpath ',path2]);
   rmpath(path2)
end
disp(' ');
