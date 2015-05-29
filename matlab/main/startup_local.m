if 1
   %%addpaths
   fxn   = @addpath;
else
   %%remove them
   fxn   = @rmpath;
end

feval(fxn,'../init');
feval(fxn,'../attenuation_youngs');
feval(fxn,'../misc');
feval(fxn,'../numerics');
feval(fxn,'../numerics/testing');
feval(fxn,'../other_deps');

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
%gitdir   = '~/GITHUB-REPOSITORIES/';
%feval(fxn,[gitdir,'matlab/Semi-Infinite-Elastic-Plate/ND_progs']);
%feval(fxn,[gitdir,'matlab/Semi-Infinite-Elastic-Plate/GEN_progs']);
%feval(fxn,[gitdir,'matlab/Semi-Infinite-Elastic-Plate/GEN_progs/OP_progs']);
