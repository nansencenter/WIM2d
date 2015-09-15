if 1
   %%addpaths
   fxn   = @addpath;
else
   %%remove them
   fxn   = @rmpath;
end

w2d   = getenv('WIM2D_PATH');
feval(fxn,[w2d,'/matlab/init']);
feval(fxn,[w2d,'/matlab/attenuation_youngs']);
feval(fxn,[w2d,'/matlab/misc']);
feval(fxn,[w2d,'/matlab/numerics']);
feval(fxn,[w2d,'/matlab/numerics/testing']);
feval(fxn,[w2d,'/matlab/other_deps']);

%%mex functions and related
fdir  = '../../fortran';
feval(fxn,[fdir,'/bin']);        %% mex functions are here (when compiled)
feval(fxn,[fdir,'/matlab_funs']);%% functions to read binary outputs (eg)

%% remove ISSM paths
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

clear;
