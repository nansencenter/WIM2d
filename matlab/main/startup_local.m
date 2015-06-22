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
