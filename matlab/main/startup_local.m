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
feval(fxn,'../other_deps');

%gitdir   = '~/GITHUB-REPOSITORIES/';
%feval(fxn,[gitdir,'matlab/Semi-Infinite-Elastic-Plate/ND_progs']);
%feval(fxn,[gitdir,'matlab/Semi-Infinite-Elastic-Plate/GEN_progs']);
%feval(fxn,[gitdir,'matlab/Semi-Infinite-Elastic-Plate/GEN_progs/OP_progs']);
