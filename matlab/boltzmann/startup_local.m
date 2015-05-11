addpath('misc');
addpath('roots');
issmdir  = getenv('ISSM_DIR');
if strcmp(issmdir,'')
   rmpath([issmdir,'/externalpackages/matlab/install/toolbox/local']);
end
