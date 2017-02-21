P     = pwd;
me    = mfilename('fullpath');
J     = strfind(me,'/');
here  = me(1:J(end)-1);
cd(here);

mexfiles = dir('*mex.c');
%mexfiles = dir('*_1d_mex.c');
for j=1:length(mexfiles)
   cmd   = ['mex -largeArrayDims ',mexfiles(j).name];
   disp(cmd);
   eval(cmd);
end
cd(P);
