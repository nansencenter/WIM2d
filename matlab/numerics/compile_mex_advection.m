P     = pwd;
me    = mfilename('fullpath');
J     = strfind(me,'/');
here  = me(1:J(end)-1);
cd(here);

mexfiles = dir('*mex.c');
for j=1:length(mexfiles)
   cmd   = ['mex COMPFLAGS=''$COMPFLAGS -Wall'' -largeArrayDims ',mexfiles(j).name];
   %cmd   = ['mex -v COMPFLAGS=''$COMPFLAGS -Wall -g3'' -largeArrayDims ',mexfiles(j).name];
   try
      disp(cmd);
      eval(cmd);
    catch ME
       disp('Compilation failed');
       disp(ME);
       cd(P);
    end
end
cd(P);
