%% /Users/twilliams/GITHUB-REPOSITORIES/WIM2d/fortran/run/run_WIM2d.m
%% Author: Timothy Williams
%% Date:   20141211, 17:49:41 CET

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make directories in which to put outputs:
dirs  = {'out',
         'out/log',
         'out/binaries',
         'out/binaries/prog'};

for j=1:length(dirs)
   dd = dirs{j};
   if ~exist(dd)
      mkdir(dd);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
   %% run non-io 2d WIM:
   disp('Running 2d WIM (no in/out)...');
   disp(' ');
   WIM2d_run_mex;
   %%
   !cat out/log/wim2d.log
   %%
   disp(' ');
   disp('Finished running 2d WIM.');
   disp(' ');
   disp('Outputs in:');
   disp(dirs);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
   %% test:
   cd test;
   startup_local;
   test_WIM2d_F;
   cd ..;
end
