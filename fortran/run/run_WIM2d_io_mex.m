%% /Users/twilliams/GITHUB-REPOSITORIES/WIM2d/fortran/run/run_WIM2d.m
%% Author: Timothy Williams
%% Date:   20141211, 17:49:41 CET

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% make directories in which to put outputs:
% dirs  = {'out',
%          'out/log',
%          'out/binaries',
%          'out/binaries/prog'};
% 
% for j=1:length(dirs)
%    dd = dirs{j};
%    if ~exist(dd)
%       mkdir(dd);
%    end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1%%testing io for mex
   x  = rand(150,50,6);
   b  = x.*x;
   y  = 0*x;
   [y(:,:,1),y(:,:,2),y(:,:,3),y(:,:,4),y(:,:,5),y(:,:,6)]  = ...
      matsq(x(:,:,1),x(:,:,2),x(:,:,3),x(:,:,4),x(:,:,5),x(:,:,6));
   for j=1:6
      tst{j}   = norm(y(:,:,j)-b(:,:,j));
   end
   disp(tst);
   y(1:8,1:8,5)
   x(1:8,1:8,5).^2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
   %% run non-io 2d WIM:
   disp('Running 2d WIM (on in/out)...');
   disp(' ');
   addpath('../bin');
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

if 0
   %% test:
   cd test;
   startup_local;
   test_WIM2d_F;
   cd ..;
end
