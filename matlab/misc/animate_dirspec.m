function [x,th_vec,Sdir_out,legend_text]   = animate_dirspec(x,th_vec,Sdir,legend_text)
%% animate_dirspec(x,th_vec,Sdir,legend_text)
%% Author: Timothy Williams
%% Date: 20160905, 11:30:20 CEST
%% Sdir  = ndir x nx matrix
%% or cell with multiple versions of Sdir
%% - in that case can use legend_text to make a legend

if nargin==0
   [out_fields,diagnostics,wave_stuff,grid_prams]  = run_WIM2d();
   figure(101);
   fn_fullscreen;
   x        = grid_prams.X(:,1)/1e3;
   th_vec   = -pi/180*(90+wave_stuff.dirs);
   %%wave_stuff.dir_spec : nx * ny * ndirs *nfreq
   Sdir  = {squeeze(wave_stuff.dir_spec(:,1,:,1))',...
            squeeze(wave_stuff.dir_spec_scattered(:,1,:,1))',...
            squeeze(wave_stuff.dir_spec(:,1,:,1))'+...
            +squeeze(wave_stuff.dir_spec_scattered(:,1,:,1))'}
   legend_text = {'normal','scattered','total'};
end

nx    = length(x);
ndir  = length(th_vec);

DO_LEG   = 0;
if exist('legend_text','var')
   DO_LEG   = iscell(Sdir);
end

if iscell(Sdir)
   Sdir_out = Sdir;
   nS       = length(Sdir);
   tmp      = zeros(ndir,nS,nx);
   leg_text = 'legend(''';
   for j=1:nS
      tmp2  = Sdir{j};%%ndir x nx
      for k=1:nx
         tmp(:,j,k)  = tmp2(:,k);
      end
      if DO_LEG
         leg_text    = [leg_text,legend_text{j},''','''];
      end
   end
   Sdir     = tmp;%% ndir x nS x nx
   leg_text = [leg_text(1:end-2),');']
   clear tmp;
else
   nS    = 1;
   tmp   = zeros(ndir,nS,nx);
   for k=1:nx
      tmp(:,1,k)  = Sdir(:,k);
   end
   Sdir     = tmp;%% ndir x nS=1 x nx
   leg_text = [];
   Sdir_out = Sdir;
   clear tmp;
end

dtheta      = th_vec(2)-th_vec(1);
norm_facs   = dtheta*sum(Sdir,1);%% 1 * nS * nx

for j=1:nx
   NF = diag(1./norm_facs(:,:,j));
   plot(180/pi*th_vec,Sdir(:,:,j)*NF)
   if DO_LEG
      eval(leg_text);
   end
   title(['x = ',num2str(x(j))])
   xlabel('\theta, degrees')
   ylabel('D(\theta)')
   hold off;

   Hs_test  = 4*sqrt(dtheta*sum(Sdir(:,:,j)))
   if max(Hs_test)<1e-4
      break;
   end

   pause(.1);
end
