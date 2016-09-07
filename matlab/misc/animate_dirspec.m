function [test_out,x,th_vec,Sdir_out,legend_text]   = animate_dirspec(x,th_vec,Sdir,legend_text,space_or_time)
%% animate_dirspec(x,th_vec,Sdir,legend_text)
%% Author: Timothy Williams
%% Date: 20160905, 11:30:20 CEST
%% Sdir  = ndir x nx matrix
%% or cell with multiple versions of Sdir
%% - in that case can use legend_text to make a legend

if nargin==0
   params_in.SCATMOD          = 3;
   params_in.DO_BREAKING      = 0;
   params_in.Dmax_init        = 100;
   params_in.PLOT_INIT        = 0;
   params_in.PLOT_PROG        = 0;
   params_in.PLOT_FINAL       = 0;
   params_in.itest            = 65;
   params_in.jtest            = 5;
   if 0
      space_or_time              = 'space';
      params_in.duration_hours   = 6;
   else
      space_or_time              = 'time';
      params_in.duration_hours   = 24;
   end

   [out_fields,diagnostics,wave_stuff,grid_prams]  = run_WIM2d(params_in);
   figure(101);
   fn_fullscreen;

   USE_EBS  = (params_in.SCATMOD==3);
   th_vec   = -pi/180*(90+wave_stuff.dirs);
   if strcmp(space_or_time,'space')
      %%animate as travelling into ice
      x  = grid_prams.X(:,1)/1e3;
      if USE_EBS
         %%wave_stuff.dir_spec : nx * ny * ndirs *nfreq
         Sdir  = {squeeze(wave_stuff.dir_spec(:,1,:,1))',...
                  squeeze(wave_stuff.dir_spec_scattered(:,1,:,1))',...
                  squeeze(wave_stuff.dir_spec(:,1,:,1))'+...
                  +squeeze(wave_stuff.dir_spec_scattered(:,1,:,1))'}
         legend_text = {'normal','scattered','total'};
      else
         Sdir        = squeeze(wave_stuff.dir_spec(:,1,:,1))';
         legend_text = [];
      end
   else
      %%animate with time - at (itest,jtest)
      x     = diagnostics.dirspec_snapshots.t;
      if USE_EBS
         Sdir  = {diagnostics.dirspec_snapshots.dir_spec,...
                  diagnostics.dirspec_snapshots.dir_spec_scattered,...
                  diagnostics.dirspec_snapshots.dir_spec+...
                  +diagnostics.dirspec_snapshots.dir_spec_scattered};
         legend_text = {'normal','scattered','total'};
      else
         Sdir        = diagnostics.dirspec_snapshots.dir_spec;
         legend_text = [];
      end
   end
end

nx    = length(x);
ndir  = length(th_vec);

DO_LEG   = 0;
if exist('legend_text','var')
   DO_LEG   = iscell(Sdir);
end
if ~exist('space_or_time','var')
   space_or_time  = 'space';
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
   NF    = diag(1./norm_facs(:,:,j));
   labs  = {'\theta, degrees','D(\theta)'};
   fn_plot1d(180/pi*th_vec,Sdir(:,:,j)*NF,labs);
   if DO_LEG
      eval(leg_text);
   end

   if strcmp(space_or_time,'space')
      title(['x = ',num2str(x(j)/1e3),'km']);
   else
      title(['t = ',num2str(x(j)/3600),'h']);
   end

   hold off;
   x0 = 180/pi*(min(th_vec)-.5*dtheta);
   x1 = 180/pi*(max(th_vec)+.5*dtheta);
   xlim([x0,x1]);
   set(gca,'Xtick',-180:90:90);

   Hmax  = 0;
   for k=1:nS
      if DO_LEG
         test_out(j,k).name   = legend_text{k};
      end

      I0                = norm_facs(1,k,j);
      test_out(j,k).Hs  = 4*sqrt(I0);
      Hmax              = max(Hmax,test_out(j,k).Hs);
      if test_out(j,k).Hs>1e-4
         S_ = Sdir(:,k,j);
         Ic = dtheta*sum(cos(th_vec).*S_)/I0;
         Is = dtheta*sum(sin(th_vec).*S_)/I0;
         r1 = sqrt(Ic^2+Is^2);

         test_out(j,k).spreading = sqrt(2*(1-r1));
         disp(test_out(j,k));
         %if k==nS
         %   pause;
         %end
      else
         disp(test_out(j,k));
      end
   end

   if strcmp(space_or_time,'space')
      if Hmax<1e-4
         break;
      end
   end

   pause(.1);
end
