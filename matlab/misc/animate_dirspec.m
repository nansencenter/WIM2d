function animate_dirspec(x,th_vec,Sdir,legend_text)
%% animate_dirspec.m
%% Author: Timothy Williams
%% Date: 20160905, 11:30:20 CEST
%% Sdir  = ndir x nx matrix

nx    = length(x);
ndir  = length(th_vec);

if iscell(Sdir)
   nS    = length(Sdir);
   tmp   = zeros(ndir,nS,nx);
   leg_text = 'legend(''';
   for j=1:nS
      tmp2  = Sdir{j};%%ndir x nx
      for k=1:nx
         tmp(:,j,k)  = tmp2(:,k);
      end
      leg_text    = [leg_text,legend_text{j},''','''];
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
   clear tmp;
end


for j=1:nx
   plot(180/pi*th_vec,Sdir(:,:,j))
   %hold on;
   %plot(180/pi*th_vec,S_out((1:ndir),j)+S_out((1:ndir)+ndir,j),'-b')
   if ~isempty(leg_text)
      eval(leg_text);
   end
   title(['x = ',num2str(x(j))])
   hold off;

   dtheta   = th_vec(2)-th_vec(1);
   Hs_test  = 4*sqrt(dtheta*sum(Sdir(:,:,j)))

   pause(.1);
end
