function S_out = EBS_steady(S_in,x,inputs)
%% y'=M_bolt*y

DO_TEST  = 0;
if nargin==0
   ndir  = 2^4;
   qs    = 1;
   M_K   = qs*ones(ndir,ndir)/ndir;
   Z     = zeros(ndir,ndir);
   %%
   dtheta         = 2*pi/ndir;
   inputs.th_vec  = dtheta*(0:ndir-1)'+.5*dtheta;%%dodge \pm\pi/2 - is this necessary?
   wavdir         = 90-180/pi*inputs.th_vec;
   M_filter       = EBS_filter_step(wavdir);%% 1 -> scattered SDF; 0 -> normal SDF
   inputs.M_bolt  = [  (1-M_filter).*M_K - qs*eye(ndir) , Z ;
                           M_filter.*M_K                , Z ];
   %%
   np = 25;
   x  = linspace(0,np,np+1);
   %x = [0,.5]';

   S_in           = zeros(2*ndir,1);
   S_in([1,ndir]) = .5/dtheta;
   DO_TEST        = 1;
   itest          = 9;
end

flds  = fieldnames(inputs);
Nf    = length(flds);
for j=1:Nf
   v     = flds{j};
   cmd   = [v,' = inputs.',v,';'];
   eval(cmd);
end
clear inputs;
%sum(M_bolt),pause

ndir  = round(size(M_bolt,1)/2);
nx    = length(x);

Dc    = diag(cos(th_vec));
Lmat  = [Dc,0*Dc;0*Dc,Dc];

if 1
   %% if we dodge \pm\pi/2:
   Dc       = diag(cos(th_vec));
   Lmat     = [Dc,0*Dc;0*Dc,Dc];
   [U,D]    = eig(M_bolt,Lmat);
   %Dc,Lmat\M_bolt,pause

   %%sort
   e_vals      = diag(D);
   [dummy,idx] = sort(abs(e_vals));
   e_vals      = e_vals(idx);
   U           = U(:,idx);

   jm    = find(e_vals<0);     %%decaying e-vals
   Jfwd  = find(cos(th_vec)>0);%%only match fwd-going waves at x=0
   cc    = U(Jfwd,jm)\S_in(Jfwd);
   S_out = zeros(2*ndir,nx);
   if DO_TEST
      S0 = U(:,jm)*cc
      %sum(U)
      %sum(U(:,jm))
      %pause
   end
   for j=1:nx
      S_out(:,j)  = U(:,jm)*diag(exp(e_vals(jm)*x(j)))*cc;
   end

else
%  %% doesn't work
%  [U_,D_,V_]  = svd(Lmat)%%Lmat=U_*D_*V_'
%  d_          = diag(D_)
%  %[d_,idx] = sort(d_,'ascend'),pause
%  JZ          = find(abs(d_)<1e-12);
%  JP          = find(abs(d_)>=1e-12);
%  d_i         = d_;
%  d_i(JP)     = 1./d_(JP);
%  
%  %%solve D_*(V_'*y)=U_'*M_bolt*V_*(V_'*y)
%  L11   = D_(JP,JP),inv(L11)
%  L12   = D_(JP,JZ)
%  %L21   = Lmat_(JZ,JP);%= 0
%  %L22   = Lmat_(JZ,JZ);%= 0
%  
%  Rmat  = U_'*M_bolt*V_
%  R11   = Rmat(JP,JP)
%  R12   = Rmat(JP,JZ)%=0
%  R21   = Rmat(JZ,JP)
%  R22   = Rmat(JZ,JZ)%=0
%  null21   = null(R21)
%  
%  [U,D] = eig(R11,L11);
%  e_vals   = diag(D),U,R21*U,pause;
%  
%  if ~exist('e_vals','var')
%     [U,D]       = eig( M_bolt,Lmat );
%     e_vals      = diag(D);
%     e_vals,U
%     clear D;
%  else
%     U  = e_vecs;
%     clear e_vecs;
%  end
%  
%  e_tol = 1e-12;
%  if ~isempty(find(e_vals>e_tol))
%     disp('Eigen-values:')
%     disp(e_vals)
%     error('\nShouldn''t be any positive eigenvalues')
%  else
%     e_vals   = min(e_vals,0);
%  end

end


if DO_TEST

   % =================================================================
   %% numerical solution
   %tic
   if ndir<=16
      soln  = ode45(@(x,y)rhs(x,y,Lmat\M_bolt),[0,max(x)],S0);
   end
   %toc

   if length(x)<10
      ls = '.b';
   else
      ls = '-b';
   end

   %%plot one of the "incident" waves
   figure(1);
   subplot(1,2,1);
   plot(x,S_out(itest,:),ls);
   hold on;
   yl = get(gca,'ylim');
   ttl  = ['\theta=',num2str(wavdir(itest))];
   title(ttl);

   plot(soln.x,soln.y(itest,:),'--r');
   legend('spectral','numeric');
   ylim(yl);
   xlabel('x')
   ylabel('SDF (normal)')
   hold off

   %%plot one of the "scattered" waves
   subplot(1,2,2);
   itest = ndir+itest;
   plot(x,S_out(itest,:),ls);
   yl = get(gca,'ylim');
   hold on;
   xlabel('x')
   ylabel('SDF (scattered)')
   title(ttl);

   plot(soln.x,soln.y(itest,:),'--r');
   legend('spectral','numeric');
   ylim(yl);
   hold off
   % =================================================================


   % =================================================================
   %%plot animation of dir spec 
   figure(2);
   subplot(1,2,1);
   %% keep x=0
   fac   = 0;
   plot(180/pi*th_vec,S0((1:ndir))+fac*S0((1:ndir)+ndir),'-k');
   hold on;
   plot(180/pi*th_vec,S0((1:ndir))+S0((1:ndir)+ndir),'-b');
   plot(180/pi*th_vec,S_in((1:ndir))+S_in((1:ndir)+ndir),'-m');
   hold off;
   legend('normal','total','B Con''s');
   title('x = 0')

   subplot(1,2,2);
   %% loops over all x
   Sdir_ = {S_out((1:ndir),:)+fac*S_out((1:ndir)+ndir,:),...
            S_out((1:ndir),:)+S_out((1:ndir)+ndir,:)};
   legend_text = {'normal','total'};
   animate_dirspec(x,th_vec,Sdir_,legend_text);
   % =================================================================
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dy=rhs(x,y,M_bolt)
dy = M_bolt*y;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,n_term] = mat_exp(M,n_term)

N  = size(M,1);
t  = M;
y  = eye(N)+t;

if exist('n_term','var');
   for n=2:n_term
      t  = M*t/n;%%M^n/(n!)=M/n*M^(n-1)/(n-1)!
      y  = y+t;
   end
else
   nM       = norm(M);
   n_term   = 200;
   for n=2:n_term
      t  = M*t/n;%%M^n/(n!)=M/n*M^(n-1)/(n-1)!
      y  = y+t;
      if norm(t)<1e-12*nM
         %disp(['stopping at n = ',num2str(n)]); pause;
         n_term   = n;
         break;
      end
   end
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
