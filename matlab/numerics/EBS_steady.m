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
   %M_filter       = EBS_filter(wavdir,'step')%% 1 -> scattered SDF; 0 -> normal SDF
   M_filter       = EBS_filter(wavdir,'delta')%% 1 -> scattered SDF; 0 -> normal SDF
   %M_filter       = EBS_filter(wavdir,'cos')%% 1 -> scattered SDF; 0 -> normal SDF
   inputs.M_bolt  = [  (1-M_filter).*M_K - qs*eye(ndir) , Z ;
                           M_filter.*M_K                , Z ];
   %%
   np = 12;
   x  = linspace(0,np,np+1)';
   %x = [0,.5]';

   if 0
      %% delta fxn
      S_in     = zeros(2*ndir,1);
      S_in(1)  = 1/dtheta;
   elseif 0
      %% wider delta fxn
      S_in           = zeros(2*ndir,1);
      S_in([1,ndir]) = .5/dtheta;
   else
      %%cos^2
      S_in           = zeros(2*ndir,1);
      S_in(1:ndir)   = 2/pi*max(0,cos(inputs.th_vec)).^2;
   end
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

ndir  = round(size(M_bolt,1)/2);
nx    = length(x);
J1    = (1:ndir);
J2    = J1+ndir;

Dc    = diag(cos(th_vec));
Lmat  = [Dc,0*Dc;0*Dc,Dc];

if 0
   %% look at e-vals/e-vecs of orig system
   [u,d] = eig(M_K-qs*eye(ndir),Dc);
   %sum(M_K-qs*eye(ndir))
   E_vals      = diag(d);
   [dummy,idx] = sort(abs(E_vals));
   E_vals      = E_vals(idx)
   u           = u(:,idx)
   pause
end

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
   %e_vals,U,pause

   jm    = find(e_vals<0);     %%decaying e-vals
   Jfwd  = find(cos(th_vec)>0);%%only match fwd-going waves at x=0
   if 0
      %%only use the "real" energy to fit the inc wave spec
      lw1   = 1.5; %linewidth-real
      lw2   = 1.25;%linewidth-total
      cc    = U(Jfwd,jm)\S_in(Jfwd);
   else
      %%use the total energy to fit the inc wave spec
      lw1   = 1.25;%linewidth-real
      lw2   = 1.5; %linewidth-total
      cc    = ( U(Jfwd,jm) + U(Jfwd+ndir,jm) )\S_in(Jfwd);
   end
   %U,sum(U(:,jm)),pause
   S0       = U(:,jm)*cc;
   cU       = U(:,jm)*diag(cc);
   cU_main  = cU(J1,1)+cU(J2,1);
   if find(cU_main<1e-8)
      disp('Main e-vec:');
      disp(cU_main);
      error('Main eigenvector has negative energy');
   end
   %cU,e_vals(jm),qs,pause
   S_out = zeros(2*ndir,nx);
   if 0*DO_TEST
      T        = 0*S_in;
      T(Jfwd)  = 1;
      [S0,S_in,T]
      pause
      %sum(U)
      %sum(U(:,jm))
      %pause
   end


   for j=1:nx
      S_out(:,j)  = cU*exp(e_vals(jm)*x(j));

      if ~isempty(find(S_out(J1,j)+S_out(J2,j)<-1e-8))
         %%check total energy >0
         disp(['x=',num2str(x(j))]);
         disp(' ');
         disp('Dir spec:');
         disp(S_out(:,j));
         disp('Total dir spec:');
         disp(S_out(J1,j)+S_out(J2,j));
         disp('Total energy < 0');
         S0(J1)+S0(J2)
         pause;
         %error('Total energy < 0');
      end
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

   if 0
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
   end


   if 1
      % =================================================================
      %%plot animation of dir spec 
      figure(2);
      fn_fullscreen;
      subplot(1,2,1);
      %% keep x=0
      fac   = 0;
      plot(180/pi*th_vec,S0((1:ndir))+fac*S0((1:ndir)+ndir),'-k');
      labs        = {'\theta, degrees','D(\theta)'};
      Y_inc       = nan*zeros(ndir,1);
      Y_inc(Jfwd) = S_in(Jfwd);
      Y  = [S0((1:ndir))+    S0((1:ndir)+ndir),...
            S0((1:ndir))+fac*S0((1:ndir)+ndir),...
                             S0((1:ndir)+ndir),...
            cU_main,...
            Y_inc]

      %% normalise to give D(\theta)
      nn       = 4;
      I0       = dtheta*sum(Y(:,1:nn));
      Y(1:nn)  = Y(1:nn)*diag(1./I0);
      Y(nn+1)  = Y(nn+1)/(dtheta*sum(S_in));

      H  = fn_plot1d(180/pi*th_vec,Y,labs);
      set(H(1),'color','k','linewidth',lw2);
      set(H(2),'color','b','linewidth',lw1);
      set(H(3),'color','c');
      set(H(4),'color','y');
      set(H(5),'color','m','linestyle','--','linewidth',1.05);
      lg = legend('total','normal','scattered','Mode 1','B Con''s');
      set(lg,'location','EastOutside');
      title('x = 0km')
      xlim(180/pi*[(th_vec(1)-.5*dtheta),th_vec(end)+.5*dtheta]);
      pause;

      subplot(1,2,2);
      %% loops over all x
                                     {cU_main,1+0*x'}
      Sdir_ = {S_out((1:ndir),:)+    S_out((1:ndir)+ndir,:),...
               S_out((1:ndir),:)+fac*S_out((1:ndir)+ndir,:),...
                                     S_out((1:ndir)+ndir,:),...
                                     cU_main*(1+0*x')}
      legend_text = {'total','normal','scattered','Mode 1'};
      animate_dirspec(1e3*x,th_vec,Sdir_,legend_text);
      % =================================================================
   end
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dy=rhs(x,y,M_bolt)
dy = M_bolt*y;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
