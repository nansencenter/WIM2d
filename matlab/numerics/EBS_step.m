function S_out = EBS_step(S_in,x,inputs,method)
%% y'=M_bolt*y

DO_TEST  = 0;
if nargin==0
   ndir  = 16;
   qs    = 1;
   M_K   = qs*ones(ndir,ndir)/ndir;
   Z     = zeros(ndir,ndir);
   %%
   dtheta         = 360/ndir;
   wavdir         = 90-dtheta*(0:ndir-1)';
   M_filter       = EBS_filter(wavdir,'delta')
   inputs.M_bolt  = [     M_filter.*M_K - qs*eye(ndir) , Z ;
                      (1-M_filter).*M_K                , Z ];
   %%
   %x        = linspace(0,5,50);
   x        = 5
   S_in     = eye(2*ndir,1);
   DO_TEST  = 1;
   itest    = 1;

   method   = 'spectral'
end

flds  = fieldnames(inputs);
Nf    = length(flds);
for j=1:Nf
   v     = flds{j};
   cmd   = [v,' = inputs.',v,';'];
   eval(cmd);
end
clear inputs;

if ~exist('method','var')
   method   = 'matrix_exponential';
   %method   = 'spectral';
end

ndir  = round(size(M_bolt,1)/2);
nx    = length(x);

if ~exist('e_vals','var')
   [U,D]       = eig(M_bolt);
   e_vals      = diag(D);
   clear D;
else
   U  = e_vecs;
   clear e_vecs;
end

e_tol = 1e-12;
if ~isempty(find(e_vals>e_tol))
   disp('Eigen-values:')
   disp(e_vals)
   error('\nShouldn''t be any positive eigenvalues')
else
   e_vals   = min(real(e_vals),0);
end


if strcmp(method,'matrix_exponential')
   %tic
   S_out          = zeros(2*ndir,nx);
   [X,J]          = sort(x,'descend');
   [eM,nterm]     = mat_exp(X(1)*M_bolt);
   S_out(:,J(1))  = eM*S_in;
   for j=2:nx
      S_out(:,J(j))  = mat_exp(X(j)*M_bolt,nterm)*S_in;
   end
   %toc
elseif strcmp(method,'spectral')

   %tic

   %% reorder
   [e_vals,II] = sort(e_vals,'descend');
   U           = U(:,II);
   %e_vals,U,sum(M_bolt),sum(U)

   %%expand in terms of e-vecs
   %% - assume all e-vecs are lin indep (det(U)~=0)
   %% - can then diagonalise sys and solve
   cc = U\S_in;
   nx = length(x);

   S_out    = zeros(2*ndir,nx);
   for j=1:nx
      S_out(:,j)  = U*diag(exp(e_vals*x(j)))*cc;
   end
   %toc
else
   error('''method'' should be ''spectral'' or ''matrix_exponential''');
end

if DO_TEST
   %% numerical solution
   %tic
   soln  = ode45(@(x,y)rhs(x,y,M_bolt),[0,max(x)],S_in);
   %toc

   if length(x)==1
      ls = '.b';
   else
      ls = '-b';
   end

   %%plot one of the "incident" waves
   subplot(1,2,1);
   plot(x,S_out(itest,:),ls);
   hold on;

   plot(soln.x,soln.y(itest,:),'--r');
   hold off

   %%plot one of the "scattered" waves
   subplot(1,2,2);
   itest = ndir+itest;
   plot(x,S_out(itest,:),ls);
   hold on;

   plot(soln.x,soln.y(itest,:),'--r');
   hold off;

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
