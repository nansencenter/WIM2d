function S_out = evolve_spectral(S_in,x,inputs)
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
   M_filter       = EBS_filter_step(wavdir)
   inputs.M_bolt  = [     M_filter.*M_K - qs*eye(ndir) , Z ;
                      (1-M_filter).*M_K                , Z ];
   %%
   x        = linspace(0,5,50);
   S_in     = eye(2*ndir,1);
   DO_TEST  = 1;
   itest    = 1;
end

flds  = fieldnames(inputs);
Nf    = length(flds);
for j=1:Nf
   v     = flds{j};
   cmd   = [v,' = inputs.',v,';']
   eval(cmd);
end
clear inputs;

if ~exist('e_vals','var')
   [U,D]       = eig(M_bolt);
   e_vals      = diag(D);
   clear D;
else
   U  = e_vecs;
   clear e_vecs;
end

if ~isempty(find(e_vals>0))
   error('Shouldn''t be any positive eigenvalues')
end

ndir  = round(length(e_vals)/2);

%% reorder
[e_vals,II] = sort(e_vals,'descend');
U           = U(:,II);
e_vals

%% zero eigenvalues: don't evolve
jz    = (1:ndir)';
Uz    = U(:,jz);

e_vals(jz)  = [];
U(:,jz)     = []
e_fac       = 1;
if exist('q_scat')
   if 1
      %% scale eigenvalues so rate of decay is q_scat
      e_fac = q_scat/e_vals(1);
   else
      %% scale eigenvalues by q_scat
      e_fac = q_scat;
   end
end
e_vals   = e_fac*e_vals;
M_bolt   = e_fac*M_bolt;

%% non-repeated ones, <0
V  = [U(:,1),U(:,ndir)];%% ~1
W  = 0*V;                  %% ~t

%% repeated ones, <0
for j=1:round(ndir/2)-1
   V0 = U(:,2*j);
   ev = e_vals(2*j);
   %V1 = (M_bolt-ev*eye(2*ndir))\V0;%%least squares solution
   V1 = lsqlin(M_bolt-ev*eye(2*ndir),V0);%%least squares solution
   V  = [V,V0,V1];
   W  = [W,0*V0,V0];
end

%%expand in terms of solution
%%TODO repeated ones, =0
cc = [Uz,V]\S_in;
nx = length(x);

S_out    = zeros(2*ndir,nx);
for j=1:nx
   S_out(:,j)  = Uz*cc(jz) + (V+x(j)*W)*diag(exp(e_vals*x(j)))*cc(ndir+jz);
end

if DO_TEST
   %% numerical solution
   soln  = ode45(@(x,y)rhs(x,y,M_bolt),[min(x),max(x)],S_in);

   %%plot one of the "incident" waves
   subplot(1,2,1);
   plot(x,S_out(itest,:));
   hold on;

   plot(soln.x,soln.y(itest,:),'--r');
   hold off

   %%plot one of the "scattered" waves
   subplot(1,2,2);
   itest = ndir+itest;
   plot(x,S_out(itest,:));
   hold on;

   plot(soln.x,soln.y(itest,:),'--r');
   hold off
end

function dy=rhs(x,y,M_bolt)
dy = M_bolt*y;
