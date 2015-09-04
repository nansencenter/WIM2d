function [an,E_coeffs] = Boltzmann_steady_numeric(th_vec,alp_scat,alp_dis)

if nargin==0
   %%test inputs
   Ndir     = 2^5
   dth      = 2*pi/Ndir;
   th_vec   = (0:Ndir-1)'*dth+.5*dth;
   %%
   period   = 12;
   h        = 2;
   young    = 5.0e9;
   visc_rp  = 0;
   om       = 2*pi/period;
   %%
   kw = om^2/9.81;
   cg = om/2/kw;
   %%
   [damping_rp,kice,kwtr,int_adm,NDprams,...
      alp_nd,modT,argR,argT] =...
         RT_param_outer(h,om,young,visc_rp);
   %%
   Dmean    = 100;
   conc     = .7;
   alp_scat = conc/Dmean*alp_nd
   alp_dis  = 2*conc*damping_rp
end

Hs    = 3;
width = 150e3;

%%solve governing eqn
Nterms   = 150;
interval = [0, width];
bn       = GEN_odesys_lin(@Af_fxn,Nterms,interval,...
                           th_vec,alp_scat,alp_dis);

%% this is the general soln:
%% 1st col of      bn(j).coeffs gives particular soln;
%% rest of cols of bn(j).coeffs give  solns to homogeneous problem;

%%solve edge condtions
TnVals0  = OP_interp_chebyshev(-1,{Nterms});
TnVals1  = OP_interp_chebyshev(1,{Nterms});

Ndir     = length(th_vec);
EdgeMat  = zeros(Ndir,Ndir+1);
EdgeVec  = zeros(Ndir,1);
C        = (Hs/4)^2*2/pi;
for j=1:Ndir
   th = th_vec(j);
   if cos(th)>0
      EdgeMat(j,:)   = TnVals0*bn(j).coeffs;
      EdgeVec(j)     = C*cos(th)^2;
   else
      EdgeMat(j,:)   = TnVals1*bn(j).coeffs;
   end
end

cn = EdgeMat(:,2:end)\(EdgeVec+EdgeMat(:,1));

%% energy
E_coeffs = zeros(Nterms+1,1);
dth      = 2*pi/Ndir;
for j=1:Ndir
   an(j).coeffs   = bn(j).coeffs*[1;cn];
   E_coeffs       = E_coeffs+dth*an(j).coeffs;
end

%% plot energy
figure(101);
nx       = 500;
tvec     = linspace(-1,1,nx)';
xvec     = interval(1)+(1+tvec)/2*(interval(2)-interval(1));
TnVals   = OP_interp_chebyshev(tvec,{Nterms});
Evec     = TnVals*E_coeffs;
plot(xvec/1e3,4*sqrt(abs(Evec)));
GEN_proc_fig('\itx, \rmkm','\itH_{\rms}, \rmm');

Y0 = zeros(Ndir,1);
Y1 = zeros(Ndir,1);
for j=1:Ndir
   Y0(j) = TnVals0*an(j).coeffs;
   Y1(j) = TnVals1*an(j).coeffs;
end

%% plot directional energy at edges
figure(102);
subplot(1,2,1);
plot(th_vec/pi,Y0);
hold on;
jfwd  = find(cos(th_vec)>0);
plot(th_vec(jfwd)/pi,EdgeVec(jfwd),'--r');
Hs0   = 4*sqrt(dth*sum(Y0));
title(['LH edge, Hs=',num2str(Hs0),' m']);
hold off;

subplot(1,2,2);
plot(th_vec/pi,Y1);
hold on;
plot(th_vec/pi,0*EdgeVec,'--r');
Hs1   = 4*sqrt(dth*sum(Y1));
title(['RH edge, Hs=',num2str(Hs1),' m']);
hold off;

%% compare deriv to source term
figure(103);
[Amat,fvec] = Af_fxn(0,th_vec,alp_scat,alp_dis);
Amat        = reshape(Amat,[Ndir,Ndir]);
Y           = zeros(nx,Ndir);
for j=1:Ndir
   Y(:,j)   = TnVals*an(j).coeffs;
end

Dc    = cg*diag(cos(th_vec));
dx    = xvec(2)-xvec(1);
dY    = (Y(3:end,:)-Y(1:end-2,:))/(2*dx);
xmid  = xvec(2:end-1);
RHS   = Y*Amat.'*Dc;
LHS   = dY*Dc;
if 1
   Ntst  = 1;
   plot(xmid/1e3,LHS(:,Ntst));
   hold on;
   plot(xvec/1e3,RHS(:,Ntst),'--r');
   ttl   = title(['testing ODE at \theta = ',num2str(180/pi*th_vec(Ntst)),'^o']);
   set(ttl,'fontname','times','fontsize',16);
   hold off;
end


%% write text file
tfil  = 'out/boltzmann_steady_numeric.dat';
disp(' ');
disp('Saving results to dat file:');
disp(tfil);
disp(' ');

fid   = fopen(tfil,'w');
fprintf(fid,'%s\n','# Steady-state results (Galerkin)');
%%
fprintf(fid,'\n%s\n','# File info:');
fprintf(fid,'%s%d\n','# Nvars : ',2);
fprintf(fid,'%s%d\n','# Lvars : ',nx);
%%
fprintf(fid,'\n%s\n','# Parameters:');
fprintf(fid,'%s%d\n','# Ndir              : ',Ndir);
fprintf(fid,'%s%d\n','# Nterms            : ',Nterms);
fprintf(fid,'%s%f\n','# visc_rp (Pa.s/m)  : ',visc_rp);
fprintf(fid,'%s%f\n','# period (s)        : ',period);
fprintf(fid,'%s%f\n','# concentration     : ',conc);
fprintf(fid,'%s%f\n','# thickness (m)     : ',h);
fprintf(fid,'%s%f\n','# Dmean (m)         : ',Dmean);
fprintf(fid,'%s%e\n','# young (Pa)        : ',young);
%%
fprintf(fid,'\n%s\n','# Variables:');
fprintf(fid,'%s\n','# x  (m)');
fprintf(fid,'%s\n','# Hs (m)');
fprintf(fid,'%s\n','#####################################################');

for j=1:nx
   Hs = 4*sqrt(Evec(j));
   fprintf(fid,'\n%f   %f',xvec(j),Hs);
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Amat,fvec] = Af_fxn(x,th_vec,alp_scat,alp_dis)
%d/dx(Y)=Amat*Y+fvec

nx    = length(x);
N     = length(th_vec);
fvec  = zeros(nx,N);
Amat  = zeros(nx,N*N);
%%
Mbolt = (alp_scat/N)*ones(N,N);
qtot  = alp_scat+alp_dis;
A     = diag(1./cos(th_vec))*(-qtot*eye(N)+Mbolt);

for n=1:nx
   Amat(n,:)   = A(:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
