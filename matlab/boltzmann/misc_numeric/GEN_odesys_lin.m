function bn=GEN_odesys_lin(Af_fxn,Nterms,interval,varargin)
%% CALL: an=GEN_odesys_lin(Af_fxn,Nterms,interval,varargin)
%%  finds general solution to Y'-A*Y=f;
%%   gives solution as
%%    Y_j(x)=\sum_{n=0}^Nterms a_n^(j)T_n(t),
%%     where the a_0^(j) are arbitrary constants,
%%      t=-1+2*(x-a)/(b-a) => x=a+(b-a)*(t+1)/2;
%% interval =[a,b] - domain eqn is over;
%% include extra arguments for Af_fxn after interval

DO_TEST  = 0;
if nargin==0%% USE SOME TEST INPUTS:
   Nterms   = 5;
   DO_TEST  = 2;
   if DO_TEST==1
      Af_fxn   = @test_fn1;
   elseif DO_TEST==2
      Af_fxn   = @test_fn2;
   end
end

if ~exist('interval','var')
   a  = -1;
   b  = 1;
else
   a  = interval(1);
   b  = interval(2);
end

Nint     = max(50,2*Nterms);
alpU     = 1;
[tt,ww]  = OP_numint_gegenbauer(alpU,Nint);
xx       = a+.5*(b-a)*(1+tt);
%%
[ipU,hn] = OP_inprod_gegenbauer(tt,ww,alpU,Nterms-1);
TnVals   = OP_interp_chebyshev(tt,{Nterms});

[Amatrix,forcing] = feval(Af_fxn,xx,varargin{:});
M                 = size(forcing,2);%%size of system
%%
nvec  = (1:Nterms)';
ipU2  = diag(.5*(b-a)./nvec)*ipU;
%%
KMAT     = zeros(M*Nterms,M*(Nterms+1));
FORCING  = KMAT(:,1);
%%
for j=1:M
   j_rows          = nvec+(j-1)*Nterms;
   FORCING(j_rows) = ipU2*forcing(:,j);
   %%
   for r=1:M
      j_cols               = [1;nvec+1]+(r-1)*(Nterms+1);
      s                    = j+(r-1)*M;
      KMAT(j_rows,j_cols)  = ipU2*diag(Amatrix(:,s))*TnVals;
   end
end

j0          = 1:(Nterms+1):(M*(Nterms+1));
FORCING     = [FORCING,KMAT(:,j0)];
KMAT(:,j0)  = [];
an0         = ( eye(M*Nterms)-KMAT )\FORCING;

II = [zeros(M,1),eye(M)];
an = zeros(size([II;an0]));
for j=1:M
   Z              = zeros(Nterms+1,M+1);
   Z(1,:)         = II(j,:); 
   Z(2:end,:)     = an0(nvec+(j-1)*Nterms,:);
   bn(j).coeffs   = Z;
end

if DO_TEST==2
   for j=1:M
      subplot(1,M,j), plot(xx,TnVals*bn(j).coeffs);
   end

   %%check particular solution
   subplot(1,M,1);
   hold on;
   plot(xx,sin(xx),'.c')
   subplot(1,M,2);
   hold on;
   plot(xx,sin(xx),'.c')

   %%check general soln 1
   subplot(1,M,1);
   c1 = TnVals(1,:)*bn(1).coeffs(:,2);
   c0 = exp(-cos(xx(1)));
   plot(xx,c1/c0*exp(-cos(xx)),'.m');
   hold off;

   %%check general soln 2
   subplot(1,M,2);
   c1 = TnVals(1,:)*bn(2).coeffs(:,3);
   c0 = exp(sin(xx(1)));
   plot(xx,c1/c0*exp(sin(xx)),'.m');
   hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,f]=test_fn1(x,varargin)

A  = cos(x);
f  = 0*x;
%% general solution is c*e^{sin(x)}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,f]=test_fn2(x,varargin)

% A = [sin(x),0;
%      0,cos(x)];
% f = [cos(x)-sin(x).^2, cos(x)-cos(x).*sin(x)];

z  = 0*x;
f  = [cos(x)-sin(x).^2, cos(x)-cos(x).*sin(x)];
A  = [sin(x),z,z,cos(x)];

%% general solution is [sin(x)+c1*e^{-cos(x)};sin(x)+c2*e^{sin(x)}]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
