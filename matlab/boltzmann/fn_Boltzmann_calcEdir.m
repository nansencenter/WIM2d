function [I,incs] = fn_Boltzmann_calcEdir(eigen_info,x,Hs_inc)
%% CALL: I  = fn_Boltzmann_calcEdir(eigen_info,x,Hs_inc)
%% Inputs:
%% eigen_info   = struct('u0',NaN,...
%%                       'v0',NaN,...
%%                       'V',V,...
%%                       'D0',D0,...
%%                       'D1',D1,...
%%                       'coeffs',c,...
%%                       'width',wth,...
%%                       'absorb',absorb,...
%%                       'angles_on_pi',th_vec);
%% - structure calculated by fn_Boltzmann_Steady.m;
%% x is vector of positions to calculate I at (default: [0,width]);
%% Hs_inc is significant wave height of incident waves (default: 1m);
%%
%% Outputs:
%% I: rows are directional energy spectrum at each point of x

V        = eigen_info.V;
D0       = eigen_info.D0;
c        = eigen_info.coeffs;
wth      = eigen_info.width;
th_vec   = eigen_info.angles/pi;
absorb   = eigen_info.absorb;
if wth~='inf'
   u0 = eigen_info.u0;% one of zero eigenvectors 
   v0 = eigen_info.v0;% complementary function to u0 (satisfies same eqn but is linearly indep) 
   D1 = eigen_info.D1;% positive eigenvals
end

if ~exist('x','var')
   x  = [0;wth];
end

if ~exist('Hs_inc','var')
   Hs_inc   = 1;
end

refs  = find(or(th_vec>0.5,th_vec<-0.5)); %% indices of reflected waves
incs  = find(~or(th_vec>0.5,th_vec<-0.5));%% indices of incident  waves
Ndir  = length(th_vec);
Nx    = length(x);
I     = zeros(Nx,Ndir);

if wth=='inf'
   I0 = V*c;
   for loop_x=1:length(x)
      Ix = V*diag(exp(D0*x(loop_x)))*c;
      I(loop_x,:) = real(Ix).';
   end
else
   if absorb==0
      %% x  = 0
      I0 = [u0,V]*diag([1;1;exp(D0*0);exp(D1*(0-wth))])*c;
      for loop_x=1:length(x)
         Ix = [x(loop_x)*v0+u0,V]*diag([1;1;exp(D0*x(loop_x));exp(D1*(x(loop_x)-wth))])*c;
         %%
         I(loop_x,:) = real(Ix).';
      end
   else
      %% x  = 0
      I0 = V*diag([exp(D0*0);exp(D1*(0-wth))])*c;
      for loop_x=1:length(x)
         Ix = V*diag([exp(D0*x(loop_x));exp(D1*(x(loop_x)-wth))])*c;
         %%
         I(loop_x,:) = real(Ix).';
      end
   end
end

dtheta   = 2*pi/Ndir;
m0_inc   = sum(I0(incs))*dtheta;          %% integrate wrt theta over incident directions
I        = (.25*Hs_inc)^2/real(m0_inc)*I; %% renormalise I so inc waves have required Hs at x=0

if 0
   dtheta
   m0_p  = sum(I(1,incs),2)*dtheta
   Hs_p  = 4*sqrt(m0_p)
end
