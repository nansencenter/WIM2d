function [alp_scat,modT,argR,argT] = RT_param_fast_v2(alp_nd,h_nd,int_adm)

do_test  = 0;
if nargin==0
   do_test  = 1;
   LOW      = 1;
   OPT      = 1;
   int_adm  = 1.33;

   if LOW==1
      h_nd  = .1;
   else
      h_nd  = .23;
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   mc_alplin   = [-3.323529252398524,3.119943407349375];
   alp_lin3    = mc_alplin(2)+mc_alplin(1)*log(h_nd);
   y0_ll       = 40;
   dy_ll       = 120;
   n_ll        = 3;
   h1_ll       = .4;
   alp_lin4    = y0_ll+dy_ll*cos(h_nd/h1_ll*pi/2).^n_ll;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if OPT==1
      alp_nd   = 1e-4;
   elseif OPT==2
      alp_nd   = 1e-1;
   elseif OPT==3
      alp_nd   = 1.2;
   elseif OPT==4
      alp_nd   = 1.5+(alp_lin3-1.5)*.8;
   elseif OPT==5
      alp_nd   = (alp_lin3+alp_lin4)/2;
   end
end

if ~exist('int_adm')
   int_adm  = 1;
end

%prams = NDphyspram(0);
%E     = prams(1);
%g     = prams(2);
%rhow  = prams(3);
%rhoi  = prams(4);
%nu    = prams(5);
%rho   = rhoi/rhow;
%
%%%
%zeta_nd  = rho*h_nd;
LOW      = (h_nd<.2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% hard-coded limits for 5 regimes of alp_nd
%% (the two upper ones depend on h_nd);
alp_nd_lims  = [1e-6,0.005,0.3,1.5];

%% start of approximately linear regime
mc_alplin      = [-3.323529252398524,3.119943407349375];
alp_lin3       = mc_alplin(2)+mc_alplin(1)*log(h_nd);
alp_nd_lims(5) = alp_lin3;

%% end of calculatable results
y0_ll          = 40;
dy_ll          = 120;
n_ll           = 3;
h1_ll          = .4;
alp_lin4       = y0_ll+dy_ll*cos(h_nd/h1_ll*pi/2).^n_ll;
alp_nd_lims(6) = alp_lin4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if LOW==1
   %%limits of h_nd range
   h0    = 1e-2;
   h1    = .2;

   %%interp mode for h_nd
   LOG_H    = 1;%%vs log(h_nd)
   Amn_fxn  = 'Amn_fxn_L';

   %%Order of approximation
   Mh_vec  = [10,10,10,10,10];
   Ma_vec  = [10,10,10,10,3];

   %%mode of interpolation
   %%1 - {log(ac),angle(R),angle(T)}
   %%2 - {(ac),angle(R),angle(T)}
   %%3 - {real(R),imag(R),real(T),imag(T)}
   INTERP_MODE_vec   = [1,1,3,2,1];
   Nthings_vec       = [3,3,4,3,3];

   %%'x' variable:
   %% 0 vs (alp_nd)
   %% 1 vs log(alp_nd)
   LOG_A_vec   = [1,1,1,0,1];
else
   %%limits of h_nd range
   h0    = .2;
   h1    = .4;

   %%interp mode for h_nd
   LOG_H    = 0;%%vs (h_nd)
   Amn_fxn  = 'Amn_fxn_H';

   %%Order of approximation
   Mh_vec  = [10,10,10,10,10];
   Ma_vec  = [10,10,10,10,4];

   %%mode of interpolation
   %%1 - {log(ac),angle(R),angle(T)}
   %%2 - {(ac),angle(R),angle(T)}
   %%3 - {real(R),imag(R),real(T),imag(T)}
   INTERP_MODE_vec   = [1,1,3,2,1];
   Nthings_vec       = [3,3,4,3,3];

   %%'x' variable:
   %% 0 vs (alp_nd)
   %% 1 vs log(alp_nd)
   LOG_A_vec   = [1,1,1,0,1];
end

%%determine alp_nd regime:
EXTRAP_a = 0;
if alp_nd<alp_nd_lims(1)
   EXTRAP_a = 1;
   %alp_nd   = alp_nd_lims(1);
   OPT      = 1;
elseif alp_nd>=alp_nd_lims(6)
   EXTRAP_a = 1;
   %alp_nd   = alp_nd_lims(6);
   OPT      = 5;
else
   OPT   = find(alp_nd>=alp_nd_lims,1,'last');
end

a0 = alp_nd_lims(OPT);
a1 = alp_nd_lims(OPT+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mh          = Mh_vec(OPT);
Ma          = Ma_vec(OPT);
LOG_A       = LOG_A_vec(OPT);
INTERP_MODE = INTERP_MODE_vec(OPT);

if LOG_A==0
   tt_a2 = -1+2*(alp_nd-a0)/(a1-a0);
else
   l0    = log(a0);
   l1    = log(a1);
   tt_a2 = -1+2*(log(alp_nd)-l0)/(l1-l0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%test if alp_nd out-of-range
if EXTRAP_a*(OPT==1)
   dx    = 1e-4;
   x     = [-1;-1+dx];
   Tn_a2 = OP_interp_chebyshev(x,{Ma});
elseif EXTRAP_a*(OPT==5)
   dx    = -1e-4;
   x     = [1;1+dx];
   Tn_a2 = OP_interp_chebyshev(x,{Ma});
else
   Tn_a2 = OP_interp_chebyshev(tt_a2,{Ma});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if LOG_H==0
   tt_h2 = -1+2*(h_nd-h0)/(h1-h0);
else
   l0    = log(h0);
   l1    = log(h1);
   tt_h2 = -1+2*(log(h_nd)-l0)/(l1-l0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%test if h_nd out-of-range
EXTRAP_h = 0;
if tt_h2<-1
   EXTRAP_h = 1;
   dy       = 1e-4;
   y        = [-1;-1+dy];
   Tn_h2    = OP_interp_chebyshev(y,{Mh});
elseif tt_h2>1
   EXTRAP_h = 1;
   dy       = -1e-4;
   y        = [1;1+dy];
   Tn_h2    = OP_interp_chebyshev(y,{Mh});
else
   Tn_h2 = OP_interp_chebyshev(tt_h2,{Mh});
end
% EXTRAP_a,EXTRAP_h
% {tt_a2,tt_h2}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get hard-coded chebyshev coefficients;
func0 = [Amn_fxn,num2str(OPT)];
tmp   = RTparam_hardcoded_v2(func0,Nthings_vec(OPT));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the interpolation;
if INTERP_MODE==1
   %%log(alp_scat) & arg(R,T)
   Amn   = zeros(Ma+1,Mh+1);
   for r=1:3
      Amn(:)   = tmp(:,r);
      Z{r}     = Tn_a2*Amn*Tn_h2.';
      %Z{r}

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%Extrapolate if out-of-range
      if EXTRAP_a==1
         z0    = Z{r}(1,:);
         m     = (Z{r}(2,:)-z0)/dx;
         Z{r}  = z0+m*(tt_a2-x(1));
         %Z{r}
      end

      if EXTRAP_h==1
         z0    = Z{r}(1);
         m     = 0;%(Z{r}(2)-z0)/dy;
         Z{r}  = z0+m*(tt_h2-y(1));
         %Z{r}
      end
      %pause
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   end

   alp_scat = exp(Z{1});
   argR     = Z{2};
   argT     = Z{3};
   modT     = sqrt(exp(-alp_scat/2)/int_adm);
end

if INTERP_MODE==2
   %%(alp_scat) & arg(R,T)
   Amn   = zeros(Ma+1,Mh+1);
   for r=1:3
      Amn(:)   = tmp(:,r);
      Z{r}     = Tn_a2*Amn*Tn_h2.';

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%Extrapolate if out-of-range
      if EXTRAP_a==1
         z0    = Z{r}(1,:);
         m     = (Z{r}(2,:)-z0)/dx;
         Z{r}  = z0+m*(tt_a2-x(1));
      end

      if EXTRAP_h==1
         z0    = Z{r}(1);
         m     = (Z{r}(2)-z0)/dy;
         Z{r}  = z0+m*(tt_h2-y(1));
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   end

   alp_scat = Z{1};
   argR     = Z{2};
   argT     = Z{3};
   modT     = sqrt(exp(-alp_scat/2)/int_adm);
end

if INTERP_MODE==3
   %%(R,T)
   Amn   = zeros(Ma+1,Mh+1);
   for r=1:4
      Amn(:)   = tmp(:,r);
      Z{r}     = Tn_a2*Amn*Tn_h2.';

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%Extrapolate if out-of-range
      if EXTRAP_a==1
         z0    = Z{r}(1,:);
         m     = (Z{r}(2,:)-z0)/dx;
         Z{r}  = z0+m*(tt_a2-x(1));
      end

      if EXTRAP_h==1
         z0    = Z{r}(1);
         m     = (Z{r}(2)-z0)/dy;
         Z{r}  = z0+m*(tt_h2-y(1));
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   end
   %%
   R        = Z{1}+1i*Z{2};
   T        = Z{3}+1i*Z{4};
   argR     = angle(R);
   argT     = angle(T);
   modT     = abs(T);
   alp_scat = -2*log(1-abs(R).^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if do_test
   disp(' ');
   disp('Inputs to RT_param_fast_v2:');
   disp(['LOW  =  ',num2str(LOW)]);
   disp(['OPT  =  ',num2str(OPT)]);
   disp(['(alp_nd,h_nd) =  (',num2str(alp_nd,'%.9f'),',',num2str(h_nd,'%.9f'),')']);
   disp(['(t_a,t_h)     =  (',num2str(tt_a2,'%.9f'),',',num2str(tt_h2,'%.9f'),')']);
   disp(' ');
   %%
   disp('Outputs from RT_param_fast_v2:');
   disp(['ac      = ',num2str(alp_scat,'%.9e')]);
   disp(['|T|     = ',num2str(modT,    '%.9e')]);
   disp(['Arg[R]  = ',num2str(argR,    '%.9e')]);
   disp(['Arg[T]  = ',num2str(argT,    '%.9e')]);
   disp(' ');
end

return;
