function [damping,kice,kwtr,int_adm,NDprams,...
            alp_scat,modT,argR,argT] =...
               RT_param_outer(om,h,ice_prams,guess)

do_test  = 0;
if nargin==0
   h     = 2;
   om    = 2*pi/10;
   H_nd  = 10;%%inf depth
   %%
   ice_prams   = fn_fill_iceprams();
   if 0
      %% test R-P drag
      ice_prams.drag_rp          = 13;
      ice_prams.viscoelastic_ws  = 0;
   elseif 1
      %% test W-S viscoelastic
      ice_prams.drag_rp          = 0;
      ice_prams.viscoelastic_ws  = 1;
   elseif 1
      %% no damping
      ice_prams.drag_rp          = 0;
      ice_prams.viscoelastic_ws  = 0;
   end
   guess       = om^2/ice_prams.g;
   do_test     = 1;
end

rho      = ice_prams.rhoice/ice_prams.rhowtr;
D        = ice_prams.young*h^3/12/(1-ice_prams.poisson^2);
L        = ( D/ice_prams.rhowtr./om.^2 ).^.2;
Lc       = ( D/ice_prams.rhowtr/ice_prams.g ).^.25;
alp_nd   = om.^2/ice_prams.g.*L;
alp      = om.^2/ice_prams.g;
alp_nd   = alp.*L;
h_nd     = h./L;
zeta_nd  = rho*h_nd;
%tst_ah   = {alp_nd,h_nd}
%%
if ~exist('H_nd')
   H_nd  = 4;%%inf depth
end
NDprams  = [alp_nd,zeta_nd,L,H_nd+0*om,Lc+0*om];

if ~exist('guess')
   guess = min(alp);
end

kice        = 0*om;
kwtr        = 0*om;
int_adm     = 0*om;
coeff_del   = 0*om;
coeff_D     = 0*om;
[dum,JJ]    = sort(om);

for j_=1:length(om)
   j  = JJ(j_);

   %%get wavenumber for ice;
   varpi = 1/alp_nd(j)-zeta_nd(j);
   [ki,BG2,coeff_del(j),coeff_D(j)]   =...
      fn_root_ice(varpi,H_nd,guess*L(j));

   kice(j)  = ki/L(j);
   guess    = kice(j);%%make guess the last root found;

   %%get wavenumber for water;
   varpi    = 1/alp_nd(j);
   Hw_nd    = H_nd+zeta_nd(j);
   [kw,BG1] = fn_root_wtr(varpi,Hw_nd,alp_nd(j));
   kwtr(j)  = kw/L(j);

   %tstBG1_  = {BG1,H_nd+zeta_nd,1/alp_nd, zeta_nd,1,1,ki,L}
   %tstBG2_  = {BG2,H_nd,1/alp_nd, zeta_nd,0,0,kw,L}
   %%get intrinsic admittance;
   %%|R|^2+int_adm*|T|^2=1
   int_adm(j)  = BG1/BG2;

   if 0
      period   = 2*pi/om(j)
      h
      wn_ice   = kice(j)
      wn_wtr   = kwtr(j)
   end
end

%% get attenuation due to Robinson-Palmer drag:
%% * dimensional eqn:
%% [ (D/rhow/om^2)*k^4+ g/om^2-rhoi*h/rhow-i*drag_rp/(om*rhow) ]*k*tanh(k*H) -1 = 0
%%    >> -i*om*drag_rp/(rhow*om^2) = -i*om*drag_rp/(rhow*om)
%% * non-dimensional eqn K=k*L:
%% [ D/(rhow*om^2*L^5)*K^4+ del-i*drag_rp/(L*om*rhow) ]*K*tanh(K*H/L) -1 = 0
%% d(del)   = -i*drag_rp/(L*om*rhow) = -i*drag_rp_nd
%% d(k)     = coeff_del*d(del)=-i*drag_rp_nd*coeff_del = i*damping_nd;
drag_rp_nd  = ice_prams.drag_rp/ice_prams.rhowtr./(om.*L);

%% W&S (still non-dim):
%% d(D)  = [-i*om*h^3*rhoi*viscoelastic_ws/6/(1+nu)]/[rhow*om^2*L^5]
%%       = [-i*om*h^3*rhoi*viscoelastic_ws/6/(1+nu)]/[rhow*om^2*L^5]
%%       = -i*h^3*rho*viscoelastic_ws/(6*(1+nu)*om*L^5)
%%       = -i*vews_nd
%% d(k)  = coeff_D*d(D)=-i*coeff_D*vews_nd = i*damping_nd;
vews_nd     = h^3*rho*ice_prams.viscoelastic_ws./(6*(1+ice_prams.poisson)*om.*L.^5);

%% add both effects together
damping_nd  = -coeff_del.*drag_rp_nd - coeff_D.*vews_nd;
damping     = damping_nd./L;

%%Get other quantities via interpolation;
alp_scat = NaN*om;
modT     = NaN*om;
argR     = NaN*om;
argT     = NaN*om;

if (nargout>5) | do_test

   for j=1:length(om)
      %% faster version with no loading of files
      %% (everything is hard-coded);
      [alp_scat(j),modT(j),argR(j),argT(j)]  =...
         RT_param_fast_v2(alp_nd(j),h_nd(j),int_adm(j));
   end

end

if do_test
   disp(['damping = ',num2str(damping, '%.9e')]);
   disp(['kice    = ',num2str(kice,    '%.9e')]);
   disp(['kwtr    = ',num2str(kwtr,    '%.9e')]);
   disp(['int_adm = ',num2str(int_adm, '%.9e')]);
   disp(' ');
   disp(['ac      = ',num2str(alp_scat,'%.9e')]);
   disp(['|T|     = ',num2str(modT,    '%.9e')]);
   disp(['Arg[R]  = ',num2str(argR,    '%.9e')]);
   disp(['Arg[T]  = ',num2str(argT,    '%.9e')]);

   c5    = 1-1i*vews_nd(1);
   del   = alp_nd(1)-zeta_nd(1);
   c1    = del - 1i*drag_rp_nd;
   Q     = roots([c5,0,0,0,c1,-1])/L(1);
   Qap   = kice(1)+1i*damping;
   %%
   [dummy,JJ]  = sort(abs(Q-Qap));
   tst         = [Q(JJ(1));Qap]
   %%
   K  = kice(1)*L(1);
   [tanh(K*H_nd),(K^4+alp_nd(1))*K*tanh(K*H_nd)-1]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function y = calc_res(Z2,hr,gamma)
%  %% y=calc_res(Z2,hr,gamma)=Res(1/f(K),gamma_n),
%  %% where gamma_n is a root of the dispersion relation
%  %% f=1/K/tanh(KH)-(Dr*K^4+lam-mr*mu);
%  %% Z2={lam,mu,H}.
%  alp_nd   = Z2{1};
%  zeta_nd  = Z2{2};
%  H        = Z2{3};
%  Dr       = hr^3;
%  mr       = hr;
%  %%
%  Gam   = Dr*gamma.^4+1/alp_nd-mr*zeta_nd;
%  Gampr = Gam+4*Dr*gamma.^4;
%  denom = H*(Gam.^2.*gamma.^2-1)+Gampr;
%  y     = -gamma./denom;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ki,BG2,coeff_del,coeff_D]=fn_root_ice(del,H,guess)
%% finds the root of the ice dispersion relation nearest to 'guess'.

fac   = 1;
tol   = 1e-12;
k0    = guess;
dk    = NR_corr_term(k0,del,H,fac);
ki     = k0-dk;

while abs(dk) > tol
   k0 = ki;
   dk = NR_corr_term(k0,del,H,fac);
   ki  = k0-dk;
end

[dk,Lam,Lampr] = NR_corr_term(ki,del,H,fac);
denom          = H*(Lam.^2.*ki.^2-1)+Lampr;
res            = -ki/denom;
BG2            = Lam^2*res;

%% perturbations to get R&P drag and W&S viscoelastic effect
%% F=Lam(k,D,del)*k*tanh(k*H)-1=0:
%%  F_k   = [ (Lam*k)_k+H*( (Lam*k)^2-1 ) ]/(Lam*k);
%%  F_del = Lam_del/Lam = 1/Lam;
%%  F_D   = Lam_D/Lam   = k^4/Lam;
%% dF/dD       = -F_D/F_k;
%% dF/d(del)   = -F_del/F_k;

coeff_del   = -ki/denom;  %% dk/d(del) - non-dimensional
coeff_D     = -ki^5/denom;%% dk/d(D)   - non-dimensional

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kw,BG1]=fn_root_wtr(del,H,guess)
%% finds the root of the ice dispersion relation nearest to 'guess'.

fac   = 0;
tol   = 1e-12;
k0    = guess;
dk    = NR_corr_term(k0,del,H,fac);
kw    = k0-dk;

while abs(dk) > tol
   k0       = kw;
   dk       = NR_corr_term(k0,del,H,fac);
   kw       = k0-dk;
   %cvg_tst  = {dk,kw},pause
end

[dk,Lam,Lampr] = NR_corr_term(kw,del,H,fac);
denom          = H*(Lam.^2.*kw.^2-1)+Lampr;
res            = -kw/denom;
BG1            = Lam^2*res;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dk,Lam,Lampr] = NR_corr_term(k,del,H,fac)
%% dk=f/f_k, where f has the same zeros as of the dispersion function, 
%% is the correction term in the Newton-Rhapson method for finding zeros in f.

Lam   = fac*k.^4+del;
Lampr = 5*fac*k.^4+del;
x     = 7.5;
if real(k*H)<=x
   f  = Lam.*k.*sinh(k*H)-cosh(k*H);
   df = Lam.*(k*H).*cosh(k*H)+(Lampr-H).*sinh(k*H);
else
   f  = Lam.*k.*tanh(k*H)-1;
   df = Lam.*k*H+(Lampr-H).*tanh(k*H);
end
dk = f./df;
