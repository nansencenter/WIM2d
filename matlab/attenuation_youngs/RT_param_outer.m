function [damping,kice,kwtr,int_adm,NDprams,...
            alp_scat,modT,argR,argT] =...
               RT_param_outer(h,om,E,visc_rp,guess)

do_test  = 0;
if nargin==0
   if 0
      h  = 2;
      om = 2*pi/10;
   else
      h  = 10;
      om = 2*pi/1;
      E  = 10e9;
   end
   visc_rp  = 13;
   g        = 9.81;
   guess    = om^2/g;
   do_test  = 1;
end

prams    = NDphyspram(0);

if ~exist('E')
   E  = prams(1);
end
if ~exist('visc_rp')
   visc_rp  = 13;
end

g     = prams(2);
rhow  = prams(3);
rhoi  = prams(4);
nu    = prams(5);
rho   = rhoi/rhow;
%%
D        = E*h^3/12/(1-nu^2);
L        = ( D/rhow./om.^2 ).^.2;
Lc       = ( D/rhow/g ).^.25;
alp_nd   = om.^2/g.*L;
alp      = om.^2/g;
alp_nd   = alp.*L;
h_nd     = h./L;
zeta_nd  = rho*h_nd;
%tst_ah   = {alp_nd,h_nd}
%%
H_nd     = 4;%%inf depth
NDprams  = [alp_nd,zeta_nd,L,H_nd+0*om,Lc+0*om];

if ~exist('guess')
   guess = min(alp);
end

kice     = 0*om;
kwtr     = 0*om;
int_adm  = 0*om;
avc      = 0*om;
[dum,JJ] = sort(om);

for j_=1:length(om)
   j  = JJ(j_);

   %%get wavenumber for ice;
   varpi             = 1/alp_nd(j)-zeta_nd(j);
   [ki,BG2,avc(j)]   = fn_root_ice(varpi,H_nd,guess*L(j));
   kice(j)           = ki/L(j);
   guess             = kice(j);%%make guess the last root found;

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
end

%if 0
%   period   = 2*pi/om
%   h        = h_nd*L
%   Hw       = Hw_nd*L
%   ahl      = {alp_nd,h_nd,L}
%   kw2      = 2*pi*L/GEN_get_wtr_wavelength(period,Hw);
%   kkaa     = {kw,kw2,alp_nd,kw*tanh(kw*Hw_nd)}
%   kw/kw2
%   %%
%   kk_ice   = {ki,2*pi*L/GEN_get_ice_wavelength(h,period,H_nd*L,E)}
%   pause
%end


%%get viscous attenuation;
visc_rp_nd  = visc_rp/rhow./(om.*L);
damping_nd  = avc.*visc_rp_nd;
damping     = damping_nd./L;

%%Get other quantities via interpolation;
alp_scat = NaN*om;
modT     = NaN*om;
argR     = NaN*om;
argT     = NaN*om;

if (nargout>4) | do_test

   for j=1:length(om)
      if 0
         %% slower version using loading of files
         %% -needs to be run first to save chebyshev coefficients to datfiles;
         [alp_scat(j),modT(j),argR(j),argT(j)]  =...
            RT_param(alp_nd(j),h_nd(j),int_adm(j));
      elseif 0
         %% faster version with no loading of files
         %% (everything is hard-coded);
         [alp_scat(j),modT(j),argR(j),argT(j)]  =...
            RT_param_fast(alp_nd(j),h_nd(j),int_adm(j));
      else
         %% faster version with no loading of files
         %% (everything is hard-coded);
         [alp_scat(j),modT(j),argR(j),argT(j)]  =...
            RT_param_fast_v2(alp_nd(j),h_nd(j),int_adm(j));
      end
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
function [ki,BG2,avc]=fn_root_ice(del,H,guess)
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
%%
avc   = ki/denom;%% -1i*dk/d(visc_rp) - non-dimensional

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kw,BG1]=fn_root_wtr(del,H,guess)
%% finds the root of the ice dispersion relation nearest to 'guess'.

fac   = 0;
tol   = 1e-12;
k0    = guess;
dk    = NR_corr_term(k0,del,H,fac);
kw     = k0-dk;

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
