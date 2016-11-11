function out = floe_scaling_smooth(Dmax,prams,mom,inverse)
%% This function computes the average floe size within a grid cell as a
%% function of the maximum floe size using a pdf
%% P(d>D)   = (D_min/D)^fsd_exp, D<=Dmax
%% P(d>D)   = 0                , D>Dmax
%%
%% fsd_exp  = 2+log(fragility)/log(xi)
%%
%% We suggest to use Dmin >= 20. Below that value, there is no scattering
%% by the floes for periods larger than 6 s (it's probably viscous however).
%%
%% Dave  = average floe size
%% Dmax  = max floe size
%% prams = [structure] 
%%            xi: 2
%%     fragility: 0.900000000000000
%%          Dmin: 20
%% moment=1,2 -> <D> or <D^2>

DO_TEST  = 0;
prams0   = struct('f',.9,...
                  'xi',2,...
                  'Dmin',20);
if ~exist('prams')
   prams = prams0;
else
   if isempty(prams)
      prams = prams0;
   end
end

if ~exist('mom')
   mom   = 1;
end
if ~exist('inverse')
   inverse   = 0;
end

f        = prams.fragility;
xi       = prams.xi;
Dmin     = prams.Dmin;
fsd_exp  = 2+log(f)/log(xi);%%power law exponent: P(d>D)=(D_min/D)^fsd_exp;


n  = mom;%need to calc <D^n> or invert it to get Dmax
b  = n-fsd_exp;


Dthresh  = 200;

if inverse==0
   %%calculate <D^n> from Dmax
   out      = Dmax.^n;

   %% small floes
   jm       = find(Dmax<=Dmin);
   out(jm)  = Dmin.^n;%%assume uniform for small Dmax

   %% bigger floes
   jp       = find((Dmax>Dmin)&(Dmax<Dthresh));
   out(jp)  = do_int(Dmax(jp),Dmin,fsd_exp,n,0);
else
   %%calculate Dmax from <D^n>
   Dc       = Dmin^n;
   Dave     = Dmax;
   out      = Dave.^(1/n);
   Dthresh  = do_int(Dthresh,Dmin,fsd_exp,n,0);

   %% small floes
   jm       = find(Dmax<=Dc);
   out(jm)  = Dmin;

   %% bigger floes
   jp = find((Dave>Dc)&(Dave<Dthresh));
   Dp = Dave(jp);
   for j=1:length(jp)
      target   = Dp(j);
      interval = [target^(1/n),1000];
      Dp(j)  = GEN_findroot_bisection(@do_int,interval,Dmin,fsd_exp,n,target);
   end
   out(jp)  = Dp;
end

if DO_TEST
   if inverse==0
      y  = test_int(max(Dmax(:)),Dmin,fsd_exp,n)
   else
      x  = max(out(:));%Dmax after inversion
      j  = find(out==x);
      j  = j(1);
      [do_int(x,Dmin,fsd_exp,n,0),Dave(j)]
   end
end

%disp(out)
%disp(Dmax)
%disp(mom)
%disp(f)
%disp(xi)
%disp(Dmin)
%error('HEY!')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=test_int(Dmax,Dmin,fsd_exp,n)
b  = n-fsd_exp;
f  = inline('d.^(-1-fsd_exp+n)');
A  = 1/quad(@(d)f(d,fsd_exp,0),Dmin,Dmax);
y  = A*quad(@(d)f(d,fsd_exp,n),Dmin,Dmax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dave = do_int(Dmax,Dmin,fsd_exp,n,target)
A     = (fsd_exp*Dmin^fsd_exp*Dmax.^fsd_exp)./(Dmax.^fsd_exp-Dmin^(fsd_exp));%%
b     = n-fsd_exp;
Dave  = -(A/b).*(Dmin^b-Dmax.^b)-target;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
