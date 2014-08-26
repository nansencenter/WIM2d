function [conc,thick,periods,dirs,S]  = init_cons(X,Y,LANDMASK,ice_props,wave_props)
%% init_cons.m
%% Author: Timothy Williams
%% Date:   20140826, 03:25:49 CEST
%% CALL: [conc,thick,periods,dirs,S]  = init_cons(X,Y,LANDMASK,ice_props,wave_props)

if nargin==0
   [X,Y,scuy,scvx,scp2,scp2i,LANDMASK] = get_grid();
   [ii,jj]                             = size(X);
   xm                                  = max(X);
   ym                                  = max(Y);

   %% ice_props = {[x0,x1],[c0,h0]};
   %% ice edge is straight line from top to bottom of grid:
   %% x  = x0+m*(y-y0)
   %% c0 = constant concentration of ice
   %% h0 = constant thickness of ice
   ice_props   = {[-70e3,-30e3],[.7,2]};

   %% wave_props={[Hs,Tp],[nfreqs,ndir]};
   wave_props  = {[2,10],[1,1]};
end

%%conc,thick
x0          = ice_props{1}(1);
x1          = ice_props{1}(2);
y0          = -ym;
y1          = ym;
m           = (x1-x0)/(y1-y0);
%%
jice        = find(X<(x0+m*(Y-Y0));
c0          = ice_props{2}(1);
h0          = ice_props{2}(2);
%%
conc        = 0*X;
thick       = 0*X;
conc(jice)  = c0;
thick(jice) = h;

%%waves
Hs    = wave_props{1}(1);
Tp    = wave_props{1}(2);
mwd   = wave_props{1}(3);
nw    = wave_props{2}(1);
ndir  = wave_props{2}(2);

%%
S  = zeros(ii,jj,nw,ndir);
if nw==1
   periods  = Tp;
   S0       = Hs^2;
else%%TODO
end

if ndir==1
   dirs        = mwd;
   theta_fac   = 1;
else%%TODO
end

S  = zeros(ii,jj,nw,ndir);
for w=1:nw
   for nth=1:ndir
      S1                      = S0(w)*theta_fac(nd)+0*X;
      S1(jice)                = 0;%%0 in ice initially
      S1(find(LANDMASK==1))   = 0;%%0 on land
      S(:,:,nw,nth)           = S1;
   end
end
