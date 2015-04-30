% Param = ModParam_def(Param)
%
% Setup of modal parameters and accuracy parameters for the solution
% method.

function Param = ModParam_def(Param,Ni,Nd,extra_pts,terms,th_res)

if ~exist('Ni','var'); Ni=1; end
if ~exist('Nd','var'); Nd=10; end
if ~exist('th_res','var'); th_res=100; end

%%% Number of vertical modes (travel + evan) for free-surf interactions
Param.Nint = Ni;

%%% Number of vertical modes (travel + evan)
Param.Ndtm = Nd;

%%% Number of evanescent horizontal modes
Param.Mev = 0;

%%% Number of extra points for irregular frequencies
if ~exist('extra_pts','var')
 Param.extra_pts = 0;
else
 if isempty(extra_pts)
  Param.extra_pts = 0;
 else
  Param.extra_pts = extra_pts;
 end
end

%%% Number of points for the intergration of the Green's function
Param.res_green = 100;

%%% Max number of terms used in Green's function
if ~exist('terms','var')
 Param.terms_green = 100;
else
 if isempty(terms)
  Param.terms_green = 100;
 else
  Param.terms_green = terms;
 end
end

%%% Limit for cutoff of Greens fn series
Param.cutoff_green = 1e-1;

%%% Tolerance on resonances

Param.tolres=1e-5; %7.5e-4; %0; %

%%% Resolution of angular spectrum (Boltzmann)
Param.th_res = th_res;

return