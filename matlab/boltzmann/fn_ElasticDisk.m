% function fn_ElasticDisk(fortyp, forval, Param, GeomDisk, ...
%   bed, th_vec, RIGID, SURGE, COMM, PLT)
%
% calculate linear scattering of monochromatic wave by floating
% thin-elastic disk using eigenfunction matching method
% see Peter et al (IJOPE, 2003) or Bennetts & Williams (JFM, 2010)
%
% nb. disk centred at origin
%
% nb. defined: Phi = (g/i\omega)*phi(x,y,z)*exp(-i*omega*t) velocity potential
%              xi  = eta(x,y)*exp(-i*omega*t)  surface elevation
%     gives dynamic surface condition: phi(z=0) = eta
%
% INPUTS:
%
% fortyp   = 'freq' (in Hz) or 'wlength' (2*pi/k0) or 'waveno' (k0)
% forval   = associated value
% Param    = see ParamDef and ModParam_def(Param,1,N,0,0)
%         field: g (gravity), rho_0 (water density), rho (disk density)
%                nu (Poisson's ratio), E (Young's mod), 
%                D (flexural rigid), beta (scaled flex rigid),
%                floe_diam (floe diameter), draft, thickness,
%                bed  (depth of fluid)
%                N (vertical modes trav + ev)
% outputs  = see below
%
% PARAMETERS
%
% th_vec   = discrete angular spectrum
% ModMass  = mass of submerged portion of body / pi / Rad
% spring   = springing coeficient
% damp     = damping coefficient
% DimG_vec = [DimG, Gegenbauer yes=1/no=0]
%            DimG = The number of Gegenbauer Modes used (>0)
%            Gegenbauer yes/no = use the Gegenbauer polys at vertical interface (0=no)
% al       =  disk mass
% be       =  flexural rigid
% s_Modes  = The limt of the Fourier expansion in horizontal plane (any nat num)
%            total number = 2*s_Modes+1
% C_mat    = matrix that controls motion in ice
% kk,k0    = wvnos in ice (kk) and water (k0)
% wt,wt0   = weight (normalising fn) attached to each vert mode
% mu_0,mu_1= apx complex roots (Dimension dependent)
%
% FLAGS:
%
% INC   = incident wave type
%         1 = plane wave, unit amplitude, travelling in x-dir
% SURGE = include surge motion
% RIGID = rigid disk (inf) or elastic disk (rigidity=10^RIGID)
% FOU   = calculate Fourier basis rep of disk motion XXX DEFUNCT XXX
% RAO   = calculate the response amplitude operators XXX DEFUNCT XXX
% PLT   = plot a contour of constant angle
% COMM  = comments on/off (1/0)
% BESNM = normalisation used for Bessel fns
% RAO_TYP = normalise RAOs by intended value (0) or actual value (1)
%           (difference due to series expansion in azimuthal direction)
%
% OUTPUTS:
%
% 'Energy' = [E0, E] where
%            E0 = int_{0}^{2\pi} E dtheta where ...
%            E  = energy radiated in unit time per unit angle 
%                (Meylan & Squire, JGR `96; Meylan et al JGR `97)
% 'En-Fou' = Fourier basis representation of energy
% 'DTM'    = diffraction transfer matrix
% 'RAOs'   = [Heave, Pitch, Surge]
%
% L Bennetts May 2013 / Adelaide
%
% REVISION HISTORY:
% - modified from directory Arb_Floe_Pool/fn_Circ_Floe_Nov09 (Otago, 2009)

%%
function out = fn_ElasticDisk(fortyp, forval, Param, outputs, ...
 th_vec, RIGID, SURGE, COMM, PLT, col)

Tol_vec(1) = 1e-16; % - Real root error - %
Tol_vec(2) = 1e-16; % - Imag root error (init) - %
Tol_vec(3) = 1e-4;  % - Tol on const soln - %
Tol_vec(4) = 1e-1;  % - Az_Dim tol - %
Tol_vec(5) = 1e-3;  % - Energy error - %

if ~exist('COMM','var'); COMM=1; end

if ~exist('INC','var'); INC=1; end

if ~exist('SURGE','var'); SURGE=0; end
% if ~exist('RIGID','var'); RIGID=5; end

if ~exist('RAO_TYP','var'); RAO_TYP=0; end

% if ~exist('FOU','var'); FOU=0; end
% if ~exist('RAO','var'); RAO=1; end

if ~exist('BESNM','var'); BESNM=2; end

if ~exist('PLT','var'); PLT=0; end

if ~exist('N','var'); N=1; end

if ~exist('outputs','var'); outputs='Energy'; end

if ~exist('col','var'); col='r'; end

%% Set up problem

if ~exist('Param','var'); Param = ParamDef_Default;
 Param = ModParam_def(Param,1,N,0,0); end

radius=Param.floe_diam/2;
draught = Param.draft;
Vert_Modes = Param.Ndtm;
bed = Param.bed;

if ~exist('fortyp','var'); fortyp = 'freq'; end
if ~exist('forval','var'); forval = 1/10; end

if ~exist('th_vec','var'); th_vec=linspace(-pi,pi,101); end

if ~exist('DimG_vec','var'); DimG_vec=[Vert_Modes,0]; end

Forcing = Force_def(Param.g,bed,fortyp,forval);

parameter_vector = [Param.rho_0, Param.rho, Param.nu, Param.E, Param.g];

if COMM
 cprintf(0.4*[1,1,1],'>>> Elastic disk problem\n');
 cprintf(0.4*[1,1,1],['>>> Radius = ' num2str(radius) '\n'])
 cprintf(0.4*[1,1,1],['>>> Draught = ' num2str(draught) '\n'])
 if SURGE; cprintf(0.4*[1,1,1],'>>> with surge\n'); 
 else cprintf(0.4*[1,1,1],'>>> without surge\n'); end
 cprintf(0.4*[1,1,1],['>>> Period = ' num2str(1/Forcing.f) '\n'])
 cprintf(0.4*[1,1,1],['>>> wavenumber/length = ',...
  num2str(2*pi/Forcing.lam0) '/' num2str(Forcing.lam0) '\n'])
 cprintf(0.4*[1,1,1],['>>> ' int2str(Vert_Modes) ' vertical mode(s)\n'])
end

% Normalise

Norm_fac = 1;

thickness = Param.thickness/Norm_fac; draught = draught/Norm_fac;
radius = radius/Norm_fac; bed = bed/Norm_fac; freq = 2*pi*Forcing.f/sqrt(Norm_fac);

al = draught;
be = Param.beta;

parameter_vector = [parameter_vector, ...
 (freq^2/Param.g)/Norm_fac, bed, draught, thickness, bed-draught, al, be];

%-------------------------------------------------------------------------%
% Surge
%
% Kinematic condition
%
% phi = kappa*u*cos(theta)
%
% Dynamic condition
%
% (S-kappa*M-i*A)*u =
%     rho_0*R*int_{-pi}^{pi}int_{-d}^{0}phi(R,theta,z)*cos(theta) dz dtheta
%   = pi*R*rho_0*int_{-d}^{0}{phi_{-1}(R,z)+phi_{1}(R,z)}dz

if SURGE
 modMass = radius*draught;
 damp = 0;
 spring = 0;
 modDamp = 2*Forcing.f*damp/Param.g/Param.rho_0;
 modSpring = spring/Param.rho_0/pi;
 clear damp spring
 kappa = parameter_vector(6);
end
%-------------------------------------------------------------------------%

%% Define matricies

kk = zeros(Vert_Modes,1); k0 = zeros(Vert_Modes,1);
wt = zeros(Vert_Modes,1); wt_0 = zeros(Vert_Modes,1);

k0(1)   = fn_RealRoot([bed*parameter_vector(6)], ...
          'fn_ReDispRel_water', 'fn_UppLimReal_water', 1e-16)/bed; 
wt_0(1) = weight_0_PWC(bed, k0(1));

for loop_Dim = 2:Vert_Modes
 k0(loop_Dim)   = 1i*fn_ImagRoot_water(loop_Dim-1, bed, ...
                                       parameter_vector(6), 1e-16); 
 wt_0(loop_Dim) = weight_0_PWC(bed, k0(loop_Dim));
end

if isfield(Param,'azi')
 s_Modes=Param.azi;
else
 Az_Dim_vec = 2*Tol_vec(4)*ones(3,1);
 Az_Dim = 0;
 count = 0;
 while max(Az_Dim_vec)>Tol_vec(4)
  Az_Dim_vec(count+1) = abs(besselj(Az_Dim,k0(1)*radius));
  Az_Dim=Az_Dim+1; count = count+1; count = mod(count,length(Az_Dim_vec));
 end

 clear Az_Dim_vec count

 Az_Dim=Az_Dim-2; s_Modes = Az_Dim; clear Az_Dim
end

if COMM
 cprintf(0.4*[1,1,1],['>>> ' int2str(s_Modes) ' azimuthal mode(s)\n'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELASTIC DISK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk(1) = fn_RealRoot([1-parameter_vector(6)*al, be/((bed-draught)^4),...
 (bed-draught)*parameter_vector(6)],...
 'fn_ReDispRel_ice', 'fn_UppLimReal_ice', 1e-16)/(bed-draught);
wt(1) = weight_PWC(parameter_vector, kk(1));

for loop_Dim = 2:Vert_Modes
 kk(loop_Dim) = 1i*fn_ImagRoot_ice(loop_Dim-1, ...
  al, be, bed-draught, parameter_vector(6), 1e-16*[1,1]);
 wt(loop_Dim) = weight_PWC(parameter_vector, kk(loop_Dim));
end

if COMM
 cprintf(0.4*[1,1,1],['>>> wavenumber/length (plate) = ',...
  num2str(kk(1)) '/' num2str(2*pi/kk(1)) '\n'])
end

%  for loop_Dim=1:Vert_Modes
%   %  if or(loop_Dim==1,~RIGID)
%   %   kk(loop_Dim) = GetRootsMMA_PWC(parameter_vector, loop_Dim, Tol_vec);
%   %  else
%   %   kk(loop_Dim) = (loop_Dim-1)*1i*pi/(bed-draught);
%   %  end
%   kk(loop_Dim) = GetRootsMMA_PWC(parameter_vector, loop_Dim, Tol_vec);
%   wt(loop_Dim) = weight_PWC(parameter_vector, kk(loop_Dim));
%  end

[mu_0, mu_1] = mu_new_PWC(parameter_vector, Vert_Modes, kk, wt);

[mat_A0, mat_A, mat_W0, mat_W] = fn_JumpMats(parameter_vector,...
 Vert_Modes, DimG_vec(1), k0, kk, wt_0, wt,DimG_vec(2));

mat_A0 = diag(mat_A0);

if SURGE
 % for use in kinematic condition
 % int_{-d}^{0} xi_m(z) dz
 % where xi(z) = wt_0*cosh{k*(z+h)}
 mat_Q0 = fn_JumpMatSg(parameter_vector, Vert_Modes, k0, wt_0);
end

C_mat = mat_C_dr_PWC(parameter_vector, length(kk), kk, mu_0, mu_1,...
 wt, parameter_vector(12));

mat_W0_T = conj(mat_W0'); mat_W_T = conj(mat_W');

%%%%%%%%%%%%%%%%
%% Solve via EMM
%%%%%%%%%%%%%%%%

Psi = zeros(Vert_Modes, 2*s_Modes+1); Phi = Psi;
if SURGE; Psi_sg = Psi; Phi_sg = Psi; end

% -- Version 1 -- %

Avec = zeros(Vert_Modes, 2*s_Modes+1); Bvec = Avec;
if SURGE; Avec_sg = Avec; Bvec_sg = Bvec; end
Ivec = zeros(Vert_Modes, 2*s_Modes+1);

v_01 = C_mat(1:Vert_Modes, Vert_Modes+1:Vert_Modes+2);

% -- Version 1.1 -- %

for loop_Az=-s_Modes:s_Modes
 
 %% Incident wave
 
 if INC==1
  % Set unit amplitude incident wave:
  % (1) divide by wt_0(1) to remove the weight 
  % (2) divide by cosh(k0(1)*bed) for correct scaling
  % nb. usually wt_0=sech(k0(1)*bed), so that (1)+(2) cancel
  Inc = (1i^loop_Az)*exp(abs(imag(k0(1)*radius))).*eye(Vert_Modes,1)...
   /wt_0(1)/cosh(k0(1)*bed);
  Ivec(:,loop_Az+s_Modes+1) = (1i^loop_Az)*eye(Vert_Modes,1);
 else
  cprintf('green','code another incident wave\n')
 end
 
 %% Bend & Shear
 % NB. scale besselj by exp(-abs(imag(z))) & besselh by exp(-iz)
 
 if BESNM==1
  [bend_mu0, shear_mu0] = ...
   BendShear3dAxi_fn( parameter_vector, [1, loop_Az, 1], mu_0*radius );
  [bend_mu1, shear_mu1] = ...
   BendShear3dAxi_fn( parameter_vector, [1, loop_Az, 1], mu_1*radius );
  [bend_K, shear_K] = ...
   BendShear3dAxi_fn( parameter_vector, [1, loop_Az, 1], diag([kk*radius]) );
 elseif BESNM==2
  [bend_mu0, shear_mu0] = ...
   fn_BendShear3dAxi( parameter_vector, [1, loop_Az], mu_0*radius );
  [bend_mu1, shear_mu1] = ...
   fn_BendShear3dAxi( parameter_vector, [1, loop_Az], mu_1*radius );
  [bend_K, shear_K] = ...
   fn_BendShear3dAxi( parameter_vector, [1, loop_Az], kk*radius );
  bend_K = diag(bend_K); shear_K = diag(shear_K);
 end
 
 % Bessel fns
 
 Amp_mu(1,1:Vert_Modes) = -C_mat(Vert_Modes+1, 1:Vert_Modes)*...
  (bend_mu1*shear_K - shear_mu1*bend_K);
 Amp_mu(2,1:Vert_Modes) = C_mat(Vert_Modes+1, 1:Vert_Modes)*...
  (bend_mu0*shear_K - shear_mu0*bend_K);
 
 Amp_mu = -Amp_mu./(bend_mu0*shear_mu1 - bend_mu1*shear_mu0);
 
 %Amp_mu_Mat(:,:,loop_Az+s_Modes+1) = diag(exp(-imag([mu_0 mu_1])*radius))*Amp_mu;
 Amp_mu_Mat(:,:,loop_Az+s_Modes+1) = Amp_mu;
 
 %% Components of full solution: Bessel & Hankel fns
 
 if BESNM==1
  J_ice = besselj(loop_Az, kk.*radius, 1);
  %Jz_ice = Bessel_dz(@besselj, loop_Az, kk.*radius, 1);
  invars = struct('name',{'BesselFn'; 'N'; 'z'; 'scale'}, ...
   'value',{@besselj; loop_Az; kk.*radius; 1 });
  Jz_ice = fn_Bessel_dz(invars);
  
  Jmu_ice = besselj(loop_Az, [ mu_0 mu_1 ].*radius, 1);
  %Jzmu_ice = Bessel_dz(@besselj, loop_Az, [ mu_0 mu_1 ].*radius, 1);
  invars = struct('name',{'BesselFn'; 'N'; 'z'; 'scale'}, ...
   'value',{@besselj; loop_Az; [ mu_0 mu_1 ].*radius; 1 });
  Jzmu_ice = fn_Bessel_dz(invars);
 elseif BESNM==2
  J_ice = ones(Vert_Modes,1);
  invars = struct('name',{'BesselFn'; 'N'; 'z'}, ...
   'value',{@besselj; loop_Az; kk.*radius});
  Jz_ice = fn_Bessel_dz(invars)./besselj(loop_Az, kk.*radius);
  
  Jmu_ice = [1;1];
  invars = struct('name',{'BesselFn'; 'N'; 'z'}, ...
   'value',{@besselj; loop_Az; [ mu_0 mu_1 ].*radius });
  Jzmu_ice = fn_Bessel_dz(invars)./besselj(loop_Az, [ mu_0 mu_1 ].*radius);
 end
 
 if 1
  J_water = besselj(loop_Az, k0.*radius, 1);
  %Jz_water = Bessel_dz(@besselj, loop_Az, k0.*radius, 1);
  invars = struct('name',{'BesselFn'; 'N'; 'z'; 'scale'}, ...
   'value',{@besselj; loop_Az; k0.*radius; 1 });
  Jz_water = fn_Bessel_dz(invars);
 elseif BESNM==2
  J_water = ones(Vert_Modes,1);
  invars = struct('name',{'BesselFn'; 'N'; 'z'; 'scale'}, ...
   'value',{@besselj; loop_Az; k0.*radius; 1 });
  Jz_water = fn_Bessel_dz(invars)./besselj(loop_Az, k0.*radius, 1);
 end
 
 if BESNM==1
  H_water = besselh(loop_Az, 1, k0.*radius, 1);
  Hz_water = Bessel_dz(@besselh, loop_Az, k0.*radius, 1, 1);
  %  invars = struct('name',{'BesselFn'; 'N'; 'z'; 'scale'; 'hank'}, ...
  %   'value',{@besselh; loop_Az; k0.*radius; 1; 1 });
  %  Hz_water = fn_Bessel_dz(invars);
 elseif BESNM==2
  H_water(1) = exp(-1i*k0(1)*radius)*...
   besselh(loop_Az, 1, k0(1).*radius);
  invars = struct('name',{'BesselFn'; 'N'; 'z'; 'hank'}, ...
   'value',{@besselh; loop_Az; k0(1).*radius; 1 });
  Hz_water(1) = exp(-1i*k0(1)*radius)*...
   fn_Bessel_dz(invars);
  if Vert_Modes>1
   invars = struct('name',{'BesselFn'; 'N'; 'z'; 'hank'}, ...
    'value',{@besselh; loop_Az; k0(2:Vert_Modes).*radius; 1 });
   H_water(2:Vert_Modes,1) = ones(Vert_Modes-1,1);
   Hz_water(2:Vert_Modes,1) = ...
    fn_Bessel_dz(invars)./...
    besselh(loop_Az, 1, k0(2:Vert_Modes).*radius);
  end
 end
 
 %% Set up soln in water:
 
 Phi_In = diag(J_water);
 Phiz_In = diag(k0.*Jz_water);
 
 Phi_Out = diag(H_water);
 Phiz_Out = diag(k0.*Hz_water);
 
 %% Set up soln from ice:
 
 % - fn psi = %
 
 Psi_S = diag(J_ice) + v_01*diag(Jmu_ice)*Amp_mu;
 Psiz_S = diag(kk.*Jz_ice) + v_01*diag([mu_0, mu_1].*Jzmu_ice)*Amp_mu;
 
 w1 = 1:Vert_Modes; w2 = Vert_Modes+1:Vert_Modes+2;
 
 % w_S = C_mat(Vert_Modes+1,w1).*J_ice.' +...
 %     C_mat(Vert_Modes+1,w2)*diag(Jmu_ice)*Amp_mu;
 %
 % w_In = C_mat(Vert_Modes+1,w1).*J_ice.' +...
 %     C_mat(Vert_Modes+1,w2)*diag(Hmu_ice)*Amp_mu_J;
 % w_Out = C_mat(Vert_Modes+1,w1).*H_ice.' +...
 %     C_mat(Vert_Modes+1,w2)*diag(Hmu_ice)*Amp_mu_H;
 % wht_In = C_mat(Vert_Modes+2,w1).*J_ice.' +...
 %     C_mat(Vert_Modes+2,w2)*diag(Hmu_ice)*Amp_mu_J;
 % wht_Out = C_mat(Vert_Modes+2,w1).*H_ice.' +...
 %     C_mat(Vert_Modes+2,w2)*diag(Hmu_ice)*Amp_mu_H;
 %
 % wz_In = C_mat(Vert_Modes+1,w1).*(kk.*Jz_ice).' +...
 %     C_mat(Vert_Modes+1,w2)*diag([mu_0, mu_1].*Hzmu_ice)*Amp_mu_J;
 % wz_Out = C_mat(Vert_Modes+1,w1).*(kk.*Hz_ice).' +...
 %     C_mat(Vert_Modes+1,w2)*diag([mu_0, mu_1].*Hzmu_ice)*Amp_mu_H;
 % wzht_In = C_mat(Vert_Modes+2,w1).*(kk.*Jz_ice).' +...
 %     C_mat(Vert_Modes+2,w2)*diag([mu_0, mu_1].*Hzmu_ice)*Amp_mu_J;
 % wzht_Out = C_mat(Vert_Modes+2,w1).*(kk.*Hz_ice).' +...
 %     C_mat(Vert_Modes+2,w2)*diag([mu_0, mu_1].*Hzmu_ice)*Amp_mu_H;
 
 %% Solve the system
 
 %%% Cty of radial velocity & submerged edge condition
 
 % following are Dim x DimG:
 
 MatB_u0 = (Phiz_Out\(diag(mat_A0)\mat_W0));
 MatB_u = (Psiz_S\(mat_A\mat_W));
 
 % following is Dim x 1:
 
 if and(SURGE,abs(loop_Az)==1)
  MatB_sg0 = 0.5*kappa*(Phiz_Out\(diag(mat_A0)\mat_Q0));
 end
 
 % following is Dim x Dim:
 
 MatB_I0 = -Phiz_Out\Phiz_In;
 
 %%% Cty of pressure
 
 % following are DimG x Dim:
 
 MatA_I0 = mat_W0_T*(Phi_In + Phi_Out*MatB_I0);
 
 % following are DimG x DimG:
 
 MatA_u0 = mat_W0_T*Phi_Out*MatB_u0;
 MatA_u = mat_W_T*Psi_S*MatB_u;
 
 % following are DimG x 1:
 
 if and(SURGE,abs(loop_Az)==1)
  MatA_sg0 = mat_W0_T*Phi_Out*MatB_sg0;
 end
 
 % find u in terms of the inc amps (DimG x Dim)
 
 MatU_I0 = (MatA_u - MatA_u0)\MatA_I0;
 
 % find u in terms of surge amplitude (DimG x 1)
 if and(SURGE,abs(loop_Az)==1)
  MatU_sg0 = (MatA_u - MatA_u0)\MatA_sg0;
 end
 
 u_vec = MatU_I0*Inc;
 
 %%% Solve
 
 % convert into scattered amps
 % 1. in the ice (floe)
 Avec(:, s_Modes+loop_Az+1) = MatB_u*u_vec;
 if and(SURGE,abs(loop_Az)==1)
  Avec_sg(:, s_Modes+loop_Az+1) = MatB_u*MatU_sg0;
 end
 % 2. in the (surrounding) water
 Bvec(:, s_Modes+loop_Az+1) = MatB_u0*u_vec + MatB_I0*Inc;
 if and(SURGE,abs(loop_Az)==1)
  Bvec_sg(:, s_Modes+loop_Az+1) = MatB_u0*MatU_sg0 + MatB_sg0;
 end
 
 % Put into vel pots and disp
 
 Phi(:, s_Modes+loop_Az+1) = Phi_In*Inc + Phi_Out*Bvec(:, s_Modes+loop_Az+1);
 Psi(:, s_Modes+loop_Az+1) = Psi_S*Avec(:, s_Modes+loop_Az+1);
 
 if and(SURGE,abs(loop_Az)==1)
  Phi_sg(:, s_Modes+loop_Az+1) = Phi_Out*Bvec_sg(:, s_Modes+loop_Az+1);
  Psi_sg(:, s_Modes+loop_Az+1) = Psi_S*Avec_sg(:, s_Modes+loop_Az+1);
 end
 
 % w(:, s_Modes+loop_Az+1) = w_S*Avec(:, s_Modes+loop_Az+1);
 % wn(:, s_Modes+loop_Az+1) = wz_In*Inc + wz_Out*Bvec(:, s_Modes+loop_Az+1);
 %
 % what(:, s_Modes+loop_Az+1) = wht_In*Inc + wht_Out*Bvec(:, s_Modes+loop_Az+1);
 % wnhat(:, s_Modes+loop_Az+1) = wzht_In*Inc + wzht_Out*Bvec(:, s_Modes+loop_Az+1);
 
end % end loop_Az

if and(SURGE,s_Modes>0)
 X_sg= ...
  ( mat_Q0.'*(Phi(:,s_Modes)+Phi(:,s_Modes+2)) ) / ...
  ( -(modSpring-kappa*modMass-1i*modDamp) - ...
  (mat_Q0.'*(Phi_sg(:,s_Modes)+Phi_sg(:,s_Modes+2))) );
 Avec(:, s_Modes)  =Avec(:, s_Modes)  +X_sg*Avec_sg(:, s_Modes);
 Avec(:, s_Modes+2)=Avec(:, s_Modes+2)+X_sg*Avec_sg(:, s_Modes+2);
 Bvec(:, s_Modes)  =Bvec(:, s_Modes)  +X_sg*Bvec_sg(:, s_Modes);
 Bvec(:, s_Modes+2)=Bvec(:, s_Modes+2)+X_sg*Bvec_sg(:, s_Modes+2);
elseif SURGE
 X_sg=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END DISK SOLUTION MAIN CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_str = ' ''dummy'' '; out_val = ' 0 ';

%% Energy check

if INC~=1; cprintf('green','Energy checks will not work!!!\n'); end

% In far-field:
% \phi=\sqrt{2/pi/k0/r}exp(i*k0*r)\sum_{n=-N}^{N}Bvec(1,n)*exp(i*n*(theta-\pi/2))
%     =\sqrt{2/pi/k0/r}exp(i*k0*r)*A(theta)
% where
% A(theta) = \sum_{n=-N}^{N}Bvec(1,n)*exp(i*n*(theta-\pi/2))

% Amplitudes of angular components of prop wave:
Bvec0 = Bvec(1,:);
% Eliminate in-built scaling of Hankel functions (see: help besselh)
Bvec0 = exp(-1i*k0(1)*radius)*Bvec0;

% Useful coefficients:
Ip = exp(1i*[-s_Modes:s_Modes]*pi/2);

% 1. E0 = int_{0}^{2\pi}|A(theta)|^2 dtheta
%    using orthogonality = ...
E0 = sum(abs(Bvec0.*conj(Ip)).^2);

% 2. A(in_angle) = A(0) = \sum_{n=-N}^{N}Bvec(1,n)*exp(n*\pi/2i))
%                       = \sum_{n=-N}^{N}Bvec(1,n)*Ip(n)
% => Real{A(0)} = E1 = ...
E1 = real(sum(conj(Ip).*Bvec0))/wt_0(1)/cosh(k0(1)*bed);

% int_{0}^{2\pi}|A(theta)|^2 dtheta = -Real{A(inc_angle)}
if abs(E0+E1) > Tol_vec(5)
 cprintf('red',['>>> energy error (' num2str(E0+E1) '): ' ...
  'period=' num2str(1/Forcing.f) '\n'])
end

%%% For circular geometry following should hold:

if max( abs( abs(Bvec0.*conj(Ip)).^2 + ...
  real(conj(Ip).*Bvec0)/wt_0(1)/cosh(k0(1)*bed) ) ) > Tol_vec(5)
 cprintf('red',['>>> energy error: ' ...
  'period=' num2str(1/Forcing.f) '\n'])
end

%%% Output the diffraction transfer matrix/vector

if strfind(outputs,'DTM')

 DTM = Bvec0.*conj(Ip)*wt_0(1)*cosh(k0(1)*bed);

 out_str = [out_str '; ''DTM'' '];
 out_val = [out_val '; DTM '];
 
end

if strfind(outputs,'wavenumber')
 out_str = [out_str '; ''wavenumber'' '];
 out_val = [out_val '; k0(1) '];
end

%% Energy

if strfind(outputs,'Energy')
 
 % FF_Amp = Far-field amplitude (function of theta)
 % E ~ |A|^2 -> scattered energy (function of energy)
 % E0 = int_{0}^{2\pi} E dtheta
 
 % nb. incang not required as axisymm problem
 % incang = 0;
 % Normalise inc wave energy to unity with incamp = sqrt(2/rho0*g)
 % incamp = sqrt(2/Param.rho_0*Param.g);
 % but this simply cancels with below, so leave !!!!
 % (CC Mei book Chap 1)
 
 % A(theta) = \sum_{n=-N}^{N}Bvec(1,n)*exp(i*n*(theta-\pi/2)) as above
 FF_Amp = (Bvec0.*conj(Ip))*exp(1i*[-s_Modes:s_Modes].'*th_vec);
 
 % Nb. Include the asymptotic factor \sqrt{2/k/pi}
 E = 2*(abs(FF_Amp).^2)/k0(1)/pi;
 
 out_str = [out_str '; ''E'' '];
 out_val = [out_val '; E '];
 
 E0 = 4*sum(abs(Bvec0.*conj(Ip)).^2)/k0(1);
 
 out_str = [out_str '; ''E0'' '];
 out_val = [out_val '; E0 '];
 
end

% FOURIER BASIS REPRESENTATION
% Array E(i,j) such that E = sum_{p,q=-N}^{N}E(p,q)*exp(i*(p-q)*theta)

if strfind(outputs,'En-Fou')

 Bvec0 = (Bvec0.*conj(Ip));
 
 E_Fourier = zeros(2*s_Modes+1,2*s_Modes+1);
 
 for loop_s1=1:2*s_Modes+1
  for loop_s2=1:2*s_Modes+1
   E_Fourier(loop_s1,loop_s2)=Bvec0(loop_s1)*conj(Bvec0(loop_s2));
  end
 end
 
 E_Fourier = 2*E_Fourier./k0(1)./pi;
 
 out_str = [out_str '; ''E_Fourier'' '];
 out_val = [out_val '; E_Fourier '];
 
end % end if FOU

%% Natural modes of vibration

if strfind(outputs,'RAOs')
 
 Jmax=5;
 
 % Displacement
 
 res = 1001;
 displ=zeros(2*s_Modes+1,res); displ_inc = displ;
 r_vec = linspace(0,radius,res);
 
 for loop_r=1:length(r_vec)
  count=1;
  for loop_Az=-s_Modes:s_Modes
   J0 = besselj(loop_Az, k0(1)*r_vec(loop_r));
   if BESNM==1
    JJ = besselj(loop_Az, kk*r_vec(loop_r)).*exp(-abs(imag(kk*radius)));
    JJ_mu = besselj(loop_Az, [mu_0;mu_1]*r_vec(loop_r)).*...
     exp(-abs(imag([mu_0;mu_1]*radius)));
   elseif BESNM==2
    JJ = besselj(loop_Az, kk*r_vec(loop_r))./besselj(loop_Az, kk*radius);
    JJ_mu = besselj(loop_Az, [mu_0;mu_1]*r_vec(loop_r))./...
     besselj(loop_Az, [mu_0;mu_1]*radius);
   end
   
   dum = C_mat(Vert_Modes+1,w1)*diag(JJ) + ...
    C_mat(Vert_Modes+1,w2(1))*JJ_mu(1)*Amp_mu_Mat(1,:,s_Modes+1+loop_Az) + ...
    C_mat(Vert_Modes+1,w2(2))*JJ_mu(2)*Amp_mu_Mat(2,:,s_Modes+1+loop_Az);
   
   displ(count,loop_r) = dum*Avec(:,s_Modes+1+loop_Az);
   
   displ_inc(count,loop_r) = diag(J0)*Ivec(1,s_Modes+1+loop_Az);
   
   count = count+1;
   
  end
 end
 clear count dum
 
 % TEST
 if 0
  % HEAVE
  displ(s_Modes+1,:)= 1 + 0*r_vec;
 elseif 0
  % PITCH
  displ(s_Modes,:) = r_vec/2;
  displ(s_Modes+2,:) = r_vec/2;
 end
 
 if PLT
  th=0;
  
  evec = exp(1i*[-s_Modes:s_Modes]*th);
  evec = reshape(evec,1,2*s_Modes+1);
  full_displ = zeros(1,length(r_vec));
  full_displ_inc = full_displ;
  
  for loop_r=1:length(r_vec)
   full_displ(loop_r)     = evec*displ(:,loop_r);
   full_displ_inc(loop_r) = evec*displ_inc(:,loop_r);
  end
  
  figure(PLT); plot(r_vec,real(full_displ),col);
  if 1
   hold on; plot(r_vec,real(exp(1i*k0(1)*r_vec)),'k:');
   plot(r_vec,real(full_displ_inc),'r-.');
  end
  clear evec full_displ
 end
 
 % Inner-products
 
 A_even=zeros(s_Modes+1,Jmax+1); A_even_inc=zeros(s_Modes+1,Jmax+1);
 A_odd =zeros(s_Modes+1,Jmax+1); A_odd_inc =zeros(s_Modes+1,Jmax+1);
 
 for loop_Az=0:s_Modes
  [wnj,Nnj] = fn_NatModes(loop_Az,Jmax,Param.nu,r_vec,radius);
  for loop_j=1:Jmax+1
   if loop_Az~=0
    A_even(loop_Az+1,loop_j) = ...
     pi*fn_Trap(r_vec.*( displ(s_Modes+1+loop_Az,:) +...
     displ(s_Modes+1-loop_Az,:) ),wnj(loop_j,:),r_vec)/Nnj(loop_j);
    A_odd(loop_Az+1,loop_j) = ...
     pi*fn_Trap(r_vec.*( displ(s_Modes+1+loop_Az,:) -...
     displ(s_Modes+1-loop_Az,:) ),wnj(loop_j,:),r_vec)/Nnj(loop_j);
    %
    A_even_inc(loop_Az+1,loop_j) = ...
     pi*fn_Trap(r_vec.*( displ_inc(s_Modes+1+loop_Az,:) +...
     displ_inc(s_Modes+1-loop_Az,:) ),wnj(loop_j,:),r_vec)/Nnj(loop_j);
    A_odd_inc(loop_Az+1,loop_j) = ...
     pi*fn_Trap(r_vec.*( displ_inc(s_Modes+1+loop_Az,:) -...
     displ_inc(s_Modes+1-loop_Az,:) ),wnj(loop_j,:),r_vec)/Nnj(loop_j);
   else
    A_even(loop_Az+1,loop_j) = ...
     2*pi*fn_Trap(r_vec.*displ(s_Modes+1,:),wnj(loop_j,:),r_vec)/Nnj(loop_j);
    A_even_inc(loop_Az+1,loop_j) = ...
     2*pi*fn_Trap(r_vec.*displ_inc(s_Modes+1,:),wnj(loop_j,:),r_vec)/Nnj(loop_j);    
   end
  end
 end
 
 % Ensure even solution
 if max(max(abs(A_odd)))>1e-5
  disp('A_odd = ')
  disp(abs(A_odd))
 end
 
 % Ensure rigid plate (if specified)
 if and(RIGID>1,COMM)
  dum_A=A_even; dum_A(1,1)=0; dum_A(2,1)=0;
  if max(max(abs(dum_A)))>1e-3
   if 1
    cprintf('green',['>>> elastic modes excited: max=' num2str(max(max(abs(dum_A)))) '\n'])
    cprintf('green',['>>> ' fortyp '=' num2str(forval) '\n'])
   else
    disp('A_even = ')
    disp(abs(A_even))
   end
  end
  clear dum_A
 end
 
 %% RAOs
 
 % Heave
 if RAO_TYP==0
  H_inc = 1;
  if COMM; cprintf('blue',['>>> intended RAO normalisation\n']); end
 else
  H_inc = abs(A_even_inc(1,1));
  if COMM; cprintf('blue',['>>> actual RAO normalisation\n']); end
 end
 H_dsk = abs(A_even(1,1));
 
 H_RAO = H_dsk/H_inc;
 
 out_str = [out_str '; ''RAO-heave'' '];
 out_val = [out_val '; H_RAO '];
 
 if COMM
  cprintf('blue',['>>> RAO-heave=' num2str(H_RAO) ...
   ' (inc wave: ' num2str(abs(A_even_inc(1,1)/H_inc))  ') \n'])
 end
 
 % Pitch
 P_dsk = abs(A_even(2,1));
 if RAO_TYP==0
  P_inc = 1; P_RAO = P_dsk/P_inc; %k0(1);
  Pang_RAO = 180*P_RAO/pi;
 else
  P_inc = abs(A_even_inc(2,1)); P_RAO = P_dsk/P_inc;
  Pang_RAO = 180*abs(A_even(2,1))/pi;
 end
 
 if RAO_TYP==0
  out_str = [out_str '; ''RAO-pitch'' '];
  out_val = [out_val '; P_RAO '];
 
  out_str = [out_str '; ''RAO-pitch/k'' '];
  out_val = [out_val '; P_RAO/k0(1) '];
 else
  out_str = [out_str '; ''RAO-pitch/k'' '];
  out_val = [out_val '; P_RAO '];
 end
 
 if COMM
  if RAO_TYP==0
   cprintf('blue',['>>> RAO-pitch/k=' num2str(P_RAO/k0(1)) ...
   ' (inc wave: ' num2str(abs(A_even_inc(2,1))/P_inc/k0(1))  ') \n'])
  else
   cprintf('blue',['>>> RAO-pitch/k=' num2str(P_RAO) ...
   ' (inc wave: ' num2str(abs(A_even_inc(2,1))/P_inc)  ') \n'])
  end
 end
 
 % Pitch angle / inc amp [deg/m] (denote p/a)
 %
 % p = atan(a[m]*P_RAO/1[m])=a*P_RA0 - O(a^3)
 %
 % => p/a = P_RAO  (nb / 1[m])
 
 out_str = [out_str '; ''ANG-pitch'' '];
 out_val = [out_val '; Pang_RAO '];
 
 if SURGE
  % Surge
  if RAO_TYP==0
   S_inc = 1;
  else
   S_inc = exp(1i*[-s_Modes:s_Modes]*0)*displ_inc(:,1);
  end
   
  S_dsk = abs(X_sg);
  
  S_RAO = S_dsk/S_inc;
  
  out_str = [out_str '; ''RAO-surge-nonorm'' '];
  out_val = [out_val '; S_RAO '];
  out_str = [out_str '; ''RAO-surge-norm'' '];
  out_val = [out_val '; S_RAO*tanh(k0(1)*bed) '];
  if COMM
   cprintf('blue',['>>> RAO-surge=' num2str(S_RAO) ...
    ' (eccentricity=' num2str(coth(k0(1)*bed)) ')' '\n'])
  end
 end
 
end % end RAO

%% Form structured output

eval(['out=struct( ''name'', {' out_str ...
 '}, ''value'', {' out_val '});'])
out(1)=[];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% derivatives of Bessel functs. using recurrence formulae

function Besdz = fn_Bessel_dz(invars)

for loop_out=1:length(invars)
 if strcmp(invars(loop_out).name,'BesselFn')
  BesselFn=invars(loop_out).value;
 elseif strcmp(invars(loop_out).name,'N')
  N=invars(loop_out).value;
 elseif strcmp(invars(loop_out).name,'z')
  z=invars(loop_out).value;
 elseif strcmp(invars(loop_out).name,'scale')
  scale=invars(loop_out).value;
 elseif strcmp(invars(loop_out).name,'hank')
  hank=invars(loop_out).value;
 end
end % end loop_out

if ~exist('BesselFn','var')
 cprintf('red','error - need BesselFn')
elseif ~exist('N','var')
 cprintf('red','error - need N')
elseif ~exist('z','var')
 cprintf('red','error - need z')
end

if and(~exist('hank','var'),~exist('scale','var'))
 
 if N == 0
  
  Besdz = - feval( BesselFn, 1, z);
  
 else
  
  Besdz = feval( BesselFn, N-1, z) - feval( BesselFn, N+1, z);
  
  Besdz = Besdz / 2;
  
 end
 
elseif ~exist('hank','var')
 
 if N == 0
  
  Besdz = - feval( BesselFn, 1, z, scale);
  
 else
  
  Besdz = feval( BesselFn, N-1, z, scale) - feval( BesselFn, N+1, z, scale);
  
  Besdz = Besdz / 2;
  
 end
 
elseif ~exist('scale','var')
 
 if N == 0
  
  Besdz = - feval( BesselFn, 1, hank, z);
  
 else
  
  Besdz = feval( BesselFn, N-1, hank, z) - feval( BesselFn, N+1, hank, z);
  
  Besdz = Besdz / 2;
  
 end
 
else
 
 if N == 0
  
  Besdz = - feval( BesselFn, 1, hank, z, scale);
  
 else
  
  Besdz = feval( BesselFn, N-1, hank, z, scale) - feval( BesselFn, N+1, hank, z, scale);
  
  Besdz = Besdz / 2;
  
 end
 
end

return

%

function [ xi, eta ] = fn_BendShear3dAxi( parameter_vector, Bes, arg )

p = 1-parameter_vector(3);

% - Bes = [ fn, N, scale ]

% - fn: 1 = J; 2 = Y; 3 = H1; 4 = I

if Bes(1) == 1
 
 B = ones(size(arg));
 invars = struct('name',{'BesselFn'; 'N'; 'z'}, ...
  'value',{@besselj; Bes(2); arg});
 Bz = fn_Bessel_dz(invars)./besselj(Bes(2), arg);
 
 xi = ( arg.^2 - p*(Bes(2).^2) ).*B + p.*arg.*Bz;
 eta = ( arg.^2 + p*(Bes(2).^2) ).*arg.*Bz - (Bes(2).^2).*p.*B;
 
elseif Bes(1) == 3 % - Hankel of 1st kind !
 
 B = ones(size(arg));
 Bz = Bessel_dz(@besselh, Bes(2), arg)./besselh(Bes(2), arg);
 
 xi = ( arg.^2 - p.*(Bes(2).^2) ).*B + p.*arg.*Bz;
 eta = ( arg.^2 + p.*(Bes(2).^2) ).*arg.*Bz - (Bes(2).^2).*p.*B;
 
elseif Bes(1) == 4
 
 B = ones(size(arg));
 
 Bz = Bessel_Idz(@besseli, Bes(2), arg)./besseli(Bes(2), arg);
 
 xi = ( arg.^2 + p.*(Bes(2).^2) ).*B - p.*arg.*Bz;
 eta = ( arg.^2 - p.*(Bes(2).^2) ).*arg.*Bz + (Bes(2).^2).*p.*B;
 
end

return

%

function [mat_A0,mat_A,mat_W0,mat_W] = ...
 fn_JumpMats(pv,Vert_Dim,DimG,k0,kk,wt_0,wt,METH)

%% -- coeffs of the JCs at r=R -- %

mat_A = matrix_A_PWC(pv, Vert_Dim, kk, wt, pv(10));
mat_A0 = matrix_A_PWC(pv, Vert_Dim, k0, wt_0, pv(7));

if METH==0
 %%% - Old way (vertical modes) - %%%
 %[mat_W0, mat_W] = jump_W_PWC(pv, Vert_Dim, DimG, kk, k0, kk, wt_0, wt);
 [mat_W0, mat_W] = jump_W_PWC_wtd(pv, Vert_Dim, kk, k0, kk, wt, wt_0, wt);
else
 %%% - New way (Gegenbauer) - %%%
 if pv(9)==0
  [mat_W0, mat_W] = Jump_Gegen(k0, kk, wt_0, wt, pv(7), pv(8), 0.5, DimG);
 else
  [mat_W0, mat_W] = Jump_Gegen(k0, kk, wt_0, wt, pv(7), pv(8), 1/6, DimG);
 end
end

%% -- Put jumps together in big matrices - %

% Big_AW0 = zeros((2*Az_Dim+1)*Vert_Dim, (2*Az_Dim+1)*DimG); Big_AW = Big_AW0;
% Big_VT0 = Big_AW0.'; Big_VT = Big_VT0;
%
% v1 = 1:Vert_Dim;
% v2 = 1:DimG;
%
% for loop_Az = 1:2*Az_Dim+1
%     Big_AW0(v1, v2) = mat_A0\mat_W0;
%     Big_AW(v1, v2) = mat_A\mat_W;
%
%     Big_VT0(v2,v1) = mat_W0.';
%     Big_VT(v2,v1) = mat_W.';
%
%     v1 = v1+Vert_Dim; v2 = v2+DimG;
% end

return

function mat_Q = fn_JumpMatSg(pv, Vert_Dim, k, wt)

mat_Q = zeros(Vert_Dim,1);

for loop_dim=1:Vert_Dim
 mat_Q(loop_dim) = wt(loop_dim)*...
  (sinh(k(loop_dim)*pv(7))-sinh(k(loop_dim)*(pv(7)-pv(8))))/k(loop_dim);
end

return

% NATURAL MODES: PLATE IN-VACUO
%
% see Montiel 2012 (Thesis, p.150)
%
% wnj(r,th) = ( Jn(Knj*r/R)+Cnj*In(Knj*r/R) )*cos(n*th) (n>0)
% wnj(r,th) = ( Jn(Knj*r/R)+Cnj*In(Knj*r/R) )           (n=0)
% wnj(r,th) = ( Jn(Knj*r/R)+Cnj*In(Knj*r/R) )*sin(n*th) (n<0)
%
% Nnj = normalisation
% nb. Knj=K-nj, Cnj=(-1)^n*C-nj and Nnj=N-nj
% asymptotic behaviour:
% Knj = (pi/2)(n+2j)   (Leissa, 1969)
% Thus use interval [(pi/2)(n+2(j-1)), (pi/2)(n+2j)] in fzero.

function [wnj,Nnj] = fn_NatModes(n,Jmax,nu,r_vec,R)

if ~exist('COMM','var'); COMM=1; end

params(1) = n; params(2)=nu;

% int = [(pi/2)*(n+2*(j-1)), (pi/2)*(n+2*j)];

int = [(pi/2)*(n+2*(-1)), (pi/2)*(n)];

wnj = zeros(Jmax+1,length(r_vec)); Nnj = zeros(Jmax+1,1);

for j=0:Jmax
 
 if and(n==0,j==0)
  
  wnj = 1+0*r_vec; int = int + pi;
  
 elseif and(abs(n)==1,j==0)
  
  wnj = r_vec; int = int + pi;
  
 else
  
  Knj = inf;
  
  while isinf(Knj)
   Knj = fn_bisection(@fn_disprelNM, int, params);
   int = int + pi;
  end
  
  fj = (Knj.^3).*fn_dbesselj(n,Knj) + (n^2)*(1-nu)*(...
   Knj.*fn_dbesselj(n,Knj) - besselj(n,Knj) );
  
  fi = (Knj.^3).*fn_dbesseli(n,Knj) - (n^2)*(1-nu)*(...
   Knj.*fn_dbesseli(n,Knj) - besseli(n,Knj) );
  
  Cnj = fj/fi; clear fi fj
  
  % CHECK????
  if COMM
   mj = (Knj.^2).*besselj(n,Knj) + (1-nu)*(...
    Knj.*fn_dbesselj(n,Knj) - (n^2)*besselj(n,Knj) );
   
   mi = (Knj.^2).*besseli(n,Knj) - (1-nu)*(...
    Knj.*fn_dbesseli(n,Knj) - (n^2)*besseli(n,Knj) );
   
   if abs(Cnj- mj/mi)>1e-5
    cprintf('red',['error in NMV: (n,j) = (' ...
     int2str(n) ',' int2str(j) ')\n'])
   end
   clear mi mj
  end
  
  wnj(j+1,:) = besselj(n,Knj*r_vec/R) + Cnj*besseli(n,Knj*r_vec/R);
  
 end
 
 if and(n==0,j==0)
  Nnj(j+1)=pi*(R^2);
 elseif and(abs(n)==1,j==0)
  Nnj(j+1) = pi*(R^4)/4;
 else
  Nnj(j+1) = besselj(n,Knj)^2 - besselj(n-1,Knj)*besselj(n+1,Knj) + ...
   (2*Cnj/Knj)*( besselj(n,Knj)*besseli(n-1,Knj) - ...
   besselj(n-1,Knj)*besseli(n,Knj) ) + ...
   (Cnj^2)*( besseli(n,Knj)^2 - besseli(n-1,Knj)*besseli(n+1,Knj) );
  Nnj(j+1) = pi*Nnj(j+1)*(R^2)/2;
 end
end

return

%

function D=fn_disprelNM(K,params)

n = params(1); nu=params(2);

fj = (K.^3).*fn_dbesselj(n,K) + (n^2)*(1-nu)*(...
 K.*fn_dbesselj(n,K) - besselj(n,K) );

fi = (K.^3).*fn_dbesseli(n,K) - (n^2)*(1-nu)*(...
 K.*fn_dbesseli(n,K) - besseli(n,K) );

mj = (K.^2).*besselj(n,K) + (1-nu)*(...
 K.*fn_dbesselj(n,K) - (n^2)*besselj(n,K) );

mi = (K.^2).*besseli(n,K) - (1-nu)*(...
 K.*fn_dbesseli(n,K) - (n^2)*besseli(n,K) );

D = zeros(1,length(K));

for loop=1:length(K)
 D(loop)=det([fj(loop),mj(loop);fi(loop),mi(loop)]);
end

return

% Fabien's bisection method

function out = fn_bisection(fn_name,int,params)

if ~exist('COMM','var'); COMM=0; end

guess_min=int(1); guess_max=int(2); clear int

if guess_min == 0
 guess_min = guess_min+0.01;
end

if feval(fn_name,guess_min, params)*feval(fn_name,guess_max, params)>0
 if COMM
  cprintf('red','no roots found on this interval or multiple roots\n')
 end
 out = inf;
 
else
 
 ans1 = guess_min; ans2 = guess_max;
 
 while abs(ans2 - ans1) > 1e-9
  fval = feval(fn_name,(ans1+ans2)/2, params);
  if fval*feval(fn_name,ans1,params) < 0
   ans2 = (ans1+ans2)/2;
  else
   ans1 = (ans1+ans2)/2;
  end
 end
 
 out = ans2;
 
end

return

%

function I = fn_Trap(f,g,x)

dx = x(2)-x(1);

f = reshape(f,1,length(x));
g = reshape(g,1,length(x));

I = f.*g;

I(1) = I(1)/2; I(end) = I(end)/2;

I = dx*sum(I);

return

function out = fn_dbesselj(n,z,scal)
% this function calcualtes the derivative of besselj(n,z)
if nargin==3
 out = (besselj(n-1,z,scal) - besselj(n+1,z,scal))/2;
else
 out = (besselj(n-1,z) - besselj(n+1,z))/2;
end

return

function out = fn_dbesseli(n,z,scal)
% this function calcualtes the derivative of besseli(n,z)
if nargin==3
 out = (besseli(n+1,z,scal) + besseli(n-1,z,scal))/2;
else
 out = (besseli(n+1,z) + besseli(n-1,z))/2;
end

return

%

function cosh_cosh = fn_cosh_cosh(root_i, root_j, water_depth)

clear parameter_vector

% - i \neq j - %

if root_i ~= root_j
 
 sinh_i = sinh(root_i*water_depth);
 sinh_j = sinh(root_j*water_depth);
 
 cosh_i = cosh(root_i*water_depth);
 cosh_j = cosh(root_j*water_depth);
 
 quot = root_i^2 - root_j^2;
 quot = 1 / quot;
 
 %-------------------------------------------%
 
 cosh_cosh = quot*(root_i*sinh_i*cosh_j - root_j*cosh_i*sinh_j);
 
 % - i = j - %
elseif root_i == root_j
 
 if root_i == 0
  
  cosh_cosh = dum_H;
  
 else
  
  arg = 2*root_i*water_depth;
  
  quot = (4*root_i)^(-1);
  
  sinh2rH = sinh(arg);
  
  cosh2rH = cosh(arg);
  
  cosh_cosh = sinh2rH + arg;
  
  cosh_cosh = cosh_cosh * quot;
 end
end

return

%% OLD OUTPUT CODE

% - Displacement - %
%
% ang=0;
% res = 50; w_vec=zeros(1,2*res-1);
% dum_w=zeros(1,res);
% r_vec = linspace(0,radius,res);
%
% for loop_Az=-s_Modes:s_Modes
%
%     JJ = besselj(loop_Az, kk*r_vec);
%     JJ_mu = besselj(loop_Az, [mu_0;mu_1]*r_vec);
%
%     dum_vec = C_mat(Vert_Modes+1,w1)*JJ +...
%         C_mat(Vert_Modes+1,w2)*(diag(Amp_mu_Mat(:,s_Modes+1+loop_Az))*Jmu_ice);
%
%     dum_w = dum_w + ...
%         exp(1i*loop_Az*(ang+pi))*w(:,s_Modes+loop_Az+1)*dum_vec;
%
%     w_vec(res:2*res-1) = w_vec(res:2*res-1) + ...
%         exp(1i*loop_Az*ang)*w(:,s_Modes+loop_Az+1)*dum_vec;
% end
%
% loop_r = 1;
% for loop_rr=res:-1:2
%     w_vec(loop_r) = dum_w(loop_rr); loop_r=loop_r+1;
% end
%
% r_vec = linspace(-radius,radius,2*res-1);
%
% % - Edge of floe - %
%
% Avec0 = Psi(1,:);
% EI = sum(abs(Avec0).^2);
% EI = EI/2/pi/radius;
%
% % - Other stuff - %
%
% k_out(1) = k0(1); k_out(2) = kk(1);
% Avec1 = Avec(1,:);

% % - Kochin function - %
%
% cprintf('red','Dont use -> not debugged -> use Far-field amplitude instead\n')
%
% fac = Param.rho_0*(freq^3)/8/pi/Param.g;
%
% for loop_th=1:length(th_vec)
%  E(loop_th) = fn_Kochin(pi+th_vec(loop_th),Psi,Amp_mu,k0(1),kk,...
%      [mu_0,mu_1],wt,bed-draught,radius);
% end
%
% E = fac*(abs(E).^2);

%% BELOW IS NOT DEBUGGED - ALSO NOT REQUIRED

% % Kochin function
%
% function H = fn_Kochin(tau,amps,amps_mu,k0,kk,mus,wt,dep,rad)
%
% H = 0;
%
% NM = size(amps);
%
% for loopN = 1:NM(1) % vertical modes
%  kap = kk(loopN); wt_n = wt(loopN);
%  M = -(NM(2)-1)/2;
%  a_mu = amps_mu(:,loopN);
%  for loopM = 1:NM(2) % azimuthal modes
%   H = H + amps(loopN,loopM)*wt_n*(k0*cosh(kap*dep) - kap*sinh(kap*dep))*...
%       (fn_IntBess(abs(M),kap,k0,rad) + ...
%       a_mu(1)*fn_IntBess(abs(M),mus(1),k0,rad) + ...
%       a_mu(2)*fn_IntBess(abs(M),mus(2),k0,rad) ) ...
%                                             *exp(1i*M*(tau-pi/2));
%   M = M + 1;
%  end
% end
%
% H = 2*pi*H;
%
% return
%
% % int_{0}^{R}[besselj(n,a*r)*besselj(-n,b*r)*r]dr
% % using wolframalpha.com
%
% function I = fn_IntBess(n,a,b,R)
%
% if a~=b
%
%  if n~=0
%
%   r = R;
%
%   I = r*( b*besselj(n,a*r)*besselj(n-1,b*r) - ...
%       a*besselj(n,b*r)*besselj(n-1,a*r) );
%
%   r = 0;
%
%   I = I - r*( b*besselj(n,a*r)*besselj(n-1,b*r) - ...
%       a*besselj(n,b*r)*besselj(n-1,a*r) );
%
%   I = ((-1)^n)*I/(a^2 - b^2);
%
%  else
%
%   r = R;
%
%   I = -r*( b*besselj(0,a*r)*besselj(1,b*r) - ...
%       a*besselj(0,b*r)*besselj(1,a*r) );
%
%   r = 0;
%
%   I = I + r*( b*besselj(0,a*r)*besselj(1,b*r) - ...
%       a*besselj(0,b*r)*besselj(1,a*r) );
%
%   I = ((-1)^n)*I/(a^2 - b^2);
%
%  end
%
% else
%
%  if n~=0
%
%   r = R;
%
%   I = 0.5*(r^2)*( besselj(n,a*r)^2 - ...
%       besselj(n+1,a*r)*besselj(n-1,a*r) );
%
%   r = 0;
%
%   I = I - 0.5*(r^2)*( besselj(n,a*r)^2 - ...
%       besselj(n+1,a*r)*besselj(n-1,a*r) );
%
%   I = ((-1)^n)*I;
%
%  else
%
%   r = R;
%
%   I = 0.5*(r^2)*( besselj(0,a*r)^2 + besselj(1,b*r)^2 );
%
%   r = 0;
%
%   I = I - 0.5*(r^2)*( besselj(0,a*r)^2 + besselj(1,b*r)^2 );
%
%   I = ((-1)^n)*I;
%
%  end
%
% end
%
% return

