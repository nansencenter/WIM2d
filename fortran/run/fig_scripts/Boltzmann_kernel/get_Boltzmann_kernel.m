function get_Boltzmann_kernel(file_info)
%% save Fourier coefficients from elastic disc code

%%add path to CITEPH:
citeph_path = getenv('CITEPH_PATH');
if strcmp(citeph_path,'')==1
   msg   = ['Trying to add path to CITEPH code...\n',...
            '\n',...
            'Get CITEPH code from\n',...
            'https://github.com/lgbennetts/CITEPH-64-2012\n',...
            '\n',...
            'Then, add ''export CITEPH_PATH=''[location of this repository]''\n',...
            'to the ~.bash_profile file.\n',...
            'NB. need to launch matlab from terminal to access this variable.'];
   error(msg,1);
else
   addpath([citeph_path,'/Fns/AttnModels']);
   addpath([citeph_path,'/Fns/AttnModels/Roots']);
   addpath([citeph_path,'/Fns/misc']);
   addpath([citeph_path,'/Fns/misc/Weight_Standard']);
   addpath([citeph_path,'/Fns/Parametrisation']);
end

if ~exist('fileinfo','var')
   fileinfo.period      = 10;%wave period    [s]
   fileinfo.thickness   = 1;%floe thickness  [m]
   fileinfo.youngs      = 4;%Young's modulus [GPa]
   fileinfo.floe_diam   = 150;%floe diameter [m]
   fileinfo.conc        = .7;%floe diameter [m]
   fileinfo.Nangles     = 2^4;%number of angles (can't be too high as it becomes ill-conditioned)
   fileinfo.Nroots      = 200;%number of roots
end

%% inputs to fn_ElasticDisk.m:
%% out = fn_ElasticDisk(fortyp, forval, Param, outputs, ...
%%          th_vec, RIGID, SURGE, COMM, PLT, col)
in    = getInputs(fileinfo)
out   = fn_ElasticDisk(in.fortyp, in.forval, in.Param, in.outputs, ...
            in.th_vec, in.RIGID, in.SURGE, in.COMM, in.PLT, in.col);

for j=1:length(out)
   if strcmp(out(j).name,'E')
      S_th  = out(j).value;
      %% value of S at theta_vec [m]
      %% (has already got the group velocity inside it?)
   elseif strcmp(out(j).name,'E0')
      beta  = out(j).value;%%scattering cross-section [m]?
   end
end

% Boltzmann eqn: cos(theta)*dE/dx = -beta*E + c/(pi*a^2)*int_{-pi}^{pi} S(th,th')*E(th') dth'
%                                 = -beta*E + int_{-pi}^{pi} K(th,th')*E(th') dth'
% => K has units of m^{-1}
a     = fileinfo.floe_diam/2;%radius [m]
fac   = fileinfo.conc/pi/a^2;%conc/[area of floe]
K_th  = fac*S_th;% [m^{-1}]
beta  = fac*beta;% [m^{-1}]

%% Get Fourier coefficients of scattering kernel:
N     = 2*fileinfo.Nangles;
No2   = round(N/2);
jo    = 2:2:N;%odd modes
%%
th_vec      = in.th_vec;%%1st and last points are same
K_fou       = 2*pi*ifft(K_th);
K_fou(jo)   = -K_fou(jo);%%th_vec is relative to [-pi,pi], fft is relative to [0,2\pi]

%% K is symmetric, so compute cosine expansion:
%% K=K_0/2+\sum_n.K_n.cos(n\theta)
Kcos  = 2*real(K_fou(1:No2+1)).';

if 1
   %%test expansion:
   Kap   = Kcos(1)/2/(2*pi)+0*th_vec;
   for n=1:No2
      Kap   = Kap+Kcos(n+1)/(2*pi)*cos(n*th_vec);
   end
   plot(th_vec/pi,K_th);
   hold on;
   plot(th_vec/pi,Kap,'--r');
end

!mkdir -p out
filename = ['out/Kfou_h',num2str(fileinfo.thickness),...
            '_E',num2str(fileinfo.youngs),...
            '_T',num2str(fileinfo.period),...
            '_D',num2str(fileinfo.floe_diam),...
            '.dat']

writefile(filename,fileinfo,Kcos);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writefile(filename,fileinfo,Kcos)

fid   = fopen(filename,'w');

Atts  = {'period','thickness','youngs','floe_diam','Nangles','Nroots','conc'};

for n=1:length(Atts)
   vbl   = getfield(fileinfo,Atts{n});
   fprintf(fid,'%f    # %s\n',vbl,Atts{n});
end

fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','Boltzmann eqn:');
fprintf(fid,'%s\n','cos(theta)*dE/dx = -beta*E + int_{-pi}^{pi} K(th-s)*E(s) ds');
fprintf(fid,'%s\n','NB: units of K are m^{-1}');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','Cosine coefficients below:');
fprintf(fid,'%s\n','2*pi*K(theta)=K_0/2+\sum_{n=1}^N[K_n*cos(n*theta)]');
fprintf(fid,'%s\n','');

for n=1:length(Kcos)
   fprintf(fid,'%0.15e\n',Kcos(n));
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inputs = getInputs(fileinfo);

inputs.fortyp  = 'freq';
inputs.forval  = 1/fileinfo.period;
inputs.outputs = 'Energy';%scattering kernel
inputs.th_vec  = linspace(-pi,pi,1+2*fileinfo.Nangles);%angles to eval scattering kernel
inputs.th_vec  = inputs.th_vec(1:end-1);
inputs.RIGID   = fileinfo.youngs;%Young's modulus [Gpa]
inputs.SURGE   = 1;%have surge
inputs.COMM    = 1;%display comments
inputs.PLT     = 0;%no plot
inputs.col     = 'r';%plot colour (redundant)

Param          = getParam(fileinfo.thickness,fileinfo.youngs,fileinfo.floe_diam);
Param.Ndtm     = fileinfo.Nroots;
Param.azi      = fileinfo.Nangles;
inputs.Param   = Param;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Param = getParam(thickness,youngs,floe_diam,rho_rel)

%%Param is structure with attributes and default values:
%%          Np: 1                         no of disks 
%%   thickness: 0.033000000000000         thickness [m]
%%           g: 9.810000000000000         gravity [m s^{-2}]
%%       rho_0: 1025                      water density [kg m^{-3}]
%%         rho: 5.590909090909091e+02     floe density [kg m^{-3}]
%%          nu: 0.300000000000000         Poisson's ratio
%%           E: 1.000000000000000e+10     Youngs modulus [Pa]
%%       draft: 0.018000000000000         draft, rho/rho_0*thickness [Pa m^3]
%%           D: 3.290934065934066e+04     flexural rigidity
%%        beta: 3.272851561059214         scaled rigidity, D/rho/g [m^4]
%%         bed: 3.100000000000000         water depth [m]
%%  MIZ_length: 5                         length of MIZ [m]
%%   floe_diam: 0.990000000000000         floe diameter [m]

% Floes are rigid:
if ~exist('youngs','var');
   youngs   = 5.45; %Young's modulus [Gpa]
end
if ~exist('thickness','var');
   thickness   = 1; %thickness [m]
end
if ~exist('floe_diam','var');
   floe_diam   = 150; %floe diameter [m]
end
if ~exist('rho_rel','var');
   rho_rel  = .9; % rho_ice/rho_wtr
end

% Number of disks
Param.Np = 1;

% thickness
Param.thickness = thickness;

% Accelaration due to gravity (in m\,s^{-2})
Param.g = 9.81;

% Density of fluid (in kg\,m^{-3})
Param.rho_0 = 1025;

%%%%% Choose to give draught 1.8/100 %%%%%

% Densities of disks (in kg\,m^{-3})
Param.rho = rho_rel*Param.rho_0; %500*ones(Param.Np,1);

% Poisson's ratios
Param.nu = 0.3;

% Young's modulii (in Pa)
Param.E = 1e10*youngs;

% Draughts (in m)
Param.draft = Param.thickness*rho_rel;

% Flexural rigidities (in MPa\,m^3)
Param.D = Param.E.*Param.thickness.^3./(12*(1-Param.nu.^2));

% Scaled rigidity 
Param.beta = Param.D./(Param.g.*Param.rho_0);

% fluid depth (open/equilibrium)

Param.bed = 300;

% MIZ length [m]

Param.MIZ_length = 5;

% diameter of floes [m]

Param.floe_diam = .99;

return

