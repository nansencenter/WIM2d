function out = calc_dirn_integrals_1freq(wavdir,Sdir)
%% function for testing:
%% split energy between fwd and back directions
%% also test integration schemes (const panel vs Simpson's)

if nargin==0
   ndir     = 16;
   th_vec   = (0:ndir-1)'*2*pi/ndir;
   wavdir   = -180/pi*th_vec-90;
   if 0
      %% isotropic
      r1    = 0
      sig1  = sqrt(2*(1-r1))
      disp('test isotropic spectrum')
      Sdir  = ones(1,1,ndir);
   elseif 0
      %% delta function
      r1    = 1
      sig1  = sqrt(2*(1-r1))
      disp('test delta function spectrum')
      Sdir     = zeros(1,1,ndir);
      Sdir(1)  = 1;
   elseif 1
      %% cos^2
      I1    = 2/pi*quad('cos(x).^3',-pi/2,pi/2);
      r1    = sqrt(I1^2)
      sig1  = sqrt(2*(1-r1))
      disp('test cos^2 (fwd only) spectrum')
      Sdir           = zeros(1,1,ndir);
      jp             = find(cos(th_vec)>=0);
      Sdir(1,1,jp)   = cos(th_vec(jp)).^2;
   end
end

[nx,ny,ndir]   = size(Sdir);


%% integrate over directions (splitting fwd/back):
th_vec   = -pi/180*(wavdir+90);%%-90 deg (from W) -> 0 radians; 0 deg (from N) -> -pi/2 radians
dth      = th_vec(2)-th_vec(1);
%%
Mom0  = zeros(nx,ny);
mom0  = zeros(nx,ny);
mom0p = zeros(nx,ny);
mom0m = zeros(nx,ny);
MomS  = zeros(nx,ny);
MomC  = zeros(nx,ny);
MomSp = zeros(nx,ny);
MomCp = zeros(nx,ny);
MomSm = zeros(nx,ny);
MomCm = zeros(nx,ny);
for wth=1:ndir
   fwd   = (cos(th_vec(wth))>=0);
   rev   = (cos(th_vec(wth))<=0);

   %% energy
   mom0  = mom0  +     dth*Sdir(:,:,wth);
   mom0p = mom0p + fwd*dth*Sdir(:,:,wth);
   mom0m = mom0m + rev*dth*Sdir(:,:,wth);

   %% mwd & directional wave spreading
   MomS  = MomS  +     dth*Sdir(:,:,wth)*sin(th_vec(wth));
   MomC  = MomC  +     dth*Sdir(:,:,wth)*cos(th_vec(wth));
   MomSp = MomSp + fwd*dth*Sdir(:,:,wth)*sin(th_vec(wth));
   MomCp = MomCp + fwd*dth*Sdir(:,:,wth)*cos(th_vec(wth));
   MomSm = MomSm + rev*dth*Sdir(:,:,wth)*sin(th_vec(wth));
   MomCm = MomCm + rev*dth*Sdir(:,:,wth)*cos(th_vec(wth));
end

%% mwd (Hauser et al., 2005)
tmp         = atan2(MomS,MomC);%\theta_1
out.mwd     = -90 -180/pi*tmp;
tmp         = atan2(MomSp,MomCp);%\theta_1 (fwd)
out.mwd_p   = -90 -180/pi*tmp;
tmp         = atan2(MomSm,MomCm);%\theta_1 (rev)
out.mwd_m   = -90 -180/pi*tmp;

%% spreading (Hauser et al., 2005)
jp             = find(mom0>0);
out.sprd       = zeros(nx,ny);
tmp            = sqrt(MomS.^2+MomC.^2);          %r_1
out.sprd(jp)   = sqrt(2*( 1-tmp(jp)./mom0(jp) ));%\sigma_1

jp             = find(mom0p>0);
out.sprd_p     = zeros(nx,ny);
tmp            = sqrt(MomSp.^2+MomCp.^2);
out.sprd_p(jp) = sqrt(2*( 1-tmp(jp)./mom0p(jp) ));

jp             = find(mom0m>0);
out.sprd_m     = zeros(nx,ny);
tmp            = sqrt(MomSm.^2+MomCm.^2);
out.sprd_m(jp) = sqrt(2*( 1-tmp(jp)./mom0m(jp) ));

out.E_tot   = mom0;
out.E_p     = mom0p;
out.E_m     = mom0m;
