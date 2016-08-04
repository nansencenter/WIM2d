function mwd = calc_mwd_1freq(wavdir,Sdir)
%% function for testing:
%% split energy between fwd and back directions
%% also test integration schemes (const panel vs Simpson's)

[nx,ny,ndir]   = size(Sdir);


%% integrate over directions (splitting fwd/back):
th_vec   = -pi/180*(wavdir+90);%%-90 deg (from W) -> 0 radians; 0 deg (from N) -> -pi/2 radians
dth      = th_vec(2)-th_vec(1);
%%
Mom0  = zeros(nx,ny);
mom0  = zeros(nx,ny);
for wth=1:ndir
   mom0  = mom0+dth*Sdir(:,:,wth);
   Mom0  = Mom0+dth*Sdir(:,:,wth)*th_vec(wth);
end
%%
mwd      = zeros(nx,ny);
jp       = find(mom0>0);
mwd(jp)  = -90 -180/pi*sqrt(Mom0(jp)./mom0(jp));
