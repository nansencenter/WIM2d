%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_dirfrac = theta_dirfrac(th1_,dtheta_,mwd_)
%% chi=pi/180*(theta-mwd) 
%%  (convert to radians and centre about mwd, the mean wave direction)
%% chi1=pi/180*(th1-mwd)
%% chi2=pi/180*(th2-mwd)
%% D(chi)=2/pi*cos^2(chi) if chi\in[-pi/2,pi/2] 
%% D(chi)=0 if chi\in[-pi/2,pi/2] (no backwards waves)
%% theta_dirfrac  = \int_{chi1}^{chi2}D(chi)dchi
%%                = 1 if -chi1=chi2=pi/2
%%                = the fraction of the energy going in the directions in [chi1,chi2]

%%get mwd inside [th1,th1+360)
mwd   = theta_in_range(mwd_,th1_); %>th1_
th2   = th1_+dtheta_;
if ((mwd>th2)&(mwd-th2)>abs(mwd-360-th1_))
   mwd   = mwd-360;
end
th1   = max(mwd-90,th1_);
th2   = min(mwd+90,th2);
th2   = max(th1,th2);%make th2>=th1

chi1  = pi*(th1-mwd)/180.;
chi2  = pi*(th2-mwd)/180.;

theta_dirfrac  = 2*(chi2-chi1)+sin(2*chi2)-sin(2*chi1);
theta_dirfrac  = theta_dirfrac/2/pi;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
