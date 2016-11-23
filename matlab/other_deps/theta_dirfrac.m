%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = theta_dirfrac(th1_,dtheta_,mwd_)
%% chi=pi/180*(theta-mwd) 
%%  (convert to radians and centre about mwd, the mean wave direction)
%% chi1=pi/180*(th1-mwd)
%% chi2=pi/180*(th2-mwd)
%% D(chi)=2/pi*cos^2(chi) if chi\in[-pi/2,pi/2] 
%% D(chi)=0 if chi\in[-pi/2,pi/2] (no backwards waves)
%% theta_dirfrac  = \int_{chi1}^{chi2}D(chi)dchi
%%                = 1 if -chi1=chi2=pi/2
%%                = the fraction of the energy going in the directions in [chi1,chi2]

if nargin==0
   %%do a test
   mwd      = 90
   N        = 110
   dtheta   = 360/N;
   thvec    = (0:N-1)'*dtheta;
   I        = 0;
   for j=1:N
      I  = I+theta_dirfrac(thvec(j),dtheta,mwd);
   end
   tst_int  = [I,1]
   return
end

%%get mwd inside [th1,th1+360)
phi1  = theta_in_range(mwd_-90,th1_); %>=th1_
phi2  = theta_in_range(mwd_+90,th1_); %>=th1_
th2_  = th1_+dtheta_;
I     = 0;
if phi2>phi1
   %% th1,phi1,phi2, and th2
   L1    = max(th1_,phi1);
   L2    = min(th2_,phi2);
   L2    = max(L1,L2);%make L2>=L1
   chi1  = pi*(L1-mwd_)/180.;
   chi2  = pi*(L2-mwd_)/180.;
   I     = I + 2*(chi2-chi1)+sin(2*chi2)-sin(2*chi1);
else
   %% th1,phi2,phi1, and th2
   %% 1st consider (th1,phi2) interval
   L1    = th1_;
   L2    = min(th2_,phi2);
   chi1  = pi*(L1-mwd_)/180.;
   chi2  = pi*(L2-mwd_)/180.;
   I     = I + 2*(chi2-chi1)+sin(2*chi2)-sin(2*chi1);

   %% 2nd consider (phi1,th2) interval
   L1    = phi1;
   L2    = max(L1,th2_);%make L2>=L1
   chi1  = pi*(L1-mwd_)/180.;
   chi2  = pi*(L2-mwd_)/180.;
   I     = I + 2*(chi2-chi1)+sin(2*chi2)-sin(2*chi1);
end

y  = I/2/pi;
