function th=theta_in_range(th_,th1)
%% th = theta_in_range(th_,th1)
%% where th1<=th<th1+360.

th2   = th1+360.0;
if (th_<th1)
   dth   = th1-th_;
   njump = ceil(dth/360.);
   th    = th_+njump*360.;
elseif (th_>th2)
   dth   = th_-th2;
   njump = ceil(dth/360.);
   th    = th_-njump*360.;
elseif (th_==th2)
   th = th1;
else
   th = th_;
end

return
