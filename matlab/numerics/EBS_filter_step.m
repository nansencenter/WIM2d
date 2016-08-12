function M_filter = filter_step(wavdir)

ndir     = length(wavdir);
M_filter = zeros(ndir,ndir);
for i=1:ndir
for j=1:ndir
   th_i           = wavdir(i);                           %% dirn energy is going to
   th_j           = theta_in_range(wavdir(j),th_i-180);  %% dirn energy is coming from (th_i-180<th_j<th_i+180)

   %% =1 if directions are close, otherwise 0 (then it is back-scattered and goes into S_scattered)
   M_filter(i,j)  = (abs(th_i-th_j)<90);
end
end

M_filter = 1-M_filter;%%everything zero => no E goes into "already-scattered" spectrum
