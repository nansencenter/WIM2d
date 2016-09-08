function M_filter = EBS_filter(wavdir,shape)
%% NxN matrix (N = no of dirns);
%% matrix = 0 => E stays in normal sdf
%% matrix = 1 => E goes to E_S

if ~exist('shape','var')
   shape = 'step';
end

ndir     = length(wavdir);
M_filter = zeros(ndir,ndir);%%zeros correspond to staying in normal E
th_vec   = -pi/180*(90+wavdir);

for i=1:ndir
   if strcmp(shape,'delta')
      i_ = mod(i+ndir/2,ndir);
      if i_==0
         i_ = ndir;
      end
      %% i_ ~ to-dirn; i~ from-dirn
      %% - back-scattered waves go to E_S
      M_filter(i_,i) = 1;

   elseif strcmp(shape,'cos')
      th_i           = th_vec(i);
      M_filter(:,i)  = max(0,cos(th_vec-th_i-pi));%%1 at th_vec=th_i+pi

   elseif strcmp(shape,'step')
      for j=1:ndir
         %% dirn-al bin energy is going to
         th_i  = wavdir(i);

         %% dirn-al bin energy is coming from (th_i-180<th_j<th_i+180)
         th_j  = theta_in_range(wavdir(j),th_i-180);

         %% =0 if directions are close, otherwise 1 (then it is back-scattered and goes into S_scattered)
         M_filter(i,j)  = (abs(th_i-th_j)>=90);
      end

   else
      disp('Unknown option for ''shape''');
   end

end
