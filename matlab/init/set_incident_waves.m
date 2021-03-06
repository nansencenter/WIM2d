%% get_wave_spec.m
%% Author: Timothy Williams
%% Date: 20141016, 18:55:20 CEST
function wave_stuff = set_incident_waves(grid_prams,wave_fields,inc_options)

%grid_prams,wave_fields,inc_options
if ~exist('inc_options','var')
   %%default options
   nw                = 1;  %% Single frequency
   ndir              = 2^4;%% Multiple directions
   Tmin              = 2.5;
   Tmax              = 25;
   DIRSPEC_INC_OPT   = 1;  %% cos^2 spreading
else
   nw                = inc_options.nw;
   ndir              = inc_options.ndir;
   Tmin              = inc_options.Tmin;
   Tmax              = inc_options.Tmax;
   DIRSPEC_INC_OPT   = inc_options.DIRSPEC_INC_OPT;
   FRQSPEC_INC_OPT   = inc_options.FRQSPEC_INC_OPT;
end

%%frequency grid:
if nw==1%%single freq
   T     = max(wave_fields.Tp(:));
   freq  = 1/T;
else
   f      = 1/Tmax;%0.042;% min freq/ resolution
   f1     = 1/Tmin;%0.4;% max freq
   freq   = linspace(f,f1,nw)';
   om     = 2*pi*freq;
end

%%direction grid:
if ndir==1
   wavdir   = -90;
else
   wavdir   = linspace(90,-270,ndir+1)';
   %%
   %if mod(ndir,2)==0
   %   %%make symmetric
   %   wavdir   = wavdir+(wavdir(2)-wavdir(1))/2;
   %end
   wavdir(end) = [];
end
wave_stuff  = struct('nfreq',nw,...
                     'ndir',ndir,...
                     'freq',freq,...
                     'dirs',wavdir);

%% change from: waves-from (degrees, 0=north, clockwise)
%%          to: waves-to   (radians, 0=east,  anti-clockwise)
theta = -pi/180*(90+wavdir);%theta/pi
dth   = 2*pi/ndir;%%direcional resolution in radians
nx    = grid_prams.nx;
ny    = grid_prams.ny;
Sdir  = zeros(nx,ny,ndir,nw);

for i=1:nx
   for j=1:ny
      wmsk  = wave_fields.WAVE_MASK(i,j);

      if wmsk==1
         Hs    = wave_fields.Hs(i,j);
         Tp    = wave_fields.Tp(i,j);
         mwd0  = wave_fields.mwd(i,j);
         mwd   = -pi/180*(90+mwd0);

         if nw>1
            %% Bretschneider spectrum
            %% - freq spectrum
            om    = 2*pi*freq;
            if FRQSPEC_INC_OPT
             Sfreq = SDF_PM(om,{Tp,Hs});
            else
             Sfreq = SDF_Bretschneider(om,{Tp,Hs});
            end
         else
            Sfreq = (Hs/4)^2;%%Hs=4*sqrt(Sfreq*wt_om), wt_om=1
         end

         if ndir>1
            %%directions
            del         = theta-mwd;
            j0          = find(del<-pi);
            del(j0)     = del(j0)+2*pi;
            j0          = find(del>pi);
            del(j0)     = del(j0)-2*pi;

            if DIRSPEC_INC_OPT==1
               %% cos^2 spreading fxn
               %dir_fac     = (1+cos(2*del))/2/(pi/2);
               %j0          = find(abs(del)>pi/2);
               %dir_fac(j0) = 0
               dir_fac  = 0*theta;
               R2D      = 180/pi;%radians to degrees
               for wth=1:ndir
                  dir_fac(wth)   = theta_dirfrac(...
                     R2D*(theta(wth)-dth/2),R2D*dth,R2D*mwd)/dth;
               end
            else
               %% Delta function
               jmwd           = find(abs(del)==min(abs(del)));
               dir_fac        = 0*wavdir;
               dir_fac(jmwd)  = 1/dth;
            end

            % mwd
            % [180/pi*[theta,del],dir_fac]
            % dth   = abs(theta(2)-theta(1));
            % sum(dth*dir_fac)
            % GEN_pause

            for jw=1:nw
            for jt=1:ndir
               Sdir(i,j,jt,jw)   = Sfreq(jw)*dir_fac(jt);
            end
            end
         else
            Sdir(i,j,1,:)  = Sfreq;
         end
      end
   end
end

wave_stuff.dir_spec  = Sdir;
if ndir==1
   wave_stuff.dirs   = mwd0;
end
