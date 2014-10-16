%% get_wave_spec.m
%% Author: Timothy Williams
%% Date: 20141016, 18:55:20 CEST
function wave_stuff = get_wave_spec(grid_params,wave_fields)

%%frequency grid:
f      = 1/16;%0.042;% min freq/ resolution
f1     = 1/2.5;%0.4;% max freq
nw     = 21;% NB needs to be odd for Simpson's rule;
freq   = linspace(f,f1,nw)';
T      = 1./freq;

%%direction grid:
ndir        = 8;
wavdir      = linspace(-180,180,ndir+1)';
wavdir(end) = [];
wave_stuff  = struct('nfreq',nw,...
                     'ndir',ndir,...
                     'freq',freq,...
                     'dirs',wavdir);

%% change from: waves-from (degrees, 0=north, clockwise)
%%          to: waves-to   (radians, 0=east,  anti-clockwise)
theta = pi/180*(90-wavdir);
dth   = abs(theta(1)-theta(2));%%dir resolution
nx    = grid_prams.nx;
ny    = grid_prams.ny;
Sdir  = zeros(nx,ny,nw,ndir);

for i=1:nx
   for j=1:ny
      wmsk  = wave_fields.WAVE_MASK(i,j);
      Hs    = wave_fields.Hs(i,j);
      Tp    = wave_fields.Tp(i,j);
      mwd   = pi/180*(90-wave_fields.mwd(i,j));%%

      if wmsk==1
         %% Bretschneider spectrum
         %% - freq spectrum
         Sfreq    = SDF_Bretschneider(om,{Tp,Hs});

         %%directions
         del         = theta-mwd;
         j0          = find(del<-pi);
         del(j0)     = del(j0)+2*pi;
         j0          = find(del>pi);
         del(j0)     = del(j0)-2*pi;

         if SHARP_DIST==0
            %% Spreading fxn
            dir_fac     = (1+cos(del))/2/(pi/2); 
            j0          = find(abs(del)>pi/2);
            dir_fac(j0) = 0;
         else
            %% Delta function
            jmwd           = find(abs(del)==min(abs(del)));
            dir_fac        = 0*wavdir;
            dir_fac(jmwd)  = 1/dth;
         end

         for jw=1:nw
         for jt=1:ndir
            Sdir(i,j,jw,jt)   = Sfreq(jw)*dir_fac(jt);
         end
         end
      end
   end
end

wave_stuff.dir_spec  = Sdir;
