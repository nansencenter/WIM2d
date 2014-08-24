function [S2,ag_eff] =...
   FBwave_advect2D_v1(S,wave_props,conc,grid_props);

%% outputs: 
%%  S2(i,j)=advected wave spectrum in (i,j) cell;
%%  NB i increases going south;
%%     j increases going east;
%%
%% inputs:
%%  S(i,j)=initial wave spectrum in (i,j) cell;
%%  wave properties: wave_props = {ag,dirn};
%%   ag = group velocity (if scalar
%%    this is taken to be uniform over whole grid);
%%   dirn = wave direction (degrees clockwise from north;
%%    NB always scalar - this is taken to be
%%     uniform over whole grid);
%%
%% conc - scalar or matrix;
%%
%%  grid properties: grid_props = {dt,dx,dy};
%%   dt = time step;
%%   dx = grid size in x (west-east) direction;
%%   dy = grid size in y (north-south) direction;

ag    = wave_props{1};
dirn  = wave_props{2};
%%
dt = grid_props{1};
dx = grid_props{2};
if length(grid_props)==2
   dy = dx;
else
   dy = grid_props{3};
end


if length(ag)==1
   ag_eff   = ag+0*S;
else
   ag_wtr   = ag{1};
   ag_ice   = ag{2};
   ag_eff   = conc.*ag_ice+(1-conc)*ag_wtr;
end


nx = size(S,2);
ny = size(S,1);
%%
zx = zeros(nx+2,1)';
zy = zeros(ny,1);
S  = [zx;[zy,S,zy];zx];
S2 = 0*S;
%%
for i = 2:ny+1
   r  = ny+3-i;
   for j = 2:nx+1
      if conc(r-1,j-1)>=0
         [r2,j2,th] = aux_fn(dirn,r,j);
         %%
         alp1     = abs(ag_eff(r-1,j-1)*dt/dx*...
                     cos(th));
         alp2     = abs(ag_eff(r-1,j-1)*dt/dy*...
                     sin(th));
         S2(r,j)  = (1-alp1)*(1-alp2)*S(r,j)+...
                     +alp1*alp2*S(r2,j2)+...
                     +alp1*(1-alp2)*S(r,j2)+...
                     +(1-alp1)*alp2*S(r2,j);
%         else
%            disp('land'),disp([i-1 j-1]),pause,
      end
   end
end

S2([1 end],:)  = [];
S2(:,[1 end])  = [];
%% S2 = S2.*exp(-alp_dim.*ag_eff*dt);


function [r2,j2,th]=aux_fn(dirn,r,j)

%dirn,i,j,ny
if dirn<0
   dirn  = dirn+360;
end

th = pi/180*(90-dirn);

if dirn<=90
   j2 = j-1;
   r2 = r+1;
elseif dirn<=180
   j2 = j-1;
   r2 = r-1;
elseif dirn<=270
   j2 = j+1;
   r2 = r-1;
elseif dirn<=360
   j2 = j+1;
   r2 = r+1;
end
