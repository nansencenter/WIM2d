x  = (-20:20)';%km
y  = 12;%km
for i=1:length(x)
   xi          = x(i);
   [lat,lon]   = mapxy(xi,y,60,-45,'N');
   ss          = sprintf('%f   %f   %f   %f',xi,y,lon,lat);
   disp(ss);
end
