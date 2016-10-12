x  = (-20:20)';%km
y  = 12;%km

TEST_MPROJ  = 1; %% 1: test m_xy2ll from m_proj; 0: test projection used in nextsim c++ code
if TEST_MPROJ
   %% m_xy2ll,m_ll2wy use sphere of radius 1
   %% stereographic projection has true scale lat at center of projection
   R  = 6378.273;
   m_proj('Stereographic','lon',-45,'lat',90,'radius',60);
end

ss          = sprintf('%s   %s   %s   %s','x','y','lon','lat');
disp(ss);
for i=1:length(x)
   xi = x(i);
   if TEST_MPROJ==1
      [lon,lat]   = m_xy2ll(xi/R,y/R);
   else
      [lat,lon]   = mapxy(xi,y,60,-45,'N');
   end
   ss          = sprintf('%f   %f   %f   %f',xi,y,lon,lat);
   disp(ss);
end
