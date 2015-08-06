function [xx,yy] = step_1d(x,y)
% vectors for constant-panel plotting

nx       = length(x);
dx       = x(2)-x(1);
%%
xx          = zeros(2*nx,1);
xx(1:2:end) = x-.5*dx;
xx(2:2:end) = x+.5*dx;
%%
yy          = zeros(2*nx,1);
yy(1:2:end) = y;
yy(2:2:end) = y;
