function H  = fn_pcolor(x,y,z,labs)
%% fn_pcolor.m
%% Author: Timothy Williams
%% Date: 20150820, 13:32:28 CEST

if nargin==0
   %%test inputs
   z  = hadamard(12);
   x  = (1:6 )';
   y  = (1:12)';
   z  = z(x,:);
end

if size(x,2)>1
   x  = x';
end
if size(y,2)>1
   y  = y';
end

%% pcolor ignores last if shading is faceted
[nx,ny]        = size(z);
Z              = zeros(nx+1,ny+1);
Z(1:nx,1:ny)   = z;

%% want x,y to refer to edges of cells
xx = cen2edges(x);
yy = cen2edges(y);

H  = pcolor(xx,yy,Z');
set(H,'EdgeColor', 'none');

if ~exist('labs','var')
   labs  = {'\itx, \rmkm','\ity, \rmkm',[]};
end
GEN_proc_fig(labs{1},labs{2});
colorbar;
GEN_font(gca);

if ~isempty(labs{3})
   ttl   = title(labs{3});
   GEN_font(ttl);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xx = cen2edges(x);
dfx   = diff(x);
xx    = (x(1:end-1)+x(2:end))/2;
xx    = [x(1)-.5*dfx(1);xx;x(end)+.5*dfx(end)];
