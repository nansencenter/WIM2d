function H  = fn_pcolor(x,y,z,labs,zlim,remove_outer)
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

if ~exist('labs','var'); labs  = {'\itx, \rmkm','\ity, \rmkm',[]}; end
if ~exist('zlim','var'); zlim  = []; end
if ~exist('remove_outer','var'); remove_outer   = 0; end

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

if (length(x)==(nx+1))&(length(y)==(ny+1))
   %% q points already (corners)
   xx = x;
   yy = y;
elseif (length(x)==nx)&(length(y)==ny)
   %% want x,y to refer to edges of cells
   xx = cen2edges(x);
   yy = cen2edges(y);
else
   error('x,y wrong size');
end

fontsize = 20;
if remove_outer
   xx(1)    = [];
   xx(end)  = [];
   yy(1)    = [];
   yy(end)  = [];
   Z(1,:)   = [];
   Z(end,:) = [];
   Z(:,1)   = [];
   Z(:,end) = [];
end

H  = pcolor(xx,yy,Z');%% rows of arg 3 correspond to y not x
daspect([1,1,1]);
set(H,'EdgeColor', 'none');

set(gca,'xlim',[xx(1),xx(end)],'ylim',[yy(1),yy(end)])
if 0%%set ticks manually
   set(gca,'ytick',[100,200]);
   %set(gca,'ytick',[100,200],'xtick',[200,400]);
end

GEN_proc_fig(labs{1},labs{2},fontsize);
cb    = colorbar;
if ~isempty(zlim)
   caxis(zlim);
   Ytick_range = [round(zlim(1)),round(mean(zlim)),round(zlim(2))]
   if length(unique(Ytick_range))==3
      set(cb,'YTick',Ytick_range,'fontname','Times','fontsize',fontsize);
   end
end
GEN_font(gca,fontsize);

if ~isempty(labs{3})
   %ttl   = title(labs{3});
   ttl   = ylabel(cb,labs{3},'rotation',90);
   GEN_font(ttl,fontsize);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xx = cen2edges(x);
dfx   = diff(x);
xx    = (x(1:end-1)+x(2:end))/2;
xx    = [x(1)-.5*dfx(1);xx;x(end)+.5*dfx(end)];
