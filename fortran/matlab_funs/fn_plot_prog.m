function out = fn_plot_prog(grid_prams,n,outdir)

out   = fn_check_prog(outdir,n);
fn_fullscreen;

if isempty(grid_prams)
   grid_prams  = fn_get_grid(outdir);
end

X        = grid_prams.X;
Y        = grid_prams.Y;
[nx,ny]  = size(X);
%%
vbls     = {'Hs','Dmax','tau_x','tau_y'};
lab3     = {'{\itH}_{\rm s}, m','{\itD}_{\rm max}, m',...
            '{\tau}_{x}, \rmPa','{\tau}_{y}, \rmPa'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%do plots

%%fix positions so figures can be compared more easily between computers
pos{1}   = [0.130000000000000   0.583837209302326   0.334659090909091   0.341162790697674];
pos{2}   = [0.570340909090909   0.583837209302326   0.334659090909091   0.341162790697674];
pos{3}   = [0.130000000000000   0.110000000000000   0.334659090909091   0.341162790697674];
pos{4}   = [0.570340909090909   0.110000000000000   0.334659090909091   0.341162790697674];

for j=1:4
   subplot('position',pos{j});
   eval(['Z = out.',vbls{j},';']);
   if ny==1
      labs  = {'\itx, \rmkm',lab3{j}};
      fn_plot1d(X/1e3,Z,labs);
   else
      labs  = {'\itx, \rmkm','\ity, \rmkm',lab3{j}};
      fn_pcolor(X(:,1)/1e3,Y(1,:)/1e3,Z,labs);
      daspect([1,1,1]);
   end
end
