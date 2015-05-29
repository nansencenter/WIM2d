function w = fn_ImagRoot_water(no, dum_h, kappa, tol_x)

%% Search for non-dimensional root w*h:

upp = pi * (no + 0.5);
low = pi * (no - 0.5);

w = fzero('fn_ImDispRel_water', [low,upp], ...
 optimset('tolX', tol_x), kappa*dum_h);

%% Dimensionalise:

w = w/dum_h;

return