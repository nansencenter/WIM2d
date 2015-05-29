function k = fn_RealRoot(args, fn, UL_fn, tol_x)

upp_limit = feval(UL_fn, args);

interval = [ 0 upp_limit ];

k = fzero(fn, interval, optimset('tolX', tol_x), args);

return