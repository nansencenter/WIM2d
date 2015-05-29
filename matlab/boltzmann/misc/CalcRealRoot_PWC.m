function k = CalcRealRoot_PWC(args, fn, UL_fn, Tol_vec)

tol_x = Tol_vec(1);

upp_limit = feval(UL_fn, args);

interval = [ 0 upp_limit ];

k = fzero(fn, interval, optimset('tolX', tol_x), args);