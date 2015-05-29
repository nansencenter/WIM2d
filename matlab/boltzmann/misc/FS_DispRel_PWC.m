function  dr = FS_DispRel_PWC(k, args)

hh = args(1);
kappa = args(2);

dr = k* tanh(k*hh) - kappa;