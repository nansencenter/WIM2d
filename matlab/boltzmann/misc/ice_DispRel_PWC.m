function dr = ice_DispRel_PWC(k, args)

dummy_alpha = args(1);
dummy_beta = args(2);
dummy_cap_H = args(3);
kappa = args(4);            

dr = (1 - kappa*dummy_alpha + dummy_beta* (k^4) ) *...
    k * tanh(k* dummy_cap_H) - kappa;