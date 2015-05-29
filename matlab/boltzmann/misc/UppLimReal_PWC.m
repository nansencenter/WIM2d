function upp_limit = UppLimReal_PWC(args)

dummy_alpha = args(1);
dummy_beta = args(2);
dummy_cap_H = args(3);
kappa = args(4);

kappa_t = tanh(kappa);

dummy_vec = [kappa / kappa_t, kappa / dummy_cap_H, ...
    (kappa*dummy_alpha / dummy_beta)^(1/4)];

upp_limit = max( dummy_vec );