function upp_limit = fn_UppLimReal_ice(args)

al_hat = args(1);
be_hat = args(2);
kap_hat = args(3);

kappa_t = tanh(kap_hat);

dummy_vec = [kap_hat / kappa_t, kap_hat, ...
    ( (1-al_hat) / be_hat)^(1/4)];

upp_limit = max( dummy_vec );

return