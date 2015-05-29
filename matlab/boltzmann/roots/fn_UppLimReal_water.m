function UL = fn_UppLimReal_water(kap_hat)

kappa_t = tanh(kap_hat);

UL = max(kap_hat/ kappa_t, kap_hat);

return