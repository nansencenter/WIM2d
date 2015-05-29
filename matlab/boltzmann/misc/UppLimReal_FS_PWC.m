function UL = UppLimReal_FS_PWC(args)

hh = args(1);
kappa = args(2);

kappa_t = tanh(kappa);

UL = max(kappa/ kappa_t, kappa / hh);