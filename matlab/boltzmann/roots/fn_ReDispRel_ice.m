function dr = fn_ReDispRel_ice(kh, args)

al_hat  = args(1);
be_hat  = args(2);
kap_hat = args(3);            

dr = (al_hat + be_hat* (kh^4) ) *kh * tanh(kh) - kap_hat;
   
return