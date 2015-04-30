function  dr = fn_ReDispRel_water(kh, kap_hat)

dr = kh* tanh(kh) - kap_hat;

return