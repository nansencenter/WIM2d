function out = fn_ImDispRel_ice(wh,prams)

out = (prams(1)+prams(2)*(wh.^4)).*wh.*sin(wh) + prams(3)*cos(wh);

return