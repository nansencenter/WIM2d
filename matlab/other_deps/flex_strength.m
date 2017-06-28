function sigf=flex_strength(inputs)

if isfield(inputs,'brine_volume_fraction')
   vb = getfield(inputs,'brine_volume_fraction');
else
   S  = inputs.salinity
   T  = inputs.temp;
   %%
   vb       = nan*T;
   jok      = find(T>-22.9 & T<-.5);
   vb(jok)  = S.*(49.185./abs(T)+.532);
end
sigf  = 1.76e6*exp(-5.88*sqrt(vb));
