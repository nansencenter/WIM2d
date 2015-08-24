function fn_plot_conts(S,lcol,lwid,lstil)

if ~exist('lcol')
   lcol  = 'k';
end

if ~exist('lwid')
   lwid  = 2;
end

if ~exist('lstil')
   lstil = '-';
end

for j=1:length(S)
   plot(S(j).x,S(j).y,'color',lcol,'linestyle',lstil,'linewidth',lwid);
   hold on;
end
