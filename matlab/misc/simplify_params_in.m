%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params   = simplify_params_in(params_in)

%%check the minimum set of fields is given,
%% and return the minimum set;
list = check_params_in(params_in);

for j=1:length(list)
   params.(list{j})  = params_in.(list{j});
end
