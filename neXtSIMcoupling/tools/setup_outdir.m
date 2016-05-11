function setup_outdir(test_i,op)

if nargin<2
   error('Need 2 inputs to setup_outdir');
end
outdir   = ['test_',num2str(test_i),'_outputs'];

if ~(strcmp(op,'mv')|strcmp(op,'cp'))
   error('3rd arg should be ''mv'' or ''cp''.')
end

eval(['!mkdir -p ',outdir]);
%%
odir  = [outdir,'/simul_in'];
eval(['!mkdir -p ',odir]);
cmd   = ['!',op,' *','simul_in*.mat ',odir];
eval(cmd);
%%
odir  = [outdir,'/simul_out_steps_mat'];
eval(['!mkdir -p ',odir]);
cmd   = ['!',op,' *','simul_out*step*.mat ',odir];
eval(cmd);
%%
odir  = [outdir,'/diagnostics'];
eval(['!mkdir -p ',odir]);
cmd   = ['!',op,' diagnostics* ',odir];
eval(cmd);
%%
odir  = [outdir,'/wim_log'];
eval(['!mkdir -p ',odir]);
cmd   = ['!',op,' test_outputs/out_2/log/* ',odir];
eval(cmd);
%%
odir  = [outdir,'/figs'];
eval(['!mkdir -p ',odir]);
odir  = [odir,'/init_final'];
eval(['!mkdir -p ',odir]);
cmd   = ['!',op,' *','test',num2str(test_i),'*.png ',odir];
eval(cmd);

%% move profile results
eval(['!mv *profile_results ',outdir,'/profile_results']);
