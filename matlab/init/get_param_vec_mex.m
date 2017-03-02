%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params_vec = get_param_vec_mex(params_mex)
%% 
%% INPUT:
%% params_mex = 
%%            SCATMOD: 3
%%            ADV_DIM: 2
%%            ADV_OPT: 2
%%      DO_CHECK_INIT: 1
%%      DO_CHECK_PROG: 1
%%     DO_CHECK_FINAL: 1
%%             STEADY: 1
%%            BRK_OPT: 1
%%           DO_ATTEN: 1
%%              young: 5.490000000000000e+09
%%            drag_rp: 13
%%            visc_ws: 0
%%           duration: 21600
%%                CFL: 0.700000000000000
%%            FSD_OPT: 1
%%         REF_Hs_ICE: 0
%%        USE_ICE_VEL: 0
%%     TAKE_MAX_WAVES: 0
%%            Hs_init: 3
%%             T_init: 12
%%           dir_init: -90
%%          conc_init: 0.700000000000000
%%             h_init: 1
%%          Dmax_init: 300
%%               Dmin: 20
%%                 xi: 2
%%          fragility: .9
%%            Dthresh: 200
%%           cice_min: 0.05
%%          model_day: 42003
%%      model_seconds: 0
%%              itest: 25
%%              jtest: 5
%%           dumpfreq: 10
%%            MEX_OPT: 1
%%            DO_DISP: 1
%%
%% OUTPUT:
%% vector with some of these fields
%% - for ordering see:
%%   - fortran/infiles/infile_nonIO.txt
%%   - read_params_vec subroutine in fortran/src/main/mod_WIM2d_run.F

%% ================================================
%% old int_prams:
fields   = {'SCATMOD',...
            'ADV_DIM',...
            'ADV_OPT',...
            'BRK_OPT',...
            'STEADY',...
            'DO_ATTEN',...
            'DO_CHECK_INIT',...
            'DO_CHECK_PROG',...
            'DO_CHECK_FINAL',...
            'young',...
            'drag_rp',...
            'visc_ws',...
            'cohesion',...
            'friction',...
            'CFL',...
            'FSD_OPT',...
            'REF_Hs_ICE',...
            'USE_ICE_VEL',...
            'TAKE_MAX_WAVES',...
            'Hs_init',...
            'T_init',...
            'dir_init',...
            'conc_init',...
            'h_init',...
            'Dmax_init',...
            'Dmin',...
            'xi',...
            'fragility',...
            'Dthresh',...
            'cice_min',...
            'duration',...
            'model_day',... % day relative to 1900-1-1
            'model_seconds',...
            'itest',...
            'jtest',...
            'dumpfreq'};
%% ================================================

if nargin==0
   %return list of fields needed by params_mex;
   fields(end+1:end+2)  = {'DO_DISP','MEX_OPT'};
   params_vec           = fields;
else
   Ni          = length(fields);
   params_vec  = zeros(Ni,1);
   for j=1:Ni
      cmd   = ['params_vec(',num2str(j),') = params_mex.',fields{j},';'];
      eval(cmd);
      if params_mex.DO_DISP==1
         disp([cmd(1:end-1),' = ',num2str(params_mex.(fields{j}))]);
      end
   end
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%