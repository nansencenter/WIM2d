function visc_rp  = RPfit_gamma()


h     = .5;
cice  = 0.5;
d_av  = 50;
hd    = [h d_av];

[Tsm,att,err,floe_density]...
         = atten_SqMo1980(cice);
asm      = att/floe_density;
err_nd   = err/floe_density;
%%
[E,g,rho_wtr,rho_ice,nu]   = NDphyspram_v2;
prams_rp                   = [E,nu,rho_wtr,rho_ice,g];
%%
om          = 2*pi./Tsm;
a_undamped  = ALPfxn_E549_cheb2(om,h);

%% function to be minimised:
%%       Y  = fn_error(vrp2,Tsm,asm,prams_rp,hd,a_undamped)
vrp_guess   = 20;
visc_rp     = fminsearch( @(vrp)fn_error(vrp,Tsm,asm,prams_rp,hd,a_undamped),...
                          vrp_guess);

%%
if nargout==1
   return;
   %%otherwise do some plots
end

nplots   = 1;
if nplots==2
   subplot(1,2,1);
   vrp2  = linspace(10,25,30)';
   Y     = fn_error(vrp2,Tsm,asm,prams_rp,hd,a_undamped);
   plot(vrp2,Y,'k');
   hold on;
   plot(visc_rp,...
      fn_error(visc_rp,Tsm,asm,prams_rp,hd,a_undamped),...
      'ok');
   GEN_proc_fig('Gamma','Error squared, m^{-2}');
   hold off;
end
%%
subplot(1,nplots,nplots);
T           = linspace(5,15,50)';
om          = 2*pi./T;
a_undamped  = ALPfxn_E549_cheb2(om,h);
%%
[lam,q_damping]   = RPget_lam_dmpg(h,om,prams_rp,visc_rp);
a_damped          = a_undamped+2*d_av*q_damping;

[lam,q_damping]   = RPget_lam_dmpg(h,om,prams_rp,1.5*visc_rp);
a_damping         = 2*d_av*q_damping;
%%
PLOT_DIM = 1;
if PLOT_DIM
   fac   = floe_density;
   ylab  = '\alpha, m^{-1}';
   ax    = [min(T) max(T) 0 5e-4];
else
   fac   = 1;
   ylab  = '\alpha';
   ax    = [min(T) max(T) 0 0.02];
end
h1 = plot(T,fac*a_undamped,'r','linewidth',2);
hold on;
h2 = plot(T,fac*a_damped,'b','linewidth',2);
h3 = plot(T,fac*a_damping,'g','linewidth',2);
%plot(Tsm,asm,'^k');
h4 = errorbar(Tsm,fac*asm,fac*err_nd,'ok','linewidth',2,'markersize',2);
legend([h1,h2,h3],{'Scattering','Scattering + drag','Drag'});
%%
set(gca,'yscale','log');
axis(ax);
%set(gca,'xscale','log','yscale','log');
GEN_proc_fig('Period, s',ylab);
hold off;

!mkdir -p out
saveas(gcf,'out/RPfit.eps','epsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y  = fn_error(vrp2,Tsm,asm,prams_rp,hd,a_undamped)
h     = hd(1);
d_av  = hd(2);
Y     = 0*vrp2;
%%
for r=1:length(vrp2)
   visc_rp  = vrp2(r);
   %%
   for j=1:3%length(Tsm)
      om = 2*pi/Tsm(j);
      [lam,q_damping]   = RPget_lam_dmpg(h,om,prams_rp,visc_rp);
      a_damped          = a_undamped(j)+2*d_av*q_damping;
      %%
      Y(r)  = Y(r)+(a_damped-asm(j))^2;
   end
end
