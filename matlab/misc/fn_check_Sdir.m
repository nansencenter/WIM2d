function fn_check_Sdir(Sf)

if nargin==0
   w2d   = getenv('WIM2D_PATH');
   Sfil  = [w2d,'/matlab/main/m_out/Sdir.mat'];
   Sf    = load(Sfil);
end

wavdir      = Sf.wavdir;
freq_vec    = Sf.freq_vec;
Sdir        = Sf.Sdir;
grid_prams  = Sf.grid_prams;
S_inc       = Sf.S_inc;
cice        = Sf.cice;
Hs          = Sf.Hs;

w_test   = 1;  %% frequency index
j_test   = 1;  %% y index
disp(' ');
disp(['testing period (s)  : ',num2str(1/freq_vec(w_test))]);
disp(' ');

dx       = grid_prams.dx;
je       = find(cice>0,1,'first');
xe       = grid_prams.X(je,1)-.5*dx;
[xx,cc]  = step_1d(grid_prams.X(:,1)-xe,cice(:,1));
[xx,HH]  = step_1d(grid_prams.X(:,1)-xe,Hs(:,1));


%% sort directions (splitting fwd/back):
th_vec   = -pi/180*(wavdir+90);%%-90 deg (from W) -> 0 radians; 0 deg (from N) -> -pi/2 radians
dth      = th_vec(2)-th_vec(1);
jp       = find(cos(th_vec)>=0);
jm       = find(cos(th_vec)<=0);
%%
A_inc    = 1.5;
S0_exact = A_inc^2/4*(2/pi)*cos(th_vec(jp)).^2;
if 1
   figure(2);
   % subplot(2,1,1);
   % plot(xx/1e3,cc);
   % ylim([0,1]);
   % xl    = xlabel('x, km');
   % yl    = ylabel('Concentration');
   % set(xl ,'fontname','times','fontsize',20);
   % set(yl ,'fontname','times','fontsize',20);
   % set(gca,'fontname','times','fontsize',20);
   % %%
   % subplot(2,1,2);
   plot(xx/1e3,HH/(2*A_inc));
   xl    = xlabel('x, km');
   yl    = ylabel('H_s/H_{s,inc}');
   set(xl ,'fontname','times','fontsize',20);
   set(yl ,'fontname','times','fontsize',20);
   set(gca,'fontname','times','fontsize',20);
   %pause;
   pause(1);
end
%%
ndir  = length(wavdir);
nw    = length(freq_vec);

for i_test=je-20:600

   if nw>1
      S_ijw    = squeeze(Sdir(i_test,j_test,:,wtest));
      S0_ijw   = squeeze(S_inc(1,j_test,:,wtest));
   else
      S_ijw    = squeeze(Sdir(i_test,j_test,:));
      S0_ijw   = squeeze(S_inc(1,j_test,:));
   end

   S_fwd       = 0*S_ijw;
   S_fwd(jp)   = S_ijw(jp);

   figure(1);
   plot(wavdir+180,S_ijw,'linewidth',2);
   hold on;
   plot(wavdir+180,S_fwd,'--r','linewidth',2);
   plot(wavdir+180,S0_ijw,'k');
   plot(wavdir(jp)+180,S0_exact,'--c','linewidth',2);
   hold off;
   leg   = legend('numeric (total)','numeric (fwd)','incident','exact incident');
   xlim([-90,270]);
   set(gca,'xtick',-90:45:270);
   %%
   xt          = grid_prams.X(i_test,1)-xe;
   Sfreq_test  = dth*[sum(S_ijw),sum(S_fwd),sum(S0_ijw)]

   ttl   = title(['x (km) : ',num2str(xt/1e3),'; Hs (m) :',num2str(Hs(i_test,1))]);
   xl    = xlabel('Wave-to direction, degrees');
   yl    = ylabel('S(\theta), m^2s^{-1}');
   set(xl ,'fontname','times','fontsize',20);
   set(yl ,'fontname','times','fontsize',20);
   set(gca,'fontname','times','fontsize',20);
   set(leg,'fontname','times','fontsize',20);
   if cice(i_test,1)==0
      set(ttl,'fontname','times','fontsize',24,'color','b');
   else
      set(ttl,'fontname','times','fontsize',36,'color','r');
   end

   drawnow;
   pause(.01);
   %pause;

end
