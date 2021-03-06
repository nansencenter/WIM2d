% function fn_Boltzmann_Steady
%
% INPUTS:
%
% fortyp      = 'freq' or 'wlength' or 'waveno'
% lam0        = value of forcing fortyp
% GeomDisk    = [x-location y-location Radius thickess] (coords not needed!)
% file_marker = string for file identifier; use 0 for no write
% tol         = tolerance on eigenvalues being real/imag
% fig         = figure handle
% col         = linecolor
%
% FLAGS:
%
% DTYP = flag for discretisation type: point-wise (0) or Fourier (1)
% COMM = flag for comments on (1) or off (0)
% PLOT = flag for plots: off (0), spectrum (1), sig wave height (2)
% SURGE = include surge motion
% RIGID = rigid disk (inf) or elastic disk (rigidity=10^RIGID)
% PS    = for PLOT = 1 -> pause(0.2) (PS=1) or pause (PS=0)
% ISO   = isotropic scattering (0=off, otherwise wavenumber in m^{-1})
% DO_SAVE = save data (1) or not (0)
%
% L Bennetts 2013 / Adelaide
% 
% Last updated:     March 2015

function out = fn_Boltzmann_Steady(fortyp, lam0, conc, wth, ISO, absorb, ...
 Param, outputs, COMM, PLOT, fig, col, DO_SAVE)

if ~exist('PLOT','var'); PLOT=2; end
if ~exist('COMM','var'); COMM=1; end
if ~exist('PS','var');   PS  =0; end
if ~exist('ISO','var');  ISO =0; end
if ~exist('DO_SAVE','var'); DO_SAVE=1; end
if DO_SAVE; sv_file='BltzSteady_saved'; end

if ~exist('outputs','var'); outputs = 'none'; end % 'transmitted energy'; end

if ~exist('th0','var'); th0=0.025; end

if ~exist('fig','var'); fig=1; end
if ~exist('col','var'); col='k'; end

%% Define test

if ~exist('absorb','var'); absorb=0; end %-1e-2; end % 

if ~exist('wth','var'); wth='inf'; end % 5e3; end % 

if ~exist('fortyp','var'); fortyp='freq'; end
if ~exist('lam0','var'); lam0=1/12; end

if ~exist('RIGID','var'); RIGID=5; end
if ~exist('SURGE','var'); SURGE=0; end

if ~exist('Vert_Modes','var'); Vert_Modes=1; end

if ~exist('conc','var'); conc=0.75; end 

if ~exist('Param','var'); Param = ParamDef_Default(RIGID);
 Param = ModParam_def(Param,1,Vert_Modes,0,0,100); end

%if ~exist('fn_inc','var'); fn_inc = 'kron_delta(0,th_vec)'; end
if ~exist('fn_inc','var'); fn_inc = 'cos(th_vec).^2'; end
%if ~exist('fn_inc','var'); fn_inc = 'series_delta(th_vec,25)'; end

if ~exist('normalisation_type','var')
   %% TODO make this an input
   %% this is only for ISO>0:
   %% - want to choose between scaling 'ISO' by c/A_floe (floe density)
   %%   (as done in non-iso scattering)
   %%   or not
   %% - default is not
   normalisation_type = 0;
end

% if absorb~=0
%  cprintf('red',['>>> Check solution for absorb~=0 ... exiting \n'])
%  return
% end

if PLOT; 
   if PLOT==1; x_res=11; else x_res=101; end
      if strcmp(wth,'inf')
         x=linspace(0,200*1e3,x_res);
      else
         x=linspace(0,wth,x_res);  
      end
   clear x_res;
end

%% Numerics

if ~exist('tol','var'); tol=1e-5; end

if ~exist('DTYP','var'); DTYP=0; end

th_res = Param.th_res;
   %th_res an integer, so always an odd no of bins
th_vec = linspace(0,1,2*th_res);%avoids 0.5 since odd no of bins and using end points
th_vec = unique([-th_vec,th_vec]);
th_vec(1)=[]; 
 
refs = find(or(th_vec>0.5,th_vec<-0.5));
incs = find(~or(th_vec>0.5,th_vec<-0.5));
th_vec = pi*th_vec;

if COMM
   cprintf([0.3,0.3,0.3],'<-------- Boltzmann steady ------->\n')
   cprintf([0.3,0.3,0.3],['>> ' fn_inc '\n'])
   if strcmp(wth,'inf')
      cprintf([0.3,0.3,0.3],['>> semi-inf problem \n'])
   else
      cprintf([0.3,0.3,0.3],['>> width = ' num2str(wth) '\n'])
   end
   cprintf([0.3,0.3,0.3],['>> ' num2str(100*conc) ' concentration \n'])
   if ~ISO
      if strcmp(fortyp,'freq')
         cprintf([0.3,0.3,0.3],['>> period = ' num2str(1/lam0) '\n'])
      else
         cprintf([0.3,0.3,0.3],['>> ' fortyp ' = ' num2str(lam0) '\n'])
      end
      cprintf([0.3,0.3,0.3],['>> floe diameter = ' num2str(Param.floe_diam) '\n'])
      cprintf([0.3,0.3,0.3],['>> ' num2str(Param.thickness) ' thick\n'])
      cprintf([0.3,0.3,0.3],['>> rigidity = ' sprintf('%0.5g',Param.E) '\n'])
      cprintf([0.3,0.3,0.3],['>>> Vertical modes = ' int2str(Param.Ndtm) '\n'])
   else%ISO
      cprintf([0.3,0.3,0.3],['>>> Isotropic scattering = ' num2str(ISO) 'm^{-1}\n'])
   end
   cprintf([0.3,0.3,0.3],['>>> angular resolution = ' int2str(th_res) '\n'])
end

if DO_SAVE
   clockout = clock; yr=clockout(1); mt=clockout(2); dy=clockout(3); 
   hr=clockout(4); mn=clockout(5); clear clockout
   what_prb = ['<-------- Boltzmann steady ------->\n' ...
    '>> ' fn_inc '\n'];
   what_prb = [what_prb '>> ' int2str(dy) '/' int2str(mt) '/' int2str(yr) ...
    ' ' int2str(hr) ':' int2str(mn) '\n']; clear yr mn dy hr mn
   if strcmp(wth,'inf')
      what_prb = [what_prb '>> semi-inf problem \n'];
   else
      what_prb = [what_prb '>> width = ' num2str(wth) '\n'];
   end
   if ~ISO
      what_prb = [what_prb '>> ' num2str(100*conc) ' concentration \n'];
      if strcmp(fortyp,'freq')
         what_prb = [what_prb '>> period = ' num2str(1/lam0) '\n'];
      else
         what_prb = [what_prb '>> ' fortyp ' = ' num2str(lam0) '\n'];
      end
      what_prb = [what_prb '>> floe diameter = ' num2str(Param.floe_diam) '\n'];
      what_prb = [what_prb '>> ' num2str(Param.thickness) ' thick\n'];
      what_prb = [what_prb '>> rigidity = ' sprintf('%0.5g',Param.E) '\n'];
      what_prb = [what_prb '>>> Vertical modes = ' int2str(Param.Ndtm) '\n'];
   else%ISO
      what_prb = [what_prb '>>> Isotropic scattering = ' num2str(ISO) '\n'];
   end
   what_prb = [what_prb '>>> angular reslution = ' int2str(th_res) '\n'];
   save(sv_file,'what_prb')
   clear what_prb
end % END IF DO_SAVE

if ISO==0
   out = fn_ElasticDisk(fortyp, lam0, Param, 'Energy', th_vec, ...
                          RIGID, SURGE, 0);
else
   out(1).name  = 'E';
   %out(1).value = ISO/length(th_vec) + 0*th_vec;
   out(1).value = ISO/2/pi + 0*th_vec;
     % NB multiply by dtheta=2*pi/Ndir later
     % => final R_mat1 = dtheta*(ISO/2/pi)*ones(Ndir,Ndir) = (ISO/Ndir)*ones(Ndir,Ndir)
end

for loop_out=1:length(out)
   if strcmp(out(loop_out).name,'E0')
      beta=out(loop_out).value;
   elseif strcmp(out(loop_out).name,'E')
      S=out(loop_out).value;
   end
end % end loop_out

%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Discrete system %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Boltzmann eqn: cos(theta)*dI/dx = -beta*I + int_{-pi}^{pi} S(th,th')*I(th') dth'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FINITE DIFFERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NB. This is the trapezium rule (because it's periodic)
dtheta = 2*pi/length(th_vec);
L_mat = cos(th_vec);
L_mat = diag(L_mat);

if ~ISO %%%%%%%%%%%%%%%%
   %%% IF NON_ISOTROPIC %%%
   %%%%%%%%%%%%%%%%%%%%%%%%
   R_mat1 = zeros(length(th_vec));
   mk = 2*th_res-1; %+1;
   for loop_th=1:2*th_res-2
      R_mat1(loop_th,:) = [S(mk:end),S(1:mk-1)];
      mk=mk-1;
   end
   R_mat1(2*th_res-1,:) = S; %+1
   mk=length(th_vec);
   for loop_th=2*th_res+0:length(th_vec)
      R_mat1(loop_th,:) = [S(mk:end),S(1:mk-1)];
      mk=mk-1;
   end
   clear mk
 
   %%% ISOTROPIC %%%
else %%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%
   R_mat1 = zeros(length(th_vec));
   for loop_th=1:length(th_vec)
      R_mat1(loop_th,:) = S;
   end
end

R_mat1 = dtheta*R_mat1;%integrate wrt theta
R_mat0 = -sum(R_mat1,1) + absorb;
R_mat = diag(R_mat0)+R_mat1;
if ISO==0 | normalisation_type==1
   %% scale by floe density, c/(pi*radius^2)
   %% - always do this for non iso scattering
   %% - optional for iso scattering
   sc_fac = conc/pi/((Param.floe_diam/2)^2);
   if ISO>0
      disp(['ISO = ',...
            num2str(ISO)]);
      disp(['scaled ISO = ',...
            num2str(ISO*sc_fac)]);
      %ISO, sc_fac*R_mat0(1), pause
   end
else
   sc_fac = 1;
   %ISO, sc_fac*R_mat0(1), pause
end
R_mat = sc_fac*R_mat;

if DO_SAVE
   K = sc_fac*R_mat1;
   alpha = sc_fac*R_mat0(1,1);
   
   save(sv_file,'th_vec','K','alpha','-append')
   
   clear K alpha
end % END DO_SAVE
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Eigenvalues & eigenvectors %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V,D] = eig(R_mat,L_mat);
D = diag(D);%convert from diagonal matrix to a vector

% min(abs(D))
% abs(det(V))

%plot(real(D),imag(D),'bx')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Boundary conditions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~DTYP
   if strcmp(wth,'inf')
      %%%%%%%%%%%%%%%%%%% 
      %%%  Inf width  %%%
      %%%%%%%%%%%%%%%%%%%
    
      %%%%%%%%%%%%%%%
      if absorb==0 %%
      %%%%%%%%%%%%%%%
      
         [Im,~,Iz] = fn_ArrangeEvals(D,0);
         
         if length([Im;Iz])~=length(incs)
            cprintf('red',['Check evals' '\n'])
         end
         
         % Rearrange
         V0 = V(incs,[Im;Iz]);
      
      %%%%%%%
      else %%
      %%%%%%% 
      
         Im = fn_ArrangeEvals(D,1);
         
         if length(Im)~=length(incs)
          cprintf('red',['Check evals' '\n'])
         end
         
         % Rearrange
         V0 = V(incs,Im);
      
      end % END IF absorb>0
      
      eval(['I0 = ' fn_inc '; I0(refs)=[];'])
      
      I0 = reshape(I0,length(I0),1);
      
      c0 = V0\I0;
      
      if absorb==0
         Vx = V(:,[Im;Iz]);
         Dx = D([Im;Iz]);
      else
         Vx = V(:,Im);
         Dx = D([Im]);
      end
    
   %%%%%%%%%%%%%%%%%%%% 
   %%% Finite width %%%
   %%%%%%%%%%%%%%%%%%%%
   else 
       
      %%% Incident wave fields
      eval(['I0 = ' fn_inc '; I0(refs)=0;'])
      eval(['I1 = 0*' fn_inc '; I1(incs)=0;'])
      Ibc = I0+I1;
      Ibc = reshape(Ibc,length(Ibc),1);
      
      %%%%%%%%%%%%%%%
      if absorb==0 %%
      %%%%%%%%%%%%%%% 
        
         %%% S(x)=b0*(x*v0+u0) + c0*v0 + Sum c+(n)*exp(D+(n)*x)*v+ + Sum c-(n)*exp(D-(n)*(x-w))*v-
         %%%                              n                           n
         %%%
         %%% where D+=[Ip] and D-=[Im]
         %%%
         %%% v0 regular eigenvector corresponding to 0 eval (multiplicity 2)
         %%% u0 generalised evec:
         %%%                      (R_mat-0*L_mat)*u0=L_mat*v0
         %%%
         %%% SL(0)=S0 and SR(w)=S1
         %%%
         %%% where SL=S[incs] and SR=S[refs] 
         
         [Im,Ip,Iz0,Iz1] = fn_ArrangeEvals(D,0);
         
         n0=length([Im;Iz0]); n1=length([Ip;Iz1]);
         
         if or(n0~=length(incs),n1~=length(refs))
            cprintf('magenta',['>>> Check evals: ' fortyp '=' num2str(lam0) '\n'])
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%% Generalised evalue for lam=0 %%%
         [VR,DR]=eig(R_mat); DR=diag(DR);
         
         % Find zero eval
         [~,i0]=min(abs(DR)); 
         
         % Solve R_mat*u0=L_mat*v0 => DR*(VR\u0)=VR\L_mat*v0;
         %
         % Pick 1st zero eval of R_mat-lam*L_mat (arbitrary)
         
         v0 = V(:,Iz0);
         
         u0=0*v0;
         
         dum_v = VR\L_mat*v0;
         
         for lp=1:length(u0)
            if lp~=i0
               u0(lp)=dum_v(lp)/DR(lp);
            end
         end
         
         u0 = VR*u0;
         
         %plot(th_vec,R_mat*u0-L_mat*v0)
         
         clear dum_v VR DR
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         %%% Rearrange
         V = V(:,[Iz0,Im.',Ip.']);
         D0 = D([Im]); D1 = D([Ip]);
           
         Vbc = zeros(length(th_vec));
         %%% BCs for x=0
         EDbc = [1;1;exp(D0*0);exp(-D1*wth)];
         Vbc(incs,:) = [u0(incs),V(incs,:)]*diag(EDbc); clear EDbc
         
         %%% BCs for x=wth
         EDbc = [1;1;exp(D0*wth);exp(D1*0)];
         Vbc(refs,:) = [wth*v0(refs)+u0(refs),V(refs,:)]*diag(EDbc); clear EDbc
         
         c = Vbc\Ibc; clear Vbc Ibc
      
      %%%%%%%
      else %%
      %%%%%%% 
      
         %%% S(x)=Sum c+(n)*exp(D+(n)*x)*v+ + Sum c-(n)*exp(D-(n)*(x-w))*v-
         %%%       n                             n
         %%%
         %%% where D+=[Im] and D-=[Ip]
         %%%
         %%% SL(0)=S0 and SR(w)=S1
         %%%
         %%% where SL=S[incs] and SR=S[refs] 
         
         [Im,Ip] = fn_ArrangeEvals(D,1);
         
         n0=length([Im]); n1=length([Ip]);
         
         if or(n0~=length(incs),n1~=length(refs))
            cprintf('magenta',['>>> Check evals: ' fortyp '=' num2str(lam0) '\n'])
         end
         
         %%% Rearrange
         V = [V(:,[Im]),V(:,[Ip])];
         D0 = D([Im]); D1 = D([Ip]);
           
         Vbc = zeros(length(th_vec));
         %%% BCs for x=0
         EDbc = [exp(D0*0);exp(-D1*wth)];
         Vbc(incs,:) = V(incs,:)*diag(EDbc); clear EDbc
         
         %%% BCs for x=wth
         EDbc = [exp(D0*wth);exp(D1*0)];
         Vbc(refs,:) = V(refs,:)*diag(EDbc); clear EDbc
         
         c = Vbc\Ibc; clear Vbc Ibc
      
      end % END IF absorb>0
      
   end % END IF wth~=inf
else
   cprintf('red',['Not coded yet' '\n'])
end

%%%%%%%%%%%%%%
%% 4. Plot? %%
%%%%%%%%%%%%%%

if PLOT
   if ~DTYP
      
      figure(fig); hold on %h1 = subplot(1,1,1);
      %%%%%%%%%%%%%%%%%%% 
      %%%  Inf width  %%%
      %%%%%%%%%%%%%%%%%%% 
      if strcmp(wth,'inf')
       
         %%% PLOT I(theta) for increasing values of x
         if PLOT == 1
         
            %%%%%%%%%%%%%%%
            if absorb==0 %%
            %%%%%%%%%%%%%%% 
            
               for loop_x=1:length(x)
                  I = Vx*diag(exp(D([Im;Iz])*x(loop_x)))*c0;
                  if ~isempty(find(abs(imag(I))>tol))
                     cprintf('red',['Check I: ', int2str(loop_x), ', ', ...
                        num2str(max(abs(imag(I)))), '\n'])
                  end
                  I = [I(end);I];
                  plot(gca,[-pi,th_vec]/pi,real(I))
                  set(gca, 'ylim',[0,max(I0)],'xlim',[-1,1])
                  title(['x=' num2str(x(loop_x))])
                  xlabel('\theta','fontsize',16); ylabel('I','fontsize',16);
                  if PS==0
                     cprintf('m',['>> paused: hit any key to continue \n'])
                     pause
                  else
                     pause(PS)
                  end
               end
               close(gcf)
               if COMM; cprintf('blue',['>> const = ' num2str(I(1)) '\n']); end
            
            %%%%%%%
            else %% absorb~=0
            %%%%%%%
            
               for loop_x=1:length(x)
                  I = Vx*diag(exp(D(Im)*x(loop_x)))*c0;
                  if ~isempty(find(abs(imag(I))>tol))
                     cprintf('red',['Check I: ', int2str(loop_x), ', ', ...
                     num2str(max(abs(imag(I)))), '\n'])
                  end
                  I = [I(end);I];
                  plot(gca,[-pi,th_vec]/pi,real(I))
                  set(gca, 'ylim',[0,max(I0)],'xlim',[-1,1])
                  title(['x=' num2str(x(loop_x))])
                  xlabel('\theta','fontsize',16); ylabel('I','fontsize',16);
                  if PS==0
                     cprintf('m',['>> paused: hit any key to continue \n'])
                     pause
                  else
                     pause(PS)
                  end
               end
               close(gcf)
               if COMM; cprintf('blue',['>> const = ' num2str(I(1)) '\n']); end
            
            end % END IF absorb ==0
         
         %%% PLOT Hs(x)
         %%%%%%%%%%%%%%%%%%%
         elseif PLOT == 2 %%
         %%%%%%%%%%%%%%%%%%% 
         
            f_inline = inline(fn_inc);
            %m0_inc = quad(fn_inc,-pi/2,pi/2)
            dtheta = th_vec(2)-th_vec(1);
            m0_inc = sum(f_inline(th_vec(incs)))*dtheta;
            Hs_inc = 4*sqrt(m0_inc);
            %%%%%%%%%%%%%%%
            if absorb==0 %%
            %%%%%%%%%%%%%%% 
               H_vec = 0*x;
               for loop_x=1:length(x)
                  I = Vx*diag(exp(D([Im;Iz])*x(loop_x)))*c0;
                  H_vec(loop_x) = real(4*sqrt(sum(I*dtheta))/Hs_inc);
                  clear I;
               end
               plot(x/1e3,H_vec,col); set(gca,'box','on')
               xlabel('x [km]','fontsize',14); ylabel('H_{s}(x)/H_{s}(0)','fontsize',14);
               if DO_SAVE
                  save(sv_file,'x','H_vec','-append')
               end % END DO_SAVE    
            %%%%%%%
            else %% absorb~=0
            %%%%%%% 
               H_vec = 0*x;
               for loop_x=1:length(x)
                  I = Vx*diag(exp(D([Im])*x(loop_x)))*c0;
                  H_vec(loop_x) = real(4*sqrt(sum(I*dtheta))/Hs_inc);
                  clear I;
               end
               plot(x/1e3,H_vec,col);
               set(gca,'box','on')
               xlabel('x [km]','fontsize',14); ylabel('H_{s}(x)/H_{s}(0)','fontsize',14);
               if DO_SAVE
                  save(sv_file,'x','H_vec','-append')
               end % END DO_SAVE
            end % END IF absorb==0
         end % END if PLOT==
       
      %%%%%%%%%%%%%%%%%%%% 
      %%% Finite width %%%
      %%%%%%%%%%%%%%%%%%%%
      
      else
         %%% PLOT I(theta) for increasing values of x
         if PLOT==1
            %%%%%%%%%%%%%%% 
            if absorb==0 %% 
               %%%%%%%%%%%%%%%
               for loop_x=1:length(x)
                  I = [x(loop_x)*v0+u0,V]*diag([1;1;exp(D0*x(loop_x));exp(D1*(x(loop_x)-wth))])*c;
                  if ~isempty(find(abs(imag(I))>tol))
                     cprintf('red',['Check I: ', int2str(loop_x), ', ', ...
                     num2str(max(abs(imag(I)))), '\n'])
                  end
                  I = [I(end);I];
                  plot(gca,[-pi,th_vec]/pi,real(I))
                  set(gca, 'ylim',[0,max(I0)],'xlim',[-1,1])
                  title(['x=' num2str(x(loop_x))])
                  xlabel('\theta','fontsize',16); ylabel('I','fontsize',16);
                  %set(gca,'yscale','log')
                  if PS==0
                     cprintf('m',['>> paused: hit any key to continue \n'])
                     pause
                  else
                     pause(PS)
                  end
               end
            %%%%%%%
            else %%
            %%%%%%% 
               for loop_x=1:length(x)
                  I = V*diag([exp(D0*x(loop_x));exp(D1*(x(loop_x)-wth))])*c;
                  if ~isempty(find(abs(imag(I))>tol))
                     cprintf('red',['Check I: ', int2str(loop_x), ', ', ...
                     num2str(max(abs(imag(I)))), '\n'])
                  end
                  I = [I(end);I];
                  plot(gca,[-pi,th_vec]/pi,real(I))
                  set(gca, 'ylim',[0,max(I0)],'xlim',[-1,1])
                  title(['x=' num2str(x(loop_x))])
                  xlabel('\theta','fontsize',16); ylabel('I','fontsize',16);
                  %set(gca,'yscale','log')
                  if PS==0
                     cprintf('m',['>> paused: hit any key to continue \n'])
                     pause
                  else
                     pause(PS)
                  end
               end
            
            end % END IF absorb==0
            
            close(gcf)
            
            % Plot Hs = 4*sqrt(m0) where m0 = int I(theta) dtheta
            % NB. actually plot ratio Hs(x)/Hs_inc
         elseif PLOT == 2
            %%%%%%%%%%%%%%%
            f_inline = inline(fn_inc);
            %m0_inc = quad(fn_inc,-pi/2,pi/2)
            dtheta = th_vec(2)-th_vec(1);
            m0_inc = sum(f_inline(th_vec(incs)))*dtheta;
            Hs_inc = 4*sqrt(m0_inc);

            if absorb==0 %%
            %%%%%%%%%%%%%%%
               H_vec = 0*x; 
               for loop_x=1:length(x)
                  I = [x(loop_x)*v0+u0,V]*diag([1;1;exp(D0*x(loop_x));exp(D1*(x(loop_x)-wth))])*c;
                  H_vec(loop_x) = real(4*sqrt(sum(I*dtheta)))/Hs_inc;
                  clear I
               end
               plot(x/1e3,H_vec,col); set(gca,'box','on')
               xlabel('x [km]','fontsize',14); 
               %ylabel('H_{s}(x)/H_{s}(0)','fontsize',14);
               ylabel('H_{s}(x)/H_{s,inc}','fontsize',14);
            %%%%%%%
            else %%
            %%%%%%%
          
               H_vec = 0*x;
               for loop_x=1:length(x)
                  I = V*diag([exp(D0*x(loop_x));exp(D1*(x(loop_x)-wth))])*c;
                  H_vec(loop_x) = real(4*sqrt(sum(I*dtheta)))/Hs_inc;
                  clear I;
               end
               plot(x/1e3,H_vec,col); set(gca,'box','on')
               xlabel('x [km]','fontsize',14); ylabel('H_{s}(x)/H_{s}(0)','fontsize',14);
            end % END IF absorb==0

         end % END IF PLOT==
      end % ENF IF wth==inf
     
   else
      cprintf('red',['Not coded yet' '\n'])
   end  
   hold off
end%%PLOT

%%%%%%%%%%%%%
%% Outputs %%
%%%%%%%%%%%%%

jj0 = find(~or(th_vec(incs)>pi*th0,th_vec(incs)<-pi*th0));

if strcmp(wth,'inf')
   %%%%%%%%%%%%%%%%%%% 
   %%%  Inf width  %%%
   %%%%%%%%%%%%%%%%%%% 
   eigen_info   = struct('V',Vx,...
                         'D0',Dx,...
                         'coeffs',c0,...
                         'width',wth,...
                         'absorb',absorb,...
                         'angles',th_vec);

   %%%%%%%%%%%%%%%
   if absorb==0 %%
      %%%%%%%%%%%%%%%
      I0 = (V(:,[Im;Iz])*c0).';
      I  = (V(:,[Im;Iz])*diag(exp(D([Im;Iz])*Param.MIZ_length))*c0).';
      if ~isempty(find(abs(imag([I,I0]))>tol))
         cprintf('red',['Check I and I0: ', int2str(loop_x), ', ', ...
         num2str(max(abs(imag(I)))), '\n'])
      end
      j0 = find(th_vec(incs)==0);
      T0 = real(I(incs(j0)))/real(I0(incs(j0)));
      Tx = sum(cos(th_vec(incs)).*real(I(incs)))/...
            sum(cos(th_vec(incs)).*real(I0(incs)));
      Tf = sum(real(I(incs)))/sum(real(I0(incs)));
      TN = sum(real(I(incs(jj0))))/sum(real(I0(incs(jj0))));
      %%%%%%%
   else %% absorb~=0
      %%%%%%%
      I0 = (V(:,Im)*c0).';
      I  = (V(:,Im)*diag(exp(D(Im)*Param.MIZ_length))*c0).';
      if ~isempty(find(abs(imag([I,I0]))>tol))
         cprintf('red',['Check I and I0: ', int2str(loop_x), ', ', ...
         num2str(max(abs(imag(I)))), '\n'])
      end
      j0 = find(th_vec(incs)==0);
      T0 = real(I(incs(j0)))/real(I0(incs(j0)));
      Tx = sum(cos(th_vec(incs)).*real(I(incs)))/...
            sum(cos(th_vec(incs)).*real(I0(incs)));
      Tf = sum(real(I(incs)))/sum(real(I0(incs)));
      TN = sum(real(I(incs(jj0))))/sum(real(I0(incs(jj0)))); 
   end % 
else
   %%%%%%%%%%%%%%%%%%%% 
   %%% Finite width %%%
   %%%%%%%%%%%%%%%%%%%% 

   %%%%%%%%%%%%%%%
   if absorb==0 %%
      %%%%%%%%%%%%%%%
      eigen_info = struct('u0',u0,...
                          'v0',v0,...
                          'V',V,...
                          'D0',D0,...
                          'D1',D1,...
                          'coeffs',c,...
                          'width',wth,...
                          'absorb',absorb,...
                          'angles',th_vec);

      I0 = [0*v0(incs)+u0(incs),V(incs,:)]*...
            diag([1;1;exp(D0*0);exp(D1*(0-wth))])*c;
      I  = [wth*v0(incs)+u0(incs),V(incs,:)]*diag([1;1;exp(D0*wth);...
      exp(D1*(wth-wth))])*c;
      if ~isempty(find(abs(imag([I;I0]))>tol))
         cprintf('red',['Check I and I0: ', int2str(loop_x), ', ', ...
         num2str(max(abs(imag(I)))), '\n'])
      end
      j0 = find(th_vec(incs)==0);
      T0 = real(I(j0))/real(I0(j0));
      Tx = sum(cos(th_vec(incs)).*real(I.'))/...
      sum(cos(th_vec(incs)).*real(I0.'));
      Tf = sum(real(I.'))/sum(real(I0.'));
      TN = sum(real(I(jj0).'))/sum(real(I0(jj0).'));

   %%%%%%%
   else %%
      %%%%%%% 
      eigen_info   = struct('u0',NaN,...
                            'v0',NaN,...
                            'V',V,...
                            'D0',D0,...
                            'D1',D1,...
                            'coeffs',c,...
                            'width',wth,...
                            'absorb',absorb,...
                            'angles_on_pi',th_vec);

      I0 = V(incs,:)*diag([exp(D0*0);exp(D1*(0-wth))])*c;
      I  = V(incs,:)*diag([exp(D0*Param.MIZ_length);...
                           exp(D1*(Param.MIZ_length-wth))])*c;
      if ~isempty(find(abs(imag([I;I0]))>tol))
         cprintf('red',['Check I and I0: ', int2str(loop_x), ', ', ...
         num2str(max(abs(imag(I)))), '\n'])
      end
      j0 = find(th_vec(incs)==0);
      T0 = real(I(j0))/real(I0(j0));
      Tx = sum(cos(th_vec(incs)).*real(I.'))/...
            sum(cos(th_vec(incs)).*real(I0.'));
      Tf = sum(real(I.'))/sum(real(I0.'));
      TN = sum(real(I(jj0).'))/sum(real(I0(jj0).'));
   end % END IF absorb==0
end % END IF wth~='inf'

out_str = ' ''dummy'' '; out_val = ' 0 ';

if strfind(outputs,'full trans energy')
 out_str = [out_str '; ''transmitted energy full'' '];
 out_val = [out_val '; Tf'];
end

if strfind(outputs,'trans X energy')
 out_str = [out_str '; ''transmitted energy x cos(\theta)'' '];
 out_val = [out_val '; Tx'];
end

if strfind(outputs,'transmitted energy')
 out_str = [out_str '; ''transmitted energy'' '];
 out_val = [out_val '; T0'];
end

if strfind(outputs,'int trans energy')
 out_str = [out_str '; ''transmitted energy (|\theta|<' num2str(th0) '\pi )'' '];
 out_val = [out_val '; TN'];
end

if strfind(outputs,'eigen-info')
 out_str = [out_str '; ''eigen-info: structure with eigenvalues, eigenvectors and coefficients of eigenfunction expansion'' '];
 out_val = [out_val '; eigen_info'];
end

eval(['out=struct( ''name'', {' out_str ...
 '}, ''value'', {' out_val '});'])
out(1)=[];

%  incs = find(~or(th_vec>0.5*pi,th_vec<0));
%
%  I = V(:,[Im;Iz])*diag(exp(D([Im;Iz])*Param.MIZ_length))*c0;
%  I = I(incs);
%
%  if max(abs(imag(I)))>tol
%   cprintf('r',['check solution: ' num2str(max(abs(imag(I)))) '\n'])
%  end
%
%  I = real(I);
%
%  th_vec = th_vec(incs);

if COMM
 cprintf([0.3,0.3,0.3],'<----- End: Boltzmann steady ----->\n')
end
 
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [Ip,Ii,Iz] = fn_ArrangeEvals(D,tol)
%
% OUTPUT
%
% Ip = positive real evals
% Im = positive real evals
% Ii = purely imaginary evals (arranged into conj pairs)
% Iz = zero evals

function [Im,Ip,IzA,IzB] = fn_ArrangeEvals(D,vers)

% Iz = find(abs(D)<tol);
% Im = find(and(~abs(D)<tol,real(D)<-tol));
% Ip = find(and(~abs(D)<tol,real(D)>tol));
% 
% if length(Iz) ~= 2
%  cprintf('magenta',['>>> Check zero evals' '\n'])
% end
% 
% IzA = Iz(1); IzB = Iz(2);
% 
% if ~isempty(find(abs(imag(D))>tol))
%  cprintf('blue',['>>> nb. imag component to evals: ' ...
%   num2str(max(abs(imag(D)))) '\n'])
% end

if ~exist('vers','var'); vers=0; end

if vers==0 % When repeated 0 eval expected
 D0=D;
 [~,IzA] = min(abs(D0));
 D0(IzA) = nan;
 [~,IzB] = min(abs(D0));
 D0(IzB) = nan;
 Ip = find(real(D0)>-imag(D0));
 Im = find(real(D0)<-imag(D0));
else % no 0 evals expected
 D0=D;
 Ip = find(real(D0)>-imag(D0));
 Im = find(real(D0)<-imag(D0));
 IzA=nan; IzB=nan;
end
return

% Kronecker Delta

function out=kron_delta(u,v)

out=zeros(1,length(v));

for loop=1:length(v)
    
 if u==v(loop)
  out(loop)=1;
 end
 
end

% sum_{n=-N}^{N} exp{i*n*theta}
% ->delta(theta-0) as N->infty

function out=series_delta(theta,N)

out=ones(1,length(theta));

for lp=1:N
    
 out = out + 2*cos(lp*theta);
 
end

return

%

% function I = fn_TrapRule(xx,yy)
% 
% xtil=diff(xx);
% ytil=0.5*(yy(1:end-1)+yy(2:end));
% 
% xtil = reshape(xtil,1,length(xtil));
% ytil = reshape(ytil,1,length(ytil));
% 
% I = sum(xtil.*ytil);
% 
% return
