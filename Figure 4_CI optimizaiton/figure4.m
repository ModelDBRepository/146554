%Reproduces Figure 4 in Sedigh-Sarvestani, Schiff, Gluckman,
%'Reconstructing Mammalian Sleep Dynamics with Data Assimiliation',
%PLoS Comp Biol, 2012 (In press). DOI:10.1371/journal.pcbi.1002788

%Reconstructs unobserved dynamics from DB model of sleep.
%DB model from: Diniz Behn and Booth, J Neurophysiol 103:1937-1953, 2010.

%Then calculates the error in reconstruction of hidden variables, when each
%model variable is used as observable, then represents these errors as a
%matrix (the EOC).

%Then determines the optimal covariance inflation (CI) parameter for
%variables of the DB model, that optimizes reconstruction performance.

%Usage: Running this .m file will produce a figure similar to Figure 4 of
%Sedigh-Sarvestani et al. 
%Make sure CD is '...\Figure Code\Figure 4_CI optimization'

%Options: You can use previously generated model data, from which noisy
%observations are obtained. This is stored as 'data_DB_output.mat' in
%'\Figure Code\Generate Data\' and is used throughout the paper.
%Alternatively, you can re-generate model data from scratch.

%Madineh Sedigh-Sarvestani, Penn State, Oct 2012
%m.sedigh.sarvestani@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]= figure4
cd('../')%go up one folder
addpath(genpath(cd)); %add path
load CI_eps.mat %this is the default CI (covariance inflation) multiplier for all variables used throughout the figures
cd([cd '/Figure 4_CI optimizaiton']) %reset folder back 

if exist('data_DB_output.mat','file')
    load data_DB_output.mat %load already generated data
else
    [Times,x,y,state,dT,P,Rs]= data_DB; %alternatively re-generate data
end

%first optimize CI for delta
optimize_CI_delta(x,y,dT,CI_eps,Rs);%fig4D
%then optimize CI for F_R, given optimized CI for delta
optimize_CI_R(x,y,dT,CI_eps,Rs);%fig4E
%then compute EOC for default CI
Full_EOC_default_CI(x,y,dT,CI_eps,Rs); %fig4A
%then compute EOC for optimized CI_delta
Full_EOC_optimized_CI_delta(x,y,dT,CI_eps,Rs)%fig4B
%then compute EOC for optimized CI_delta and CI_F_R
Full_EOC_optimized_CI_delta_CI_rem(x,y,dT,CI_eps,Rs);%fig4C

%then plot all together
plot_figure;
return

%this function optimizes CI_delta when observable is F_LC
function[]=optimize_CI_delta(x,y,dT,CI_eps,Rs)
%outputs
%inputs:
%x(matrix of DB model data)
%y(matrix of DB model observations-noise added version of x)
%dT(sampling,RK4 integration time)
%CI_eps (multiplier (by variance of each variable) for covariance inflation
%matrix)
%Rs (measurement noise covariance)

%setup default CI matrix
vars=var(x');
CI=CI_eps.*blkdiag(vars(1),vars(2),vars(3),vars(4),vars(5), ...
        vars(6),vars(7),vars(8),vars(9),vars(10),vars(11),vars(12));

CI_delta_opt=CI;

%choose observable
varn=1; %observe F_LC
y=y(1,:); %observe F_LC
dy=length(varn);
R=(Rs([varn])'*(ones(1,dy)))'.*eye(dy,dy); %measurement  covariance matrix

%get a series of CI multiplier coefficients to try
tt=[-3 0 3 4 6 7 9 10 10.5 11 12 12.5 13 14 16];
CI_coeff=[2.^tt]; 

%preallocate xhat
[dx,N]=size(x);
xhat=zeros(dx,N,length(CI_coeff));

%go through and reconstruct hidden variables using different CI's defined
%above
for j=1:length(CI_coeff)
    j
    CI_delta_opt(12,12)=CI(12,12)*CI_coeff(j);
    [xhat(:,:,j)]=UKF(y,R,dT,varn,CI_delta_opt); %reconstruct
end

%calculate EOC
for j=1:length(CI_coeff)
    for i=[1:dx]
        chisq_all(i,j)=mean((x(i,N/2:end)-xhat(i,N/2:end,j)).^2,2)/vars(i); %MSE
    end
    EOC_delta(j,:)=1./(chisq_all(:,j)+1); %EOC
end

save figure4D.mat
return

function []=optimize_CI_R(x,y,dT,CI_eps,Rs)
%outputs
%inputs:
%x(matrix of DB model data)
%y(matrix of DB model observations-noise added version of x)
%dT(sampling,RK4 integration time)
%CI_eps (multiplier (by variance of each variable) for covariance inflation
%matrix)
%Rs (measurement noise covariance)

vars=var(x');
CI=CI_eps.*blkdiag(vars(1),vars(2),vars(3),vars(4),vars(5), ...
        vars(6),vars(7),vars(8),vars(9),vars(10),vars(11),vars(12)*4000);

CI_rem_opt=CI;

%     
 varn=1; %observe F_LC
 y=y(1,:); %observe F_LC
 dy=length(varn);
 R=(Rs([varn])'*(ones(1,dy)))'.*eye(dy,dy); %measurement  covariance matrix
% 
%get a series of CI multipliers
tt=[-3 0 3 4 6 7 9 10 10.5 11 11.5 12 12.5 13 14 16];
CI_coeff=[2.^tt]; 
CI_coeff=[CI_coeff];

% %preallocate xhat
[dx,N]=size(x);
xhat=zeros(dx,N,length(CI_coeff));

%go through and reconstruct with all different CIs
for j=1:length(CI_coeff)
    CI_rem_opt(4,4)=CI(4,4)*CI_coeff(j);
    [xhat(:,:,j)]=UKF(y,R,dT,varn,CI_rem_opt);
end

%calculate EOC
for j=1:length(CI_coeff)
    for i=[1:dx]
        chisq_all(i,j)=mean((x(i,N/2:end)-xhat(i,N/2:end,j)).^2,2)/vars(i); %MSE
    end
    EOC_REM(j,:)=1./(chisq_all(:,j)+1) %EOC
end
save figure4E.mat

%now calculate for CI=0 (not included in above set since we plot on log
%scale)
xhat=zeros(dx,N);
CI_rem_opt(4,4)=0;
[xhat_0(:,:)]=UKF(y,R,dT,varn,CI_rem_opt);
mean_err_sq=mean((x(:,N/2:end)-xhat_0(:,N/2:end)).^2,2)./vars';
EOC_0=(1./(mean_err_sq+1));
save figure4E_q_rem_0.mat
return


function[]=Full_EOC_default_CI(x,y,dT,CI_eps,Rs)
%outputs
%inputs:
%x(matrix of DB model data)
%y(matrix of DB model observations-noise added version of x)
%dT(sampling,RK4 integration time)
%CI_eps (multiplier (by variance of each variable) for covariance inflation
%matrix)
%Rs (measurement noise covariance)


%setup Covariance inflation (CI) matrix
vars=var(x'); %variance of all model variables
%CI is a diagonal matrix where the diag entries=CI_eps*variance(variable)
CI=CI_eps.*blkdiag(vars(1),vars(2),vars(3),vars(4),vars(5), ...
        vars(6),vars(7),vars(8),vars(9),vars(10),vars(11),vars(12));

[dx,N]=size(x);
xhat=zeros(dx,N,dx); %preallocate
EOC=zeros(dx,dx);%preallocate

%track data using every variable as an observable
dobs=1; %dimension of obsevable

for i=1:dx %go through each state variable as observable
    varn=i
    R=(Rs([varn])'*(ones(1,dobs)))'.*eye(dobs,dobs); %measurement  covariance matrix
    [xhat(:,:,i)]=UKF(y(i,:),R,dT,varn,CI);
end

%calculate EOC
for j=[1:dx] %observed
    mean_err_sq=mean((x(:,N/2:end)-xhat(:,N/2:end,j)).^2,2)./vars';
    EOC(:,j)=(1./(mean_err_sq+1));
end

save figure4A.mat
return

function[]=Full_EOC_optimized_CI_delta(x,y,dT,CI_eps,Rs)
%outputs
%inputs:
%x(matrix of DB model data)
%y(matrix of DB model observations-noise added version of x)
%dT(sampling,RK4 integration time)
%CI_eps (multiplier (by variance of each variable) for covariance inflation
%matrix)
%Rs (measurement noise covariance)


vars=var(x');
CI=CI_eps.*blkdiag(vars(1),vars(2),vars(3),vars(4),vars(5), ...
        vars(6),vars(7),vars(8),vars(9),vars(10),vars(11),vars(12));

CI_optimized=CI;
CI_optimized(12,12)=CI(12,12)*4000; 

[dx,N]=size(x);
xhat=zeros(dx,N,dx); %preallocate
EOC=zeros(dx,dx);%preallocate

%track data using every variable as an observable
dobs=1; %dimension of obsevable

for i=1:dx %go through each state variable as observable
    varn=i
    R=(Rs([varn])'*(ones(1,dobs)))'.*eye(dobs,dobs); %measurement  covariance matrix
    [xhat(:,:,i)]=UKF(y(i,:),R,dT,varn,CI_optimized);
end

%calculate EOC
for j=[1:dx] %observed
    mean_err_sq=mean((x(:,N/2:end)-xhat(:,N/2:end,j)).^2,2)./vars';
    EOC(:,j)=(1./(mean_err_sq+1));
end

save figure4B.mat
return


function[]=Full_EOC_optimized_CI_delta_CI_rem(x,y,dT,CI_eps,Rs)
%outputs
%inputs:
%x(matrix of DB model data)
%y(matrix of DB model observations-noise added version of x)
%dT(sampling,RK4 integration time)
%CI_eps (multiplier (by variance of each variable) for covariance inflation
%matrix)
%Rs (measurement noise covariance)

vars=var(x');
CI=CI_eps.*blkdiag(vars(1),vars(2),vars(3),vars(4),vars(5), ...
        vars(6),vars(7),vars(8),vars(9),vars(10),vars(11),vars(12));

CI_optimized=CI;
CI_optimized(12,12)=CI(12,12)*4000;%delta
CI_optimized(4,4)=0; %F_R

[dx,N]=size(x);
xhat=zeros(dx,N,dx); %preallocate
EOC=zeros(dx,dx);%preallocate

%track data using every variable as an observable
dobs=1; %dimension of obsevable

for i=1:dx %go through each state variable as observable
    varn=i
    R=(Rs([varn])'*(ones(1,dobs)))'.*eye(dobs,dobs); %measurement  covariance matrix
    [xhat(:,:,i)]=UKF(y(i,:),R,dT,varn,CI_optimized);
end

%calculate EOC
for j=[1:dx] %observed
    mean_err_sq=mean((x(:,N/2:end)-xhat(:,N/2:end,j)).^2,2)./vars';
    EOC(:,j)=(1./(mean_err_sq+1));
end

save figure4C.mat
return

%this is the UKF function, embedded functions carry out filter-model
%integration, unscented transform, etc.
function [xhat]=UKF(y,R,dT,varn,CI)
%outputs:
%xhat (12 dimensional vector of reconstructioned estimates of x)

%inputs:
%y(12 dimensional vector of observables, noise-added version of x)
%R (measurement noise covariance matrix)
%dT (sampling, RK4 integration time)
%varn (index of observed variables)
%CI (covariance inflation matrix)


global dt nn %sampling time step, integration time stemp
dx=12; %dimension of aug state vector (DB2 has 12 variables)
dy=length(varn); % Dimensions of observation
fct='UKF_model'; % this is the model function F(x) used in filtering
obsfct='UKF_obsfct'; % this is the observation function G(x)

%RK4 integration setup
dt=dT; %dt is integration time, which can be smaller than dT, the sampling time
N=length(y);
nn=fix(dT/dt);  % the integration time step can be smaller than dT
Times = dT*(1:N);

%preallocate matrices
xhat=zeros(dx,N); % Preallocate estimated x
Pxx=zeros(dx,dx,N); % Prallocate Covariance in x
Ks=zeros(dx,dy,N); % Preallocate Kalman gains

%Store which variable was observed and dT
P.varn=varn;
P.dT=dT;

%initialize xhat and Pxx
xhat(:,1)=[3 3 2.5 2.5 2.5 0.5 0.5 0.5 0.5 0.5 0.25 0]; %arbitrary
Pxx(:,:,1)=CI;
Pxx(4,4,1)=1e-10; %for no covariance inflation, still have to initialize non-zero

%%%%%%%%%%%%%%%%%%%%%%%%%model parameters
P=OriginalDBParams;

%Main loop for recursive UKF estimation
for k=2:N
    [xhat(:,k),Pxx(:,:,k),Ks(:,:,k)]=UKF_ut(xhat(:,k-1),Pxx(:,:,k-1),y(:,k),fct,obsfct,dx,R,CI,varn,P);
end
return


%this function selects the index of variables in G(x), which is a noise
%added version of F(x)
function observ=UKF_obsfct(y,varn)
%outputs:
%observ: noisy observations vector, whose number of rows=length(varn)

%inputs:
%y(12 dimensional vector of observables, noise-added version of x)
%varn (index of observed variables)

observ=[];
for i=1:length(varn)
observ=[observ; y(varn(i),:)];
end
return


%this is the Unscented transform, the heart of the UKF
function [xhat_new,Pxx_new,K]=UKF_ut(xhat,Pxx,y,fct,obsfct,dx,R,CI,obs_var,P) 
% Unscented transformation from Voss et al 2004
% This Function has been modified by S. Schiff and T. Sauer 2008
% and by Madineh Sarvestani 2012

%outputs:
%xhat_new:new reconstruction at time t=k
%Pxx_new: updated model covariance (model uncertainty) matrix, at time t=k
%K: the Kalman gain

%inputs:
%xhat: reconstructed estimate from the last iteration, at time t=k-1
%Pxx: model covariance (model uncertainty) matrix, at time t=k-1
%y: observation matrix
%fct: the UKF filter model that iterates sigma points
%obsfct: the observation function that selects index of observed variables
%dx: dimension of state-space
%R: measurement covariance matrix
%CI: model covariance inflation matrix
%obs_var: index of observed variables
%P: struct which holds parameter values


N=2*dx; %Number of Sigma Points
Pxx=(Pxx+Pxx')/2;
xsigma=chol((dx*Pxx ))'; % Cholesky decomposition - note that Pxx=chol'*chol
Xa=xhat*ones(1,N)+[xsigma, -xsigma]; %Generate Sigma Points around xhat
X=feval(fct,Xa,P); %Iterate sigma points through model all at once
X(find(isnan(X)==1))=0;%remove bad points
xtilde=sum(X')'/N; %Mean of X's
X1=X-xtilde*ones(1,size(X,2)); % subtract mean from X columns
Pxx=X1*X1'/(N)+CI; %Pxx covariance once iterated through model, plus inflation
Pxx=(Pxx+Pxx')/2; 
Y=feval(obsfct,X,obs_var); %the select columns corresponding to observed variables
ytilde=sum(Y')'/(N); %mean of the sigma points for observed variables
Y1=Y-ytilde*ones(1,size(Y,2)); % subtract mean from Y columns
Pyy=Y1*Y1'/(N) + R; %Pyy covariance calculation

Pxy=X1*Y1'/(N); %cross-covariance calculation
K=Pxy*inv(Pyy); %Kalman gain calculation
xhat_new=xtilde+K*(y-ytilde); %update xhat after observation is made
Pxx=Pxx-K*Pxy'; Pxx_new=(Pxx+Pxx')/(2); %update Pxx after observation is made
  
%set maxes and minimums for xhat estimate
min=[0 0 0 0 0 0 0 0 0 0 0 0]';
max=[6.5 6.5 5 10 5 1 1 1 1 1 1 20]';
%
xhat_new(1:dx)=xhat_new(1:dx).*(xhat_new(1:dx)>min)+(min.*(xhat_new(1:dx)<min));
xhat_new(1:dx)=xhat_new(1:dx).*(xhat_new(1:dx)<max)+(max.*(xhat_new(1:dx)>max));
return


%this function integrates the filter model for the UKF to iterate sigma
%points Xa
function [X]=UKF_model(x,P)
%outputs:
%X(iterated sigma points X)

%intputs:
%x(sigma point matrix Xa: 2dx+1 by dx)
%P(struct of parameter values)


global dt nn
xnl=x(1:size(x(:,1)),:);

[Rows,Columns]=size(xnl);
xn2=zeros(Rows,Columns);

for NumCol=1:Columns
    xn2(:,NumCol)=x(1:size(x(:,1)),NumCol);
    for n=1:nn %if dt!=dT, then nn would be >1
        %RK4
        k1=dt.*DB_SS(xn2(:,NumCol),P);
        k2=dt.*DB_SS(xn2(:,NumCol)+k1/2,P);
        k3=dt.*DB_SS(xn2(:,NumCol)+k2/2,P);
        k4=dt.*DB_SS(xn2(:,NumCol)+k3,P);
        xn2(:,NumCol)=xn2(:,NumCol)+k1./6+k2./3+k3./3+k4./6;
    end
end
X=[ xn2];
return

%this function holds the DB model state-space
function [x_dot]=DB_SS(x,P)
%output:
%x_dot: 12 dimensional vector of instantaneous derivatives for firing rates, transmitter concentrations,
%h and delta at time t=k+1

%inputs:
%x ( 12-dimensional data vector at time t=k)
%P (parameter values stored in struct format)


%what comes in:
%[F_LC,F_DR,F_VLPO,F_R,F_WR,C_N,C_S,C_G,C_AR,C_AWR,h]
F_LC=x(1);
F_DR=x(2);
F_VLPO=x(3);
F_R=x(4);
F_WR=x(5);
C_N=x(6);
C_S=x(7);
C_G=x(8);
C_AR=x(9);
C_AWR=x(10);
h=x(11);
C_A=C_AR+C_AWR;
deltaLC=x(12);

F=[x(1:5)]'; %firing rates in vector format
C=[x(6:10)]';%neurotransmitters in vector format


%steady state neurotransmitter concentrations 
C_ss=tanh(F./P.cgamma).*1; 
C_Dot=(C_ss-C)./P.ctau;

%homeostatic sleep drive
heavarg1=(F_LC+F_DR)-P.thetaW;
heavarg2=P.thetaW-(F_LC+F_DR);
hDot=((heavarg1>=0)*((1-h)/P.tauhw))-((heavarg2>=0)*(h/P.tauhs));

%delta noise input to LC and DR
deltaLCDot=-(0.1*deltaLC);

%summed neurotransmmitter input
cLC=P.gALC*C_A-P.gNLC*C_N-P.gGLC*C_G+deltaLC; 
cDR=P.gADR*C_A-P.gSDR*C_S-P.gGDR*C_G+deltaLC;
cVLPO=-P.gNVLPO*C_N-P.gSVLPO*C_S-P.gGVLPO*C_G;
cR=P.gAR*C_A-P.gNR*C_N-P.gSR*C_S-P.gGR*C_G;
cWR=P.gAWR*C_A-P.gGWR*C_G;
Cin=[cLC cDR cVLPO cR cWR];

%firing rate param
Fbeta=[P.betaLC P.betaDR -7*h P.betaR P.betaWR];

%firing rate
F_ss=(P.Fmax.*(0.5.*(1+tanh((Cin-Fbeta)./P.Falpha))));
F_Dot=(F_ss-F)./P.Ftau;

%[F_LC,F_DR,F_VLPO,F_R,F_WR,C_N,C_S,C_G,C_AR,C_AWR,h,deltaLC]
x_dot=[F_Dot(:); C_Dot(:); hDot(:);deltaLCDot];
return

function []=plot_figure


%% plot various pieces

cd('../')%go up one foler
addpath(cd); %add path
cd([cd '/Figure 4_CI optimizaiton']) %reset folder back 

size_label=18;
size_tick=18;

labels={'F_{LC}','F_{DR}','F_{VLPO}','F_{R}','F_{W/R}','C_N',...
    'C_S','C_G','C_{A(R)}','C_{A(WR)}','h','\delta'};

no_labels={'','','','','','',...
    '','','','','',''};


%%%%%%%%%%%%% Figure 4A default EOC
figure;
subplot(3,2,1);
load figure4A.mat

%add empty column and row necessary for plotting
temp=flipud(EOC);
[N2,M]=size(temp);
plot_dist=[[temp,zeros(N2,1)];[zeros(1,M+1)]];

pcolor(plot_dist);
colormap(flipud(gray));
%ylabel('Reconstructed''fontsize', size_label,'fontweight','bold');
xlabel('Observed','fontsize', size_label,'fontweight','bold');
zlabel('ln(Distance)')
set(gca,'XTick',(1:12)+0.5); 
set(gca,'YTick',(1:12)+0.5); 
set(gca, 'fontsize', 14,'LineWidth',2)
set(gca,'XAxisLocation','top');
[hx,hy] = format_ticks(gca,labels',fliplr(labels)',[],[],90,0);



%%%%%%%%%%%%%Figure 4B EOC with IC_delta optimized
subplot(3,2,3);
%goes EOC with delta optimized
load figure4B.mat

%add empty column and row necessary for plotting
temp=flipud(EOC);
[N2,M]=size(temp);
plot_dist=[[temp,zeros(N2,1)];[zeros(1,M+1)]];

pcolor(plot_dist);
colormap(flipud(gray));
ylabel('Reconstructed','fontsize', size_label,'fontweight','bold');
zlabel('ln(Distance)')
set(gca,'XTick',(1:12)+0.5);
set(gca,'YTick',(1:12)+0.5);
set(gca, 'fontsize', 14,'LineWidth',2)
[hx,hy] = format_ticks(gca,labels',fliplr(labels)',[],[],90,0);
colorbar;
set(gca, 'fontsize', 14,'LineWidth',2)



%%%%%%%%%%%%%%Figure 4C EOC with IC optimized for delta and F_R
load figure4C.mat
subplot(3,2,5)

%add empty column and row necessary for plotting
temp=flipud(EOC);
[N2,M]=size(temp);
plot_dist=[[temp,zeros(N2,1)];[zeros(1,M+1)]];

pcolor(plot_dist);
colormap(flipud(gray));
%ylabel('Reconstructed');
%xlabel('Observed');
zlabel('ln(Distance)')
set(gca,'XTick',(1:12)+0.5);
set(gca,'YTick',(1:12)+0.5);
set(gca, 'fontsize', 14,'LineWidth',2)
[hx,hy] = format_ticks(gca,labels',fliplr(labels)',[],[],90,0);


%%%%%%%%%%%Figure 4D EOC(delta,F_LC) as a function of CI_delta
load figure4D.mat
subplot(2,2,2);

Q_plot=log10((CI(12,12).*CI_coeff)./(vars(12))); %this
plot(Q_plot,EOC_delta(:,12),'b*--'); hold on;
xlim([-6 0]);
xlabel('log_{10}(qdelta/var)');
ylabel('EOC_{\delta,F_{LC}}','fontsize', size_label,'fontweight','bold');
set(gca, 'fontsize', 14,'LineWidth',2)


%%%%%%%%%%%Figure 4D EOC(F_R,F_LC) as a function of CI_F_R
load figure4E.mat
subplot(2,2,4);

Q_plot=log10((CI(4,4).*CI_coeff)./(vars(4))) %this
plot(Q_plot(1:end-1),EOC_REM(1:end-1,4),'b*--'); hold on; %can't take log(0)
xlim([-6 0]);
xlabel('log_{10}(qrem/var)');
ylabel('EOC_{F_{R},F_{LC}}','fontsize', size_label,'fontweight','bold');
set(gca, 'fontsize', 14,'LineWidth',2)

%now draw black line
load figure4E_q_rem_0.mat EOC_0
xaxis=Q_plot(1:end-1);
yaxis=EOC_0(4,1).*ones(size(xaxis));

hold on; plot(xaxis,yaxis,'k','LineWidth',2);

return

