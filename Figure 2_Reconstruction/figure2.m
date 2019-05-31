%Reproduces Figure 2 in Sedigh-Sarvestani, Schiff, Gluckman,
%'Reconstructing Mammalian Sleep Dynamics with Data Assimiliation',
%PLoS Comp Biol, 2012 (In press). DOI:10.1371/journal.pcbi.1002788

%Reconstructs unobserved dynamics from DB model of sleep.
%DB model from: Diniz Behn and Booth, J Neurophysiol 103:1937-1953, 2010.

%Usage: Running this .m file will produce a figure similar to Figure 2 of
%Sedigh-Sarvestani et al. 
%Make sure CD is '...\Figure Code\Figure2_Reconstruction'

%Options: You can use previously generated model data, from which noisy
%observations are obtained. This is stored as 'data_DB_output.mat' in
%'\Figure Code\Generate Data\' and is used throughout the paper.
%Alternatively, you can re-generate model data from scratch.

%Madineh Sedigh-Sarvestani, Penn State, Oct 2012
%m.sedigh.sarvestani@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[]= figure2
cd('../')%go up one foler
addpath(genpath(cd)); %add path

load CI_eps.mat %this is the default CI (covariance inflation) multiplier for all variables used throughout the figures

cd([cd '/Figure 2_Reconstruction']) %reset folder back 

if exist('data_DB_output.mat','file')
    load data_DB_output.mat %load already generated data
else
    [Times,x,y,state,dT,P,Rs]= data_DB; %alternatively re-generate data
end

%now setup Covariance inflation (CI) matrix
vars=var(x'); %variance of all model variables
%CI is a diagonal matrix where the diag entries=CI_eps*variance(variable)
CI=CI_eps.*blkdiag(vars(1),vars(2),vars(3),vars(4),vars(5), ...
        vars(6),vars(7),vars(8),vars(9),vars(10),vars(11),vars(12));

%optimized covariance inflation (see Figure 4)
CI_optimized=CI;
CI_optimized(4,4)=0; %F_R
CI_optimized(12,12)=CI(12,12)*4000; %delta

%We'll reconstruct 2 hours of data
%dT=0.5 seconds is how often we'll sample 
totaltime=3600*2/dT;
x=x(:,1:totaltime);
y=y(:,1:totaltime);
Times=dT:dT:totaltime*dT;

%We'll observe firing rate of the wake-active LC group
varn=1; %F_LC index
y=y(1,:); %observe F_LC
dy=length(varn);
R=(Rs([varn])'*(ones(1,dy)))'.*eye(dy,dy); %measurement  covariance matrix 


%now run UKF with the settings above for CI,R,y and dT
%x is the mostly hidden true data set
[xhat,y,P]=UKF(x,y,R,dT,varn,CI); %unoptimized CI
save ('figure2A.mat')

[xhat,y,P]=UKF(x,y,R,dT,varn,CI_optimized); %optimzed CI
save ('figure2B.mat')
plot_figure;
return;


%this is the UKF function, embedded functions carry out filter-model
%integration, unscented transform, etc.
function [xhat,y,P]=UKF(x,y,R,dT,varn,CI)
%outputs:
%xhat (12 dimensional vector of reconstructioned estimates of x)
%y (12 dimensional vector of observables-noisy versions of x-output for plotting)
%P (struct that holds parameter value)

%inputs:
%x(12 dimensional vector of DB model output)
%y(12 dimensional vector of observables, noise-added version of x)
%R (measurement noise covariance matrix)
%dT (sampling, RK4 integration time)
%varn (index of observed variables)
%CI (covariance inflation matrix)

global dt nn %sampling time step, integration time stemp
dx=12; %dimension of state vector (DB2 has 12 variables)
dy=length(varn); % Dimensions of observation
fct='UKF_model'; % this is the model function F(x) used in filtering
obsfct='UKF_obsfct'; % this is the observation function G(x)

%RK4 integraiton setup
dt=dT; %dt is integration time, which can be smaller than dT, the sampling time
N=length(x);%number of data samples in sample time steps
nn=fix(dT/dt);  %number of cycles of model intergration before we sample

%preallocate matrices
xhat=zeros(dx,N); % Preallocate estimated x
Pxx=zeros(dx,dx,N); % Prallocate Covariance in x
errors=zeros(dx,N); % Preallocate errors
Ks=zeros(dx,dy,N); % Preallocate Kalman gains

%Store which variable was observed and dT
P.varn=varn;
P.dT=dT;

%initialize xhat and Pxx
xhat(:,1)=[3 3 2.5 2.5 2.5 0.5 0.5 0.5 0.5 0.5 0.25 0]; %initialize xhat abritrarily
xhat(varn,1)=x(varn,1); %initialize value of observed variables correctly
Pxx(:,:,1)=CI;
Pxx(4,4,1)=1e-10; %for no covariance inflation, still have to initialize non-zero

%%%%%%%%%%%%%%%%%%%%%%%%%model parameters
P=OriginalDBParams;

%Main loop for recursive UKF estimation
for k=2:N
    [xhat(:,k),Pxx(:,:,k),Ks(:,:,k)]=UKF_ut(xhat(:,k-1),Pxx(:,:,k-1),y(:,k),fct,obsfct,dx,R,CI,varn,P);
    errors(:,k)=sqrt(diag(Pxx(:,:,k))); %model uncertainty after each assimilation step
end;
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
X=[xn2];
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

%this function produces Figure 2
function[]=plot_figure

clear all;

size_label=18;
size_tick=18;

load figure2A.mat %UKF estimation with unoptimized CI

figure;
subplot(4,2,1);
plot(Times,y,'b*','MarkerSize',3); hold on;
plot(Times,x(1,:),'k','LineWidth',4);hold on;
plot(Times,xhat(1,:),'r','LineWidth',2); 
ylabel('F_{LC}','fontsize', size_label,'fontweight','bold');
xlim([0  3600]); ylim([0 7]); 
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'','',''});
set(gca,'YTick',[0 5]);set(gca,'YTickLabel',{'0','5'});
set(gca, 'fontsize', size_tick,'LineWidth',2)


subplot(4,2,5); 
plot(Times,x(4,:),'k','LineWidth',4); hold on;
plot(Times,xhat(4,:),'r','LineWidth',2); 
ylabel('F_{R}','fontsize', size_label,'fontweight','bold');
xlim([0  3600]); ylim([0 7]);
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'','',''});
set(gca,'YTick',[0 5]);set(gca,'YTickLabel',{'0','5'});
set(gca, 'fontsize', size_tick,'LineWidth',2)

subplot(4,2,3); 
plot(Times,x(3,:),'k','LineWidth',4);hold on;
plot(Times,xhat(3,:),'r','LineWidth',2); 
ylabel('F_{VLPO}','fontsize', size_label,'fontweight','bold');
xlim([0  3600]); ylim([0 7]);
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'','',''});
set(gca,'YTick',[0 5]);set(gca,'YTickLabel',{'0','5'});
set(gca, 'fontsize', size_tick,'LineWidth',2)

subplot(4,2,7);
plot(Times,x(12,:),'k','LineWidth',4);
hold on; plot(Times,xhat(12,:),'r','LineWidth',2);
ylabel('\delta','fontsize', size_label,'fontweight','bold'); 
xlabel('Time (Hours)','fontsize', size_label,'fontweight','bold');
xlim([0  3600]); ylim([0 10]);
set(gca,'YTick',[0 10]);set(gca,'YTickLabel',{'0','10'});
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'0','0.5','1'});
set(gca, 'fontsize', size_tick,'LineWidth',2)



load figure2B.mat %UKF estimation with optimized CI

subplot(4,2,2);
plot(Times,y,'b*','MarkerSize',3); hold on;
plot(Times,x(1,:),'k','LineWidth',4);hold on;
plot(Times,xhat(1,:),'r','LineWidth',2); %ylabel('F_{LC}');
xlim([0  3600]);ylim([0 7]); 
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'','',''});
set(gca,'YTick',[0 5]);set(gca,'YTickLabel',{'',''});
set(gca, 'fontsize', size_tick,'LineWidth',2)

subplot(4,2,6);
plot(Times,x(4,:),'k','LineWidth',4); hold on;
plot(Times,xhat(4,:),'r','LineWidth',2); %ylabel('F_{R}');
xlim([0  3600]); ylim([0 7]);
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'','',''});
set(gca,'YTick',[0 5]);set(gca,'YTickLabel',{'',''});
set(gca, 'fontsize', size_tick,'LineWidth',2)

subplot(4,2,4);
plot(Times,x(3,:),'k','LineWidth',4);hold on; 
plot(Times,xhat(3,:),'r','LineWidth',2); %ylabel('F_{VLPO}');
xlim([0  3600]); ylim([0 7]);
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'','',''});
set(gca,'YTick',[0 5]);set(gca,'YTickLabel',{'',''});
% set(gca,'Yaxislocation','right');
set(gca, 'fontsize', size_tick,'LineWidth',2)


subplot(4,2,8); 
plot(Times,x(12,:),'k','LineWidth',4); hold on; 
plot(Times,xhat(12,:),'r','LineWidth',2); 
xlabel('Time (Hours)','fontsize', size_label,'fontweight','bold');
xlim([0  3600]); ylim([0 10]);
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'0','0.5','1'});
set(gca,'YTick',[0 5]);set(gca,'YTickLabel',{'',''});
set(gca, 'fontsize', size_tick,'LineWidth',2)

return