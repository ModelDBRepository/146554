%Reproduces Figure 7 in Sedigh-Sarvestani, Schiff, Gluckman,
%'Reconstructing Mammalian Sleep Dynamics with Data Assimiliation',
%PLoS Comp Biol, 2012 (In press). DOI:10.1371/journal.pcbi.1002788

%Reconstructs unobserved dynamics from DB model of sleep.
%DB model from: Diniz Behn and Booth, J Neurophysiol 103:1937-1953, 2010.

%Observed variables used in the UKF are derived from the State-of-vigilance (see map in code below)
%rather than directly from the variables of the DB model.

%Usage: Running this .m file will produce a figure similar to Figure 7 of
%Sedigh-Sarvestani et al. 
%Make sure CD is '...\Figure Code\Figure 7_Hypno Tracking'

%Options: You can use previously generated model data, from which noisy
%observations are obtained. This is stored as 'data_DB_output.mat' in
%'\Figure Code\Generate Data\' and is used throughout the paper.
%Alternatively, you can re-generate model data from scratch.

%Madineh Sedigh-Sarvestani, Penn State, Oct 2012
%m.sedigh.sarvestani@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[]=figure7

cd('../')%go up one folder
addpath(genpath(cd)); %add path
load CI_eps.mat %default multiplier for Covariance Inflation (CI)
cd([cd '/Figure 7_Hypno Tracking']) %reset folder back

if exist('data_DB_output.mat','file')
    load data_DB_output.mat %load already generated data
else
    [Times,x,y,state,dT,P,Rs]=data_DB; %alternatively re-generate data
end

%now setup Covariance inflation (CI) matrix
vars=var(x'); %variance of all model variables
%CI is a diagonal matrix where the diag entries=CI_eps*variance(variable)
CI=CI_eps.*blkdiag(vars(1),vars(2),vars(3),vars(4),vars(5), ...
        vars(6),vars(7),vars(8),vars(9),vars(10),vars(11),vars(12));

CI_optimized=CI;
CI_optimized(4,4)=0; %F_R
CI_optimized(12,12)=CI(12,12)*1000; %delta
CI_optimized(1,1)=CI(1,1)*30; %F_LC
CI_optimized(2,2)=CI(2,2)*30; %F_DR

%% now extract observables from sleep-state of DB data
%state is generated in data_DB.m based on a rule that thresholds the firing
%rates of LC,DR,R,and VLPO
wake_on=x(:,state==1); %x values when in wake
rem_on=x(:,state==3); %x values when in REM
nrem_on=x(:,state==2); %xvalues when in NREM

%get medians and standard deviations
mean_wake=mean(wake_on,2); med_wake=median(wake_on,2); s_wake=std(wake_on')';
mean_nrem=mean(nrem_on,2); med_nrem=median(nrem_on,2); s_nrem=std(nrem_on')';
mean_rem=mean(rem_on,2); med_rem=median(rem_on,2); s_rem=std(rem_on')';

%wid^2 = var+(mean-median)^2
%now correct width of distribution since using median
var_wake=(s_wake.^2)+(mean_wake-med_wake).^2;
var_nrem=(s_nrem.^2)+(mean_nrem-med_nrem).^2;
var_rem=(s_rem.^2)+(mean_rem-med_rem).^2;

varn=[1 2 3 4]; %choose which variables we'll create from mapping SOV: F_LC,F_DR,F_VLPO,F_R
y=zeros(length(varn),length(x)); %preallocate
R=zeros(length(varn),length(varn),length(x)); %preallocate R

%create observable set y-derived from sleep-state
for k=1:length(x);
    if state(k)==1 %during wake
        y(varn,k)=med_wake(varn); %firing rates are set to medians
        Rnot=(1/1).*(var_wake(varn));%R is set to variance
        R(:,:,k)=diag(Rnot,0);
    elseif state(k)==2 %during nrem
        y(varn,k)=med_nrem(varn);
        Rnot=(1/1).*(var_nrem(varn));
        R(:,:,k)=diag(Rnot,0);
    else
        y(varn,k)=med_rem(varn); %during rem
        Rnot=(1/1).*(var_rem(varn));
        R(:,:,k)=diag(Rnot,0);
    end
end

y=y.*(y>0); %make sure all observations are >0
y(2,:)=y(1,:); %make sure LC and DR observations match (as they do in DB model)

%this is needed in case we want to map SOV onto variables not in order
for i=1:length(varn)
    y(i,:)=y(varn(i),:);
end
yinput=y(1:length(varn),:);

%now reconstruct unobserved variables using these mapped observables
%reminder: mapped observables are F_LC,F_DR,F_VLPO,F_R
ic=ones(12,1); %initial condition
[xhat]=UKF(varn,yinput,dT,R,CI_optimized,ic);

save figure7.mat
plot_figure;
return

%this is the UKF function, embedded functions carry out filter-model
%integration, unscented transform, etc.
function [xhat]=UKF(varn,y,dT,R,CI,ic)
%outputs:
%xhat (12 dimensional vector of reconstructioned estimates of x)

%inputs:
%varn (index of observed variables)
%y(12 dimensional vector of observables, noise-added version of x)
%dT (sampling, RK4 integration time)
%R (measurement noise covariance matrix)
%CI (covariance inflation matrix)
%parameter (current value for parameter of interest)
%ic (model initial condition)


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
Params.varn=varn;
Params.dT=dT;

%initialize xhat and Pxx
xhat(:,1)=ic;
Pxx(:,:,1)=CI;
Pxx(4,4,1)=1e-10; %for no covariance inflation, still have to initialize non-zeroize non-zero

%%%%%%%%%%%%%%%%%%%%%%%%%model parameters
Params=OriginalDBParams;

%Main loop for recursive UKF estimation
for k=2:N       
[xhat(:,k),Pxx(:,:,k),Ks(:,:,k)]=UKF_ut(xhat(:,k-1),Pxx(:,:,k-1),y(:,k),fct,obsfct,dx,R(:,:,k),CI,varn,Params);
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
max=[6.5 6.5 5 10 5 1 1 1 1 1 1 10]';
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


function plot_figure

size_label=18;
size_tick=18
load figure7.mat
Times=Times(1:length(x));
[dx,N]=size(x);

%plot hypnogram, or the sleep-state for 1 hour
figure;
subplot(6,1,2);
plot(Times,state);
plot(Times,state(1,:),'k','LineWidth',3); 
ylim([0.8 3.2]);set(gca,'YTick',[1 2 3]);
set(gca,'YTickLabel',{'Wake','NREM','REM'},'fontsize',size_label,'fontweight','bold');
xlim([0 3600]); set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'','',''});
%xlabel('Time (Hours)','fontsize', size_label,'fontweight','bold');
set(gca, 'fontsize', size_tick,'LineWidth',2)
hold on;
%now plot color
symbols={'g*','r*','b*'};
for i=1:3600*2;
    plot(Times(i),state(i),symbols{state(i)});
end
% now get distributions of firing rates for LC,VLPO,and R during different
% sleep states and plot

% wake_on=x(:,state==1); %x values when in wake
% rem_on=x(:,state==3); %x values when in REM
% nrem_on=x(:,state==2); %xvalues when in NREM

%then get histograms
%plot F_LC during different SOV
edges=0:0.5:6.5;
subplot(6,3,1)
[N1,xout1] = hist(wake_on(1,:),edges);
[N2,xout2] = hist(nrem_on(1,:),xout1);
[N3,xout3] = hist(rem_on(1,:),xout1);
ybar=[N1/length(wake_on(1,:));N2/length(nrem_on(1,:));N3/length(rem_on(1,:))];
stairs(xout1,ybar(1,:)','g','LineWidth',3);hold on;
stairs(xout2,ybar(2,:),'r','LineWidth',3);
stairs(xout3,ybar(3,:),'b','LineWidth',3);
title('F_{LC}','fontsize', size_label,'fontweight','bold');
set(gca,'YTickLabel',[]);
xlim([0 7]);xlabel('Firing Rate','fontsize', 16,'fontweight','bold');
set(gca, 'fontsize', 12,'LineWidth',2);

%plot F_VLPO dring different SOV
subplot(6,3,2);
[N1,xout1] = hist(wake_on(3,:),edges);
[N2,xout2] = hist(nrem_on(3,:),xout1);
[N3,xout3] = hist(rem_on(3,:),xout1);
ybar=[N1/length(wake_on(3,:));N2/length(nrem_on(3,:));N3/length(rem_on(3,:))];
bar1=bar(xout1',ybar(1,:)','FaceColor', 'k', 'EdgeColor', 'k');
stairs(xout1,ybar(1,:)','g','LineWidth',3);hold on;
stairs(xout2,ybar(2,:),'r','LineWidth',3);
stairs(xout3,ybar(3,:),'b','LineWidth',3);
title('F_{VLPO}','fontsize', size_label,'fontweight','bold');
set(gca,'YTickLabel',[]);
xlim([0 6.5]);xlabel('Firing Rate','fontsize', 16,'fontweight','bold');
set(gca, 'fontsize', 12,'LineWidth',2);

%plot F_R during different SOV
subplot(6,3,3);
[N1,xout1] = hist(wake_on(4,:),edges);
[N2,xout2] = hist(nrem_on(4,:),xout1);
[N3,xout3] = hist(rem_on(4,:),xout1);
ybar=[N1/length(wake_on(4,:));N2/length(nrem_on(4,:));N3/length(rem_on(4,:))];
bar1=bar(xout1',ybar(1,:)','FaceColor', 'k', 'EdgeColor', 'k');
stairs(xout1,ybar(1,:)','g','LineWidth',3);hold on;
stairs(xout2,ybar(2,:),'r','LineWidth',3);
stairs(xout3,ybar(3,:),'b','LineWidth',3);
title('F_{R}','fontsize', size_label,'fontweight','bold');
set(gca,'YTickLabel',[]);
xlim([0 6.5]);xlabel('Firing Rate','fontsize', 16,'fontweight','bold');
set(gca, 'fontsize', 12,'LineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%now plot reconstruction results
subplot(6,1,3);
%yinput is the observables generated from SOV
%xhat is the reconstruction
%x is true underlying value
j=1;%F_LC
plot(Times, yinput(j,:),'b.'); hold on; plot(Times,xhat(j,:),'r'); plot(Times,x(j,:),'k'); 
ylabel('F_{LC}','fontsize', size_label,'fontweight','bold');
xlim([0  3600]); ylim([0 7]); 
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'','',''});
set(gca, 'fontsize', size_tick,'LineWidth',2)


subplot(6,1,4);
j=3;%F_VLPO
plot(Times,xhat(j,:),'r'); hold  on;  plot(Times,x(j,:),'k'); 
ylabel('F_{VLPO}','fontsize', size_label,'fontweight','bold');
xlim([0  3600]); ylim([0 7]); 
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'','',''});
set(gca, 'fontsize', size_tick,'LineWidth',2)


subplot(6,1,5);
j=4;%F_R
plot(Times,xhat(j,:),'r');hold on; plot(Times,x(j,:),'k');
ylabel('F_{R}','fontsize', size_label,'fontweight','bold');
xlim([0  3600]); ylim([0 7]); 
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'','',''});
set(gca, 'fontsize', size_tick,'LineWidth',2)


subplot(6,1,6);
j=11; %h
plot(Times,xhat(j,:),'r');hold on; plot(Times,x(j,:),'k');
ylabel('h','fontsize', size_label,'fontweight','bold');
xlabel('Time (Hours)','fontsize', size_label,'fontweight','bold');
xlim([0  3600]); ylim([0 1]); 
set(gca,'XTick',[0 1800 3600]);set(gca,'XTickLabel',{'0','0.5','1'});
set(gca, 'fontsize', size_tick,'LineWidth',2)


return