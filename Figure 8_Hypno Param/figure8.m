%Reproduces Figure 8 in Sedigh-Sarvestani, Schiff, Gluckman,
%'Reconstructing Mammalian Sleep Dynamics with Data Assimiliation',
%PLoS Comp Biol, 2012 (In press). DOI:10.1371/journal.pcbi.1002788

%Reconstructs unobserved dynamics from DB model of sleep.
%DB model from: Diniz Behn and Booth, J Neurophysiol 103:1937-1953, 2010.

%Estimates unobserved parameter gALC, using observables that are 
%extracted from overall SOV of the animal and not direct observations of
%model variables

%Usage: Running this .m file will produce a figure similar to Figure 8 of
%Sedigh-Sarvestani et al. 
%Make sure CD is '...\Figure Code\Figure 8_Hypno Param'

%note: The input observables to this UKF/Param estimation function is
%generated in figure7.m. You can either upload that saved matrix here
%(figure7.mat) or you can run figure7.m, which will produce this matrix
%from scratch, then run this file.

%Madineh Sedigh-Sarvestani, Penn State, Oct 2012
%m.sedigh.sarvestani@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[]=figure8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%get data
cd('../')%go up one foler
addpath(genpath(cd));
cd([cd '/Figure 8_Hypno Param']) %reset folder back 

if exist('figre7.mat','file')
load figure7.mat yinput varn dT R CI_optimized ic Times %get generated data
else
    figure7; %alternatively re-generate data
    load figure7.mat yinput varn dT R CI_optimized ic Times
end

y_obs=yinput; %match terminology 
N=length(yinput); %number of data points
dx=12;%dimension of DB state space


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%parameter estimation+UKF
Mvec=zeros(dx,1);
Mvec(varn)=1;
M=diag(Mvec',0); %creates a diagonal matrix to weigh cost-function by observed varibles

%now setup param estimation windows Twin, and overlap
pTimeWindow=30;%param estimate window length in minutes
pWindow = floor(pTimeWindow*60/dT);
%overlap is pTimeWindow-pWindowOverlap;
pWindowOverlap = 6; %minutes (how much the window is moved)
pOverlap = floor(pWindowOverlap*60/dT);
i1 = 1:pOverlap:(N-pWindow);
i2 = i1+pWindow-1;
divd=length(i1); %number of separate windows we'll estimate over

%how often we force sync trajectories onto xhat
% set up offsets - assume i2(n)-i1(n)=const
ForceFix_step = 120; %seconds

%setup how much parameter is allowed to vary on each update
DeltaMax = 0.1; %(2*pi/(24*3600))*pWindow;

%preallocate matrices
y=zeros(dx,N); %test trajectories
xhat=zeros(dx,N); %UKF reconstruction

% initial value for parameter we'll estimate
p_initial=[7];
p=zeros(1,1+divd);
dp=zeros(1,divd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%begin
for n=1:divd
    if n<=1 %for first iteration use initial values
        p_input=p_initial;
        ic=[3 3 2.5 2.5 2.5 0.5 0.5 0.5 0.5 0.5 0.25 0];
    else
        p_input=p(n-1);%update parameter
        ic=xhat(:,i1(n-1)+pOverlap); %update UKF ic
    end
    
    %Iterate between UKF estimation and parameter estimation
    %UKF estimation using p_input as paramater value for gALC
    [xhat(:,i1(n):i2(n))]=UKF(y_obs(:,i1(n):i2(n)),R(:,:,i1(n):i2(n)),dT,varn,CI_optimized,p_input,ic);
    
    %setup force sync (xhat to trajectory) values
    DivLength = i2(n)-i1(n);
    RangeStart = 0:ForceFix_step:DivLength-5;
    RangeEnd = circshift(RangeStart',-1)';
    RangeEnd(end)=DivLength;
    NSteps = length(RangeStart);
    
    p_new_plus = p_input+DeltaMax; %current param +delta
    p_new_minus = p_input-DeltaMax; %current param -delta

    
    %generate the FIRST test trajectory with param=current param (use force
    %sync)
    for step = NSteps:-1:1
        iS = i1(n)+RangeStart(step);
        iE = i1(n)+RangeEnd(step);
        [yy(:,iS:iE)]...
            =test_traj(xhat(:,iS:iE),p_input,{'gALC'});
    end
    %calculate error between this trajectory and xhat
    error_n=xhat(:,i1(n):i2(n))-yy(:,i1(n):i2(n));
    error_n = M*error_n;
    costs_old = mean(mean(error_n.*error_n,2));

    
    %generate the SECOND test trajectory with param=current param+delta (use force
    %sync)
    for step = NSteps:-1:1
        iS = i1(n)+RangeStart(step);
        iE = i1(n)+RangeEnd(step);
        [yyplus(:,iS:iE)]...
            =test_traj(xhat(:,iS:iE),p_new_plus,{'gALC'});
    end
    %calculate error between this trajectory and xhat
    error_n=xhat(:,i1(n):i2(n))-yyplus(:,i1(n):i2(n));
    error_n = M*error_n;
    costs_plus = mean(mean(error_n.*error_n,2));

    %generate the THIRD test trajectory with param=current param-delta (use force
    %sync)      
    for step = NSteps:-1:1
        iS = i1(n)+RangeStart(step);
        iE = i1(n)+RangeEnd(step);
        [yyminus(:,iS:iE)]...
            =test_traj(xhat(:,iS:iE),p_new_minus,{'gALC'});
    end
    %calculate error between this trajectory and xhat
    error_n=xhat(:,i1(n):i2(n))-yyminus(:,i1(n):i2(n));
    error_n = M*error_n;
    costs_minus = mean(mean(error_n.*error_n,2));
    
    
    
    %now determine which test trajectory minimizes the cost_function
    if (costs_old==min([costs_old, costs_plus,costs_minus]))
        dp(n) = 0; %if current parameter minimizes, don't change parameter values
    else
        if (costs_plus<costs_minus)
            dp(n) = DeltaMax;%otherwise update parameter by +delta
            y(:,iS:iE)=yyplus(:,iS:iE);%and keep corresponding trajectory
        else
            dp(n) = -DeltaMax;%otherwise update parameter by -delta
            y(:,iS:iE)=yyminus(:,iS:iE);%and keep corresponding trajectory
        end
    end
    
    p_input = p_input+dp(n);%update current parameter for next UKF step
    p(n) = p_input;% store in a separate vector

end

%now create a full vector of true gALC values for plotting
TrueP=OriginalDBParams;
p=[p_initial,p]; %estimated values
for n=1:divd
    p_full(i1(n):i2(n))=p(n); %our estimates vector
    p_best(i1(n):i2(n))=TrueP.gALC; %true value vector
end

save figure8.mat;

plot_figure;
return

%this is the UKF function, embedded functions carry out filter-model
%integration, unscented transform, etc.
function [xhat]=UKF(y,R,dT,varn,CI, parameter,ic)
%outputs:
%xhat (12 dimensional vector of reconstructioned estimates of x)

%inputs:
%y(12 dimensional vector of observables, noise-added version of x)
%R (measurement noise covariance matrix)
%dT (sampling, RK4 integration time)
%varn (index of observed variables)
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
Pxx(4,4,1)=1e-10; %for no covariance inflation, still have to initialize non-zero

%%%%%%%%%%%%%%%%%%%%%%%%%model parameters
Params=OriginalDBParams;
Params.gALC=parameter; %this is the parameter we're estimating, its value will be input

%Main loop for recursive estimation
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
max=[6.5 6.5 5 10 5 1 1 1 1 1 1 20]';
%
xhat_new(1:dx)=xhat_new(1:dx).*(xhat_new(1:dx)>min)+(min.*(xhat_new(1:dx)<min));
xhat_new(1:dx)=xhat_new(1:dx).*(xhat_new(1:dx)<max)+(max.*(xhat_new(1:dx)>max));
return


%this is the tracker model (when not doing param estimation)

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
    for n=1:nn%if dt!=dT, then nn would be >1
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

%this is the function that produces a system of test-trajectories using 
%variable parameter values
function[y]=test_traj(xhat,params,paramnames)
%outputs:
%y (matrix of points for test-trajectory)

%inputs:
%xhat(matrix of reconstructed estimate)
%params(parameter of interest that is being estimated)
%paranmanes(name of parameter of interest)


global nn dt

dx=12; %dimension of DB
N=length(xhat(1,:)); %num data points

y=zeros(size(xhat));%preallocate
y(:,1)=[xhat(1:dx,1)]; %initialize test-trajectories
dqsync=length(params); %number of params we'll estimate


P=OriginalDBParams; %fix rest of parameters
%set the value for gALC, the parameter we're estimating
names=fieldnames(P);
for i=1:dqsync
    index=strmatch(paramnames{i},names);
    P=setfield(P,names{index},params(i));
end


%main loop for model (trajectory) integration
for n=1:N-1;
    
    yy=y(:,n);
    %error blowup protection (put boundaries on F and C values)
    min=[0 0 0 0 0 0 0 0 0 0 0 0 ]';
    max=[6.5 6.5 6 6 6 1 1 1 1 1 1 20]';
    
    yy(1:dx)=yy(1:dx).*(yy(1:dx)>min)+(min.*(yy(1:dx)<min));
    yy(1:dx)=yy(1:dx).*(yy(1:dx)<max)+(max.*(yy(1:dx)>max));
    
    %RK4
    for i=1:nn
        k1=dt*DB_SS(yy,P);
        k2=dt*DB_SS(yy,P);
        k3=dt*DB_SS(yy,P);
        k4=dt*DB_SS(yy,P);
        yy=yy+k1/6+k2/3+k3/3+k4/6;
    end;
    y(:,n+1)=yy; %add last data point to vector
end;

return;


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

load figure7.mat 
size_label=18;
size_tick=18

%compute group VARS from means
% this approximates the error estimate of the 
% median-estimated SOV dependent states
%
% recall that if  a=const
% then:
% <(X-a)^2>= var(X) + (<x>-a)^2
% first compute fractional time in each state
P_wake = sum(state==1)/sum(state>0);
P_nrem = sum(state==2)/sum(state>0);
P_rem  = sum(state==3)/sum(state>0);
eps_base = (P_wake*(var_wake + (mean_wake-med_wake).^2) + ...
    P_nrem*(var_nrem + (mean_nrem-med_nrem).^2) + ...
    P_rem*(var_rem + (mean_rem-med_rem).^2))'./vars;
    

%% now open figure
fig8 = figure();%'PaperUnits','inches','PaperPosition',[2, 2, 4, 4]);

%% now plot estimation
load figure8.mat
Times=Times(1:N);

%calculate mean errors
for n=1:divd
    error=(x(:,i1(n):i2(n))-xhat(:,i1(n):i2(n))).^2;
    mean_err(:,n)=mean(error,2)./vars';
end

%plot test trajectories overlayed with x(true) and xhat(UKF) trajectories
%for different time points of parameter estimation
Time=Times/3600;
subplot(5,3,4);
plot(Time,x(1,:),'k.','LineWidth',4); hold on;
plot(Time,xhat(1,:),'r.','LineWidth',3); hold on;
plot(Time(1:length(yy)),yy(1,:),'m.','LineWidth',3); xlim([0.15 0.35]); ylim([0 8]);
set(gca,'YTick',[0 5]);set(gca,'YTickLabel',{'0','5'});
set(gca,'XTick',[0.25]);set(gca,'XTickLabel',{'0.25'});
set(gca, 'fontsize', size_tick,'LineWidth',2)
ylabel('F_{LC}','fontsize', size_label,'fontweight','bold');

subplot(5,3,5);
plot(Time,x(1,:),'k.','LineWidth',4); hold on;
plot(Time,xhat(1,:),'r.','LineWidth',3); hold on;
plot(Time(1:length(yy)),yy(1,:),'m.','LineWidth',3); xlim([4.6 4.8]);ylim([0 8]);
set(gca,'YTick',[0 5]);set(gca,'YTickLabel',{});
set(gca,'XTick',[4.7]);set(gca,'XTickLabel',{'4.7'});
set(gca, 'fontsize', size_tick,'LineWidth',2)

subplot(5,3,6);
plot(Time,x(1,:),'k.','LineWidth',4); hold on;
plot(Time,xhat(1,:),'r.','LineWidth',3); hold on;
plot(Time(1:length(yy)),yy(1,:),'m.','LineWidth',3); xlim([11 11.2]);ylim([0 8]);
set(gca,'YTick',[0 5]);set(gca,'YTickLabel',{});
set(gca,'XTick',[11.1]);set(gca,'XTickLabel',{'11.1'});
set(gca, 'fontsize', size_tick,'LineWidth',2)

%plot parameter estimate (p_full) and true parameter value (p_best)
subplot(5,1,1);
plot(Time(1:length(p_full)),p_best,'k','LineWidth',4); hold on;
plot(Time(1:length(p_best)),p_full,'m','LineWidth',3); hold on; 
set(gca,'XTick',[0:2:12]);set(gca,'XTickLabel',{[]});
set(gca,'YTick',[4 6]);set(gca,'YTickLabel',{'4','6'});
xlim([0 12]); ylim([2 8]);
set(gca, 'fontsize', size_tick,'LineWidth',2)
ylabel('g_{ALC}','fontsize', size_label,'fontweight','bold');

%plot mean errors in reconstruction overlayed with base error calculated
%above
subplot(5,1,3);
errbaseplot = ones(size(i1))*eps_base(1); 
plot(Time(i1),mean_err(1,:),'m',...
       Time(i1),errbaseplot,'--k',...
       'LineWidth',3); 
set(gca,'XTick',[0:2:12]);set(gca,'XTickLabel',{[]});
xlim([0 12]);
set(gca, 'fontsize', size_tick,'LineWidth',2)
ylabel('\epsilon^2 _{F_{LC }}','fontsize', size_label,'fontweight','bold');


subplot(5,1,4);
errbaseplot = ones(size(i1))*eps_base(5); 
plot(Time(i1),mean_err(5,:),'m',...
       Time(i1),errbaseplot,'--k',...
    'LineWidth',3);
ylabel('\epsilon^2_{F_{W/R }}','fontsize', size_label,'fontweight','bold');
set(gca,'XTick',[0:2:12]);set(gca,'XTickLabel',{[]});
xlim([0 12]);
set(gca, 'fontsize', size_tick,'LineWidth',2)

subplot(5,1,5);
errbaseplot = ones(size(i1))*eps_base(11); 
plot(Time(i1),mean_err(11,:),'m',...
       Time(i1),errbaseplot,'--k',...
    'LineWidth',3); 
set(gca,'XTick',[0:2:12]);set(gca,'XTickLabel',{0:2:12});
xlim([0 12]);
ylabel('\epsilon^2 _{h }','fontsize', size_label,'fontweight','bold');
xlabel('Time (Hours)','fontsize', size_label,'fontweight','bold');
set(gca, 'fontsize', size_tick,'LineWidth',2)

return
