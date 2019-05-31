%Reproduces Figure 6 in Sedigh-Sarvestani, Schiff, Gluckman,
%'Reconstructing Mammalian Sleep Dynamics with Data Assimiliation',
%PLoS Comp Biol, 2012 (In press). DOI:10.1371/journal.pcbi.1002788

%Reconstructs unobserved dynamics from FBFD model of sleep.
%FBFD model from: Fleshner, Booth, Forger, Diniz Behn, Philos Transact A
%Math Phys Eng Sci. 2011 Oct 13;369(1952):3855-83.

%Uses as UKF model the DB model which does not consist SCN components:
%DB model from: Diniz Behn and Booth, J Neurophysiol 103:1937-1953, 2010.

%Estimates unobserved (and unmodeled) parameter gSCN

%Usage: Running this .m file will produce a figure similar to Figure 6 of
%Sedigh-Sarvestani et al. 
%Make sure CD is '...\Figure Code\Figure 6_Parameter Estimation SCN'

%Options: You can use previously generated model data, from which noisy
%observations are obtained. This is stored as 'data_FBFD_output.mat' in
%'\Figure Code\Generate Data\' and is used throughout the paper.
%Alternatively, you can re-generate model data from scratch.

%Note: this function takes a while to run as 72 hours of data at dT=0.5
%seconds sampling is reconstructed.

%Madineh Sedigh-Sarvestani, Penn State, Oct 2012
%m.sedigh.sarvestani@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=figure6
%%%%%%%%%%%%%%%%%%%%%%%%%%get data
cd('../')%go up one foler
addpath(genpath(cd)); %add path

load CI_eps.mat

cd([cd '/Figure 6_Parameter Estimation SCN']) %reset folder back 

%if exist('data_FBFD_output.mat','file')
  % load data_FBFD_output.mat %load already generated data
%else
    [Times,x,y,state,dT,P,Rs]= data_FBFD; %alternatively re-generate data
%end

%% create new system for tracking (match 14D FBFD variables to 12D DB varialbes)
xnew(1:5,:)=x(1:5,:); %match firing rates
xnew(6:10,:)=x(7:11,:); %match neurotransmitters
xnew(11,:)=x(13,:); %match h
xnew(12,:)=x(14,:); %match delta

%generate measurements (a noisy version of each variable)
[dy,N]=size(xnew);
Rs=0.2^2*var(xnew');
R=(Rs'*(ones(1,dy)))'.*eye(dy,dy);
y=xnew+sqrtm(R)*randn(dy,N); % noisy data

%enforce >=0 constraint on all observations
for i=1:dy
    temp=find(y(i,:)<0);
    y(i,temp)=0;
    clear temp
end

%now setup Covariance inflation (CI) matrix

vars=var(xnew');%variance of all model variables
%CI is a diagonal matrix where the diag entries=CI_eps*variance(variable)
CI=CI_eps.*blkdiag(vars(1),vars(2),vars(3),vars(4),vars(5), ...
        vars(6),vars(7),vars(8),vars(9),vars(10),vars(11),vars(12));

%optimized covariance inflation (see Figure 4)
CI_optimized=CI;
CI_optimized(4,4)=0; %F_R
CI_optimized(12,12)=CI(12,12)*4000; %delta

%pick observed variables
varn=[1 4]; %observe F_LC and F_R
y_obs=y(varn,:); %observe F_LC and F_R
dy=length(varn);
R=(Rs([varn])'*(ones(1,dy)))'.*eye(dy,dy); %measurement  covariance matrix 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%parameter estimation+UKF
dx=12; %12 DB
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

divd=divd;
%how often we force sync trajectories onto xhat
% set up offsets - assume i2(n)-i1(n)=const
ForceFix_step = 120; %seconds

%setup how much parameter is allowed to vary on each update
DeltaMax = 0.02; %(2*pi/(24*3600))*pWindow;

%preallocate matrices
y=zeros(dx,N); %test trajectories
xhat=zeros(dx,N); %UKF reconstruction

% initial value for parameter we'll estimate
p_initial=[1];
p=zeros(1,1+divd);
dp=zeros(1,divd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%begin
for n=1:divd %go through all windows
    if n<=1 %for first iteration use initial values
        p_input=p_initial;
        ic=[3 3 2.5 2.5 2.5 0.5 0.5 0.5 0.5 0.5 0.25 0];
    else
        p_input=p(n-1); %update parameter
        ic=xhat(:,i1(n-1)+pOverlap); %update UKF ic
    end
    
    
    %Iterate between UKF estimation and parameter estimation
    %UKF estimation using p_input as paramater value for gALC
    [xhat(:,i1(n):i2(n))]=UKF(y_obs(:,i1(n):i2(n)),R,dT,varn,CI_optimized,p_input,ic);
    
    %setup force sync (xhat to trajectory) values
    DivLength = i2(n)-i1(n);
    RangeStart = 0:ForceFix_step:DivLength-5;
    RangeEnd = circshift(RangeStart',-1)';
    RangeEnd(end)=DivLength;
    NSteps = length(RangeStart);
    
    p_new_plus = p_input+DeltaMax; %current param +delta
    p_new_minus = p_input-DeltaMax; %current param -delta

    %generate the FIRST test trajectory with param=current param (use force sync)
    for step = NSteps:-1:1
        iS = i1(n)+RangeStart(step);
        iE = i1(n)+RangeEnd(step);
        [yy(:,iS:iE)]...
            =test_traj(xhat(:,iS:iE),p_input,{'pSCN'});
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
            =test_traj(xhat(:,iS:iE),p_new_plus,{'pSCN'});
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
            =test_traj(xhat(:,iS:iE),p_new_minus,{'pSCN'});
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
            dp(n) = DeltaMax; %otherwise update parameter by +delta
            y(:,iS:iE)=yyplus(:,iS:iE); %and keep corresponding trajectory
        else
            dp(n) = -DeltaMax; %otherwise update parameter by -delta
            y(:,iS:iE)=yyminus(:,iS:iE);%and keep corresponding trajectory
        end
    end
    
    p_input = p_input+dp(n); %update current parameter for next UKF step     
    p(n) = p_input;% store in a separate vector

    true_p(n)=mean(x(12,i1(n):i2(n))); %get average SCN output to compare to estimated value
end


%now create a full vector of true pSCN values for plotting
p=[p_initial,p]; %estimated values
for n=1:divd
    p_full(i1(n):i2(n))=p(n); %our estimates vector
end

% calculate average trajectory 
p_best=x(12,1:length(p_full)); %the underlying true GABA SCN output
mean_true=mean(p_best,1); %average over window
mean_estimate=mean(p_full,1);%mean of estimate over window

%calculate mean errors in reconstruction of variables
for n=1:divd
    error=(xnew(:,i1(n):i2(n))-xhat(:,i1(n):i2(n))).^2; %MSE for each Twin
    mean_err(:,n)=mean(error,2); %MSE averaged over Twin
    var_err2(:,n)=var(xnew');%variance of each variable
    mean_err_n2(:,n)=mean_err(:,n)./var_err2(:,n); %MSE averaged over Twin normalized by var
end

%clearvars -except state state_hat xhat y q_full q initial varn pWindow pOverlap DeltaMax RangeStep M turn savedState x state Time
save ('figure6.mat')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%model parameters
Params=OriginalDBParams_SCN;
Params.pSCN=parameter; %this is the parameter we'll estimate

%Main loop for recursive estimation
for k=2:N
    [xhat(:,k),Pxx(:,:,k),Ks(:,:,k)]=UKF_ut(xhat(:,k-1),Pxx(:,:,k-1),y(:,k),fct,obsfct,dx,R,CI,varn,Params);
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
  
%set max and means to constraint ouput
min=[0 0 0 0 0 0 0 0 0 0 0 0 ]';
max=[6.5 6.5 5 5 5 1 1 1 1 1 1 30]';


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
    for n=1:nn%if dt!=dT, then nn would be >1
        %RK4
        k1=dt.*DB_SS_mod(xn2(:,NumCol),P);
        k2=dt.*DB_SS_mod(xn2(:,NumCol)+k1/2,P);
        k3=dt.*DB_SS_mod(xn2(:,NumCol)+k2/2,P);
        k4=dt.*DB_SS_mod(xn2(:,NumCol)+k3,P);
        xn2(:,NumCol)=xn2(:,NumCol)+k1./6+k2./3+k3./3+k4./6;
    end
end
X=[ xn2];
return

%this function holds the DB model state-space
function [x_dot]=DB_SS_mod(x,Params)
%output:
%x_dot: 12 dimensional vector of instantaneous derivatives for firing rates, transmitter concentrations,
%h and delta at time t=k+1

%inputs:
%x ( 12-dimensional data vector at time t=k)
%Params (parameter values stored in struct format)


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
    
    
    F=[x(1:5)]';
    C=[x(6:10)]';
    
    %steady state neurotransmitter concentrations
    C_ss=tanh(F./Params.cgamma).*1;
    C_Dot=(C_ss-C)./Params.ctau;
    
    %homeostatic sleep drive
    heavarg1=(F_LC+F_DR)-Params.thetaW;
    heavarg2=Params.thetaW-(F_LC+F_DR);
    hDot=((heavarg1>=0)*((1.1-h)/Params.tauhw))-((heavarg2>=0)*(h/Params.tauhs));
    
    
    %delta noise input to LC and DR
    deltaLCDot=-(0.1*deltaLC);

    %firing rate param
    Fbeta=[Params.betaLC Params.betaDR (Params.beta_h)*h Params.betaR Params.betaWR];
    
    %summed neurotransmmitter input
    cLC=Params.gALC*C_A-Params.gNLC*C_N-Params.gGLC*C_G-Params.gGSCNLC*Params.pSCN-deltaLC ; %took out noise terms
    cDR=Params.gADR*C_A-Params.gSDR*C_S-Params.gGDR*C_G-Params.gGSCNDR*Params.pSCN-deltaLC;
    cVLPO=-Params.gNVLPO*C_N-Params.gSVLPO*C_S-Params.gGVLPO*C_G+Params.gGSCNVLPO*Params.pSCN;
    cR=Params.gAR*C_A-Params.gNR*C_N-Params.gSR*C_S-Params.gGR*C_G-Params.gGSCNR*Params.pSCN;
    cWR=Params.gAWR*C_A-Params.gGWR*C_G;
    Cin=[cLC cDR cVLPO cR cWR];
    
    %firing rate
    F_ss=(Params.Fmax.*(0.5.*(1+tanh((Cin-Fbeta)./Params.Falpha))));
    F_Dot=(F_ss-F)./Params.Ftau;
        
    x_dot=[F_Dot(:); C_Dot(:); hDot(:);deltaLCDot];
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


P=OriginalDBParams_SCN; %fix rest of parameters
%set the value for pSCN, the parameter we're estimating
names=fieldnames(P);
for i=1:dqsync
    index=strmatch(paramnames{i},names);
    P=setfield(P,names{index},params(i));
end

%main loop for model (trajectory) integration
for n=1:N-1;
    yy=y(:,n);
    %error blowup protection (put boundaries on F and C values)
    min=[0 0 0 0 0 0 0 0 0 0 0 0]';
    max=[6.5 6.5 6 6 6 1 1 1 1 1 1 30]';
    
    yy(1:dx)=yy(1:dx).*(yy(1:dx)>min)+(min.*(yy(1:dx)<min));
    yy(1:dx)=yy(1:dx).*(yy(1:dx)<max)+(max.*(yy(1:dx)>max));
    
    %RK4
    for i=1:nn
        k1=dt*DB_SS_mod(yy,P);
        k2=dt*DB_SS_mod(yy,P);
        k3=dt*DB_SS_mod(yy,P);
        k4=dt*DB_SS_mod(yy,P);
        yy=yy+k1/6+k2/3+k3/3+k4/6;
    end;
    y(:,n+1)=yy; %add last data point to vector
end;

return;

function[]=plot_figure

load figure6.mat
size_label=18;
size_tick=18;

%now plot
Time_plot=Times./(3600);
N=length(mean_true);
ranges=1:N;

%plot color
circ=2+(2.*sin((2*pi)*(1/(3600*24)).*Times));
WayBig = 10;
WaySmall = -10;
delta = round(length(ranges)/1000);
rangesArea = [ranges(1:delta:end) ranges(end)];

LightTime = WayBig*(circ>2);
DarkTime  = WayBig*(circ<2);
Yellow = [1 1 0];
LightGrey = [.90 .90 .90];

figure;
%plot trajectories of observable and reconstructed
a7=subplot(6,3,1);
area(Time_plot(rangesArea),LightTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',Yellow);hold on;
area(Time_plot(rangesArea),DarkTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',LightGrey);

plot(Time_plot(1:N),y_obs(1,1:N),'b.');  hold on;
plot(Time_plot(1:N),xhat(1,1:N),'r','LineWidth',2); 
ylim([0 7]); ylabel('F_{LC}','fontsize', size_label,'fontweight','bold');
xlim([12.5 13.5]); set(gca,'XTick',[13]);set(gca,'XTickLabel',{''});
set(gca, 'fontsize', size_tick,'LineWidth',2)
% 
a8=subplot(6,3,2);
area(Time_plot(rangesArea),LightTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',Yellow);hold on;
area(Time_plot(rangesArea),DarkTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',LightGrey);

plot(Time_plot(1:N),y_obs(1,1:N),'b.');hold on;
plot(Time_plot(1:N),xhat(1,1:N),'r','LineWidth',2); 
ylim([0 7]);set(gca,'YTickLabel',{''});
xlim([29.5 30.5]);set(gca,'XTick',[30]);set(gca,'XTickLabel',{''});
set(gca, 'fontsize', size_tick,'LineWidth',2)
% 
a9=subplot(6,3,3)
area(Time_plot(rangesArea),LightTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',Yellow);hold on;
area(Time_plot(rangesArea),DarkTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',LightGrey);

plot(Time_plot(1:N),y_obs(1,1:N),'b.');hold on;
plot(Time_plot(1:N),xhat(1,1:N),'r','LineWidth',2); 
xlim([45.5 46.5]);set(gca,'XTick',[46]);set(gca,'XTickLabel',{''});
ylim([0 7]);set(gca,'YTickLabel',{''});
set(gca, 'fontsize', size_tick,'LineWidth',2)
% 
% %plot trajectories of true and UKF reconstructed
a1=subplot(6,3,4)
area(Time_plot(rangesArea),LightTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',Yellow);hold on;
area(Time_plot(rangesArea),DarkTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',LightGrey);

plot(Time_plot(1:N),xhat(3,1:N),'r.');  hold on;
plot(Time_plot(1:N),xnew(3,1:N),'k','LineWidth',2); 
ylim([0 7]);ylabel('F_{VLPO}','fontsize', size_label,'fontweight','bold');
xlim([12.5 13.5]); set(gca,'XTick',[13]);set(gca,'XTickLabel',{'19'});
set(gca, 'fontsize', size_tick,'LineWidth',2)
% 
a2=subplot(6,3,5)
area(Time_plot(rangesArea),LightTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',Yellow);hold on;
area(Time_plot(rangesArea),DarkTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',LightGrey);

plot(Time_plot(1:N),xhat(3,1:N),'r.');hold on;
plot(Time_plot(1:N),xnew(3,1:N),'k','LineWidth',2); 
xlim([29.5 30.5]);set(gca,'XTick',[30]);set(gca,'XTickLabel',{'12'});
ylim([0 7]);set(gca,'YTickLabel',{''});
set(gca, 'fontsize', size_tick,'LineWidth',2)
% 
a3=subplot(6,3,6)
area(Time_plot(rangesArea),LightTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',Yellow);hold on;
area(Time_plot(rangesArea),DarkTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',LightGrey);

plot(Time_plot(1:N),xhat(3,1:N),'r.');hold on;
plot(Time_plot(1:N),xnew(3,1:N),'k','LineWidth',2); 
xlim([45.5 46.5]);set(gca,'XTick',[46]);set(gca,'XTickLabel',{'3'});
ylim([0 7]);set(gca,'YTickLabel',{''});
set(gca, 'fontsize', size_tick,'LineWidth',2)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%now plot parameter estimates
subplot(6,1,3);
area(Time_plot(rangesArea),LightTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',Yellow);hold on;
area(Time_plot(rangesArea),DarkTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',LightGrey);

plot(Time_plot(1:N),mean_true,'k','LineWidth',3); hold on;
plot(Time_plot(1:N),mean_estimate,'m','LineWidth',3)
set(gca,'XTick',[12:6:48]);set(gca,'XTickLabel',{''});
ylabel('C_{G,SCN}','fontsize', size_label,'fontweight','bold');
xlim([12 48]);
set(gca, 'fontsize', size_tick,'LineWidth',2);
ylim([0 1]);set(gca,'YTick',[0 0.5 1]);set(gca,'YTickLabel',{[0 0.5 1]});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%now plot errors in reconstruction
subplot(6,1,4);
area(Time_plot(rangesArea),LightTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',Yellow);hold on;
area(Time_plot(rangesArea),DarkTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',LightGrey);

plot(Time_plot(i1),mean_err_n2(1,:),'m','LineWidth',3);
line([0 Time_plot(i1(end))],[R(1,1) R(1,1)],'Color','b','LineWidth',4);
set(gca,'XTick',[12:6:48]);set(gca,'XTickLabel',{''});
set(gca, 'fontsize', size_tick,'LineWidth',2)
ylabel('\epsilon^2 _{F_{LC }}','fontsize', size_label,'fontweight','bold');
xlim([12 48]);
ylim([0 0.3]);set(gca,'YTick',[0 0.3]);set(gca,'YTickLabel',{[0 0.3]});

subplot(6,1,5);
area(Time_plot(rangesArea),LightTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',Yellow);hold on;
area(Time_plot(rangesArea),DarkTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',LightGrey);

plot(Time_plot(i1),mean_err_n2(5,:),'m','LineWidth',3); 
set(gca,'XTick',[12:6:48]);set(gca,'XTickLabel',{''});
set(gca, 'fontsize', size_tick,'LineWidth',2)
ylabel('\epsilon^2 _{F_{W/R }}','fontsize',size_label,'fontweight','bold');
xlim([12 48]);
ylim([0 0.3]);set(gca,'YTick',[0 0.3]);set(gca,'YTickLabel',{[0 0.3]});

subplot(6,1,6);
area(Time_plot(rangesArea),LightTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',Yellow);hold on;
area(Time_plot(rangesArea),DarkTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',LightGrey);

plot(Time_plot(i1),mean_err_n2(11,:),'m','LineWidth',3);
set(gca,'XTick',[12:6:48]);set(gca,'XTickLabel',{[18 0 6 12 18 0 6]});
set(gca, 'fontsize', size_tick,'LineWidth',2)
ylabel('\epsilon^2 _{h }','fontsize', size_label,'fontweight','bold');
xlim([12 48]);
ylim([0 0.3]);set(gca,'YTick',[0 0.3]);set(gca,'YTickLabel',{[0 0.3]});

linkaxes([a1 a7 ],'x')
linkaxes([a2  a8],'x')
linkaxes([a3  a9],'x')

return
