%Reproduces Figure 3 in Sedigh-Sarvestani, Schiff, Gluckman,
%'Reconstructing Mammalian Sleep Dynamics with Data Assimiliation',
%PLoS Comp Biol, 2012 (In press). DOI:10.1371/journal.pcbi.1002788

%Reconstructs unobserved dynamics from DB model of sleep.
%DB model from: Diniz Behn and Booth, J Neurophysiol 103:1937-1953, 2010.

%Then calculates the error in reconstruction of hidden variables, when each
%model variable is used as observable, then represents these errors as a
%matrix (the EOC).

%Usage: Running this .m file will produce a figure similar to Figure 3 of
%Sedigh-Sarvestani et al. 
%Make sure CD is '...\Figure Code\Figure 3_ EOC'

%Madineh Sedigh-Sarvestani, Penn State, Oct 2012
%m.sedigh.sarvestani@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=figure3
cd('../')%go up one foler
addpath(genpath(cd)); %add path
cd([cd '/Figure 3_ EOC']) %reset folder back 


% Generate DB data (without delta noise, see paper for explanation)
dT=0.5; totaltime=12*3600; %Generate 12 hours of data
%set initial condition
ic=[2 2 6 4 1 ... %firing rates
    0 0 1 1 1 ... %transmitter concentrations
    0.4]; %h 
[Times,x,P]=DB_Generate_nodelta(totaltime,dT,ic); %we won't use previously generated
%data (data_DB_output.mat) since that set has delta noise included


%generate measurements (a noisy version of each variable)
[dx,N]=size(x); [dy,N]=size(x);
vars=var(x');
Rs=0.2^2*vars;
R=(Rs'*(ones(1,dy)))'.*eye(dy,dy);
y=x+sqrtm(R)*randn(dy,N); % noisy data
%enforce >=0 constraint on all observations
for i=1:dy
    temp=find(y(i,:)<0);
    y(i,temp)=0;
    clear temp
end

%generate covariance inflation matrix CI
load CI_eps.mat
CI=CI_eps.*blkdiag(vars(1),vars(2),vars(3),vars(4),vars(5), ...
        vars(6),vars(7),vars(8),vars(9),vars(10),vars(11));


%track data using every variable as an observable
dobs=1; %dimension of obsevable

for i=1:dx %go through each state variable as observable
    varn=i
    R=(Rs([varn])'*(ones(1,dobs)))'.*eye(dobs,dobs); %measurement covariance matrix
    [xhat(:,:,i)]=UKF(x,y(i,:),R,dT,varn,CI); %reconstruct
end

%calculate EOC
for j=[1:dy] %observed
    mean_err_sq=mean((x(:,N/2:end)-xhat(:,N/2:end,j)).^2,2)./vars'; %calculate reconstruction error
    EOC(:,j)=(1./(mean_err_sq+1)); %calculate EOC
end

save figure3.mat
%now plot
plot_figure;

return

%This function generates data from DB model without delta noise
function [Times,x,P]=DB_Generate_nodelta(totaltime,dT,ic)
%outputs: 
%Times (sampling times)
%x (11 dimensional matrix output of model: firing rates, neurotransmitter
%concentrations, homeostatic drive)
%state (sleep-state for each data point, derived from threshholding the
%firing rates)
%P (struct which holds parameter values)

%inputs 
%totaltime (total simulation time)
%dT (sampling, RK4 integration time)
%ic (model initial condition)


global dt nn %sampling time step, integration time stemp
dx=11; %number of variables in DB without delta
dt=dT; %integration time

N=totaltime/dT;%number of data samples in sample tim steps (being fed by GUI)
nn=fix(dT/dt);  % the integration time step can be smaller than dT
Times = dT*(1:N);
x=zeros(dx,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%initial value conditions
xx(:,1)=ic;
x(:,1)=xx;

%%%%%%%%%%%%%%%%%%%%%%%%%model parameters
P=OriginalDBParams;

%%%%%%%%%%%%%%%%%%%%%%%%%neurotransmitter noise
dProbNoise = dT * 10;
Cnow=1;
state(1)=0;
%%%%%%%%%%%%%%%%%%%%%%%RK4 solver
for n=1:N-1;
    
    %set noise on neurotrans concentrations
    randomnum=rand(1,5);
    for i=1:5
        if (randomnum(i)<dProbNoise)
            CNow(i) = 1 + 0.1*randn(1,1);
        else
            CNow(i)=1;
        end
    end
    P.cnoise=[CNow];
    %for deterministic C:
    P.cnoise=[1, 1,1,1,1]; %don't set this to zero ever! If you want 
    %deterministic neurotrans, set all to 1.
    
    %normal noise
    %xx(12)= xx(12)+((8+ 0.1*randn(1,1)).*(poissrnd(0.003*dT))); %noise for LC,D 
     
    %RK4 integration
    for i=1:nn 
         [k1]=DB_SS_nodelta(xx,P);
         [k2]=DB_SS_nodelta(xx+dt*k1/2,P);
         [k3]=DB_SS_nodelta(xx+dt*k2/2,P);
         [k4]=DB_SS_nodelta(xx+dt*k3,P);
         xx=xx+dt*(k1+(k2+k3)*2+k4)/6;
     end 

    x(:,n+1)=xx; %add last point to data vector
end 
return


%this function holds the DB model state-space
function [x_dot]=DB_SS_nodelta(x,P)
%output:
%x_dot: 11 dimensional vector of instantaneous derivatives for firing rates, transmitter concentrations,
%h at time t=k+1

%inputs:
%x ( 11-dimensional data vector at time t=k)
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
%deltaLC=x(12);

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
%deltaLCDot=-(0.1*deltaLC);

%summed neurotransmmitter input
cLC=P.gALC*C_A-P.gNLC*C_N-P.gGLC*C_G; 
cDR=P.gADR*C_A-P.gSDR*C_S-P.gGDR*C_G;
cVLPO=-P.gNVLPO*C_N-P.gSVLPO*C_S-P.gGVLPO*C_G;
cR=P.gAR*C_A-P.gNR*C_N-P.gSR*C_S-P.gGR*C_G;
cWR=P.gAWR*C_A-P.gGWR*C_G;
Cin=[cLC cDR cVLPO cR cWR];

%firing rate param
Fbeta=[P.betaLC P.betaDR -7*h P.betaR P.betaWR];

%firing rate
F_ss=(P.Fmax.*(0.5.*(1+tanh((Cin-Fbeta)./P.Falpha))));
F_Dot=(F_ss-F)./P.Ftau;

%[F_LC,F_DR,F_VLPO,F_R,F_WR,C_N,C_S,C_G,C_AR,C_AWR,h]
x_dot=[F_Dot(:); C_Dot(:); hDot(:)];
return


%this is the UKF function, embedded functions carry out filter-model
%integration, unscented transform, etc.
function [xhat]=UKF(x,y,R,dT,varn,CI)
%outputs:
%xhat (11 dimensional vector of reconstructioned estimates of x)

%inputs:
%x(11 dimensional vector of DB model output)
%y(11 dimensional vector of observables, noise-added version of x)
%R (measurement noise covariance matrix)
%dT (sampling, RK4 integration time)
%varn (index of observed variables)
%CI (model covariance inflation matrix)


global dt nn %sampling time step, integration time stemp
dx=11; %dimension of DB state space without delta
dy=length(varn); % Dimensions of observation
fct='UKF_model'; % this is the model function F(x) used in filtering
obsfct='UKF_obsfct'; % this is the observation function G(x)

%RK4 integration setup
totaltime=length(x)*dT;
dt=dT; %dt is integration time, which can be smaller than dT, the sampling time
N=totaltime/dT;%number of data samples in sample time steps (being fed by GUI)
nn=fix(dT/dt);  %number of cycles of model intergration before we sample
Times = dT*(1:N);

%preallocate matrices
xhat=zeros(dx,N); % Preallocate estimated x
Pxx=zeros(dx,dx,N); % Prallocate Covariance in x
errors=zeros(dx,N); % Preallocate errors
Ks=zeros(dx,dy,N); % Preallocate Kalman gains

%Store which variable was observed and dT
P.varn=varn;
P.dT=dT;

%initialize xhat and Pxx
xhat(:,1)=[3 3 2.5 2.5 2.5 0.5 0.5 0.5 0.5 0.5 0.25];
xhat(varn,1)=x(varn,1); %initialize value of observed variables correctly
Pxx(:,:,1)=CI;

%%%%%%%%%%%%%%%%%%%%%%%%%model parameters
P=OriginalDBParams;

%Main loop for recursive UKF estimation
for k=2:N
    [xhat(:,k),Pxx(:,:,k),Ks(:,:,k)]=UKF_ut(xhat(:,k-1),Pxx(:,:,k-1),y(:,k),fct,obsfct,dx,R,CI,varn,P);
    errors(:,k)=sqrt(diag(Pxx(:,:,k)));
end; 
return

%this funciton selects the index of variables in G(x), which is a noise
%added version of F(x)
function observ=UKF_obsfct(y,varn)
%outputs:
%observ: noisy observations vector, whose number of rows=length(varn)

%inputs:
%y(11 dimensional vector of observables, noise-added version of x)
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
min=[0 0 0 0 0 0 0 0 0 0 0]';
max=[6.5 6.5 5 10 5 1 1 1 1 1 1]';
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
        k1=dt.*DB_SS_nodelta(xn2(:,NumCol),P);
        k2=dt.*DB_SS_nodelta(xn2(:,NumCol)+k1/2,P);
        k3=dt.*DB_SS_nodelta(xn2(:,NumCol)+k2/2,P);
        k4=dt.*DB_SS_nodelta(xn2(:,NumCol)+k3,P);
        xn2(:,NumCol)=xn2(:,NumCol)+k1./6+k2./3+k3./3+k4./6;
    end
end
X=[xn2];
return

function[]= plot_figure
%figure3; 
%% plot EOC matrix
load figure3.mat EOC
size_label=18;
size_tick=18;

cd('../')%go up one foler
addpath(genpath(cd)); %add path
cd([cd '/Figure 3_ EOC']) %reset folder back 

figure;
labels={'F_{LC}','F_{DR}','F_{VLPO}','F_{R}','F_{W/REM}','C_{N}',...
    'C_{S}','C_{G}','C_{A(R)}','C_{A(W/R})','h'};


%add empty column and row necessary for plotting
temp=flipud(EOC);
[N2,M]=size(temp);
plot_dist=[[temp,zeros(N2,1)];[zeros(1,M+1)]];

pcolor(plot_dist);
colormap(flipud(gray));
ylabel('Reconstructed','fontsize', size_label,'fontweight','bold');
xlabel('Observed','fontsize', size_label,'fontweight','bold');
zlabel('ln(Distance)')
set(gca,'XTick',(1:11)+0.5);
set(gca,'YTick',(1:11)+0.5);
set(gca,'XTickLabel',labels');
set(gca,'YTickLabel',fliplr(labels)');
set(gca, 'fontsize', 14,'LineWidth',2)
[hx,hy] = format_ticks(gca,labels',fliplr(labels)',[],[],90,0);
set(gca,'XAxisLocation','top');
colorbar;
set(gca, 'fontsize', 14,'LineWidth',2)
return
