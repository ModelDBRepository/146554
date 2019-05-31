%This function generates 12 hours of data, sampled at 0.5 seconds, from the
%DB model using an RK4 integrator. The stochastic variable delta is included 

%DB model from: Diniz Behn and Booth, J Neurophysiol 103:1937-1953, 2010.

%Usage: Running this .m file will produce 'data_DB.mat'. The rest of the
%figure files will check use this generated data.

%Sleep-state for each 0.5 seconds is scored and included in the set.

%A noisy observation function is used to create a measurement set (12
%variables)

%Madineh Sedigh-Sarvestani, Penn State, Oct 2012
%m.sedigh.sarvestani@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Times,x,y,state,dT,P,Rs]= data_DB
%outputs: 
%Times (sampling times)
%x (12 dimensional matrix output of model: firing rates, neurotransmitter
%concentrations, homeostatic drive, and delta)
%y (12 dimensional noise-added version of x)
%state (vector of sleep-state for each data point, derived from threshholding the
%firing rates)
%dT (sampling time)
%P (struct which holds parameter values)
% Rs (measurement noise covariance)

dT=0.5; %integration sampling time in seconds
%arbitrary initial conditions
ic=[2 2 6 4 1 ... %firing rates
    0 0 1 1 1 ... %transmitter concentrations
    0.4 0]; %h and noise

totaltime=12*3600; %time in seconds

%generate data
[Times,x,state,P]=DB_Generate(totaltime,dT,ic) %generate data

%generate measurements ( a noisy version of each variable)
[dy,N]=size(x);
Rs=0.2^2*var(x'); %20% of variance of each variable
R=(Rs'*(ones(1,12)))'.*eye(dy,dy);
y=x+sqrtm(R)*randn(dy,N); % noisy data


%enforce >=0 constraint on all observations
for i=1:dy
    temp=find(y(i,:)<0);
    y(i,temp)=0;
    clear temp
end

clear R dy N
save data_DB_output.mat
return


%this function holds the state-space DB model
%uses and RK4 to integrate the equations and generate output
function [Times,x,state,P]=DB_Generate(totaltime,dT,ic)
%outputs: 
%Times (sampling times)
%x (12 dimensional matrix output of model: firing rates, neurotransmitter
%concentrations, homeostatic drive, and delta)
%state (sleep-state for each data point, derived from threshholding the
%firing rates)
%P (struct which holds parameter values)

%inputs 
%totaltime (total simulation time)
%dT (sampling, RK4 integration time)
%ic (model initial condition)

global nn %sampling time step, integration time stemp
dx=12; %number of variables in the model
% sampling time step (global variable) in seconds
dt=dT; %integration time=sampling time

N=totaltime/dT;%number of data samples in sample time steps (being fed by GUI)
nn=fix(dT/dt);  % the integration time step can be smaller than dT
Times = dT*(1:N); %setup Times vector
x=zeros(dx,N); %preallocate x
state=zeros(1,N);%preallocate state
%%%%%%%%%%%%%%%%%%%%%%%%%model parameters
P=OriginalDBParams;

%%%%%%%%%%%%%%%%%%%%%%%%%%initial value conditions
xx(:,1)=ic;
x(:,1)=ic;
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
     xx(12)= xx(12)+((8+ 0.1*randn(1,1)).*(poissrnd(0.003*dT))); %noise for LC,D

    %Do both integrations at once
    for i=1:nn 
         [k1]=DB_SS(xx,P);
         [k2]=DB_SS(xx+dt*k1/2,P);
         [k3]=DB_SS(xx+dt*k2/2,P);
         [k4]=DB_SS(xx+dt*k3,P);
         xx=xx+dt*(k1+(k2+k3)*2+k4)/6;
     end 

    x(:,n+1)=xx; %add last point to data vector
       
end;
%get state of vigilance
for n=1:length(x);
if x(1,n)>0.5 && x(4,n)<0.5 %if F_LC is high and F_R is low
    state(n)=1; %state=wake
else if x(1,n)>0.5 && x(4,n)>0.5 %if F_LC is high and F_R is high
        state(n)=3; %state=rem
    else
        state(n)=2; %state=nrem
    end
end
end


return


%DB state-space equations
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
deltaLC=x(12);
C_A=C_AR+C_AWR;

F=[x(1:5)]'; %firing rates in vector format
C=[x(6:10)]';%neurotransmitters in vector format

%steady state neurotransmitter concentrations 
C_ss=tanh(F./P.cgamma).*P.cnoise;
C_Dot=(C_ss-C)./P.ctau;

%homeostatic sleep drive
heavarg1=(F_LC+F_DR)-P.thetaW;
heavarg2=P.thetaW-(F_LC+F_DR);
hDot=((heavarg1>=0)*((1-h)/P.tauhw))-((heavarg2>=0)*(h/P.tauhs));

%noise input for LC and DR
deltaLCDot=-(0.1*deltaLC);


%summed neurotransmmitter input
cLC=P.gALC*C_A-P.gNLC*C_N-P.gGLC*C_G + deltaLC;
cDR=P.gADR*C_A-P.gSDR*C_S-P.gGDR*C_G+ deltaLC;
cVLPO=-P.gNVLPO*C_N-P.gSVLPO*C_S-P.gGVLPO*C_G;
cR=P.gAR*C_A-P.gNR*C_N-P.gSR*C_S-P.gGR*C_G;
cWR=P.gAWR*C_A-P.gGWR*C_G;
Cin=[cLC cDR cVLPO cR cWR];

%firing rate param
Fbeta=[P.betaLC P.betaDR -7*h P.betaR P.betaWR];

%firing rate
F_ss=P.Fmax.*(0.5.*(1+tanh((Cin-Fbeta)./P.Falpha)));
F_Dot=(F_ss-F)./P.Ftau;

%[F_LC,F_DR,F_VLPO,F_R,F_WR,C_N,C_S,C_G,C_AR,C_AWR,h,deltaLC]
x_dot=[F_Dot(:); C_Dot(:);...
     hDot(:); deltaLCDot(:)];
 return