%This function generates 72 hours of data, sampled at 0.5 seconds, from the
%FBFD model using an RK4 integrator. The stochastic variable delta is included 

%FBFD model from: Fleshner, Booth, Forger, Diniz Behn, Philos Transact A
%Math Phys Eng Sci. 2011 Oct 13;369(1952):3855-83.

%Usage: Running this .m file will produce 'data_FBFD.mat'. The rest of the
%figure files will use this generated data.

%Sleep-state for each 0.5 seconds is scored and included in the set.

%A noisy observation function is used to create a measurement set (12
%variables)

%Madineh Sedigh-Sarvestani, Penn State, Oct 2012
%m.sedigh.sarvestani@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Times,x,y,state,dT,P,Rs]=data_FBFD

%outputs: 
%Times (sampling times)
%x (14 dimensional matrix output of model: firing rates, neurotransmitter
%concentrations, homeostatic drive, and delta)
%y (14 dimensional noise-added version of x)
%state (sleep-state for each data point, derived from threshholding the
%firing rates)
%dT (sampling time)
%P (struct which holds parameter values)
% Rs (measurement noise pre-multiplier (gets multiplied by variance of each
% variable))


dT=0.5; %sampling time step, integration time stemp
ic=[2 2 6 4 1 1 ... %firing rates
    0 0 1 1 1 1 ... %transmitter concentrations
    0.4 0]; %h and noise

totaltime=72*3600; %time in seconds

%generate data
[Times,x,state,P]=FBFD_Generate(totaltime,dT,ic) %generate data

%generate measurements ( a noisy version of each variable)
[dy,N]=size(x);
Rs=0.2^2*var(x');
R=(Rs'*(ones(1,14)))'.*eye(dy,dy);
y=x+sqrtm(R)*randn(dy,N); % noisy data


%enforce >=0 constraint on all observations
for i=1:dy
    temp=find(y(i,:)<0);
    y(i,temp)=0;
    clear temp
end

clear R dy N
save data_FBFD_output.mat

return


%this function holds the state-space FBFD model
%uses and RK4 to integrate the equations and generate output
function [Times,x,state,P]=FBFD_Generate(totaltime,dT,ic) %generate data
%outputs: 
%Times (sampling times)
%x (14 dimensional matrix output of model: firing rates, neurotransmitter
%concentrations, homeostatic drive, and delta)
%state (sleep-state for each data point, derived from threshholding the
%firing rates)
%P (struct which holds parameter values)

%inputs 
%totaltime (total simulation time)
%dT (sampling, RK4 integration time)
%ic (model initial condition)

global nn
dx=14; %number of variables in the model
% sampling time step (global variable) in seconds
dt=dT; %integration time

N=totaltime/dT;%number of data samples in sample time steps (being fed by GUI)
nn=fix(dT/dt);  % the integration time step can be smaller than dT
x=zeros(dx,N); %preallocate x
Times = dT*(1:N);%setup Times vector
state=zeros(1,N); %preallocate state
%%%%%%%%%%%%%constants from paper%%%%%%%%%%%%%%%%%
P=OriginalFBFDParams;
%%%%%%%%%%%%%%%%%%%%%%%%%%initial value conditions
xx(:,1)=ic;
x(:,1)=ic;


%%%%%%%%%%%%%%%%%%%%%%%%%neurotransmitter noise
dProbNoise = dT * 10;
Cnow=1;

state(1)=0;
%%%%%%%%%%%%%%%%%%%%%%%RK4 solver
for n=1:N-1;
    
    %use below for variable neurotransmitter release (see Diniz Behn 2010)
    %     randomnum=rand(1,6);
    %         for i=1:6
    %             if (randomnum(i)<dProbNoise)
    %                 CNow(i) = 1 + 0.1*randn(1,1);
    %             else
    %                 CNow(i)=1;
    %             end
    %         end
    %         P.cnoise=[CNow];
    TimeNow = Times(n);
    period=60*60*24;
    P.CIRC=(sin(2*pi*(1/(period))*TimeNow)); %period is 3600
    
    P.cnoise=[1, 1, 1, 1, 1,1];
    
    
    xx(14)= xx(14)+(((16.5-12*P.CIRC)+ 0.1*randn(1,1)).*(poissrnd(0.0004*dT))); %noise for LC,D
    
    %Do both integrations at once
    for i=1:nn
        [k1]=FBFDIntegrate(xx,P);
        [k2]=FBFDIntegrate(xx+dt*k1/2,P);
        [k3]=FBFDIntegrate(xx+dt*k2/2,P);
        [k4]=FBFDIntegrate(xx+dt*k3,P);
        xx=xx+dt*(k1+(k2+k3)*2+k4)/6;
    end
    
    x(:,n+1)=xx; %add last point to data vector
end   

    %get state (rule is different than DB)
    for n=1:length(x);
        if x(1,n)>4 %if F_LC is high
            state(n)=1; %wake
        else if  x(4,n)>4 %if F_R is high
                state(n)=3; %REM
            else
                state(n)=2; %NREM
            end
        end
    end
return
    
%DB state-space equations    
function [x_dot]=FBFDIntegrate(x,P)
%output:
%x_dot: 14 dimensional vector of instantaneous derivatives for firing rates, transmitter concentrations,
%h and delta at time t=k+1

%inputs:
%x (14-dimensional data vector at time t=k)
%P (parameter values stored in struct format)


%what comes in:
%[F_LC,F_DR,F_VLPO,F_R,F_WR,C_N,C_S,C_G,C_AR,C_AWR,h]
F_LC=x(1);
F_DR=x(2);
F_VLPO=x(3);
F_R=x(4);
F_WR=x(5);
F_SCN=x(6);

C_N=x(7);
C_S=x(8);
C_G=x(9);
C_AR=x(10);
C_AWR=x(11);
C_GSCN=x(12);
h=x(13);
deltaLC=x(14);
C_A=C_AR+C_AWR;

SYN=P.gASCN*(C_A)+P.gSSCN*C_S;

%steady state neurotransmitter concentrations (there should be a noise term
%here)
F=[x(1:6)]';
C=[x(7:12)]';

C_ss=tanh(F./P.cgamma).*P.cnoise;
C_Dot=(C_ss-C)./P.ctau;

%homeostatic sleep drive
heavarg1=(F_LC+F_DR)-P.thetaW;
heavarg2=P.thetaW-(F_LC+F_DR);
hDot=((heavarg1>=0)*((P.Hmax-h)/P.tauhw))-((heavarg2>=0)*(h/P.tauhs));

%firing rate param
Fbeta=[P.betaLC P.betaDR -P.k*h P.betaR P.betaWR P.betaSCN];

% %noise input for LC and DR
deltaLCDot=-((1/(10))*deltaLC);

%summed neurotransmmitter input
cLC=P.gALC*C_A-P.gNLC*C_N-P.gGLC*C_G -P.gGSCNLC*C_GSCN-deltaLC;
cDR=P.gADR*C_A-P.gSDR*C_S-P.gGDR*C_G- P.gGSCNDR*C_GSCN-deltaLC ;
cVLPO=-P.gNVLPO*C_N-P.gSVLPO*C_S-P.gGVLPO*C_G+P.gGSCNVLPO*C_GSCN;
cR=P.gAR*C_A-P.gNR*C_N-P.gSR*C_S-P.gGR*C_G-P.gGSCNR*C_GSCN;
cWR=P.gAWR*C_A-P.gGWR*C_G;
cSCN=P.CIRC+SYN;
Cin=[cLC cDR cVLPO cR cWR cSCN];

%firing rate
F_ss=P.Fmax.*(0.5.*(1+tanh((Cin-Fbeta)./P.Falpha)));
F_Dot=(F_ss-F)./P.Ftau;


x_dot=[F_Dot(:); C_Dot(:);...
    hDot(:); deltaLCDot];

    return;
    
