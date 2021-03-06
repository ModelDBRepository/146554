%this function holds parameters of the DB model of sleep, duplicated from:
%Diniz Behn and Booth, J Neurophysiol 103:1937-1953, 2010.
function P= OriginalDBParams

P.cgamma=[5 5 4 3 3]; %
P.ctau=[25 25 10 10 10];

P.gALC=3.5; P.gNLC=1.5; P.gGLC=1.5;
P.gADR=3.5; P.gSDR=1.5; P.gGDR=1.5;
P.gAR=2.5; P.gNR=3.5; P.gSR=3.5; P.gGR=1.25;
P.gAWR=1; P.gGWR=1.7;
P.gNVLPO=2; P.gSVLPO=2; P.gGVLPO=0.5;

%firing rate parameters (LC,DR,VLPO,R,WR)
P.Ftau=[25 25 10 1 10];

P.Fmax=[6.5 6.5 5 5 5];
P.Falpha=[0.75 0.75 0.25 0.25 0.25];
P.betaLC=2; P.betaDR=2; P.betaR=-0.5; P.betaWR=-0.2;
%we've left betaVLPO out because it is dependent on h

%homeostatic sleep constants
P.thetaW=3; P.tauhs=320; P.tauhw=600;
P.tauhs_i=1/320; P.tauhw_i=1/600;
return;
