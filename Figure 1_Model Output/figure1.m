%Reproduces Figure 1 in Sedigh-Sarvestani, Schiff, Gluckman,
%'Reconstructing Mammalian Sleep Dynamics with Data Assimiliation',
%PLoS Comp Biol, 2012 (In press). DOI:10.1371/journal.pcbi.1002788

%Fig 1B) Reconstructs unobserved dynamics from DB model of sleep.
%DB model from: Diniz Behn and Booth, J Neurophysiol 103:1937-1953, 2010.

%Fig 1D) Reconstructs unobserved dynamics from FBFD model of sleep.
%FBFD model from: Fleshner, Booth, Forger, Diniz Behn, Philos Transact A
%Math Phys Eng Sci. 2011 Oct 13;369(1952):3855-83.

%Usage: Running this .m file will produce a figure similar to Figure 1B,D of
%Sedigh-Sarvestani et al. 
%Make sure CD is '...\Figure Code\Figure 1_Model Output''

%the .mfiles to generate model data are stored in '\Figure Code\Generate
%Data\'

%Madineh Sedigh-Sarvestani, Penn State, Oct 2012
%m.sedigh.sarvestani@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot Figure 1D first
cd('../')%go up one foler
addpath(genpath(cd)); %add path
cd([cd '/Figure 1_Model Output']) %reset folder back \

if exist('data_FBFD_output.mat','file')
    load data_FBFD_output.mat %load already generated data
else
    [Times,x,y,state,dT,P,Rs]=data_FBFD; %alternatively re-generate data
end

size_label=18; %font for labels
size_tick=18; %font for ticks

%%now plot FBFD output
figure
hold on;
subplot(2,1,2);

%select range for plotting, the data contains 72 hours, but we only want to
%plot 36 hours
ranges=[12*3600/dT:48*3600/dT];
circ=2+(1.*sin((2*pi)*(1/(3600*24)).*Times)); %this is the 24 hour periodic factor CIRC from FBFD 2012 paper

%Setup lights-on/off background coloring 
WayBig = 3.2;
WaySmall = 0.8;

delta = round(length(ranges)/1000);
rangesArea = [ranges(1:delta:end) ranges(end)];

LightTime = WayBig*(circ>2);
DarkTime  = WayBig*(circ<2);

Yellow = [1 1 0];
LightGrey = [.90 .90 .90];
area(Times(rangesArea),LightTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',Yellow);
hold on;
area(Times(rangesArea),DarkTime(rangesArea),WaySmall,'LineStyle','none','FaceColor',LightGrey);
xlabel('Time (Hours)','fontsize', size_label,'fontweight','bold');
set(gca,'YTick',[1 2 3]);
set(gca,'YTickLabel',{'Wake','NREM','REM'},'fontsize', size_label,'fontweight','bold');

%now plot States of vigilance (see paper or data_FBDF.m for how this is
%derived from firing rates)
plot(Times(ranges),state(1,(ranges)),'k.','LineWidth',3);  hold on;
plot(Times(ranges),circ(ranges),'m--','LineWidth',2);
xlabel('Time (Hours)','fontsize', size_label,'fontweight','bold');
 ylim([0.8 3.2]);
 set(gca,'YTick',[1 2 3]);
 set(gca,'YTickLabel',{'Wake','NREM','REM'},'fontsize', size_label,'fontweight','bold');
 set(gca,'XTick',([12 18 24 30 36 42 48])*3600);
 set(gca,'XTickLabel',[18 0 6 12 18 0 6]);
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% now plot DB model (1 hour)
cd('../')%go up one folder
addpath(cd); %add path
cd([cd '/Figure 1_Model Output']) %reset folder back 


if exist('data_DB_output.mat','file')
    load data_DB_output.mat %load already generated data
else
    [Times,x,y,state,dT,P,Rs]= data_DB; %alternatively re-generate data
end


subplot(8,1,1);
plot(Times,x(1,:),'g','LineWidth',3); 
ylabel('F_{LC}','fontsize', size_label,'fontweight','bold');
ylim([0 7]); xlim([0 3600]);
set(gca,'YTick',[0 5]);
set(gca,'YTickLabel',[0 5]);
set(gca,'XTick',[0:1:2]*1800);
set(gca,'XTickLabel',{});
set(gca, 'fontsize', size_tick,'LineWidth',2)

subplot(8,1,2);
plot(Times,x(3,:),'r','LineWidth',3); 
ylabel('F_{VLPO}','fontsize', size_label,'fontweight','bold');
ylim([0 7]); xlim([0 3600]);
set(gca,'YTick',[0 5]);
set(gca,'YTickLabel',[0 5]);
set(gca,'XTick',[0:1:2]*1800);
set(gca,'XTickLabel',{});
set(gca, 'fontsize', size_tick,'LineWidth',2)

subplot(8,1,3);
plot(Times,x(4,:),'b','LineWidth',3); 
ylabel('F_{R}','fontsize', size_label,'fontweight','bold');
ylim([0 7]); xlim([0 3600]);
set(gca,'YTick',[0 5]);
set(gca,'YTickLabel',[0 5]);
set(gca,'XTick',[0:1:2]*1800);
set(gca,'XTickLabel',{});
set(gca, 'fontsize', size_tick,'LineWidth',2)

subplot(8,1,4);
plot(Times,state(1,:),'k','LineWidth',3); 
xlabel('Time (Hours)','fontsize', size_label,'fontweight','bold');
ylim([0.8 3.2]); xlim([0 3600]);
set(gca,'YTick',[1 2 3]);
set(gca,'YTickLabel',{'Wake','NREM','REM'},'fontsize',size_label,'fontweight','bold');
set(gca,'XTick',[0:1:2]*1800);
set(gca,'XTickLabel',[0:0.5:2]);
set(gca, 'fontsize', size_tick,'LineWidth',2)


