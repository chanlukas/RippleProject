% Script For Semi Automated Ripple Detection for Recordings with linear
% Silicon Probes, Jan Klee 23.10.16
   

%% lets establish the physical channel order tip to base, so that we do not get confused at later stages

clear all

%ProbeBase2TipOmnetics=flip([20,35,21,34,22,33,23,32,24,31,25,30,26,29,27,28,2,17,3,16,4,15,5,14,9,10,8,11,6,13,7,12]);%100micron
ProbeBase2TipOmnetics=([20,17,2,35,21,16,3,34,22,15,4,33,23,14,5,32,24,13,6,31,25,12,7,30,26,11,8,29,9,28,27,10]);%50micron
%OmneticsToIntan=[nan,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,nan,nan,24,25,26,27,28,29,30,31,0,1,2,3,4,5,6,7,nan];
OmneticsToIntan=[nan,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,nan,nan,7,6,5,4,3,2,1,0,31,30,29,28,27,26,25,24,nan];
%IntanToOmnetics=[27,28,29,30,31,32,33,34,35,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,20,21,22,23,24,25,26];

for i=1:32;   
    channelOrder(i)=OmneticsToIntan(ProbeBase2TipOmnetics(i))+1;
end


%% load LFP 

%define the filter, in this case a 3D butterworth low-pass filter
%(at 500hz)

Cut=500; %cut off Frequency
Fn=30000; %Sampling Rate Raw Data
[b,a]   = butter(3,Cut/Fn,'low');

samplerate=2000;   % desired LFP Sampling Rate

LFPchannels=1:32;  % which channels to include, avoid ACC channels here

for i = 1:length(LFPchannels) %loops through all Channels
    
filename=['100_CH',num2str(channelOrder(i)),'.continuous']; %creates the filename for each channel

[Data] = load_open_ephys(filename); % open ephys function reads individual channel voltage trace data into vector "Data"

DataDown=double(Data(1:(30000/samplerate):end)); % downsample 'Data' to 2000hz by taking every 15th sample and convert to 'double' data format (because that is needed for filtering)

LFPData(i,:)=filter(b,a,DataDown); % filtering downsapled data, with filter settings above, and saving it in as "LFPData"
end

figure(1)         
plotmat(LFPData(:,400000:402000)); 

%% load and plot small strech of MUAto manually look for cells 
% and determine cell layer

%define the filter, in this case a 3D butterworth High-pass filter
%(at 500hz)

Cut=500; %cut off Frequency
Fn=30000; %Sampling Rate Raw Data
[d,c]   = butter(3,Cut/Fn,'high');

MUAchannels=1:32;  % which channels to include, avoid ACC channels here

for i = 1:length(MUAchannels) %loops through all Channels
    
filename=['100_CH',num2str(channelOrder(i)),'.continuous']; %creates the filename for each channel

[Data] = load_open_ephys(filename); % open ephys function reads individual channel voltage trace data into vector "Data"

MUAData(i,:)=filter(d,c,double(Data(6000000:6030000))); % filtering downsapled data, with filter settings above, and saving it in as "LFPData"
end

VRev=repmat(mean(MUAData),32,1);
MUAData1=MUAData-VRev;

figure(2)
plotmat(MUAData1); 
%%
!! carefull also applies not for all Sessions approach!!
Ctx=1; % Manually Select Cortial Channel for TFA (use MUAplot to find out)
Hpc=8; % Manually Select Hpc Pyramidal Layer Channel for TFA (")


%% Load Sleep Scoring
[Sleep,SleepLong,EEGslow,AvgACC,Mov]=SleepScoring(Ctx,2,0);

%% Figure out best way for Ripple Detection

% 1. Approach:  Filtering Data in Ripple Frequency Band

%define the filter, in this case a 3D butterworth bandpass filter for
%Ripple Frequency band

CutLow=150; %Low Cut off
CutHigh=250; %High Cut Off   
Fn=30000; % Sampling Rate Raw Data
[b,a]   = butter(3,[CutLow CutHigh]/Fn,'bandpass');
samplerate=2000;   % desired LFP Sampling Rate

LFPchannels=1:32;  % which channels to include, avoid ACC channels here

for i =1:length(LFPchannels) %loops through all Channels
    
filename=['100_CH',num2str(channelOrder(i)),'.continuous']; %creates the filename for each channel

[Data] = load_open_ephys(filename); % open ephys function reads individual channel voltage trace data into vector "Data"

DataDown=double(Data(1:(30000/samplerate):end)); % downsample 'Data' to 2000hz by taking every 15th sample and convert to 'double' data format (because that is needed for filtering)

RipFilt(i,:)=filter(b,a,DataDown); % filtering downsapled data, with filter settings above 

end

CtxEnv=envelope(RipFilt(Ctx,:))/std(mean(RipFilt([Ctx,Hpc],:))); %get Envelop of RipFilt Data, otherwise thresholding does not work
HpcEnv=envelope(RipFilt(Hpc,:))/std(mean(RipFilt([Ctx,Hpc],:)));    %get Envelop of RipFilt Data, otherwise thresholding does not work
HpcEx=HpcEnv-CtxEnv;

Trs=5; % Set Threshold for Rip detection manually as multiple of Std 
RippleTresh=std(HpcEnv)*Trs;
% Thresholding by looping through data and asking wheter i is bigger than
% threshold and i-1 is still smaller or equal to threshold, this way we get only upward

Trs=1; % Set Threshold for Rip detection manually as multiple of Std 
RippleTresh=std(HpcEnv)*Trs;
% Thresholding by looping through data and asking wheter i is bigger than
% threshold and i-1 is still >= threshold, this way we get only upward

% threshold crossings
for ii=2:length(HpcEnv)
        if (HpcEnv(ii)>RippleTresh&&HpcEnv(ii-1)<=RippleTresh )
        RipCrossings(ii) =1  ;
        else
        RipCrossings(ii) =0   ;
        end
end


RipTS=find(RipCrossings==1); % get timestamps of RipCrossings

figure(1)
plot(HpcEnv)
hold on
tresh(1:length(RipCrossings))=0;
tresh(RipTS)=RippleTresh;
plot(tresh)
hold on
tresh1(1:length(RipCrossings))=0;

tresh1(SleepLong)=5;

tresh1(SleepLong)=50;

plot(tresh1,'k')


% % 2. Approach:  Time Frequency Analysis (TFA) of Ripple Frequency Band, for
% % time and space reasons:
% % only for Pyramidal Layer Channel and 1 Cortical Channel as control 
% 
% HpcTFA=mean(TFA(LFPData(17,:),[150:10:250],6,2000)); %Takes mean of RipTFA band Power
% CtxTFA=mean(TFA(LFPData(5,:),[150:10:250],6,2000)); %Takes mean of RipTFA band Power
% HpcExclusiveTFA=mean(HpcTFA)-mean(CtxTFA);
% 
% TrsTs=6; % Set Threshold for Rip detection manually as multiple of Std 
% RippleTreshTFA=std(HpcTFA(1:end-500))*TrsTs; %discard last 500 data points due to filtering aretfact
% 
% for ii=2:size(HpcTFA,2)
%         if (HpcTFA(ii)>RippleTreshTFA&HpcTFA(ii-1)<=RippleTreshTFA )
%         RipCrossingsTFA(ii) =1  ;
%         else
%         RipCrossingsTFA(ii) =0   ;
%         end
% end
% 
% RipTSTFA=find(RipCrossingsTFA==1);
% 
% figure(2)
% plot([HpcTFA(1:end-500),zeros(1,499)])
% hold on
% tresh2(1:length(RipCrossings))=0;
% tresh2(RipTSTFA)=RippleTreshTFA;
% plot(tresh2)
% hold on
% tresh3(1:length(RipCrossings))=0;
% tresh3(SleepLong)=RippleTreshTFA*2;
% plot(tresh3,'k')

%% Load SleepScoring

SleepRips=ismember(RipTS,SleepLong);

SleepRipsTS=RipTS(SleepRips==1);

WakeRipsTS=RipTS(SleepRips==0);

pre=.5; %in ms
pad=0.1
pre=pre+pad;
post=0.5; %in ms
SleepRipsLFPMatrix=[];  

for i=1:length(SleepRipsTS)
    if ((SleepRipsTS(i)-((pre*2000)-1))>0&&(SleepRipsTS(i)+(post*2000))<size(LFPData,2))      %filtering out events for which the pre post cut off would violate recording bounds    
        SleepRipsLFPMatrix=cat(3,SleepRipsLFPMatrix,LFPData(:,(SleepRipsTS(i)-((pre*2000)-1)):(SleepRipsTS(i)+post*(2000))));  
    end
end 

pre=.5; %in ms
pad=0.1
pre=pre+pad;
post=0.5; %in ms

WakeRipsLFPMatrix=[];  

for i=1:length(WakeRipsTS)
    if ((WakeRipsTS(i)-((pre*2000)-1))>0&&(WakeRipsTS(i)+(post*2000))<size(LFPData,2))      %filtering out events for which the pre post cut off would violate recording bounds    
        WakeRipsLFPMatrix=cat(3,WakeRipsLFPMatrix,LFPData(:,(WakeRipsTS(i)-((pre*2000)-1)):(WakeRipsTS(i)+post*(2000))));   
    end
end

%% CSD  Y = diff(f)/h

for i= 1:size(SleepRipsLFPMatrix,3)
SleepRipsCSDMatrix1(:,:,i)=diff(diff(SleepRipsLFPMatrix(:,:,i)));
end

for i= 1:size(WakeRipsLFPMatrix,3)
    
WakeRipsCSDMatrix1(:,:,i)=diff(diff(WakeRipsLFPMatrix(:,:,i)));
end


%% Plots

PSTHlfp1=squeeze(mean(SleepRipsLFPMatrix(:,:,:),3));
figure(1);
imagesc(PSTHlfp1(:,pad*2000:end))
ax = gca;
ticks=(((pre-pad+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre-pad)+(ticks/2000)*1),num2str(-(pre-pad)+(ticks/2000)*2),num2str(-(pre-pad)+(ticks/2000)*3),num2str(-(pre-pad)+(ticks/2000)*4),num2str(-(pre-pad)+(ticks/2000)*5),num2str(-(pre-pad)+(ticks/2000)*6),num2str(-(pre-pad)+(ticks/2000)*7),num2str(-(pre-pad)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Recording Site' ) 


PSTHlfp2=squeeze(mean(WakeRipsLFPMatrix(:,:,:),3));
figure(2);
imagesc(PSTHlfp2(:,pad*2000:end))
ax = gca;
ticks=(((pre-pad+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre-pad)+(ticks/2000)*1),num2str(-(pre-pad)+(ticks/2000)*2),num2str(-(pre-pad)+(ticks/2000)*3),num2str(-(pre-pad)+(ticks/2000)*4),num2str(-(pre-pad)+(ticks/2000)*5),num2str(-(pre-pad)+(ticks/2000)*6),num2str(-(pre-pad)+(ticks/2000)*7),num2str(-(pre-pad)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Recording Site' ) 



PSTHcsd1=squeeze(mean(SleepRipsCSDMatrix1(:,:,:),3));
figure(3);
imagesc(PSTHcsd1(:,pad*2000:end))
ax = gca;
ticks=(((pre-pad+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre-pad)+(ticks/2000)*1),num2str(-(pre-pad)+(ticks/2000)*2),num2str(-(pre-pad)+(ticks/2000)*3),num2str(-(pre-pad)+(ticks/2000)*4),num2str(-(pre-pad)+(ticks/2000)*5),num2str(-(pre-pad)+(ticks/2000)*6),num2str(-(pre-pad)+(ticks/2000)*7),num2str(-(pre-pad)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Recording Site (*150 equals distance from surface)' ) 


PSTHcsd1=squeeze(mean(WakeRipsCSDMatrix1(:,:,:),3));
figure(4);
imagesc(PSTHcsd1(:,pad*2000:end))
ax = gca;
ticks=(((pre-pad+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre-pad)+(ticks/2000)*1),num2str(-(pre-pad)+(ticks/2000)*2),num2str(-(pre-pad)+(ticks/2000)*3),num2str(-(pre-pad)+(ticks/2000)*4),num2str(-(pre-pad)+(ticks/2000)*5),num2str(-(pre-pad)+(ticks/2000)*6),num2str(-(pre-pad)+(ticks/2000)*7),num2str(-(pre-pad)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Recording Site (*150 equals distance from surface)' ) 

%% scratch
for i = 1:size(EventLFPMatrix,3)
    figure(1)
    plotmat(SleepRipsLFPMatrix(:,:,i));
   % imagesc(test(:,pad*2000:end,i))
%    figure(2)
%    plot(EventRipFrqMatrix(18,:,i))
%    hold on
%    plot(EventRipFrqMatrix(8,:,i),'r')
    plotmat(EventRipFrqMatrix(:,:,i));
   % imagesc(test(:,pad*2000:end,i))
   figure(2)
   plot(EventRipFrqMatrix(18,:,i))
   hold on
   plot(EventRipFrqMatrix(8,:,i),'r')
   pause
    close all
end

tresh(1:length(HpcExclusiveCrossings))=0;
tresh(HpcExclusiveTS)=HpcExclusiveTresh;

tresh1(1:length(HpcEnv))=0;
tresh1(RipTS)=RippleTresh;

SleepRipsRawMatrix=[]
for i=1:length(SleepRipsTS)
    
        if ((SleepRipsTS(i)-((pre*2000)-1))>0&&(SleepRipsTS(i)+(post*2000))<size(LFPData,2))      %filtering out events for which the pre post cut off would violate recording bounds    
    
for ii = 1:32 %loops through all Channels
    
        filename=['100_CH',num2str(channelOrder(i)),'.continuous']; %creates the filename for each channel

        [Data] = load_open_ephys(filename); % open ephys function reads individual channel voltage trace data into vector "Data"

        SleepRipsRawMatrix(ii,:,i)=Data(SleepRipsTS(i)*15-(pre*30000)-1:(SleepRipsTS(i)*15+post*30000));
end
        end
end
 
