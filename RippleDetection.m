% Script For Semi Automated Ripple Detection for Recordings with linear
% Silicon Probes, Jan Klee 23.10.16
   

%% 1. lets establish the physical channel order tip to base, so that we do not get confused at later stages

clear all

ProbeBase2TipOmnetics=([20,17,2,35,21,16,3,34,22,15,4,33,23,14,5,32,24,13,6,31,25,12,7,30,26,11,8,29,9,28,27,10]); % I got these numbers in the Probe specs
OmneticsToIntan=[nan,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,nan,nan,7,6,5,4,3,2,1,0,31,30,29,28,27,26,25,24,nan]; % and I got these numbers from the Intan headstage specs

for i=1:32;   
    channelOrder(i)=OmneticsToIntan(ProbeBase2TipOmnetics(i))+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.to be able to do Ripple detection we need to figure out which channels correspond to CA1 cell layer, and str. Radiatum. 
% (Histology provides a good starting point but swelling and movement of
% the brain has created substantial differences between sessions)
% Therefore, I load small streches of LFP and MUA activity, plot them and
% try to figure out the regions of interest by physiological
% characteristics, specifically I try to identify cell waveforms in the MUA Data to
% determine the cell layer.
%  to determine broken channels that would screw up
% the rest of the analysis.
% The result of this analysis step is saved for every session as AnatomyAndChannelInfo

%% 2.1 load LFP 

%define the filter, in this case a 3D butterworth low-pass filter
%(at 750hz)
Cut=625; %cut off Frequency
Fn=2000; %Sampling Rate Raw Data
[b,a]   = butter(3,Cut/(Fn),'low');

samplerate=2000;   % desired LFP Sampling Rate (Down-sampled)

for i = 1:32 %loops through all Channels
    
filename=['100_CH',num2str(channelOrder(i)),'.continuous']; %creates the filename for each channel
[Data] = load_open_ephys(filename); % open ephys function reads individual channel voltage trace data into vector "Data"
DataDown=resample(double(Data),1,15); % downsample 'Data' to 2000hz by taking every 15th sample and convert to 'double' data format (because that is needed for filtering
LFPData(i,:)=filter(b,a,DataDown); % filtering downsapled data, with filter settings above, and saving it in as "LFPData"

end

%% 2.2 Plot 10 seconds LFP (all channels)
figure(1)         
plotmat(LFPData(:,400000:420000)); 

%% 2.3 Load and Preprocess small strech of MUA data

Cut=500; %cut off Frequency
Fn=30000; %Sampling Rate Raw Data
[d,c]   = butter(3,Cut/(Fn/2),'high');

for i = 1:32 %loops through all Channels
    
filename=['100_CH',num2str(channelOrder(i)),'.continuous']; %creates the filename for each channel
[Data] = load_open_ephys(filename); % open ephys function reads individual channel voltage trace data into vector "Data"
MUAData(i,:)=filter(d,c,double(Data(6300000:6315000))); % Takes only .5 second of data (the first second of the LFP above)

end

VRev=repmat(mean(MUAData),32,1); % By averaging and substracting the avergage from all channels I create a
MUAData1=MUAData-VRev;           % so called virual reference this helps to get rid of some of the noise

%% 2.4 Plot 1 second of MUA (all channels)
figure(2)
plotmat(MUAData1(:,:)); 

%% 2.5 Saving Session Specific Info about the Anatomical location and quality of the channels
!! carefull!! This has to be done manually session by session !!!
% I allready did this for some of the sessions, so maybe load the variable frist and see for
% your self if it makes sense

% Ctx=32; !!
% Hpc=2; !!
% Rad=8; !!
% DG=30;  %% Thise one does not have to be DG, can be CA3, it is just a channel outside of CA1
% badchannels=[6 7 14 19 ]; !!
% 
% allHpcChannels=1:30;
% goodchannels=setdiff(allHpcChannels,badchannels);
% save('AnatomyAndChannelInfo.mat','goodchannels','Ctx','Hpc','Rad','DG')
% 
load AnatomyAndChannelInfo.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4 Ripple Detection
% 1. The idea is to filter the CA1 cell layer, Ctx and the channel outside of
% CA1 in the Ripple Range and check when the Ripple band crosses a threshold in 
% the cell layer but not outside of CA1, this gets rid of movement noise which 
% is often also in the same frq band but allways shows up on all channels.
% 2. Then I check if there is also a sharp wave (bigger than threshold)
% pressent, only if these two criteria a met, the trial is included.

%% 4.1 filter in Ripple band

CutLow=200; %Low Cut off
CutHigh=250; %High Cut Off   
Fn=2000; % Sampling Rate Raw Data
[b,a]   = butter(3,[CutLow CutHigh]/(Fn/2),'bandpass');
LFPchannels=[Hpc,Ctx,DG];  

for i =1:length(LFPchannels) %loops through only the the important channels for time reasons
    
filename=['100_CH',num2str(channelOrder(LFPchannels(i))),'.continuous']; %creates the filename for each channel
[Data] = load_open_ephys(filename); % open ephys function reads individual channel voltage trace data into vector "Data"
DataDown=resample(double(Data),1,15); % downsample 'Data' to 2000hz by taking every 15th sample and convert to 'double' data format (because that is needed for filtering)
RipFilt(i,:)=filter(b,a,DataDown); % filtering downsapled data, with filter settings above 

end

%% 4.2 Get the envelopes and substract the mean of non CA1 channels from the CA1 channel to get instances where the Ripple band is HIGH only in CA1
%(similar to virtual referencing above)
HpcEnv=envelope(RipFilt(1,:));  
CtxEnv=envelope(RipFilt(2,:));
DGEnv=envelope(RipFilt(3,:));

VRefRip=(CtxEnv+DGEnv)/2;

HpcCor=HpcEnv-VRefRip;

%% 4.2 Threshold Crossings

!!
Abs=1.5; % Set Threshold for for Ripple in Ca1 in std
Tresh=std(HpcEnv)*Abs;

!!
Cor=6; %  Threshold for how many standart deviations the CA1 trace needs to be bigger than the DG trace (to exclude noise)
TreshCor=std(HpcCor)*Cor;

% Thresholding by looping sample to sample through one channel data and asking wheter sample i is bigger than
% threshold and sample i-1 is still >= threshold,this way we get only the
% upward threshold crossings (in the noise corrected HpcData), at the same
% time I ask if the threshold in the absolute HPC ripple Range is also met

for ii=2:length(HpcCor)
        if ((HpcCor(ii)>TreshCor&&HpcCor(ii-1)<=TreshCor)&&(HpcEnv(ii)>Tresh))
        RipCrossings(ii) =1  ;
        else
        RipCrossings(ii) =0   ;
        end
end

RipTS=find(RipCrossings==1); % get timestamps of RipCrossings

%% 4.3 Cutting out periods around each deteced Threshold crossing to create one big data 

pre=.5; %in ms
post=0.5; %in ms
pad=0; % taking a little bit more on both sides is sometimes a good idea for later plotting and stuff

pre=pre+pad; 
post=post+pad;
RipsLFPMatrix=[];  
included=[];
for i=1:length(RipTS)
    if ((RipTS(i)-((pre*2000)-1))>0&&(RipTS(i)+(post*2000))<size(LFPData,2))      %filtering out events for which the pre post cut off would violate recording bounds    
        RipsLFPMatrix=cat(3,RipsLFPMatrix,LFPData(goodchannels,(RipTS(i)-((pre*2000)-1)):(RipTS(i)+post*(2000))));        
    included=[included,i];
    end    
end
RipTS=RipTS(included);

%4.3.2 cut out events in original ripple filt data for later troubleshooting
RipsRipMatrix=[];  
for i=1:length(RipTS)
    if ((RipTS(i)-((pre*2000)-1))>0&&(RipTS(i)+(post*2000))<size(LFPData,2))      %filtering out events for which the pre post cut off would violate recording bounds    
        RipsRipMatrix=cat(3,RipsRipMatrix,RipFilt(:,(RipTS(i)-((pre*2000)-1)):(RipTS(i)+post*(2000))));        
    included=[included,i];
    end    
end

%% 4.4 Checking if a Sharpe Wave (bigger then threshold) had its peak within a x s window around the detected Ripple
% save that info in SWcounter
!!
wind=15; %window for sharp wave peak around ripple in ms


[Mins1,PT] =min(RipsLFPMatrix(Rad,(size(RipsLFPMatrix,2)/2)-wind:(size(RipsLFPMatrix,2)/2)+wind,:));
%STD=std(RipsLFPMatrix(Rad,:,:),0,2);

%Mins=Mins1./STD;
SWcounter(1:size(RipsLFPMatrix,3))=0;
!!
%ThresSW=2;
SW=find(Mins1<-750);
SWcounter(SW)=1;

%% 5.1 Load Sleep Scoring

% This uses my other script to get Periods of Sleep, Sleeplong that had a
% certain minimum duration. cortical EEG, AvgACC data, and calculated
% Movement Data. 

[Sleep,SleepLong,EEGslow,AvgACC,Mov]=SleepScoring(Ctx,2,0);

SleepCounter=ismember(RipTS,SleepLong);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%% 6.1 CSD  Y = diff(f)/h    

RipsCSDMatrix1=[]  ; 
for i= 1:size(RipsLFPMatrix,3)    
RipsCSDMatrix1(:,:,i)=diff(diff(RipsLFPMatrix(:,:,i)));
end

dist=diff(goodchannels);
x=repmat(dist(2:end),size(RipsCSDMatrix1,2),1)';
y=repmat(x,1,1,size(RipsCSDMatrix1,3));

RipsCSDMatrix=RipsCSDMatrix1./y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots

%% 1. Plot CA1 Ripple band (actually CA1-(Ctx+DG/2) for entire session with epochs of sleep and Detected Ripples overlayed

figure(1)
subplot(2,1,1)
plot(HpcEnv)
hold on
tresh(1:length(RipCrossings))=0;
tresh(RipTS)=Tresh;
plot(tresh)
hold on
tresh1(1:length(RipCrossings))=0;
tresh1(SleepLong)=Tresh/2;
plot(tresh1,'k')
axis tight
title('CA1 Ripple Filtered (Envelope)')
legend ('Ripple trace','Detected Ripples','SleepScore')


subplot(2,1,2)
plot(HpcCor)
hold on
tresh(1:length(RipCrossings))=0;
tresh(RipTS)=TreshCor;
plot(tresh)
hold on
tresh1(1:length(RipCrossings))=0;
tresh1(SleepLong)=TreshCor/2;
plot(tresh1,'k')
axis tight
title('CA1 Ripple Filtered- DG RippleFiltered')
legend ('Ripple Trace-DG Trace','Detected Ripples','SleepScore')
legend boxoff 

%% 2. Go through every detected Ripple, to see if the script worked allright, change settings if not
figure (2)
x=find(SWcounter==1);
for i= 1: length(x)
   
    subplot(4,2,[1])
    plot(RipsRipMatrix(1,:,x(i)))
    SleepCounter(i)
    hold on
    plot(RipsLFPMatrix(Hpc,:,x(i)),'r')
    title('CA1 Cell Layer LFP & Ripple Filtered during Ripple')
    legend ('Ripple Filtered ','LFP cell layer')
    legend boxoff 
 ax = gca;
ticks=(((pre+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre)+(ticks/2000)*1),num2str(-(pre)+(ticks/2000)*2),num2str(-(pre)+(ticks/2000)*3),num2str(-(pre)+(ticks/2000)*4),num2str(-(pre)+(ticks/2000)*5),num2str(-(pre)+(ticks/2000)*6),num2str(-(pre)+(ticks/2000)*7),num2str(-(pre)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Voltage mV' ) 

    subplot(4,2,[2])
    plot(RipsLFPMatrix(Rad,:,x(i)))
    title('CA1 Radiatum LFP during Ripple')
    legend ('Radiatum LFP')
    legend boxoff 
   ax = gca;
ticks=(((pre+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre)+(ticks/2000)*1),num2str(-(pre)+(ticks/2000)*2),num2str(-(pre)+(ticks/2000)*3),num2str(-(pre)+(ticks/2000)*4),num2str(-(pre)+(ticks/2000)*5),num2str(-(pre)+(ticks/2000)*6),num2str(-(pre)+(ticks/2000)*7),num2str(-(pre)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Voltage mV' ) 
 
    subplot(4,2,[3,5,7])
    plotmat(RipsLFPMatrix(:,:,x(i)));
    title('LFP during Ripple (all channels)')
    ax = gca;
ticks=(((pre+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre)+(ticks/2000)*1),num2str(-(pre)+(ticks/2000)*2),num2str(-(pre)+(ticks/2000)*3),num2str(-(pre)+(ticks/2000)*4),num2str(-(pre)+(ticks/2000)*5),num2str(-(pre)+(ticks/2000)*6),num2str(-(pre)+(ticks/2000)*7),num2str(-(pre)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Recording Site' ) 

    subplot(4,2,[4,6,8])
    imagesc(imgaussfilt(RipsCSDMatrix(:,:,i),2));
    title('CSD during Ripple (all channels)')
    ax = gca;
ticks=(((pre+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre)+(ticks/2000)*1),num2str(-(pre)+(ticks/2000)*2),num2str(-(pre)+(ticks/2000)*3),num2str(-(pre)+(ticks/2000)*4),num2str(-(pre)+(ticks/2000)*5),num2str(-(pre)+(ticks/2000)*6),num2str(-(pre)+(ticks/2000)*7),num2str(-(pre)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Recording Site' ) 

pause
clf
end

%% 3. Averages for Sleep and Wake Ripple LFPs
PSTHlfp1=squeeze(mean(RipsLFPMatrix(:,:,SleepCounter==1&SWcounter==1),3));
figure(3);
imagesc(PSTHlfp1)
ax = gca;
ticks=(((pre+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre)+(ticks/2000)*1),num2str(-(pre)+(ticks/2000)*2),num2str(-(pre)+(ticks/2000)*3),num2str(-(pre)+(ticks/2000)*4),num2str(-(pre)+(ticks/2000)*5),num2str(-(pre)+(ticks/2000)*6),num2str(-(pre)+(ticks/2000)*7),num2str(-(pre)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Recording Site' ) 
title('Average Ripple LFP during Sleep')
nSleep=size(RipsLFPMatrix(:,:,SleepCounter==1&SWcounter==1),3);
text((pre+post)*1750,5,['n =',num2str(nSleep)])

PSTHlfp2=squeeze(mean(RipsLFPMatrix(:,:,SleepCounter==0&SWcounter==1),3));
figure(4);
imagesc(PSTHlfp2)
ax = gca;
ticks=(((pre+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre)+(ticks/2000)*1),num2str(-(pre)+(ticks/2000)*2),num2str(-(pre)+(ticks/2000)*3),num2str(-(pre)+(ticks/2000)*4),num2str(-(pre)+(ticks/2000)*5),num2str(-(pre)+(ticks/2000)*6),num2str(-(pre)+(ticks/2000)*7),num2str(-(pre)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Recording Site' ) 
title('Average Ripple LFP during Wake')
nWake=size(RipsLFPMatrix(:,:,SleepCounter==0&SWcounter==1),3)
text((pre+post)*1750,5,['n =',num2str(nWake)])

%% 4. Averages for Sleep and Wake Ripple CSDs
figure(5)
PSTHCSD1=squeeze(mean(RipsCSDMatrix(:,:,SleepCounter==1&SWcounter==1),3));
imagesc(imgaussfilt(PSTHCSD1,2));
ax = gca;
ticks=(((pre+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre)+(ticks/2000)*1),num2str(-(pre)+(ticks/2000)*2),num2str(-(pre)+(ticks/2000)*3),num2str(-(pre)+(ticks/2000)*4),num2str(-(pre)+(ticks/2000)*5),num2str(-(pre)+(ticks/2000)*6),num2str(-(pre)+(ticks/2000)*7),num2str(-(pre)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Recording Site' ) 
title('Average Ripple CSD during Sleep')
text((pre+post)*1750,5,['n =',num2str(nSleep)])

figure(6)
PSTHCSD2=squeeze(mean(RipsCSDMatrix(:,:,SleepCounter==0&SWcounter==1),3));
imagesc(imgaussfilt(PSTHCSD2,2));
ax = gca;
ticks=(((pre+post)*2000)/10);
ax.XTick = [ticks,ticks*2,ticks*3,ticks*4,ticks*5,ticks*6,ticks*7,ticks*8,ticks*9];
ax.XTickLabel = {num2str(-(pre)+(ticks/2000)*1),num2str(-(pre)+(ticks/2000)*2),num2str(-(pre)+(ticks/2000)*3),num2str(-(pre)+(ticks/2000)*4),num2str(-(pre)+(ticks/2000)*5),num2str(-(pre)+(ticks/2000)*6),num2str(-(pre)+(ticks/2000)*7),num2str(-(pre)+(ticks/2000)*8)};
xlabel('Time (s)')
ylabel('Recording Site' ) 
title('Average Ripple CSD during Sleep')
text((pre+post)*1750,5,['n =',num2str(nWake)])

%% scratch


