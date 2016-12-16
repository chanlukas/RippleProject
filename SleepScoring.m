%  Sleep Scoring

% Script for automated Sleep-Scoring using
% data from cortical EEG screws, ACC data (and Bonsai Video)
% recorded with the Open-Ephys system (and dataformat)
% and Accelerometer Data

%Input:
% EEGChannel= supposed to be Channel in Cortex
% Delay= scored sleep that is shorter than this (seconds) will not be
% included (this is to exclude very brief instances of sleep that are very
% difficult to sperate from artefacts
% plot= 1 if data should be plotted
% BonsaiFname= include if movement should be scored by video tracking and
% not only by Acc Data.

% Jan Klee 9.9.2016

function [Sleep,SleepLong,EEGslow,AvgACC,Mov]=SleepScoring(EEGChannel,delay,plot,BonsaiFname)

%%TroubleShooting
% EEGChannel=8;
% delay = 1;

%% Channel Order
ProbeBase2TipOmnetics=flip([20,17,2,35,21,16,3,34,22,15,4,33,23,14,5,32,24,13,6,31,25,12,7,30,26,11,8,29,9,28,27,10]);%50micron
OmneticsToIntan=[nan,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,nan,nan,24,25,26,27,28,29,30,31,0,1,2,3,4,5,6,7,nan];
%IntanToOmnetics=[27,28,29,30,31,32,33,34,35,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,20,21,22,23,24,25,26];

for i=1:32;   
    channelOrder(i)=OmneticsToIntan(ProbeBase2TipOmnetics(i))+1;
end

%% Skull screw EEG 

%define the filter, in this case a 3D butterworth low-pass filter
%(at 500hz)

Cut=500;
Fn=30000;
[b,a]   = butter(3,Cut/Fn,'low');

samplerate=2000;   % desired EEG sampling Rate

filename=['100_CH',num2str(channelOrder(EEGChannel)),'.continuous']; %creates the filename for each channel

[Data] = load_open_ephys(filename); % open ephys function reads individual channel voltage trace data into vector "DATA"

DataDown=double(Data(1:(30000/samplerate):end)); % downsample 'DATA' to 2000hz by taking every 15th sample and convert to 'double' data format (because that is needed for filtering)

EEGData=filter(b,a,DataDown); % filtering downsapled data, with filter settings above, and saving it in as "EEGData"

%% TFA time frequency analysis of cortical channel

% define parameters for TFA

Frequencys=[1:20];  %Frequncys used in TFA, for sleep scoring we are mostly intereste in low frequencys,   
cycles=6;           % number of cycles that are used for TFA, higher numbers mean better frequency resolution, lower number better temporal resolution
 
eegTFA=TFA(EEGData',Frequencys,cycles,samplerate); % runs TFA function for EEG channel

%imagesc(eegTFA)  % can be used for a quick look at the data

EEGLowFreqMean=mean(eegTFA(2:4,:)); % takes mean of low Frequency data for sleep scoring.

% use low pass filter to get smoother curve

[b,a]=butter(3,1/(samplerate*4),'low'); % filter at 0.25 hz, this is a bit arbitrary, I just checked what worked best

EEGslow=filter(b,a,EEGLowFreqMean);

%% ACC Data

ACCchannel=1:3;
[b,a]   = butter(3,0.0001,'low');

for i = ACCchannel  % made it a loop in case more EEG channels have to be included
    filename=['100_AUX',num2str(i),'.continuous']; %creates the filename for each channel
    [Data] = load_open_ephys(filename); % open ephys function reads individual channel voltage trace data into vector "DATA"
    ACCData(i,:)=double(Data(1:(30000/samplerate):end)); % downsample 'DATA' to 2000hz by taking every 15th sample and convert to 'double' data format (because that is needed for filtering)
    ACCDataFiltDiffSqrt(i,:)=sqrt(diff(filter(b,a,ACCData(i,:))).^2);
end

AvgACC=zeros(1,size(ACCData,2));
AvgACC(2:end)=mean(ACCDataFiltDiffSqrt);
%% CAM Timestamps 

[EventData, EventTimestamps, EventInfo] = load_open_ephys('all_channels.events'); % loads ttl input data; open ephys function

indexEvent3ON=find(EventInfo.eventType==3&EventInfo.eventId==1);  %% this gives you the indeces of events that are of type 3 and have the ID 1. This means only the indeces when input channel 3 was turned ON (which happens when a new Frame is recorded)

StartOfYourRecording=EventTimestamps(1); %% simply finds the first timestamp in your recording, (basically the time that passed from when you pressed 'play' to when you pressed 'Record')

CorrectedTimestamps=EventTimestamps-StartOfYourRecording; %% substracts the time from when you pressed 'play' from all other timestamps, this way your timestamps start with 0 at the time when you actually pressed record

TimestampsInSamples=floor(CorrectedTimestamps*30000); %% Transformes timestamps which were in seconds into samples, 30000 is the sampling rate of your recording system

Event3Timestamps=TimestampsInSamples(indexEvent3ON); %% gets the EventTimestamps with the right Indeces

%% Bonsai Postion and Movement
Mov=0;
if exist('BonsaiFname')==1
%BonsaiFname = uigetfile('.txt')
%BonsaiFname = 'RippleLED2016-10-18T13_38_42.txt';

[condition,x,y]=ImportBonsaiPositionFile(BonsaiFname);  %reads in raw X Y position data

IndexC = strfind(condition, 'True');
Index = find(not(cellfun('isempty', IndexC)));

% use low pass filter to get rid of sharp, noisy peaks in tracking data

[b,a]=butter(3,0.01,'low'); % filter design is a bit arbitrary here, I just tested out what works

XlowFilt=filter(b,a,x(Index));
YlowFilt=filter(b,a,y(Index));

Distances=sqrt(diff(XlowFilt).^2+diff(YlowFilt).^2); % gets the distance the center of the mouse moved between consecutive Frames
Distances(2:length(Distances)+1)=Distances;

    if length(Distances)==length(Event3Timestamps)

    Mov=zeros(1,size(Event3Timestamps))

        for i=1:size(EEGData,2)
        bigger=find(Event3Timestamps>=i);
        Mov(i)=Distances(bigger(1));
         end

    end
end
%% Sleep Scoring %%

EEGtresh=std(EEGslow)*.5;
ACCtresh=std(AvgACC)*.2;
 
Sleep=find(EEGslow>EEGtresh&AvgACC<ACCtresh);

%Sleep Long
SleepLongPre=[Sleep(2000*delay:end),zeros(1,2000*delay-1)];

SleepLongPost=[zeros(1,2000*delay-1),Sleep(1:end-(2000*delay-1))];

SleepLong=Sleep(find(Sleep+2000*delay-1==SleepLongPre&Sleep-(2000*delay-1)==SleepLongPost));

%% plot
if plot==1

time=(1/2000)/60:(1/2000)/60:(length(EEGslow)/2000)/60;

subplot(2,1,1)
plot(time,EEGslow,'b')
hold on
tresh(1:length(EEGslow))=0;
tresh(Sleep)=EEGtresh;
tresh1(1:length(EEGslow))=0;
tresh1(SleepLong)=EEGtresh;
plot(time,tresh,'r')
plot(time,tresh1,'g')
axis tight

subplot(2,1,2)
plot(time,AvgACC,'b')
axis tight
hold on
tresh(1:length(AvgACC))=0;
tresh(Sleep)=ACCtresh;
tresh1(1:length(AvgACC))=0;
tresh1(SleepLong)=ACCtresh;
plot(time,tresh,'r')
plot(time,tresh1,'g')
axis tight
end








