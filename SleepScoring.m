% Sleep Scoring

% Script for automated Sleep-Scoring using data from cortical EEG screws
% recorded with the Open-Ephys system (and dataformat)
% and Accelerometer Data
% Jan Klee 9.9.2016

%%
clear all

% lets establish the physical channel order tip to base, so that we do not get confused at later stages

ProbeBase2TipOmnetics=flip([20,35,21,34,22,33,23,32,24,31,25,30,26,29,27,28,2,17,3,16,4,15,5,14,9,10,8,11,6,13,7,12]);
OmneticsToIntan=[nan,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,nan,nan,24,25,26,27,28,29,30,31,0,1,2,3,4,5,6,7,nan];
%IntanToOmnetics=[27,28,29,30,31,32,33,34,35,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,20,21,22,23,24,25,26];

for i=1:32;   
    channelOrder32_100(i)=OmneticsToIntan(ProbeBase2TipOmnetics(i))+1;
end


%% Skull screw LFP 

%define the filter, in this case a 3D butterworth low-pass filter
%(at 500hz)

Cut=500;
Fn=30000;
[b,a]   = butter(3,Cut/Fn,'low');

samplerate=2000;   % desired LFP sampling Rate

EEGchannel=1;

for i = EEGchannel  % made it a loop in case more EEG channels have to be included

filename=['100_CH',num2str(channelOrder32_100(i)),'.continuous']; %creates the filename for each channel

[Data] = load_open_ephys(filename); % open ephys function reads individual channel voltage trace data into vector "DATA"

DataDown=double(Data(1:(30000/samplerate):end)); % downsample 'DATA' to 2000hz by taking every 15th sample and convert to 'double' data format (because that is needed for filtering)

EEGData(i,:)=filter(b,a,DataDown); % filtering downsapled data, with filter settings above, and saving it in as "EEGData"

end

%% TFA time frequency analysis of cortical channel

% define parameters for TFA

Frequencys=[1:25];  %Frequncys used in TFA, for sleep scoring we are mostly intereste in low frequencys,   
cycles=6;           % number of cycles that are used for TFA, higher numbers mean better frequency resolution, lower number better temporal resolution
 
eegTFA=TFA(EEGData,Frequencys,cycles,samplerate); % runs TFA function for EEG channel

%imagesc(eegTFA)  % can be used for a quick look at the data

EEGLowFreqMean=mean(eegTFA(2:4,:)); % takes mean of low Frequency data for sleep scoring.

% use low pass filter to get smoother curve

[b,a]=butter(3,1/(samplerate*4),'low'); % filter at 0.25 hz, this is a bit arbitrary, I just checked what worked best

EEGsleep=filter(b,a,EEGLowFreqMean);
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

%% Bonsai Postion and Movement

%BonsaiFname = uigetfile('.txt')
BonsaiFname = 'test.txt';

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

for i=1:size(Event3Timestamps)
    bigger=find(Event3Timestamps>=i);
    Mov(i)=Distances(bigger(1));
end

end

%% Sleep Scoring %%

EEGthres=median(EEGsleep)*3;
ACCthres=median(AvgACC)*3;

sleep=find(EEGsleep>EEGthres&AvgACC<ACCthres);

%% plot

time=(1/2000)/60:(1/2000)/60:(length(EEGsleep)/2000)/60;

subplot(2,1,1)
plot(time,EEGsleep,'b')
hold on
thres(1:length(EEGsleep))=0;
thres(sleep)=EEGthres;
plot(time,thres,'r')
axis tight

subplot(2,1,2)
plot(time,AvgACC,'b')
axis tight
hold on
thres(1:length(AvgACC))=0;
thres(sleep)=ACCthres;
plot(time,thres,'r')










