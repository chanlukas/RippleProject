%% Scrappy Pipeline for Place cell analysis 
% 1. Preprocessing (Vref, splitting,cat)to be fed into klusta
% 2. Extracting Timestamps and ClusterIdentity (runs on klusta Kwik files)
% 3. Extracting Space Info from Blender or Bonsai pos file
% 4. Import TTL events
% 5. Align Position and Time

%% 1.
clear all

ProbeBase2TipOmnetics=([20,17,2,35,21,16,3,34,22,15,4,33,23,14,5,32,24,13,6,31,25,12,7,30,26,11,8,29,9,28,27,10]);%50micron
OmneticsToIntan=[nan,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,nan,nan,7,6,5,4,3,2,1,0,31,30,29,28,27,26,25,24,nan];

for i=1:32;   
    channelOrder(i)=OmneticsToIntan(ProbeBase2TipOmnetics(i))+1;
end


for ii =1:10
for i = 1:32 %loops through all Channels
    
filename=['100_CH',num2str(channelOrder(i)),'.continuous']; %creates the filename for each channel

[Data] = load_open_ephys(filename); % open ephys function reads individual channel voltage trace data into vector "Data"

if ii==1
DATA(i,:)=Data(1:round(length(Data)/10));
elseif ii==10
        
DATA(i,:)=Data(round(length(Data)/10)*(ii-1)+1:end);
else
    
DATA(i,:)=Data(round(length(Data)/10)*(ii-1)+1:round(length(Data)/10)*(ii));
end
i
end
ii
VRef1=mean(DATA);

DATAVRef=double(DATA)-repmat(VRef1,32,1);

DATAVRef=int16(DATAVRef);

if ii==1
    fid = fopen('Raw4Klusta.txt', 'w');   
    fwrite(fid,DATAVRef,'int16')
    fclose(fid)
else
    fid = fopen('Raw4Klusta.txt', 'a');   
    fwrite(fid,DATAVRef,'int16')
    fclose(fid)
end

clearvars DATA

end

movefile('Raw4Klusta.txt','Raw4Klusta.dat') 


%% 2. 
%h5disp('Raw4Klusta.kwik')
spikeTS=hdf5read('Raw4Klusta.kwik', '/channel_groups/0/spikes/time_samples');
Clusters=hdf5read('Raw4Klusta.kwik', '/channel_groups/0/spikes/clusters/main');


%% 3.

% ImportBonsaiPositionFile('pos.txt') %Bonsai

pos=ImportBonsaiPositionFile('Pos.txt'); % Blender

for i= 1:size(pos,1)
    
Pos(i,:)=strsplit(pos{i},',');
for ii=1:3
    POS(i,ii)=str2num(cell2mat(Pos(i,ii)));
end
end

%% 4.

filename=['100_CH',num2str(1),'.continuous'];
[Data, Timestamps, Info] = load_open_ephys(filename);

[EventData, EventTimestamps, EventInfo] = load_open_ephys('all_channels.events')

indexEvent3ON=find(EventData==2&EventInfo.eventId==1);  %% this gives you the indeces of events that are of type 2 and have the ID 1. This means only the indeces when input channel 2 was turned ON

StartOfYourRecording=EventTimestamps(1); %% simply finds the first timestamp in your recording, (basically the time that passed from when you pressed 'play' to when you pressed 'Record')

CorrectedTimestamps=EventTimestamps-StartOfYourRecording; %% substracts the time from when you pressed 'play' from all other timestamps, this way your timestamps start with 0 at the time when you actually pressed record

TimestampsInSamples=round(CorrectedTimestamps*30000); %% Transformes timestamps which were in seconds into samples, 30000 is the sampling rate of your recording system

Event3Timestamps=TimestampsInSamples(indexEvent3ON); %% gets the EventTimestamps with the right Indeces


if size(POS,1)==length(Event3Timestamps)

Mov=zeros(1,size(Data,1));

    for i=1:size(Data,1)
    bigger=find(Event3Timestamps>=i);
        if isempty(bigger)
        Mov(i)=POS(end,2);
        else
        Mov(i)=POS(bigger(1),2);
        end
    end
end

% Plotting

SpikesCellOne=spikeTS(Clusters==53);

x(1:size(Mov,2))=1;
plot(x,Mov)
hold on
MovCellOne=Mov(SpikesCellOne);
y(1:size(MovCellOne,2))=1;
plot(y,MovCellOne,'*r')

x=magic(5)
x=int16(x)
for i = 1:5 %loops through all Channels
 
    if i==1
    fid = fopen('Raw4Klusta2.txt', 'w');
    fwrite(fid,x(i,:),'double')
    fclose(fid)
    else
    fid = fopen('Raw4Klusta2.txt', 'a');
    fwrite(fid,x(i,:),'double')
    fclose(fid)
    end

end
    fid = fopen('Raw4Klusta3.txt', 'w');
    fwrite(fid,x,'integer*1')
    fclose(fid)



