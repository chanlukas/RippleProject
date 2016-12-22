% function for plotting 2D Data matrices, mainly usefull to plot LFP or MUA
% recorded with linear silicon probes
% Jan Klee 10.10.2016

function [DataMat]=plotmat(data)

maxData=max(data);
minData=min(data)*-1;

limit=max([maxData,minData]);

for i=1:size(data,1)   
    
    if i==1
    DataMat(i,:)=data(i,:);
    else
    DataMat(i,:)=data(i,:)-(limit*i);
    end
    plot(DataMat(i,:))
    hold on
    axis tight
    
end
