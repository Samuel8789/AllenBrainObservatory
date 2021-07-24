


%number of cells



nCellVIP=zeros(length(FullDataVIP),1)
for nExp=1:length(FullDataVIP)
    nCellVIP(nExp)=size(FullDataVIP(nExp).dataset.Full.events,1);
end
  
nCellPV=zeros(length(FullDataPvalb),1)
for nExp=1:length(FullDataVIP)
    nCellPV(nExp)=size(FullDataPvalb(nExp).dataset.Full.events,1);
end
    
nCellSST=zeros(length(FullDataVIP),1)
for nExp=1:length(FullDataSst)
    nCellSST(nExp)=size(FullDataSst(nExp).dataset.Full.events,1);
end
    



PVDataset=FullDataPvalb(6).dataset;
VIPDataset=FullDataVIP(7).dataset;
SSTDataset=FullDataSst(4).dataset;

FullDatasetCombined(1).dataset=PVDataset
FullDatasetCombined(2).dataset=VIPDataset
FullDatasetCombined(3).dataset=SSTDataset





 
spo=[size(FullDatasetCombined(1).dataset.Spont.events,2);size(FullDatasetCombined(2).dataset.Spont.events,2);size(FullDatasetCombined(3).dataset.Spont.events,2)]
for i=1:length(spo)
    if spo(i)<max(spo);
        diffff=max(spo)-spo(i);
        FullDatasetCombined(i).dataset.Spont.events=[ FullDatasetCombined(i).dataset.Spont.events zeros(size(FullDatasetCombined(i).dataset.Spont.events,1),diffff)];
        FullDatasetCombined(i).dataset.Spont.speed=[ FullDatasetCombined(i).dataset.Spont.speed zeros(size(FullDatasetCombined(i).dataset.Spont.speed,1),diffff)];

    else
    continue
    
    end
end
spo=[size(FullDatasetCombined(1).dataset.Spont.events,2);size(FullDatasetCombined(2).dataset.Spont.events,2);size(FullDatasetCombined(3).dataset.Spont.events,2)]

    


dft=[size(FullDatasetCombined(1).dataset.Drift.events,2);size(FullDatasetCombined(2).dataset.Drift.events,2);size(FullDatasetCombined(3).dataset.Drift.events,2)]
for j=1:length(dft)
    if dft(j)<max(dft);
        diffff=max(dft)-dft(j);
        FullDatasetCombined(j).dataset.Drift.events=[ FullDatasetCombined(j).dataset.Drift.events zeros(size(FullDatasetCombined(j).dataset.Drift.events,1),diffff)];
        FullDatasetCombined(j).dataset.Drift.speed=[ FullDatasetCombined(j).dataset.Drift.speed zeros(size(FullDatasetCombined(j).dataset.Drift.speed,1),diffff)];

    else
    continue
    
    end
end
dft=[size(FullDatasetCombined(1).dataset.Drift.events,2);size(FullDatasetCombined(2).dataset.Drift.events,2);size(FullDatasetCombined(3).dataset.Drift.events,2)]


combined={}

combined.CombinedEventsSpo=[FullDatasetCombined(1).dataset.Spont.events;FullDatasetCombined(2).dataset.Spont.events;FullDatasetCombined(3).dataset.Spont.events];
combined.CombinedSpeedSpo=[FullDatasetCombined(1).dataset.Spont.speed;FullDatasetCombined(2).dataset.Spont.speed;FullDatasetCombined(3).dataset.Spont.speed];

combined.CombinedEventsDrift=[FullDatasetCombined(1).dataset.Drift.events;FullDatasetCombined(2).dataset.Drift.events;FullDatasetCombined(3).dataset.Drift.events];
combined.CombinedSpeedDrift=[FullDatasetCombined(1).dataset.Drift.speed;FullDatasetCombined(2).dataset.Drift.speed;FullDatasetCombined(3).dataset.Drift.speed];








