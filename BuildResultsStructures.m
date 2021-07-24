ProjectCodePath=pwd;
[Repos,ProjectName,~] = fileparts(pwd);
ProjectDataPath=strcat('G:\Dropbox\CodeData\', ProjectName);
ProjectRAWPath=strcat('F:\CodeRawData\', ProjectName);
ProjectTempRAWPath=strcat('G:\CodeTempRawData\', ProjectName);
addpath(genpath(ProjectCodePath));
addpath(genpath(ProjectDataPath));
addpath(genpath(strcat(Repos,'\Neural_Ensemble_Analysis')));
addpath('C:\Users\sp3660\Documents\GitHub\YuryEnsembles')

SelectedCreLine='Vip-IRES-Cre';
SelectedStructure='VISp';
SelectedDepth=275;
SelectedSession='three_session_A';



%%
FilesToPrcess=[];
F=FilesToPrcess;
FOlderOfFiles=[];
D = strcat(ProjectDataPath,'\',SelectedCreLine,'\',SelectedStructure,'\',string(SelectedDepth),'\',SelectedSession);
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for mm = 1:numel(N)
    T = dir(fullfile(D,N{mm},'*.mat')); % improve by specifying the file extension.
    if isempty(T)
        IdxToDele=mm
        continue
    end
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for ll = 1:numel(C)
        [FileFolder,Fr,~] = fileparts(fullfile(D,N{mm},C{ll}));
        F=[F Fr];
        FOlderOfFiles=[FOlderOfFiles FileFolder];
    end
end
N(IdxToDele)=[];
FullData={};

for i=1:numel(F)
    Dataset={};
    Dataset.name=F(i);
    FullData(i).name= Dataset.name;
    FullData(i).dataset= Dataset;
    load(fullfile(D,N{i},strcat(Dataset.name,'.mat')));
    
    Spont=StimTab{5,2}:StimTab{5,3};
    Drift1=StimTab{1,2}:StimTab{1,3};
    Drift2=StimTab{4,2}:StimTab{4,3};
    Drift3=StimTab{7,2}:StimTab{7,3};

    SortedGrats=sortrows(Grats,[1 2]);
    SpikSpont=Events(:,Spont);
    SpeedSpont=Speed(:,Spont);
    SpikDrift1=Events(:,Drift1);
    SpikDrift2=Events(:,Drift2);
    SpikDrift3=Events(:,Drift3);
    TotalDriftSpikes=[SpikDrift1 SpikDrift2 SpikDrift3];
    SpeedGrats=Speed(:,[Drift1 Drift2 Drift3]);
    %%

    TempFreqs=[1 2 4 8 15];
    Degr=[0 45 90 135 180 225 270 315];
    Comb=zeros(40,17);
    a3 = combvec(Degr,TempFreqs);
    Comb(:,1:2)=a3';
    %%
    % get indexs of grting start from grating array
    for ii=TempFreqs
        for j=Degr       
            xi=find(Grats(:,1)==ii);
            xj=find(Grats(:,2)==j);
            xixj=intersect(xi, xj);
            if length(xixj)<15
                difff=15-length(xixj);
                pad=zeros(difff);          
                xixj(end+1:15)=pad;
            end
            x=Comb(:,1)==j;
            y=Comb(:,2)==ii;
            xy=[x y];
            idx=find(xy(:,1)==1 & xy(:,2)==1);   
            Comb(idx,3:end)=xixj;
        end
    end
    %% Get indexes of starting stimuli frame in full data
    Stim_Start_Indx_FullData=Comb;
    for iii=1:40
        if Comb(iii,end)==0
            Stim_Start_Indx_FullData(iii,3:end-1)=Grats(Comb(iii,3:end-1)',end-1)';
        else
            Stim_Start_Indx_FullData(iii,3:end)=Grats(Comb(iii,3:end)',end-1)';
        end
    end

    Stim_End_Indx_FullData=Comb;
    for iiii=1:40
        if Comb(iiii,end)==0
            Stim_End_Indx_FullData(iiii,3:end-1)=Grats(Comb(iiii,3:end-1)',end)';
        else
            Stim_End_Indx_FullData(iiii,3:end)=Grats(Comb(iiii,3:end)',end)';
        end
    end


    StimuliLength=Stim_End_Indx_FullData(:,3:end)-Stim_Start_Indx_FullData(:,3:end);

    %% Shift drifting to get a single continuouns activity

    lSp1=length(SpikDrift1);
    lSp2=length(SpikDrift2);
    lSp3=length(SpikDrift3);

    framesonly=Stim_Start_Indx_FullData(:,3:end);

    %get index of each grating session
    tochn1=find(framesonly>=Drift1(1) & framesonly<=Drift1(end));
    tochn2=find(framesonly>=Drift2(1) & framesonly<=Drift2(end));
    tochn3=find(framesonly>=Drift3(1) & framesonly<=Drift3(end));


    % shift timestamps for each sesssion for all stimuli starts
    framesonly(tochn1)=framesonly(tochn1)-double(Drift1(1)-1);
    framesonly(tochn2)=framesonly(tochn2)-double(Drift2(1)-1)+lSp1;
    framesonly(tochn3)=framesonly(tochn3)-double(Drift3(1)-1)+lSp1+lSp2;


    ShiftedGratingsIdx=Stim_Start_Indx_FullData;
    ShiftedGratingsIdx(:,3:end)=framesonly;

    %% Build UDF for Dariks Code
    ShiftedGratingsIdx=ShiftedGratingsIdx';
    ShiftedGratingsIdxEnd=ShiftedGratingsIdx;
    ShiftedGratingsIdxEnd(3:end,:)=ShiftedGratingsIdx(3:end,:)+StimuliLength';


    ByOrientationEnd=[sortrows(ShiftedGratingsIdxEnd',1)]';


    BuildUDF=zeros(length(TotalDriftSpikes),length(ShiftedGratingsIdx));
    BuildUDFByOrient=zeros(length(TotalDriftSpikes),length(ShiftedGratingsIdx));
    BuildUDFByOrientColl=zeros(length(TotalDriftSpikes),length(ShiftedGratingsIdx));

    [NaRow NaCol]=find(ShiftedGratingsIdx(3:end,:)==0);
    NaCorrd=[NaRow+2 NaCol];

    for iiiiii=1:size(NaCorrd,1)
        ShiftedGratingsIdx(NaCorrd(iiiiii,1),NaCorrd(iiiiii,2))=1;
    end
    ByOrientation=[sortrows(ShiftedGratingsIdx',1)]';


    %Collapse opposite orientations
    ByOrientationCollapsed=ByOrientation;
    ByOrientationCollapsedEnd=ByOrientationEnd;
    A=[0 45 90 135];
    for iiiiiii=1:4
        ByOrientationCollapsed(1,[find(ByOrientationCollapsed(1,:)==A(iiiiiii)+180)])=A(iiiiiii);
        ByOrientationCollapsedEnd(1,[find(ByOrientationCollapsedEnd(1,:)==A(iiiiiii)+180)])=A(iiiiiii);
    end
    ByOrientationCollapsed=sortrows(ByOrientationCollapsed')';
    ByOrientationCollapsedEnd=sortrows(ByOrientationCollapsedEnd')';




    for k=1:40
        for jj=3:17
        BuildUDF(ShiftedGratingsIdx(jj,k):1:ShiftedGratingsIdxEnd(jj,k),k)=1;
        BuildUDFByOrient(ByOrientation(jj,k):1:ByOrientationEnd(jj,k),k)=1;
        BuildUDFByOrientColl(ByOrientationCollapsed(jj,k):1:ByOrientationCollapsedEnd(jj,k),k)=1;
        end
    end

    for m=1:size(NaCorrd,1)
        ShiftedGratingsIdx(NaCorrd(m,1),NaCorrd(m,2))=0;
    end
%% drift

raster=double(TotalDriftSpikes>0);

SpeedToPlot=SpeedGrats;
spine_ensembles;
smooth_window_drift=sigma;
raster_norm_drift=raster_norm;

dim_corr_drift=dim_corr;
raster_LR_drift=raster_LR;
dend_order_cell_drift=dend_order_cell;
Z_drift=Z;
raster_LR_NMF_drift=raster_LR_NMF;
SI_cells_norm_drift=SI_cells_norm;
SI_cells_LR_drift=SI_cells_LR;

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FOlderOfFiles(i), strcat(FigName,'drift', '.fig')));
end


%% spont
SpeedToPlot=SpeedSpont;
raster=double(SpikSpont>0);
spine_ensembles;
smooth_window_spont=sigma;
raster_norm_spont=raster_norm;

dim_corr_spont=dim_corr;
raster_LR_spont=raster_LR;
dend_order_cell_spont=dend_order_cell;
Z_spont=Z;
raster_LR_NMF_spont=raster_LR_NMF;
SI_cells_norm_spont=SI_cells_norm;
SI_cells_LR_spont=SI_cells_LR;

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FOlderOfFiles(i), strcat(FigName,'spont', '.fig')));
end







%%
    OrderedBy=BuildUDFByOrientColl;  
    MatrixToOrder=raster_norm_drift;
    %MatrixToOrder=raster_LR_NMF;
    
   
    Test=[OrderedBy [1:size(OrderedBy,1)]'];
    zzz=sortrows(Test, [1:40]);
    neworder=zzz(:,end);
    neworder=[[1:size(OrderedBy,1)]' neworder ];
    anothertest=sortrows(neworder,2);

    SimOrdered = sortrows([anothertest MatrixToOrder'])';
    SimOrdered=SimOrdered(3:end,:);

    SImOrderedSpeed=sortrows([anothertest SpeedGrats'])';
    SImOrderedSpeed=SImOrderedSpeed(3:end,:);


 %%   
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).events=SpikSpont;
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).spikes=double(SpikSpont>0);
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).speed=SpeedSpont;
    
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).window=smooth_window_spont;
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).normalized=raster_norm_spont;
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).dimensions=dim_corr_spont;
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).rasterSVD=raster_LR_spont;
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).clusterorder=dend_order_cell_spont;
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).cluster=Z_spont;
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).NMF=raster_LR_NMF_spont;
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).Similaritynomral=SI_cells_norm_spont;
    FullData(find([FullData.name]== Dataset.name)).dataset.Spont(1).SimilaritySVD=SI_cells_LR_spont;

    

    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).events=TotalDriftSpikes;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).spikes=double(TotalDriftSpikes>0);
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).speed=SpeedGrats;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).sortedByOrient=[];
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).sortedByOrientColl=SimOrdered;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).sortedByFreq=[];
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).sortedSpeedByOrient=[];
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).sortedSpeedByOrientColl=SImOrderedSpeed;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).sortedSpeedByFreq=[];
    
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).window=smooth_window_drift;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).normalized=raster_norm_drift;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).dimensions=dim_corr_drift;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).rasterSVD=raster_LR_drift;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).clusterorder=dend_order_cell_drift;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).cluster=Z_drift;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).NMF=raster_LR_NMF_drift;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).Similaritynomral=SI_cells_norm_drift;
    FullData(find([FullData.name]== Dataset.name)).dataset.Drift(1).SimilaritySVD=SI_cells_LR_drift;

    
    FullData(find([FullData.name]== Dataset.name)).dataset.Full(1).events=Events;
    FullData(find([FullData.name]== Dataset.name)).dataset.Full(1).spikes=double(Events>0);
    FullData(find([FullData.name]== Dataset.name)).dataset.Full(1).speed=Speed;

end
close all
save(strcat(Datapath,'\FullData','.mat'),'FullData')  
clear all