% tHIS IS FOR SESSION a EXPERIMENTS ONLY
clear all

ProjectCodePath=pwd;
[Repos,ProjectName,~] = fileparts(pwd);
ProjectDataPath=strcat('G:\Dropbox\CodeData\', ProjectName);
ProjectRAWPath=strcat('F:\CodeRawData\', ProjectName);
ProjectTempRAWPath=strcat('G:\CodeTempRawData\', ProjectName);
addpath(genpath(ProjectCodePath));
addpath(genpath(ProjectDataPath));
addpath(genpath(strcat(Repos,'\Neural_Ensemble_Analysis')));
addpath('C:\Users\sp3660\Documents\GitHub\YuryEnsembles')

SelectedCreLine='Emx1-IRES-Cre';
SelectedStructure='VISp';
SelectedDepth=275;
SelectedSession='three_session_A';

[file,path,indx] =uigetfile(strcat(ProjectDataPath,'\',SelectedCreLine,'\',SelectedStructure,'\',string(SelectedDepth),'\',SelectedSession));
load(file);


%% Separate Spikes By Stimuli
Spont=StimTab{5,2}:StimTab{5,3};
Drift1=StimTab{1,2}:StimTab{1,3};
Drift2=StimTab{4,2}:StimTab{4,3};
Drift3=StimTab{7,2}:StimTab{7,3};
NatMov1=StimTab{3,2}:StimTab{3,3};
NatMov3A=StimTab{2,2}:StimTab{2,3};
NatMov3B=StimTab{6,2}:StimTab{6,3};


SortedGrats=sortrows(Grats,[1 2]);

SpikSpont=Events(:,Spont);
SpeedSpont=Speed(:,Spont);
SpikDrift1=Events(:,Drift1);
SpikDrift2=Events(:,Drift2);
SpikDrift3=Events(:,Drift3);
SpikNatMov1=Events(:,NatMov1);
SpikNatMov3A=Events(:,NatMov3A);
SpikNatMov3B=Events(:,NatMov3B);

%%

TempFreqs=[1 2 4 8 15];
Degr=[0 45 90 135 180 225 270 315];
Comb=zeros(40,17);
a3 = combvec(Degr,TempFreqs);
Comb(:,1:2)=a3';
%%
% get indexs of grting start from grating array
for i=TempFreqs
    for j=Degr       
        xi=find(Grats(:,1)==i);
        xj=find(Grats(:,2)==j);
        xixj=intersect(xi, xj);
        if length(xixj)<15
            difff=15-length(xixj);
            pad=zeros(difff);          
            xixj(end+1:15)=pad;
        end
        x=Comb(:,1)==j;
        y=Comb(:,2)==i;
        xy=[x y];
        idx=find(xy(:,1)==1 & xy(:,2)==1);   
        Comb(idx,3:end)=xixj;
    end
end
%% Get indexes of starting stimuli frame in full data
Stim_Start_Indx_FullData=Comb;
for i=1:40
    if Comb(i,end)==0
        Stim_Start_Indx_FullData(i,3:end-1)=Grats(Comb(i,3:end-1)',end-1)';
    else
        Stim_Start_Indx_FullData(i,3:end)=Grats(Comb(i,3:end)',end-1)';
    end
end

Stim_End_Indx_FullData=Comb;
for i=1:40
    if Comb(i,end)==0
        Stim_End_Indx_FullData(i,3:end-1)=Grats(Comb(i,3:end-1)',end)';
    else
        Stim_End_Indx_FullData(i,3:end)=Grats(Comb(i,3:end)',end)';
    end
end


StimuliLength=Stim_End_Indx_FullData(:,3:end)-Stim_Start_Indx_FullData(:,3:end);

%% Shift drifting to get a single continuouns activity
TotalDriftSpikes=[SpikDrift1 SpikDrift2 SpikDrift3];
SpeedGrats=Speed(:,[Drift1 Drift2 Drift3]);


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

for i=1:size(NaCorrd,1)
    ShiftedGratingsIdx(NaCorrd(i,1),NaCorrd(i,2))=1;
end
ByOrientation=[sortrows(ShiftedGratingsIdx',1)]';


%Collapse opposite orientations
ByOrientationCollapsed=ByOrientation;
ByOrientationCollapsedEnd=ByOrientationEnd;
A=[0 45 90 135];
for i=1:4
    ByOrientationCollapsed(1,[find(ByOrientationCollapsed(1,:)==A(i)+180)])=A(i);
    ByOrientationCollapsedEnd(1,[find(ByOrientationCollapsedEnd(1,:)==A(i)+180)])=A(i);
end
ByOrientationCollapsed=sortrows(ByOrientationCollapsed')';
ByOrientationCollapsedEnd=sortrows(ByOrientationCollapsedEnd')';




for i=1:40
    for j=3:17
    BuildUDF(ShiftedGratingsIdx(j,i):1:ShiftedGratingsIdxEnd(j,i),i)=1;
    BuildUDFByOrient(ByOrientation(j,i):1:ByOrientationEnd(j,i),i)=1;
    BuildUDFByOrientColl(ByOrientationCollapsed(j,i):1:ByOrientationCollapsedEnd(j,i),i)=1;
    end
end

for i=1:size(NaCorrd,1)
    ShiftedGratingsIdx(NaCorrd(i,1),NaCorrd(i,2))=0;
end


%%

spine_ensembles
%% Sorting by Stimuli
%OrderedBy=BuildUDFByOrient;    
%OrderedBy=BuildUDF ;  
OrderedBy=BuildUDFByOrientColl;  
MatrixToOrder=raster_norm;
%MatrixToOrder=raster_LR_NMF;
%%    
Test=[OrderedBy [1:size(OrderedBy,1)]'];
zzz=sortrows(Test, [1:40]);
neworder=zzz(:,end);
neworder=[[1:size(OrderedBy,1)]' neworder ];
anothertest=sortrows(neworder,2);

SimOrdered = sortrows([anothertest MatrixToOrder'])';
SimOrdered=SimOrdered(3:end,:);

SImOrderedSpeed=sortrows([anothertest SpeedGrats'])';
SImOrderedSpeed=SImOrderedSpeed(3:end,:);




f_vh_plot_raster(SimOrdered(dend_order_cell,:), frame_times);
hold on
for i=1:8
    coordins=find(sum(zzz(:,5*(i-1)+1:5*(i-1)+5),2)==1);
  
yellows=[1 1 0;1 1 0.2;1 1 0.4;1 1 0.6;1 1 0.8];
red=[1 0 0;1 0.2 0.2;1 0.4 0.4;1 0.6 0.6;1 0.8 0.8];
green=[0 1 0;0.2 1 0.2;0.4 1 0.4;0.6 1 0.6;0.8 1 0.8];
blues= [0 1 1; 0.2 1 1; 0.4 1 1; 0.6 1 1; 0.8 1 1];

colors=[ red; yellows ;green ;blues;red; yellows ;green ;blues  ];
  c = colors(5*(i-1)+1,:);
xxxx=[coordins(1) coordins(end) coordins(end) coordins(1)];
yyyy=[numel(dend_order_cell)+1 numel(dend_order_cell)+1 0 0];
v=[coordins(1)+50 numel(dend_order_cell)+1;coordins(end)-50 numel(dend_order_cell)+1 ;coordins(end)-50 0;coordins(1)+50 0];
f=[1 2 3 4];
patch('Faces',f,'Vertices',v,...
    'EdgeColor',c,'FaceColor','none','LineWidth',2);




end



%% get coordinates

B = permute(all_roi_masks,[2 3 1]);


centroids=zeros(size(B,3),2);

for i=1:size(B,3)
props=regionprops(B(:,:,i), 'Centroid');
centroids(i,:)=props.Centroid;

end
%%
data=TotalDriftSpikes;   %frames x cells
data=double(data>0)';
dataJesus=data';
[filepath,name,ext] = fileparts(file);
filename=strcat(name,'Processed',ext); %string.mat
filename2=strcat(path,filename);
coords=centroids;  %cells x coords
SelectedUDF=BuildUDF;

Neural_Ensemble_Analysis

UDF=SelectedUDF;
%% Dimensionality Reduced By Vectors
if 0

UDFWithPeaks=[BuildUDFByOrient dataJesus_analysis.Peaks.Indices ];
VectorIdx=find(UDFWithPeaks(:,end)~=0);
VectorOnlyUDF=UDFWithPeaks(VectorIdx,:);
VectorOnlyUDF=VectorOnlyUDF(:,1:end-1);
data=data(VectorIdx,:);
UDF=VectorOnlyUDF;
UDF(:,find(sum(UDF)==0))=[];







if 0
   
save(filename2,'data','coords','filename','UDF');
end