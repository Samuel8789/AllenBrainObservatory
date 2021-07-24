pth=pwd;
[~,Project]=fileparts(pth);
TempDataDir=strcat('G:\CodeTempRawData\', Project);
PermDataDir=strcat('F:\CodeRawData\', Project);
CodeDataDir=strcat('G:\Dropbox\CodeData', Project);

d=dir(TempDataDir);
dfolders = d([d(:).isdir]);
b = extractfield(dfolders,'name');
strrr = string(b);
L = cellfun(@length,strrr);
strrr(L<3)=[];
already=[];
for el=1:length(strrr)
    already(el) = str2num(strrr{el});
end

videos = dir(TempDataDir + "\"+'*.h5');
a = extractfield(videos,'name');
str = string(a);
ExpIDS=[];
for el=1:length(str);
    ExpIDS(el) = str2num(str{el}(end-11:end-3));
end
ExpIDS=setdiff(ExpIDS,already);

%Params
for lol=1:length(ExpIDS)
    mov=h5read(fullfile(TempDataDir,strcat('ophys_experiment_', num2str(ExpIDS(lol)), '.h5')) , '/data');
    size=length(mov);
    divis=divisors(size);
    d=size./divis;
    Fact=find(d<7500);
    F=divis((Fact(1)));
    offset=0
    while F>20
        offset=offset+1
        sizecorr=size-offset
        divis=divisors(sizecorr);
        d=sizecorr./divis;
        Fact=find(d<7500);
        F=divis((Fact(1)));
    end
    div=1/F;
    mkdir(fullfile(TempDataDir,string(ExpIDS(lol))));
    Total=waitbar(0, 'Please wait...');
    size=sizecorr;
%% Saving Tiffs

    for i=0:(div^-1)-1
        Indiv=waitbar(0, 'Writing...');
        pos_Total=get(Total,'position');
        pos_Indiv=[pos_Total(1) pos_Total(2)+pos_Total(4) pos_Total(3) pos_Total(4)];
        set(Indiv,'position',pos_Indiv,'doublebuffer','on')
        init=(size*i*div)+2;
        finish=size*(i*div+div);
        if i==(div^-1)-1
            finish=finish+offset
        else
        end
        name=string(init-1)+'_'+string(finish);
        fullname=strcat(TempDataDir,'\',string(ExpIDS(lol)),'\' ,name);
        imwrite(mov(:,:,init-1), fullname+'.tif');
        count=0;
        for j=init:finish
            count=count+1;
            imwrite(mov(:,:,j), fullname+'.tif','WriteMode', 'append');
            waitbar(count/(size*div),Indiv);      
        end      
        waitbar(i/((div^-1)-1),Total)
        close(Indiv)
    end
close(Total)
clear mov
end