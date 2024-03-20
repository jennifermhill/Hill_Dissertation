%this is the version for cropping ROIs of STORM image using coordinates of
%centers from other methods (uses variable 'coords' from HistTxt2D,
%otherwise explicity declare 'coords' as [X;Y] Nx2 matrix before running)
%also crops and aligns one or two associated .dax files
% function CropROIHR_measureFI_JH_231011
clearvars -except coords LastFolder;

AverageMode=0;
ConvMatchMode=1;
ConvThresh=1000;
%for avg. and stdev of Z
SaveData=1;
%ROIWidth should be even

ROIWidth=4;
ROIHalf=ROIWidth/2;
offset=10;
scale=10;

ROIWidthScale=ROIWidth*scale;
%-0.5 to account for pixel vs matrix indexing mismatch
offsetScale=((offset-0.5)*scale-ROIWidthScale/2);

% map
if exist('LastFolder','var')
    GetFileName=sprintf('%s/*.bin',LastFolder);
else
    GetFileName='*.bin';
end 

getSTORM='Select the STORM bin file to crop';
if ~ispc; menu(getSTORM,'OK'); end
[FileNameL,PathNameL] = uigetfile(GetFileName,getSTORM);
LastFolder=PathNameL;

getdax='Select dax files to crop';
if ~ispc; menu(getdax,'OK'); end
GetFileName=sprintf('%s/*.dax',LastFolder);
[FileNameDax,PathNameDax] = uigetfile(GetFileName, getdax, MultiSelect='on');
%%
% check number of dax files selected and load each file
if isfloat(FileNameDax), error('No dax files selected');  
else
    if isstr(FileNameDax), FileNameDax={FileNameDax}; end
end

for i=1:length(FileNameDax)
    dax{i}=char(string(sprintf('%s%s',PathNameDax,string(FileNameDax(i)))));
    InfName{i}=dax{i}(1:end-4);
    inf{i}=sprintf('%s.inf',InfName{i});
    movie{i}=ReadDax(dax{i});
    info{i}=ReadInfoFile(inf{i});
    movieScale{i}=imresize(movie{i},scale);
end
%% 


% GetFileName=sprintf('%s/*.inf',LastFolder);
% GetFileName=sprintf('%s/*.inf',LastFolder);
% [FileNameInf,PathNameInf] = uigetfile(GetFileName,'Select inf file');

% load bin file (molecule list)
LeftFile =sprintf('%s%s',PathNameL,FileNameL);
if LeftFile(end-2)=='b'
    list_left = readbinfileNXcYcZc(LeftFile);
else
    list_left=readbintext(LeftFile);
end

% Set Lx and Ly as X and Y drift corrected values from molecule list
filehead = LeftFile(1:end-4);
Lx=list_left.xc;
Ly=list_left.yc;
LxN=Lx;
LyN=Ly;
Frame=list_left.frame;

%make lists of x and y coordinates to crop around (from coords)
Rx=coords(:,1);
Ry=coords(:,2);


%% initialize data matrices for indexing
if ~~dax{1}
    for i=1:length(movie)
        MovSize{i}=size(movie{i});
        NewMov{i}=[];
        NewMov_match{i}=[];
        MovFrame{i}=[];
    end
end

m=StructToMat(list_left);
data_final=[];
data_final_match=[];
frame_count=1;
frame_count_match=1;
match_list_ind = [];
avg_intensity_final_match=[];
avg_intensity_final=[];
coords_final_match=[];
coords_final=[];
%%
for i=1:numel(Rx)
    
    fprintf('Cropping %d of %d\n',i,numel(Rx))

    % Set ROI X/Y upper/lower boundaries based on the coordinate and size of
    % ROI
    XMax=Rx(i)+ROIHalf;
    XMin=Rx(i)-ROIHalf;
    YMax=Ry(i)+ROIHalf;
    YMin=Ry(i)-ROIHalf;
    % Find index of all molecules that fall within ROI
    ROIInd=find(Lx>XMin&Lx<XMax&Ly>YMin&Ly<YMax);
    % data_now is data from molecule list for all molecules that fall
    % within ROI
    data_now = m(ROIInd,:);
    
    if ConvMatchMode==1

        % Recenter molecule positions around 0+offset
        xCenter=Rx(i);
        yCenter=Ry(i);
        
        LxNow=Lx(ROIInd)-Rx(i);
        LyNow=Ly(ROIInd)-Ry(i);
%         
        data_now(:,2)=data_now(:,2)-xCenter+offset;
        data_now(:,3)=data_now(:,3)-yCenter+offset;     
        data_now(:,5)=data_now(:,5)-xCenter+offset;
        data_now(:,6)=data_now(:,6)-yCenter+offset; 
        CatNow=data_now(:,1);
        
        % Find the same ROI in dax files
        daxCenter=[Ry(i),Rx(i)];
        daxCenterScale=round(daxCenter*scale);
        
        avg_intensity_now=0;
        if ~~dax{1}
            daxLoc=[daxCenterScale+ROIWidthScale/2; daxCenterScale-ROIWidthScale/2];
            for j=1:length(dax)
                movieROI{j}=movieScale{j}(daxLoc(2,1):daxLoc(1,1),daxLoc(2,2):daxLoc(1,2),1);
                %find average intensity of conventional image in ROI
                avg_intensity_now=mean(movieROI{1},"all");
                MovFrame{j}=zeros(2560,2560);
                MovFrame{j}(offsetScale:offsetScale+ROIWidthScale,offsetScale:offsetScale+ROIWidthScale)=movieROI{j};
                MovFrameI{j}=imresize(MovFrame{j},1/scale);
            end
        end
        
        % Assign current ROI to matchlist or nonmatch list
        if AverageMode==1 && avg_intensity_now>ConvThresh
            data_now(:,15)=frame_count_match;
            coords_final_match(frame_count_match,1) = Rx(i);
            coords_final_match(frame_count_match,2) = Ry(i);
            frame_count_match=frame_count_match+1;
            data_final_match = [data_final_match; data_now];
            avg_intensity_final_match=[avg_intensity_final_match; avg_intensity_now];
            match_list_ind_now = i;
            match_list_ind = [match_list_ind; match_list_ind_now];
            for j=1:length(dax), NewMov_match{j} = cat(3,NewMov_match{j},MovFrameI{j}); end
        else
            data_now(:,15)=frame_count;
            coords_final(frame_count,1) = Rx(i);
            coords_final(frame_count,2) = Ry(i);
            avg_intensity_final=[avg_intensity_final; avg_intensity_now];
            frame_count=frame_count+1;
            data_final = [data_final; data_now];
            if ~~dax{1}
                for j=1:length(dax), NewMov{j}=cat(3,NewMov{j},MovFrameI{j}); end
            end
        end
        
        % Assign ROIInd frame to current structure number (i)
        Frame(ROIInd)=i;
        
    end

end
fprintf('Cropping done! \rMatched ROIs: %d \rNon-matched ROIs: %d \rWriting output... \r',frame_count_match-1, frame_count-1)

%% 

% Write new molecule list for match and nonmatches
if AverageMode==1
    outfile=sprintf('%s-CropROIs-average.bin',filehead);
    outfile_match=sprintf('%s-CropROIs-average-match.bin',filehead);
else
    outfile=sprintf('%s-CropROIs-combined.bin',filehead);
end


list_final = MatToStruct(data_final);
if AverageMode==1
    if ~isempty(data_final_match)
        list_final_match = MatToStruct(data_final_match);
        LxN_match = list_final_match.xc;
        LyN_match = list_final_match.yc;
        WriteMolBinNXcYcZc(list_final_match,outfile_match);
    end
end
WriteMolBinNXcYcZc(list_final,outfile);
% LxN = list_final.xc;
% LyN = list_final.yc;



%Write corresponding matched and nonmatched cropped dax files and info
%files
if ~~dax{1}
    for i=1:length(dax)
        if AverageMode==1
            MovSize{i}=size(NewMov{i});
            MovSize_match{i}=size(NewMov_match{i});
            MovSize{i}=MovSize{i}(3);
            %in case there's a single matched ROI
            try
                MovSize_match{i}=MovSize_match{i}(3);
            catch
                MovSize_match{i}=1;
            end
            NewMov{i}=int16(NewMov{i});
            NewMov_match{i}=int16(NewMov_match{i});
            NewMov{i}=abs(NewMov{i});
            NewMov_match{i}=abs(NewMov_match{i});
    
            fileheadDax{i} = dax{i}(1:end-4);
            DaxName{i}=char(string(sprintf('%s-CropROIs-nonmatch.dax',string(fileheadDax{i}(1:end-4)))));
            DaxName_match{i}=char(string(sprintf('%s-CropROIs-match.dax',string(fileheadDax{i}(1:end-4)))));
            FileNameInf{i}=sprintf('%s.inf',FileNameDax{i}(1:end-4));
            % InfName=sprintf('%s-CropROIs.inf',fileheadDax);
            % InfName_match=sprintf('%s-CropROIs-match.inf',fileheadDax);
            FileNameInfNew{i}=char(string(sprintf('%s-CropROIs-nonmatch.inf',string(FileNameInf{i}(1:end-4)))));
            FileNameInfNew_match{i}=char(string(sprintf('%s-CropROIs-match.inf',string(FileNameInf{i}(1:end-4)))));
            info{i}.number_of_frames=MovSize{i};
            info{i}.file=DaxName{i};
            info{i}.notes=[];
            info{i}.localName=FileNameInfNew{i};
            WriteDAXFiles(NewMov{i},info{i});
        
            info{i}.number_of_frames=MovSize_match{i};
            info{i}.file=DaxName_match{i};
            info{i}.notes=[];
            info{i}.localName=FileNameInfNew_match{i};
            WriteDAXFiles(NewMov_match{i},info{i});
        else
            MovSize{i}=size(NewMov{i});
            MovSize{i}=MovSize{i}(3);
            NewMov{i}=int16(NewMov{i});
            NewMov{i}=abs(NewMov{i});
    
            fileheadDax{i} = dax{i}(1:end-4);
            DaxName{i}=char(string(sprintf('%s-CropROIs.dax',string(fileheadDax{i}(1:end-4)))));
            FileNameInf{i}=sprintf('%s.inf',FileNameDax{i}(1:end-4));
            FileNameInfNew{i}=char(string(sprintf('%s-CropROIs.inf',string(FileNameInf{i}(1:end-4)))));
            info{i}.number_of_frames=MovSize{i};
            info{i}.file=DaxName{i};
            info{i}.notes=[];
            info{i}.localName=FileNameInfNew{i};
            WriteDAXFiles(NewMov{i},info{i});
        end
    end
end


% Write csv files containing intensities and coords
if AverageMode==1
    csvwrite(sprintf('%sintensity_list_matched.csv', LastFolder), avg_intensity_final_match);
    csvwrite(sprintf('%scoords_list_matched.csv', LastFolder), coords_final_match);
end
csvwrite(sprintf('%sintensity_list.csv', LastFolder), avg_intensity_final);
csvwrite(sprintf('%scoords_list.csv', LastFolder), coords_final);


