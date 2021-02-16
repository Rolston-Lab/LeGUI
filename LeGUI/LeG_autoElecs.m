function [WC,T] = LeG_autoElecs(app)

Img = app.CTImg;

% ElecRad = str2double(get(app.ElecRadEditH,'Value'));
ElecRad = 2; %total radius of standard ad-tech grid electrode (4mm dia, 2.3mm exposed)
XYZScale = app.XYZScale; %need to swap x/y if comparing with connected components
MRInfo = app.MRInfo;
ProjSurfRaw = app.ProjSurfRaw;

TMax = (app.CTRng(4)-app.CTRng(1))/(app.CTRng(2)-app.CTRng(1)); %max (CTRng(4)) normalized with respect to 1st (CTRng(1)) and 99th (CTRng(2)) percentiles
TMin = 1; %99th percentile
Thresh = linspace(TMax,TMin,21); Thresh(1) = []; %start with 1 step below max and progressively lower threshold to search for electrodes
ThreshHU = Thresh*(app.CTRng(2)-app.CTRng(1))+app.CTRng(1); %threholds in hounsfield units

ElecVolVox = ceil(4/3*pi*mean(ElecRad./XYZScale).^3); %approximate volume of an electrode (standard ecog - typically largest intracranial electrode) in number of voxels
ElecVolRng = [4,ElecVolVox]; %4 voxels as minimum since this is the minimum to form a spherical volume (1 voxel is too small and could be artifact)

StartTime = tic;
CCList = cell(length(Thresh),1);
WCList = cell(length(Thresh),1);
for k=1:length(Thresh)
    app.WaitH.Value = k/length(Thresh);
    
    ImgBin = Img>Thresh(k);
    
    CC = bwconncomp(ImgBin,6);
    CCSize = cellfun(@length,CC.PixelIdxList);
    
    idx = CCSize<ElecVolRng(1)|CCSize>ElecVolRng(2);
    CC.PixelIdxList(idx) = [];
    CC.NumObjects = length(CC.PixelIdxList);
    
    if CC.NumObjects>0
        s = regionprops3(CC,Img.*ImgBin,'weightedcentroid','meanintensity');
                
        WC = s.WeightedCentroid;
        WC(:,[1,2]) = WC(:,[2,1]);
        
        WCmm = [WC,ones(size(WC,1),1)]*MRInfo.mat'; WCmm(:,4) = [];
        
        idx = LeG_intriangulation(ProjSurfRaw.vertices,ProjSurfRaw.faces,WCmm);
        
        s(~idx,:) = [];
        WC(~idx,:) = [];
        CC.PixelIdxList(~idx) = [];
        CC.NumObjects = length(CC.PixelIdxList);
        
        [~,sidx] = sort([s.MeanIntensity],'descend');
        
        WC = round(WC(sidx,:));
        WCList(k) = {WC};
        CCList(k) = {CC};
    end
end
NumObj = cellfun(@(x)size(x,1),WCList);

cc = bwconncomp(abs(diff(NumObj))<=5 & NumObj(1:end-1)>10);%change in number of detected electrodes (as threshold decreases) should be less than 5 and total number greater than 10
skipflag = true;
if cc.NumObjects>0
    ccsize = cellfun(@length,cc.PixelIdxList);
    cc.PixelIdxList(ccsize<2) = [];
    cc.NumObjects = length(cc.PixelIdxList);
    if cc.NumObjects>0
        ccval = cellfun(@(x)mean(NumObj(x)),cc.PixelIdxList);
        [~,midx] = max(ccval); %find the largest continuous cluster that meet the criteria above
        idx = cc.PixelIdxList{midx};
        [~,midx] = max(NumObj(idx)); %find the index within the chosen cluster that has the most electrodes
        idx = idx(midx);
        skipflag = false;
    end
end
if skipflag
    idx = find(NumObj>10 & NumObj<150);
    if ~isempty(idx)
        idx = idx(round(end/2));
    else
        idx = round(length(Thresh)/2);
    end
end

figure('position',[50,50,400,400],'name',app.PatientIDStr); 
plot(Thresh,NumObj); 
hold on; 
plot(Thresh(idx),NumObj(idx),'or');
EndTime = toc(StartTime);

% CC = CCList{idx};
WC = WCList{idx};
THU = ThreshHU(idx);
T = Thresh(idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s = regionprops3(CC,Img.*ImgBin,'solidity');
% pd = pdist2(WC,WC,'euclidean','smallest',2);
% [TF,L] = isoutlier(pd(end,:));

% minpts = size(WC,2)+1;
% epsln = floor(40/mean(XYZScale)); %maximum search radius (voxels) for nearest minpts given standard ecog
% idx = dbscan(WC,epsln,minpts); %-1 is an outlier

% WC(idx==-1|s.Solidity<0.7,:) = []; %remove cluster outliers or electrodes that are not spherical in shape
% WC(idx==-1,:) = []; %remove cluster outliers
% WC(s.Solidity<0.7,:) = []; %remove cluster outliers

% b = TF&pd(end,:)<L;
% WC(b,:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% title(sprintf('%0.1f (%0.0f)',EndTime,sum(b)))
title(sprintf('%0.1f (%0.0f)',EndTime,THU))

WC(151:end,:) = []; %remove if more than 150



