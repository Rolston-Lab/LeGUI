function [WC,Thresh] = LeG_autoElecs(handles)

Img = handles.DispImgSub;

ElecRad = str2double(get(handles.ElecRadEditH,'Value'));
XYZScale = handles.XYZScale;
MRInfo = handles.MRInfo;
ProjSurfRaw = handles.ProjSurfRaw;

Thresh = 1.5:-0.05:1;
WCList = cell(length(Thresh),1);
for k=1:length(Thresh)
    ImgBin = Img>Thresh(k);
    
    CC = bwconncomp(ImgBin,6);
    CCSize = cellfun(@length,CC.PixelIdxList);
    ElecVoxCnt = prod(ceil(2*ElecRad./XYZScale)); %approximate volume of an electrode in number of voxels
    ElecVoxRng = [ElecVoxCnt/20,ElecVoxCnt*2]; %range for removing false positives
    
    idx = CCSize<ElecVoxRng(1)|CCSize>ElecVoxRng(2);
    CC.PixelIdxList(idx) = [];
    CC.NumObjects = length(CC.PixelIdxList);
    
    if CC.NumObjects>0
        s = regionprops3(CC,Img.*ImgBin,'weightedcentroid','meanintensity','volume','solidity');
        s(s.Solidity<0.7,:) = [];
        
        WC = s.WeightedCentroid;
        WC(:,[1,2]) = WC(:,[2,1]);
        
        WCmm = [WC,ones(size(WC,1),1)]*MRInfo.mat'; WCmm(:,4) = [];
        
        idx = LeG_intriangulation(ProjSurfRaw.vertices,ProjSurfRaw.faces,WCmm);
        
        WC(~idx,:) = [];
        WC = round(WC);
        
        s(~idx,:) = [];
        [~,sidx] = sort([s.Volume].*[s.MeanIntensity],'descend');
        
        WC = WC(sidx,:);
        WCList(k) = {WC};
    end
end
NumObj = cellfun(@(x)size(x,1),WCList);

idx = find(abs(diff(NumObj))<3 & NumObj(1:end-1)>10,1); %change in number of detected electrodes (as threshold decreases) should be less than 3 and total number greater than 10
if isempty(idx)
    idx = find(NumObj<=100);
    [~,midx] = max(NumObj(idx));
    idx = idx(midx);
end

WC = WCList{idx};
WC(101:end,:) = []; %remove if more than 100

Thresh = Thresh(idx);


