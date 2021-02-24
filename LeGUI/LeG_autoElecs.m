function [WC,T] = LeG_autoElecs(app)

Img = app.CTImg;

ElecRad = 2; %total radius of standard ad-tech grid electrode (4mm dia, 2.3mm exposed)
XYZScale = app.XYZScale; %need to swap x/y if comparing with connected components
MRInfo = app.MRInfo;
ProjSurfRaw = app.ProjSurfRaw;

TMax = (app.CTRng(4)-app.CTRng(1))/(app.CTRng(2)-app.CTRng(1)); %max (CTRng(4)) normalized with respect to 1st (CTRng(1)) and 99th (CTRng(2)) percentiles
TMin = 1; %99th percentile
Thresh = linspace(TMax,TMin,21); Thresh(1) = []; %start with 1 step below max and progressively lower threshold to search for electrodes
ThreshHU = Thresh*(app.CTRng(2)-app.CTRng(1))+app.CTRng(1); %threholds in hounsfield units

ElecVolVox = ceil(4/3*pi*mean(ElecRad./XYZScale).^3); %approximate volume of an electrode (standard ecog - typically largest intracranial electrode) in number of voxels
ElecVolRng = [6,ElecVolVox]; %6 voxels as minimum seems to work well for a wide range of electrode types

% StartTime = tic;
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
    cc.PixelIdxList(ccsize<2) = []; %need at least 3 (2 diffs) stable thresholds where number of detected electrodes does not change by more than 5
    cc.NumObjects = length(cc.PixelIdxList);
    if cc.NumObjects>0
        ccval = cellfun(@(x)mean(NumObj(x)),cc.PixelIdxList);
        [~,midx] = max(ccval); %find the largest continuous cluster that meet the criteria above
        idx = cc.PixelIdxList{midx};
%         [~,midx] = max(NumObj(idx)); %find the index within the chosen cluster that has the most electrodes
%         idx = idx(midx);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tidx = 1:length(Thresh); %threshold indices from max to stop point
%         wc = WCList;
%         ccsize = [];
%         for k=1:length(CCList)
%             cc = CCList{k};
%             if ~isempty(cc)
%                 ccsize = [ccsize;reshape(cellfun(@length,cc.PixelIdxList),[],1)];
%             end
%         end
%         trep = cell2mat(cellfun(@(x,y)repmat(y,size(x,1),1),wc,num2cell(tidx)','UniformOutput',false));
% 
%         wc = cell2mat(wc);
%         pd = pdist2(wc,wc)<1/mean(XYZScale); %find centroids that are <1mm from each other
%         ridx = sum(pd,2)<3; %each potential electrode should be detected across at least 3 thresholds
%         pd(ridx,:) = [];
%         pd(:,ridx) = [];
%         trep(ridx) = [];
%         wc(ridx,:) = [];
%         ccsize(ridx) = [];
%         
%         WC = [];
%         GrpSize = [];
%         GrpVar = [];
%         B = false(length(trep),1);
%         for k=1:length(trep)
%             P = pd(pd(k,:),:); %potential electrode
%             [~,c] = find(P);
%             p = unique(c);
%             t = trep(p); %matching threshold indices
%             Grps = findMonGrps(t(:)');
%             if ~isempty(Grps)
%                 [mval,midx] = max([Grps{:,3}]);
%                 if ~(mval<3 || any(B(p))) %electrode is detected for at least 3 continuous thresholds
%                     pidx = Grps{midx,2};
% %                     cs = ccsize(p(pidx)); %electrode size as threshold decreases
%                     WC = [WC;wc(p(pidx(round(end/2))),:)];
%                     GrpSize = [GrpSize;mval];
%                     GrpVar = [GrpVar;sum(sum(diff(wc(p(pidx),:))).^2)];
%                 end
%             end
%             B(p) = true;
%         end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        idx = idx(round(end/2)); %choose the middle index of cluster
        skipflag = false;
    end
end
if skipflag %if no stable clusters are found, do this
    idx = find(NumObj>10 & NumObj<250);
    if ~isempty(idx)
        idx = idx(round(end/2));
    else
        idx = round(length(Thresh)/2);
    end
end

% fH = figure('position',[50,50,400,400],'name',app.PatientIDStr); 
% aH = axes('parent',fH);
% plot(aH,Thresh,NumObj); 
% hold(aH,'on');
% plot(aH,Thresh(idx),NumObj(idx),'or');
% EndTime = toc(StartTime);

% CC = CCList{idx};
WC = WCList{idx};
% THU = ThreshHU(idx);
T = Thresh(idx);

%%%%%%%%%%%%%%%% Outlier removal %%%%%%%%%%%%%%%%%%%%%%%%%%%
% s = regionprops3(CC,Img.*ImgBin,'solidity');
pd = pdist2(WC,WC,'euclidean','smallest',2); %find closest electrode to each detected electrode
% [TF,L] = isoutlier(pd(end,:));

% minpts = size(WC,2)+1;
% epsln = floor(40/mean(XYZScale)); %maximum search radius (voxels) for nearest minpts given standard ecog
% idx = dbscan(WC,epsln,minpts); %-1 is an outlier

% WC(idx==-1|s.Solidity<0.7,:) = []; %remove cluster outliers or electrodes that are not spherical in shape
% WC(idx==-1,:) = []; %remove cluster outliers
% WC(s.Solidity<0.7,:) = []; %remove cluster outliers
% WC(TF&pd(end,:)<L,:) = [];

WC(pd(end,:)*mean(XYZScale)<1,:) = []; %remove detections that are closer than 1mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% title(aH,sprintf('%0.1f (%0.0f)',EndTime,THU))
% print(fH,fullfile(app.SaveDir,'AutoElec.png'),'-dpng','-r300')

WC(251:end,:) = []; %remove if more than 250 detections


% function Grps = findMonGrps(A)
% % Find duplicates
% % A = [1,2,3,8,9,9,9,10,11,12,13,13,16,17,18,18,18,18,19,20];
% 
% [uA,~,uIdx] = unique(A);
% dup = false(1,numel(A));
% for k=1:length(uA)
%     idx = (uIdx==k);
%     if sum(idx)>1
%         dup(idx) = true;
%     end
% end
% dup = find(dup);
% 
% %Find monotonically increasing groups excluding duplicates. However, group
% %length is calculated using one of the duplicates.
% mon = (diff(A)==0|diff(A)==1);
% mongrp = reshape(find(diff([false,mon,false])),2,[]);
% Grps = cell(size(mongrp,2),3);
% for k=1:size(mongrp,2)
%     grp = mongrp(1,k):mongrp(2,k);
%     ugrp = unique(A(grp));
%     grp = setdiff(grp,dup);
%     Grps(k,:) = {A(grp),grp,length(ugrp)};
% end




