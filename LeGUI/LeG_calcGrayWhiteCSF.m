function [ElecType,ElecIdxLin] = LeG_calcGrayWhiteCSF(app,ElecFullIdx)
%Calculates gray/white/csf/unknown class for an electrode given the index
%location "ElecFullIdx". ElecFullIdx is a matrix of index values specifying
%the volume of the electrode. If there is any ambiguity, gray is selected
%over white, white is selected over csf, and csf is selected over unknown.
%
%Tyler Davis
%20240523

Thresh = app.graywhite.threshold; %0.1
% DefaultType = app.graywhite.defaulttype; %'Gray'

ElecIdxLin = sub2ind(size(app.GrayImg),ElecFullIdx(:,1),ElecFullIdx(:,2),ElecFullIdx(:,3));
GrayVals = app.GrayImg(ElecIdxLin);
WhiteVals = app.WhiteImg(ElecIdxLin);
CSFVals = app.CSFImg(ElecIdxLin);
UnknownVals = 1-(GrayVals+WhiteVals+CSFVals);
mGrayVal = mean(GrayVals);
mWhiteVal = mean(WhiteVals);
mCSFVal = mean(CSFVals);
mUnknownVal = mean(UnknownVals);
if any([mGrayVal,mWhiteVal,mCSFVal,mUnknownVal]>Thresh)
    [~,~,stats] = kruskalwallis([GrayVals,WhiteVals,CSFVals,UnknownVals],1:4,'off');
    c = multcompare(stats,'alpha',0.05,'display','off');
    [~,midx] = max(stats.meanranks); %if ranks are the type will be gray
    switch midx
        case 1 %gray
            ElecType = 'Gray';
        case 2 %white
            idx = find(c(1,6)>0.05,1); %If gray is not different, choose gray. If csf is not different, still choose white.
            if isempty(idx)
                ElecType = 'White';
            else
                ElecType = 'Gray';
            end
        case 3 %csf
            idx = find(c([2,4],6)>0.05,1); %If gray|white are not different, choose gray|white (i.e. gray will be default if multiple are >0.05)
            if isempty(idx)
                ElecType = 'CSF';
            else
                switch idx
                    case 1 %gray
                        ElecType = 'Gray';
                    case 2 %white
                        ElecType = 'White';
                end
            end
        case 4 %unknown
            idx = find(c([3,5,6],6)>0.05,1); %If gray|white|csf are not different, choose gray|white|csf (i.e. gray will be default if multiple are >0.05)
            if isempty(idx)
                ElecType = 'Unknown';
            else
                switch idx
                    case 1 %gray
                        ElecType = 'Gray';
                    case 2 %white
                        ElecType = 'White';
                    case 3 %csf
                        ElecType = 'CSF';
                end
            end
    end
else
    ElecType = 'Unknown';
end %if >thresh