function [ElecType,ElecIdxLin] = LeG_calcGrayWhite(GrayImg,WhiteImg,ElecFullIdx)

%Calculates gray/white class for and electrode given the index location
%"ElecFullIdx". ElecFullIdx is a 3d matrix of index values specifying the
%volume of the electrode.
%
%Tyler Davis
%20190629

ElecIdxLin = sub2ind(size(GrayImg),ElecFullIdx(:,1),ElecFullIdx(:,2),ElecFullIdx(:,3));
GrayVals = GrayImg(ElecIdxLin);
WhiteVals = WhiteImg(ElecIdxLin);
mGrayVal = mean(GrayVals);
mWhiteVal = mean(WhiteVals);
if mGrayVal>0.1 || mWhiteVal>0.1
    [~,h,stats] = ranksum(GrayVals,WhiteVals,'alpha',0.001);
    if h
        if isfield(stats,'zval')
            if stats.zval>0
                ElecType = 'Gray';
            else
                ElecType = 'White';
            end
        else
            if mGrayVal>mWhiteVal
                ElecType = 'Gray';
            else
                ElecType = 'White';
            end
        end
    else
        switch stats.ranksum
            case 1 %if only one voxel provided (i.e. click on 2D away from electrode)
                ElecType = 'White';
            case 2
                ElecType = 'Gray';
            otherwise
                ElecType = 'Gray';
        end
    end
else
    ElecType = 'Unknown';
end