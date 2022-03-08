function NiiFileOut = LeG_rotateImg2Standard(NiiFile,varargin)
%Performs 90 deg rotations to place image in standard space with dim1 (x ->
%saggital plane), dim2 (y -> coronal plane), dim3 (z -> axial plane). Finds
%matrix to transform nifti affine matrix to identity matrix with 1st dim
%mirrored ([-1 0 0; 0 1 0; 0 0 1]). Saves with "sd" prefix.
%
%Tyler Davis
%20190731

imgInfo = spm_vol(NiiFile);
img = spm_read_vols(imgInfo);
imgSize = size(img);
imgCenter = round(imgSize/2);
imgScale = sqrt(sum(imgInfo.mat(1:3,1:3).^2));
% imgCenterMM = imgCenter.*imgScale;

[path,file,ext] = fileparts(imgInfo.fname);
if isempty(regexp(file,'^sd','once'))
    imgInfo.fname = fullfile(path,['sd',file,ext]);
end

%transform can be provided as input instead of getting from nifti file
%header
if nargin>1
    imgInfo.mat = varargin{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = imgInfo.mat(1:3,1); [~,xi] = max(abs(x)); xs = sign(x(xi));
% y = imgInfo.mat(1:3,2); [~,yi] = max(abs(y)); ys = sign(y(yi));
% z = imgInfo.mat(1:3,3); [~,zi] = max(abs(z)); zs = sign(z(zi));
% 
% omat = imgInfo.mat; %original transform
% omat(1,4) = imgCenterMM(1)*sign(omat(1,4));
% omat(2,4) = imgCenterMM(2)*sign(omat(2,4));
% omat(3,4) = imgCenterMM(3)*sign(omat(3,4));
% 
% trans = [imgCenterMM(xi),imgCenterMM(yi),imgCenterMM(zi)]; trans(2:3) = -trans(2:3);
% scl = [imgScale(xi),imgScale(yi),imgScale(zi)]; scl(1) = -scl(1);
% 
% nmat = eye(4); %new (desired) transform
% nmat(1:3,1:3) = diag(scl);
% nmat(1:3,4) = trans;
% 
% rmat = omat\nmat;
% rmat(1:3,4) = zeros(3,1);
% 
% tform = affine3d(rmat);
% sameAsInput = affineOutputView(imgSize([2,1,3]),tform,'BoundsStyle','SameAsInput');

% imgR = imwarp(permute(img,[2,1,3]),tform,'OutputView',sameAsInput);
% imgR = imwarp(img,tform);
% imgR = imwarp(permute(img,[2,1,3]),tform);
% imgR = ipermute(imgR,[2,1,3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Performing 90 deg rotations to put in standard space
eyemat = eye(4); 
eyemat(1,1) = -1;
rmat = imgInfo.mat\eyemat; %imgInfo.mat*rmat = eyemat
b = spm_imatrix(rmat);
rvec = b(4:6);

imgR = img;
nrot = round(rvec/(pi/2)); %number of rotations
rvec = nrot*(pi/2); %rounded to nearest 90 deg rotation
for m=1:3 %spm order (shear, scale, rotation, translation)
    pmat = [setdiff(1:3,m),m];
    if nrot(m)~=0
        imgR = permute(imgR,pmat);
        imgR = rot90(imgR,nrot(m));
        imgR = ipermute(imgR,pmat);
    end
end

b0 = zeros(1,12);
b0(4:6) = rvec;
b0(7:9) = 1;
rmat0 = spm_matrix(b0); %this is now the rotation matrix rounded to nearest 90 deg

trans = abs((imgCenter.*imgScale)*rmat0(1:3,1:3)); trans(2:3) = -trans(2:3);
scl = abs(imgScale*rmat0(1:3,1:3)); scl(1) = -scl(1);

eyemat = eye(4);
eyemat(1:3,1:3) = diag(scl);
eyemat(1:3,4) = trans;

%%%%%%%%%%%%%%%%%%%%% flip l/r for RAS orientation %%%%%%%%%%%%%%%%%%%%%
% x = imgInfo.mat(1:3,1); [~,xi] = max(abs(x)); xs = sign(x(xi));
% y = imgInfo.mat(1:3,2); [~,yi] = max(abs(y)); ys = sign(y(yi));
% z = imgInfo.mat(1:3,3); [~,zi] = max(abs(z)); zs = sign(z(zi));
% if xs==1 %need to perform a left/right flip (check this!!)
%     p = [2,1,3];
%     imgR = permute(imgR,p);
%     imgR = fliplr(imgR);
%     imgR = ipermute(imgR,p);
%     disp('Left/right flip was performed!');
% end
% 
% if ~all([sum(x==0),sum(y==0),sum(z==0)]==2)
%     disp(num2str(imgInfo.mat));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imgInfo.mat = eyemat;
% imgInfo.mat = nmat;
imgInfo.dim = size(imgR);
spm_write_vol(imgInfo,imgR);

NiiFileOut = imgInfo.fname;