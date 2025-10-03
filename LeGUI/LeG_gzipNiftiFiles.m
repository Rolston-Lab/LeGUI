function LeG_gzipNiftiFiles(SaveDir,FigH)
%Finds all .nii files in SaveDir and compresses them using gzip. Second
%argument should be a handle to the main figure window in LeGUI.
%
% Tyler 20241107

NiftiFiles = dir(fullfile(SaveDir,'*.nii'));

if ~isempty(NiftiFiles)
    if strcmp(uiconfirm(FigH,'Uncompressed nifti files have been detected in the Registered folder. Compress now to save disk space?',''),'Cancel')
        return;
    end
else
    return;
end

WaitH = uiprogressdlg(FigH,'Indeterminate','off','Message','Compressing nifti files...');
for k=1:length(NiftiFiles)
    if ~isempty(WaitH)
        WaitH.Value = k/length(NiftiFiles);
    end
    src = fullfile(NiftiFiles(k).folder,NiftiFiles(k).name);
    dest = regexprep(src,'\.nii$','.nii.gz');
    gzip(src);
    if isfile(dest)
        delete(src);
    end
end
close(WaitH);
pause(0.1);
