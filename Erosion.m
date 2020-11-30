function Erosion(input_filename)

% Erosion of white matter or CSF segmentations.
%
% Input argument (input_filename): Specify the name of the file to be
% eroded, e.g.:
%       input_filename = 'F:\Classification of Amygdala Reactivity (CAR)\Data\NIFTI_MARS\xm13101101\xm13101101_4_1\c2xm13101101_4_1-0001.nii'
%       input_filename = 'F:\Classification of Amygdala Reactivity (CAR)\Data\NIFTI_MARS\xm13101101\xm13101101_4_1\c3xm13101101_4_1-0001.nii'
% Input argument (reslice_mask): Indicate whether the nuisance variable map
% Subfunctions: -


% ----- Obtain (meta)data of the nuisance mask itself ----- %
H_Scan = spm_vol(input_filename);
D_Scan = spm_read_vols(H_Scan);
D_Scan(D_Scan > 0.8) = 1;
D_Scan(D_Scan <= 0.8) = 0;

% ----- Create empty matrices with size identical to original data ----- %
horizontal_matrix = zeros(size(D_Scan));
sagittal_matrix = zeros(size(D_Scan));
coronal_matrix = zeros(size(D_Scan));

% ----- Extract row, column and slice subscripts of the data ----- %
[nRows, nCols, nSlices] = size(D_Scan);

% ----- Erosion in the horizontal plane ----- %
for iHorizontalSlice = 2:(nSlices-1)
    this_horizontal_slice = D_Scan(:,:,iHorizontalSlice);
    
    % ----- Obtain slices above and below current slice ----- %
    slice_below = D_Scan(:,:,(iHorizontalSlice - 1));
    slice_above = D_Scan(:,:,(iHorizontalSlice + 1));
    
    % ----- Multiply with slices above and below current slice ----- %
    this_horizontal_slice = this_horizontal_slice .* slice_below .* slice_above;
    
    % ----- Replace current slice with resampled slice ----- %
    horizontal_matrix(:,:,iHorizontalSlice) = this_horizontal_slice;
end

% ----- Erosion in the coronal plane ----- %
for iCoronalSlice = 2:(nCols-1)
    this_coronal_slice = D_Scan(:,iCoronalSlice,:);
    
    % ----- Reshape current slice ----- %
    this_coronal_slice = reshape(this_coronal_slice,[nRows,nSlices]);
    
    % ----- Obtain (and reshape) slices behind current slice ----- %
    slice_behind = D_Scan(:,(iCoronalSlice - 1),:);
    slice_behind = reshape(slice_behind,[nRows,nSlices]);
    
    % ----- Obtain (and reshape) slices in front of current slice ----- %
    slice_front = D_Scan(:,(iCoronalSlice + 1),:);
    slice_front = reshape(slice_front,[nRows,nSlices]);
    
    % ----- Multiply current slice with slice in front and behind ----- %
    this_coronal_slice = this_coronal_slice .* slice_behind .* slice_front;
    
    % ----- Replace current slice with resampled slice ----- %
    coronal_matrix(:,iCoronalSlice,:) = this_coronal_slice;
end

% ----- Erosion in the sagittal plane ----- %
for iSagittalSlice = 2:(nRows-1)
    this_sagittal_slice = D_Scan(iSagittalSlice,:,:);

    % ----- Reshape current slice ----- %
    this_sagittal_slice = reshape(this_sagittal_slice,[nCols,nSlices]);
    
    % ----- Obtain (and reshape) slices left from current slice ----- %
    slice_left = D_Scan((iSagittalSlice - 1),:,:);
    slice_left = reshape(slice_left,[nCols,nSlices]);
    
    % ----- Obtain (and reshape) slices right from current slice ----- %
    slice_right = D_Scan((iSagittalSlice + 1),:,:);
    slice_right = reshape(slice_right,[nCols,nSlices]);
    
    % ----- Multiply current slice with adjacent (l/r) slices ----- %
    this_sagittal_slice = this_sagittal_slice .* slice_left .* slice_right;
    
    % ----- Replace current slice with resampled slice ----- %
    sagittal_matrix(iSagittalSlice,:,:) = this_sagittal_slice;
end

% ----- Multiply horizontal, coronal, and sagittal matrices ----- %
eroded_matrix = coronal_matrix .* sagittal_matrix .* horizontal_matrix;

% ----- Save eroded mask to (.nii) file ----- %
filename_parts = strsplit(input_filename, '\');
H_Scan.fname = [];
for iParts = 1:length(filename_parts) - 1
    H_Scan.fname = [H_Scan.fname filename_parts{iParts} '\'];
end
H_Scan.fname = [H_Scan.fname 'e' filename_parts{end}];
% H.pinfo(1) = 1;
spm_write_vol(H_Scan, eroded_matrix);

end
