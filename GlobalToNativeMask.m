function GlobalToNativeMask(working_dir, input_native_img, input_global_img)

% Transform NIFTI image from global (MNI) to native space.
%
% Input arguments (examples):
%       working_dir = 'O:\MGGZ\WO\BETER_fMRI\Amygdala Prediction Project\Analysis';
%       input_native_img = 'F:\Classification of Amygdala Reactivity (CAR)\Data\NIFTI_MARS\xm13101101\xm13101101_4_1\cxm13101101_4_1-0001.nii';
%       input_global_img = 'Amygdala_total_probability_map.nii';      
% Subfunctions: GlobalToNative_job.m
%
% List of open inputs
%       Normalise: Estimate & Write: Image to Align - cfg_files
%       Normalise: Estimate & Write: Images to Write - cfg_files
%       Coregister: Reslice: Images to Reslice - cfg_files
%       Deformations: Image to base inverse on - cfg_files
%       Normalise: Write: Bounding box - cfg_entry
%       Normalise: Write: Voxel sizes - cfg_entry


% ----- Specify base directory ----- %
fname_pts = strsplit(input_native_img, '\');
if length(fname_pts) > 1
    subjects_scan_dir = strjoin(fname_pts(1:end-1), '\');
    cd(subjects_scan_dir);
end

% ----- Obtain bounding box and voxel size from native scan ----- %
H_Raw = spm_vol(input_native_img);
[BB, vx] = spm_get_bbox(H_Raw, 0);

% ----- Define input for spm jobfile ----- %
inputs = {};
inputs{1} = {[input_native_img ',1']};
inputs{2} = inputs{1};
inputs{3} = {[working_dir '\' input_global_img ',1']};
inputs{4} = {[input_native_img]};
inputs{5} = BB;
inputs{6} = abs(vx);
inputs{7} = inputs{1};

% ----- Run spm job ----- %
jobfile = {[working_dir '\' 'GlobalToNativeMask_job.m']};
spm_jobman('run', jobfile, inputs{:});

% ----- Move output native mask to subject directory ----- %
movefile([working_dir '\c' input_global_img], [subjects_scan_dir '\c' input_global_img])
movefile([working_dir '\ic' input_global_img], [subjects_scan_dir '\ic' input_global_img])
movefile([working_dir '\cic' input_global_img], [subjects_scan_dir '\cic' input_global_img])

% ----- Return to working directory ----- %
cd(working_dir)
