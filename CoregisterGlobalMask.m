function CoregisterGlobalMask(working_dir, input_reference_img, input_global_img)

% Transform NIFTI image from global (MNI) to native space.
%
% Input arguments (examples):
%       working_dir = 'E:\AMYGDALA_RECON\Analysis';
%       input_reference_img = 'E:\AMYGDALA_RECON\Data\NIFTI_MARS_EMO\xm13101101\xm13101101_4_1\nxm13101101_4_1-0001.nii';
%       input_global_img = 'Amygdala_total_probability_map.nii';
% Subfunctions: CoregisterGlobalMask_job.m
%
% List of open inputs
%       Coregister: Reslice: Image Defining Space - cfg_files
%       Coregister: Reslice: Images to Reslice - cfg_files


% ----- Specify subject directory ----- %
fname_pts = strsplit(input_reference_img, '\');
if length(fname_pts) > 1
    subjects_scan_dir = strjoin(fname_pts(1:end-1), '\');
    cd(subjects_scan_dir);
end

% ----- Define inputs for spm jobfile ----- %
inputs = {};
inputs{1} = {[input_reference_img ',1']};
inputs{2} = {[working_dir '\' input_global_img ',1']};

% ----- Run spm job ----- %
jobfile = {[working_dir '\' 'CoregisterGlobalMask_job.m']};
spm_jobman('run', jobfile, inputs{:});

% ----- Move output coregistered mask to subject directory ----- %
movefile([working_dir '\r' input_global_img], [subjects_scan_dir '\r' input_global_img]);

% ----- Return to working directory ----- %
cd(working_dir);

end
