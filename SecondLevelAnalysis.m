function SecondLevelAnalysis(working_dir, study_ID, contrast_map, container_file, mask_fn)

% Univariate second level analysis of the task fMRI data.
%       working_dir = 'E:\AMYGDALA_RECON\Analysis';
%       study_ID = 'be';
%       contrast_map = 'con_0001';
%       contrast_map = 'connectivity_map_lh_';
%       container_file = 'E:\AMYGDALA_RECON\Analysis\scan_filenames.txt';
%       mask_fn = 'GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii';
% Subfunctions: ParametricSecondLevelAnalysis_job
%
% List of open inputs
%       Factorial design specification: Directory - cfg_files
%       Factorial design specification: Scans - cfg_files


inputs = {};

% ----- Print progress to command window ----- %
fprintf('\nConducting 2nd level analysis\n');

% ----- Specify the (1st level) output directory ----- %
inputs{1} = {[working_dir '\' study_ID '_second_level_analysis_' contrast_map]};
if exist(inputs{1}{1}, 'dir')
    rmdir(inputs{1}{1}, 's')
end
if ~exist(inputs{1}{1}, 'dir')
    mkdir(inputs{1}{1});
end

% ----- Specify the contrast input images ----- %
fid = fopen(container_file);

iScan = 1;
while ~feof(fid)
    inputs{2}{iScan, 1} = [fgetl(fid) ',1'];
    
    iScan = iScan + 1;
end
fclose(fid);

% ----- Specify the explicit (grey matter) mask ----- %
if ~isempty(mask_fn)
    inputs{3} = {[working_dir '\' mask_fn ',1']};
else
    inputs{3} = {''};
end

% ----- Run first level analysis ----- %
jobfile = {[working_dir '\' 'SecondLevelAnalysis_job.m']};
spm_jobman('run', jobfile, inputs{:});

end
