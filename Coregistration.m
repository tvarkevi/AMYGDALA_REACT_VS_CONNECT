function log_fname = Coregistration(iSubject, subject_ID, t1_dir, t2_dir, scans_dir, prefix, study_ID, working_dir, log_fname)

% Preprocessing of the fMRI data: Coregistration
%
% Input arguments:
%       iSubject = 85;                          (subject iterator)
%       subject_ID = 'be1l080213';              (subject identifier)
%       t1_dir = '_4_1';                        (anatomical folder)
%       t2_dir = '_7_1';                        (functional folder)
%       scans_dir = 'F:\Classification of Amygdala Reactivity
%           (CAR)\Data\NIFTI_BETER';            (scan data directory)
%       prefix = 'r';                           (realigned images)
%       study_ID = 'be';                        (study identifier)
%       working_dir = 'F:\Classification of Amygdala Reactivity
%           (CAR)\Analysis';                    (working directory)
%       log_fname = 'log.txt';                  (log filename)
% Subfunctions: -
%
% List of open inputs:
%       inputs{1} = Coregister: Estimate & Reslice: Reference Image - cfg_files
%       inputs{2} = Coregister: Estimate & Reslice: Source Image - cfg_files
%       inputs{3} = Coregister: Estimate & Reslice: Other Images - cfg_files


inputs = {};

% ----- Print progress to command window ----- %
fprintf(['\n' num2str(iSubject) '\tCoregistration for subject: \t' subject_ID '\n']);

% ----- Specify subfolders ----- %
anatomical_folder = [subject_ID t1_dir];
functional_folder = [subject_ID t2_dir];

% ----- Define inputs{1}: reference image - anatomical scan ----- %
this_scan = dir([scans_dir '\' subject_ID '\' anatomical_folder '\' study_ID '*.nii']);
for iScan = 1:length(this_scan)
    inputs{1} = {[scans_dir '\' subject_ID '\' anatomical_folder '\' this_scan(iScan).name ',1']};
end

% ----- Define inputs{2}: source image - mean image ----- %
this_scan = dir([scans_dir '\' subject_ID '\' functional_folder '\mean*' study_ID '*.nii']);
for iScan = 1:length(this_scan)
    inputs{2} = {[scans_dir '\' subject_ID '\' functional_folder '\' this_scan(iScan).name ',1']};
end

% ----- Define inputs{3}: other imgages - functional images ----- %
all_scans = dir([scans_dir '\' subject_ID '\' functional_folder '\' prefix study_ID '*.nii']);
for iScan = 1:length(all_scans)
    inputs{3}{iScan, 1} = [scans_dir '\' subject_ID '\' functional_folder '\' all_scans(iScan).name ',1'];
end

% ----- Run preprocessing ----- %
jobfile = {[working_dir '\Coregistration_job.m']};
try
    spm_jobman('run', jobfile, inputs{:});
    % ----- Write progress to log file ----- %
    fileID = fopen([working_dir '\' log_fname], 'a');
    fprintf(fileID, ['\t\tMatlab code executed with no errors\n']);
    fclose(fileID);
catch
    % ----- Write failure to log file ----- %
    fileID = fopen([working_dir '\' log_fname], 'a');
    fprintf(fileID, ['\t\tError at subject: ' subject_ID '\n']);
    fclose(fileID);
end

end
