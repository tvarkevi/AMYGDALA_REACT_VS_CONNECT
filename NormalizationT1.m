function log_fname = NormalizationT1(iSubject, subject_ID, t1_dir, scans_dir, prefix, study_ID, working_dir, log_fname)

% Preprocessing of the fMRI data: Normalization
%
% Input arguments:
%       iSubject = 85;                          (subject iterator)
%       subject_ID = 'be1l080213';              (subject identifier)
%       t1_dir = '_4_1';                        (anatomical folder)
%       scans_dir = 'F:\Classification of Amygdala Reactivity
%           (CAR)\Data\NIFTI_BETER';            (scan data directory)
%       prefix = 'beta';                        (beta maps)
%       study_ID = 'be';                        (study identifier)
%       working_dir = 'F:\Classification of Amygdala Reactivity
%           (CAR)\Analysis';                    (working directory)
% Subfunctions: -
%
% List of open inputs:
%       inputs{1} = Normalise: Estimate & Write: Image to Align - cfg_files
%       inputs{2} = Normalise: Estimate & Write: Images to Write - cfg_files
%       inputs{3} = Normalise: Location of spm12\tpm\TPM.nii


inputs = {};

% ----- Print progress to command window ----- %
fprintf(['\n' num2str(iSubject) '\tNormalization of T1 scan for subject: \t' subject_ID '\n']);

% ----- Specify subfolders ----- %
anatomical_folder = [subject_ID t1_dir];

% ----- Define inputs{1}: anatomical scan ----- %
this_scan = dir([scans_dir '\' subject_ID '\' anatomical_folder '\' study_ID '*.nii']);
for iScan = 1:length(this_scan)
    inputs{1} = {[scans_dir '\' subject_ID '\' anatomical_folder '\' this_scan(iScan).name ',1']};
end

% ----- Define inputs{2}: anatomical scan ----- %
this_scan = dir([scans_dir '\' subject_ID '\' anatomical_folder '\' prefix study_ID '*.nii']);
for iScan = 1:length(this_scan)
    inputs{2} = {[scans_dir '\' subject_ID '\' anatomical_folder '\' this_scan(iScan).name ',1']};
end

% ----- Define inputs{3}: tissue probability map ----- %
inputs{3} = {[working_dir '\Programs\spm12\tpm\TPM.nii']};

% ----- Run preprocessing ----- %
jobfile = {[working_dir '\NormalizationT1_job.m']};
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