function log_fname = Realignment(iSubject, subject_ID, ~, t2_dir, scans_dir, prefix, study_ID, working_dir, log_fname)

% Preprocessing of the fMRI data: Realignment
%
% Input arguments:
%       iSubject = 85;                          (subject iterator)
%       subject_ID = 'be1l080213';              (subject identifier)
%       t2_dir = '_7_1';                        (functional folder)
%       scans_dir = 'F:\Classification of Amygdala Reactivity
%           (CAR)\Data\NIFTI_BETER';            (scan data directory)
%       study_ID = 'be';                        (study identifier)
%       working_dir = 'F:\Classification of Amygdala Reactivity
%           (CAR)\Analysis';                    (working directory)
%       log_fname = 'log.txt';                  (log filename)
% Subfunctions: -
%
% List of open inputs:
%       inputs{1} = Realign: Estimate: Data - cfg_repeat
%       inputs{2} = Realign: Reslice: Images - cfg_files


inputs = {};

% ----- Print progress to command window ----- %
fprintf(['\n' num2str(iSubject) '\tRealignment for subject: \t' subject_ID '\n']);

% ----- Specify subfolders ----- %
functional_folder = [subject_ID t2_dir];

% ----- Define inputs{1}: raw functional scans ----- %
all_scans = dir([scans_dir '\' subject_ID '\' functional_folder '\' study_ID '*.nii']);
for iScan = 1:length(all_scans)
    inputs{1}{iScan, 1} = [scans_dir '\' subject_ID '\' functional_folder '\' all_scans(iScan).name ',1'];
end

% ----- Define inputs{2}: (slice-time corrected) functional scans ----- %
all_scans = dir([scans_dir '\' subject_ID '\' functional_folder '\' prefix study_ID '*.nii']);
for iScan = 1:length(all_scans)
    inputs{2}{iScan, 1} = [scans_dir '\' subject_ID '\' functional_folder '\' all_scans(iScan).name ',1'];
end

% ----- Run preprocessing ----- %
jobfile = {[working_dir '\Realignment_job.m']};
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
