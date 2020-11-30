function FirstLevelAnalysis(iSubject, subject_ID, working_dir, t1_dir, t2_dir, prefix, mask_fn, scans_dir, onsets_dir, log_fname)

% Run first level analysis in SPM12.
%
% Input arguments:
%       iSubject = 0;                           (subject iterator)
%       subject_ID = 'be1a031010';              (subject identifier)
%       working_dir = 'F:\Classification of Amygdala Reactivity
%           (CAR)\Analysis';                    (working directory)
%       t2_dir = '_6_1';                        (functional folder)
%       t1_dir = '_4_1';                        (anatomical folder)
%       prefix = 'r';                           (realigned images)
%       mask_fn = '';                           (explicit mask)
%       scans_dir = 'F:\Classification of Amygdala Reactivity
%           (CAR)\Data\NIFTI_BETER';            (scan data directory)
%       onsets_dir = 'F:\Classification of Amygdala Reactivity
%           (CAR)\Analysis\Onsets';             (onsets data directory)
%       log_fname = 'log.txt';                  (log filename)
% Subfunctions: FirstLevelAnalysis_job
%
% List of open inputs
%       fMRI model specification: Directory - cfg_files
%       fMRI model specification: Scans - cfg_files
%       fMRI model specification: Onsets - cfg_entry
%       fMRI model specification: Durations - cfg_entry
%       fMRI model specification: Onsets - cfg_entry
%       fMRI model specification: Durations - cfg_entry
%       fMRI model specification: Onsets - cfg_entry
%       fMRI model specification: Durations - cfg_entry
%       fMRI model specification: Multipe regressors - cfg_entry
%       fMRI model specification: Explicit mask - cfg_files


inputs = {};

% ----- Print progress to command window ----- %
fprintf(['\n' num2str(iSubject) '\tConducting 1st level analysis for subject: \t' subject_ID '\n']);

% ----- Specify (scans) subfolder ----- %
subject_t2_directory = [subject_ID t2_dir];
subject_t1_directory = [subject_ID t1_dir];

% ----- Specify the (1st level) output directory ----- %
inputs{1} = {[scans_dir '\' subject_ID '\' subject_t2_directory '\first_level_analysis']};
if ~exist(inputs{1}{1}, 'dir')
    mkdir(inputs{1}{1});
end

% ----- Delete pre-existing SPM(.mat) files ----- %
preexisting_file = [inputs{1}{1} '\SPM.mat'];
if exist(preexisting_file, 'file')
    delete(preexisting_file);
end

% ----- Specify the input images ----- %
all_scans = dir([scans_dir '\' subject_ID '\' subject_t2_directory '\' prefix subject_ID '*.nii']);
for iScan = 1:length(all_scans)
    inputs{2}{iScan, 1} = [all_scans(iScan).folder '\' all_scans(iScan).name ',1'];
end

% ----- Find onsets data in subject directory ----- %
onsets_file = dir([onsets_dir '\' subject_ID '*.csv']);

% ----- Pre-allocate the onsets inputs ----- %
inputs{3} = []; % onsets neutral
inputs{5} = []; % onsets positive
inputs{7} = []; % onsets negative

% ----- Extract logdata from file ----- %
this_line = 1;
filename = [onsets_dir '\' onsets_file.name];
try
    fileID = fopen(filename, 'r');
    while ~feof(fileID)
        this_line_data = fgetl(fileID);
        
        % ----- Read column headers ----- %
        if this_line == 1
            column_headers = textscan(this_line_data, '%s', 'delimiter', ',');
            column_headers = column_headers{1};
            % ----- Find column indices ----- %
            onset_column = strcmp(column_headers, 'onsets');
            % condition_column = strcmp(column_headers, 'conditions');
            % response_column = strcmp(column_headers, 'response');
            congruent_column = strcmp(column_headers, 'congruent');
            % ----- Proceed to next line ----- %
            this_line = this_line + 1;
            continue % continue to next line
        end
        
        % ----- Read data entry ----- %
        this_data_entry = textscan(this_line_data, '%s', 'delimiter', ',');
        this_data_entry = this_data_entry{1};
        
        % ----- Extract onsets ----- %
        this_onset = str2double(this_data_entry{onset_column});
        % this_condition = str2double(this_data_entry{condition_column});
        % this_condition = str2double(this_data_entry{response_column});
        this_condition = str2double(this_data_entry{congruent_column}); % Use only congruent trials (condition == response)
        if this_condition == 1 % onsets neutral
            inputs{3} = [inputs{3}; this_onset];
        end
        if this_condition == 2 % onsets positive
            inputs{5} = [inputs{5}; this_onset];
        end
        if this_condition == 3 % onsets negative
            inputs{7} = [inputs{7}; this_onset];
        end
        
        % ----- Increment line enumerator ----- %
        this_line = this_line + 1;
    end
    fclose(fileID);
    
    % ----- Re-scale onsets ----- %
    inputs{3} = inputs{3} ./ 10000;
    inputs{5} = inputs{5} ./ 10000;
    inputs{7} = inputs{7} ./ 10000;
    
    % ----- Warn if subject has too little trials in any category  ----- %
    cut_off_value = 10; % corresponds to 30%
    if length(inputs{3}) < cut_off_value
        % ----- Write progress to log file ----- %
        fileID = fopen([working_dir '\' log_fname], 'a');
        fprintf(fileID, ['\t\tWarning: less than ' num2str(cut_off_value) ' neutral trials remaining for subject: ' subject_ID '\n']);
        fclose(fileID);
    end
    if length(inputs{5}) < cut_off_value
        % ----- Write progress to log file ----- %
        fileID = fopen([working_dir '\' log_fname], 'a');
        fprintf(fileID, ['\t\tWarning: less than ' num2str(cut_off_value) ' positive trials remaining for subject: ' subject_ID '\n']);
        fclose(fileID);
    end
    if length(inputs{7}) < cut_off_value
        % ----- Write progress to log file ----- %
        fileID = fopen([working_dir '\' log_fname], 'a');
        fprintf(fileID, ['\t\tWarning: less than ' num2str(cut_off_value) ' negative trials remaining for subject: ' subject_ID '\n']);
        fclose(fileID);
    end
    
    
    % ----- Specify durations ----- %
    durations = 2;
    inputs{4} = zeros(size(inputs{3})) + durations;
    inputs{6} = zeros(size(inputs{5})) + durations;
    inputs{8} = zeros(size(inputs{7})) + durations;
    
    % ----- Specify the nuisance variables (motion parameters) ----- %
    confound_regressors_file = dir([scans_dir '\' subject_ID '\' subject_t2_directory '\' subject_ID '*confound_regressors.csv']);
    confound_regressors = importdata([confound_regressors_file.folder '\' confound_regressors_file.name]);
    
    % spike_regressors_file = dir([scans_dir '\' subject_ID '\' subject_t2_directory '\' subject_ID '*spike_regressors.csv']);
    % spike_regressors = importdata([spike_regressors_file.folder '\' spike_regressors_file.name]);
    
    R = [confound_regressors];
    % R = [confound_regressors spike_regressors];
    savename = [scans_dir '\' subject_ID '\' subject_t2_directory '\R.mat'];
    save(savename, 'R');
    
    inputs{9} = {savename};
    
    % ----- Specify the explicit (grey matter) mask ----- %
    if ~isempty(mask_fn)
        inputs{10} = {[scans_dir '\' subject_ID '\' subject_t1_directory '\' mask_fn ',1']};
    else
        inputs{10} = {''};
    end
    
    % ----- Run first level analysis ----- %
    jobfile = {[working_dir '\' 'FirstLevelAnalysis_job.m']};
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
