function high_res_normalization_T1(input_workdir, input_xlsfile, input_datadir)

% High-resolutation normalization (and reslicing) of T1 scans.
%
% Input argument(input_workdir): Specify the working directory, e.g.:
%       input_workdir = 'E:\AMYGDALA_RECON\Analysis';
% Input argument(input_xlsfile): Specify the full directory of the input
% xls file, e.g.:
%       input_xlsfile = 'E:\AMYGDALA_RECON\Analysis\Subjects_BETER.xlsx';
%       input_xlsfile = 'E:\AMYGDALA_RECON\Analysis\Subjects_MARS.xlsx';
% Input argument(input_datadir): Specify the data directory, e.g.:
%       input_datadir = 'E:\AMYGDALA_RECON\Data\NIFTI_BETER_EMO';
%       input_datadir = 'E:\AMYGDALA_RECON\Data\NIFTI_BETER_REST';
%       input_datadir = 'E:\AMYGDALA_RECON\Data\NIFTI_MARS_EMO';
%       input_datadir = 'E:\AMYGDALA_RECON\Data\NIFTI_MARS_REST';
% Subfunctions: high_res_normalization_T1_job


% ----- Read the input xls file ----- %
[~, ~, xls_data] = xlsread(input_xlsfile);
subjects_list = xls_data(2:end, contains(xls_data(1, :), 'SubjT0ID'));
subjects_include = xls_data(2:end, contains(xls_data(1, :), 'Include'));
t1_dirs = xls_data(2:end, contains(xls_data(1, :), 'T1'));

% ----- Define study ID and output filename ----- %
scans_type = input('Enter REST or EMO: ', 's');
if contains(input_xlsfile, 'BETER')
    log_fname = ['BETER_HighRes_T1_Normalization_' scans_type '_log.txt'];
end
if contains(input_xlsfile, 'MARS')
    log_fname = ['MARS_HighRes_T1_Normalization_' scans_type '_log.txt'];
end

% ----- Delete pre-existing logfile ----- %
if exist(log_fname, 'file')
    delete([input_workdir '\' log_fname]);
end

% ----- Loop over subjects ----- %
for iSubject = 1:length(subjects_list)
    this_subject = subjects_list{iSubject};
    if contains(this_subject, 'missing999')
        continue % continue to next subject
    end
    
    % ----- Print progress to command window ----- %
    fprintf(['\n' num2str(iSubject) '\tHigh-res T1 normalization for subject: \t' this_subject '\n']);
    
    % ----- Print progress to logfile ----- %
    fileID = fopen([input_workdir '\' log_fname], 'a');
    fprintf(fileID, [num2str(iSubject) '\tHigh-res T1 normalization for subject: \t' this_subject '\n']);
    fclose(fileID);
    
    % ----- Skip excluded subjects ----- %
    this_subject_include = subjects_include{iSubject};
    if this_subject_include == 0
        fprintf(['\t\t' 'Warning: This subject has an inclusion index of zero.' '\n\n']);
        
        % ----- Write progress to logfile ----- %
        fileID = fopen([input_workdir '\' log_fname], 'a');
        fprintf(fileID, ['\t\t' 'Warning: This subject has an inclusion index of zero' '\n']);
        fclose(fileID);
        
        continue % continue to next subject
    end
    
    % ----- Define subject t1 directory ----- %
    this_subject_t1_dir = [this_subject t1_dirs{iSubject}];
    if contains(this_subject_t1_dir, 'miss')
        fprintf(['\t\t' 'Error: missing anatomical data for this subject' '\n']);
        
        % ----- Write progress to logfile ----- %
        fileID = fopen([input_workdir '\' log_fname], 'a');
        fprintf(fileID, ['\t\t' 'Error: missing anatomical data for this subject' '\n']);
        fclose(fileID);
        
        continue  % continue to next subject
    end
    
    % ----- Rename pre-existing deformation field ----- %
    this_subject_df_file = dir([input_datadir '\' this_subject '\' this_subject_t1_dir '\y_' this_subject_t1_dir '-0001.nii']);
    if isempty(this_subject_df_file)
        fprintf(['\t\t' 'Error: missing deformation field for this subject' '\n']);
        
        % ----- Write progress to logfile ----- %
        fileID = fopen([input_workdir '\' log_fname], 'a');
        fprintf(fileID, ['\t\t' 'Error: missing deformation field for this subject' '\n']);
        fclose(fileID);
        
        continue % continue to next subject
    end
    
    src = [input_datadir '\' this_subject '\' this_subject_t1_dir '\' this_subject_df_file(1).name];
    dst = [input_datadir '\' this_subject '\' this_subject_t1_dir '\' this_subject_df_file(1).name(1:end-4) '_standard_res.nii'];
    
    movefile(src, dst);
    
    % ----- Pre-allocate inputs variable ----- %
    inputs = {};
    
    % ----- Define inputs{1} and inputs{2}: anatomical scan ----- %
    this_scan = dir([input_datadir '\' this_subject '\' this_subject_t1_dir '\' this_subject_t1_dir '*.nii']);
    if isempty(this_scan)
        
        % ----- Write progress to logfile ----- %
        fileID = fopen([input_workdir '\' log_fname], 'a');
        fprintf(fileID, ['\t\t' 'Error: missing data for this subject' '\n']);
        fclose(fileID);
        
        continue  % continue to next subject
    end
    inputs{1} = {[input_datadir '\' this_subject '\' this_subject_t1_dir '\' this_scan(1).name ',1']};
    
    inputs{2} = inputs{1};
    
    % ----- Define inputs{3}: TPM map ----- %
    inputs{3} = {[input_workdir '\Programs\spm12\tpm\TPM.nii']};
    
    % ----- Run preprocessing ----- %
    jobfile = {[input_workdir '\high_res_normalization_T1_job.m']};
    try
        spm_jobman('run', jobfile, inputs{:});
        
        % ----- Write progress to logfile ----- %
        fileID = fopen([input_workdir '\' log_fname], 'a');
        fprintf(fileID, ['\t\t' 'Matlab code executed with no errors' '\n']);
        fclose(fileID);
    catch
        % ----- Write failure to logfile ----- %
        fileID = fopen([input_workdir '\' log_fname], 'a');
        fprintf(fileID, ['\t\t' 'Error: code failed for this subject' '\n']);
        fclose(fileID);
    end
    
    % ----- Rename output deformation field ----- %
    this_subject_df_file = dir([input_datadir '\' this_subject '\' this_subject_t1_dir '\y_' this_subject_t1_dir '-0001.nii']);
    if isempty(this_subject_df_file)
        fprintf(['\t\t' 'Error: missing deformation field for this subject' '\n']);
        
        % ----- Write progress to logfile ----- %
        fileID = fopen([input_workdir '\' log_fname], 'a');
        fprintf(fileID, ['\t\t' 'Error: missing deformation field for this subject' '\n']);
        fclose(fileID);
        
        continue % continue to next subject
    end
    
    src = [input_datadir '\' this_subject '\' this_subject_t1_dir '\' this_subject_df_file(1).name];
    dst = [input_datadir '\' this_subject '\' this_subject_t1_dir '\' this_subject_df_file(1).name(1:end-4) '_high_res.nii'];
    
    movefile(src, dst);
end

end
