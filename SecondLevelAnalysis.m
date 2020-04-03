function SecondLevelAnalysis(working_dir, analysis_dir, scans_dir, prefix, predictor_filename, covariate_filename, explicit_mask)

% Introudction
%
% Input arguments:
%       working_dir = 'E:\AMYGDALA_REACT_VS_CONNECT\Analysis';
%                                               (working directory)
%       analysis_dir = 'E:\AMYGDALA_REACT_VS_CONNECT\Analysis\SnPM';
%                                               (output directory)
%       scans_dir = 'E:\AMYGDALA_REACT_VS_CONNECT\Data\SnPM';
%                                               (scan data directory)
%       prefix = 'nb_map_lh_';                  (normalized left beta maps)
%       predictor_filename = 'snpm_predictor_vector_lhemi.txt';
%                                               (left amygdala reactivity)
%       covariate_filename = '';                (covariate data)
%       explicit_mask = 'GRAND_Inclusive_GM_FOV_Mask_REST.nii';
%                                               (inlusive FOV-weighted grey matter mask)
% Subfunctions: -
%
% List of open inputs
%       MultiSub: Simple Regression; 1 covariate of interest: Analysis Directory - cfg_files
%       MultiSub: Simple Regression; 1 covariate of interest: Images to analyze - cfg_files
%       MultiSub: Simple Regression; 1 covariate of interest: Covariate - cfg_entry
%       MultiSub: Simple Regression; 1 covariate of interest: Vector - cfg_entry
%       MultiSub: Simple Regression; 1 covariate of interest: Name - cfg_entry


inputs = {};

% ----- Specify the (2nd level) output directory ----- %
inputs{1} = {analysis_dir};

% ----- Specify the input (beta map) images ----- %
all_scans = dir([scans_dir '\' prefix '*.nii']);
for iScan = 1:length(all_scans)
    inputs{2}{iScan, 1} = [all_scans(iScan).folder '\' all_scans(iScan).name ',1'];
end

% ----- Extract predictor-of-interest data from file ----- %
inputs{3} = importdata(predictor_filename);

% ----- Extract the covariate-of-no-interest data from file ----- %
if ~isempty(covariate_filename)
    inputs{4} = struct('c', {importdata(covariate_filename)}, 'cname', {covariate_filename(1:end-4)});
else
    inputs{4} = struct('c', {}, 'cname', {});
end

% ----- Specify the explicit (grey matter) mask ----- %
if ~isempty(explicit_mask)
    inputs{5} = {[explicit_mask ',1']};
else
    inputs{5} = {''};
end

% ----- Run first level analysis ----- %
jobfile = {[working_dir '\' 'SecondLevelAnalysis_job.m']};
spm_jobman('run', jobfile, inputs{:});

end