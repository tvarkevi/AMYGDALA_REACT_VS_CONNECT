function ROIBasedRegressionAnalysis(workdir, rs_dir, connectivty_prefix, reactivity_fn, reactivity_hemi, mask_fn, fwhm, n_perms, n_tests)

% Conduct ROI-based linear regression analysis.
%
% Input arguments:
%       workdir = 'E:\AMYGDALA_RECON\Analysis';
%       rs_dir = 'E:\AMYGDALA_RECON\Data\NIFTI_GRAND_REST';
%       connectivty_prefix = 'connectivity_b_map_lh_';
%       reactivity_fn = 'con_0001.nii';
%       reactivity_hemi = '1'
%       reactivity_hemi = '2'
%       mask_fn = 'GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii';
%       fwhm = 6;
%       n_perms = 2000;
%       n_tests = 30;


warning('off','MATLAB:nearlySingularMatrix');

% ----- Define hemisphere of ROI ----- %
reactivity_hemi = str2num(reactivity_hemi);
if (reactivity_hemi ~= 1) && (reactivity_hemi ~= 2)
    error('Error: This is not a valid input value for reactivity_hemi. Enter 1 (left hemisphere) or 2 (right hemisphere) ');
end

% ----- Define the number of iterations ----- %
R.nIts = n_perms;
R.nTests = n_tests; % FWER-correction factor

% ----- Read text files with subject ID's ----- %
fcontents = importdata([workdir '\ROI_' reactivity_fn(1:end-4) '_based_regression_analysis.txt']);
subjects_list = fcontents.textdata;
subjects_roi_vals = fcontents.data;

% ----- Loop over subjects ----- %
S = struct;

S.mask_fn = [workdir '\' mask_fn];
S.smoothing_fwhm = fwhm;

for iSubject = 1:length(subjects_list)
    this_subject = subjects_list{iSubject};
    
    S.subject_names{iSubject} = this_subject;
    
    % ----- Extract functional connectivity map for this subject ----- %
    this_subject_connect_map_ID = dir([rs_dir '\' this_subject '\' this_subject '*\' connectivty_prefix '*.nii']);
    
    H = spm_vol([this_subject_connect_map_ID.folder '\' this_subject_connect_map_ID.name]);
    D = spm_read_vols(H);
    
    spm_smooth(D, D, S.smoothing_fwhm);
    
    S.connectivity_maps{iSubject} = D;
    
    % ----- Extract task-based ROI values for this subject ----- %
    S.reactivity_vals{iSubject} = subjects_roi_vals(iSubject, reactivity_hemi);
end

% ----- Landscape-based cluster analysis: permutation T-maps ----- %
[R.nBlobs, R.maxBlobPowers, R.maxBlobPower95] = inner_permutation_analysis(R, S);

% ----- Landscape-based cluster analysis: real T-map ----- %
p = 1:length(S.subject_names);
[R.nBlobsReal, R.maxBlobPowerReal, R.info, R.sumTReal, R.L, R.sumT, R.Tcorr, R.tcrit, R.Traw] = inner_ROI_based_regression_clustering(S, p);

% ----- Find significant clusters in real data and print results ----- %
R.f_sig_blobs = find(R.sumTReal > R.maxBlobPower95);
R.nBlobsSig = length(R.f_sig_blobs);

fprintf(['\tMax T sum 95 = ' num2str(R.maxBlobPower95) '\n']);
fprintf(['\tReal data: nBlobs = ' num2str(R.nBlobsReal) ', max T sum = ' num2str(max(R.maxBlobPowerReal)) ', nBlobs sig = ' num2str(R.nBlobsSig) '.\n']);

% ----- Clean up cluster label map ----- %
R.Lsig = R.L;
f_nonsig_blobs = find(R.sumTReal <= R.maxBlobPower95);
for n = 1:length(f_nonsig_blobs)
    f = find(R.L == R.info.label(f_nonsig_blobs(n)));
    R.Lsig(f) = 0;
end

% ----- Create output directory ----- %
if reactivity_hemi == 1 % Left hemisphere
    output_dir = [workdir '\roi_based_regression_' connectivty_prefix 'vs_' reactivity_fn(1:end-4) '_lh'];
end
if reactivity_hemi == 2 % Right hemisphere
    output_dir = [workdir '\roi_based_regression_' connectivty_prefix 'vs_' reactivity_fn(1:end-4) '_rh'];
end
if exist(output_dir, 'dir')
    rmdir(output_dir, 's');
end
mkdir(output_dir)

% ----- Save results to output directory ----- %
save([output_dir '\R.mat'], 'R');

% ----- Write T-maps and L-maps to output directory ----- %
H.fname = [output_dir '\L.nii'];
spm_write_vol(H, R.L);

H.fname = [output_dir '\sumT.nii'];
spm_write_vol(H, R.sumT);

H.fname = [output_dir '\Tcorr.nii'];
spm_write_vol(H, R.Tcorr);

H.fname = [output_dir '\Lsig.nii'];
spm_write_vol(H, R.Lsig);

H.fname = [output_dir '\Traw.nii'];
spm_write_vol(H, R.Traw);

end

function [nBlobs, maxBlobPowers, maxBlobPower95] = inner_permutation_analysis(R, S)
% Conduct permutation-based landscape cluster analysis (Gladwin, Vink,
% Mars, 2016; MethodsX, 3)

nBlobs = [];
maxBlobPowers = [];
for iIt = 1:R.nIts
    p = randperm(length(S.subject_names));
    
    % ----- Conduct ROI-based regression clustering ----- %
    [nBlobs0, maxBlobPower0, ~, ~, ~, ~, ~, ~, ~] = inner_ROI_based_regression_clustering(S, p);
    
    nBlobs = [nBlobs; nBlobs0];
    maxBlobPowers = [maxBlobPowers; maxBlobPower0];
    
    % ----- Print results of permutation to command window ----- %
    fprintf(['Iteration ' num2str(iIt) ': nBlobs = ' num2str(nBlobs0) ', max T sum = ' num2str(mean(maxBlobPower0)) '\n']);
end

% ----- Calculate threshold of significance ----- %
sorted = sort(maxBlobPowers, 'ascend');
a = ceil((1 - (0.05 / R.nTests)) * length(sorted)); % Bonferonni FWER-correction: p / n_tests

maxBlobPower95 = sorted(a);

end

function [nBlobs, maxBlobPower0, info, sumtv, L, sumT, T, tcrit, Traw] = inner_ROI_based_regression_clustering(S, p)
% Conduct ROI-based regression and apply landscape-based cluster
% analysis of resultant T-map (Gladwin, Vink, Mars, 2016; MethodsX, 3)

% ----- Load mask ----- %
H_Mask = spm_vol(S.mask_fn);
D_Mask = spm_read_vols(H_Mask);

vox_ind = find(D_Mask > 0);

% ----- Supress matrix is singular warning ----- %
warning('off', 'MATLAB:singularMatrix');

% ----- Loop over voxels ----- %
T = zeros(size(D_Mask));
for iVoxel = 1:length(vox_ind)
    this_voxel_index = vox_ind(iVoxel);
    
    % ----- Loop over subjects ----- %
    X = [ones(length(S.subject_names), 1) zeros(length(S.subject_names), 1)];
    y = zeros(length(S.subject_names), 1);
    for iSubject = 1:length(S.subject_names)
        
        % ----- Concatenate variables ----- %
        X(iSubject, end) = S.connectivity_maps{iSubject}(this_voxel_index);
        y(iSubject) = S.reactivity_vals{p(iSubject)};
    end
    
    % ----- Run linear regression ----- %
    try
        B = (X' * X)^(-1) * X' * y;
        
        MSE = sum((y - (X * B)).^2) / (size(X, 1) - size(X, 2)); % Mean square error
        VAR = MSE * (X' * X)^(-1); % Variance-covariance matrix
        
        T(this_voxel_index) = B(end) / sqrt(VAR(end, end));
    catch
        T(this_voxel_index) = 0;
    end
end 

Traw = T;

% ----- Set NaNs to zero ----- %
T(isnan(T)) = 0;  % Set NaNs to zero
T = abs(T);  % Convert to absolute T-values
tcrit = tinv(1 - 0.05, length(S.subject_names) - (length(B) - 1) - 1);  
    % Calculate cutoff t-value (from p-value using df = n - p -1, where p 
    % is the number of predictors not counting the intercept). The default 
    % is to set the cutoff at a t-value corresponding to a p-value of 0.05
    % (Gladwin, Vink, & Mars,2016). Note: this should be considered as a
    % form of pre-processing, not as a means of FWER-adjustment.
% tcrit = 0;
T(T <= tcrit) = 0;

% ----- Landscape-based cluster analysis ----- %
[L, info] = islander_recursive2016(T);

% ----- Extract cluster analysis specifications ----- %
nBlobs = info.N;

% ----- Loop over clusters: specifications ----- %
u = info.label;

sumtv = [];
sumT = zeros(size(T));
if ~isempty(u)
    for iu = 1:length(u)
        f = find(L == u(iu));
        
        % ----- Compute sum of T ----- %
        sumt0 = sum(T(f));
        sumtv = [sumtv; sumt0];
        sumT(f) = sumt0;
    end
    
    % ----- Extract maximum sum of T ----- %
    saveval = max(sumtv);
    if isnan(saveval)
        saveval = 0;
    end
    maxBlobPower0 = saveval;
else
    maxBlobPower0 = 0;
end

end
