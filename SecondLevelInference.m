function SecondLevelInference(analysis_dir, analysis_type, effect_sign)

% Run inference of SnPM output data.
%
% Input argument:
%       analysis_dir = 'E:\AMYGDALA_REACT_VS_CONNECT\Analysis\SnPM'
%                                               (SnPM output directory)
%       analysis_type = 1                       (voxel-wise, uncorrected)
%       analysis_type = 2                       (voxel-wise, FDR-corrected)
%       analysis_type = 3                       (voxel-wise, FWE-corrected)
%       analysis_type = 4                       (cluster-wise, FWE-corrected)
%       effect_sign = 1                         (positive effects)
%       effect_sign = -1                        (negative effects)
% Subfunctions: -


cd(analysis_dir);

% ----- Specify the (2nd level) output directory ----- %
matlabbatch{1}.spm.tools.snpm.inference.SnPMmat = {[analysis_dir '\SnPM.mat']};

% ----- Specify inference) type and output filename ----- %
if analysis_type == 1 % Voxel-wise: p < 0.001, uncorrected
    matlabbatch{1}.spm.tools.snpm.inference.Thr.Vox.VoxSig.Pth = 0.001;
end
if analysis_type == 2 % Voxel-wise: p < 0.025, FDR-corrected
    matlabbatch{1}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FDRth = 0.025;
end
if analysis_type == 3 % Voxel-wise: p < 0.025, FWE-corrected
    matlabbatch{1}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FWEth = 0.025;
end
if analysis_type == 4 % Cluster-wise: p < 0.025, FWE-corrected
    matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = NaN;
    matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.FWEthC = 0.025;
end

% ----- Specify the sign of the analysis/inference ----- %
matlabbatch{1}.spm.tools.snpm.inference.Tsign = effect_sign;

% ----- Specify the output filename ----- %
if analysis_type == 1 % Voxel-wise: p < 0.001, uncorrected
    if effect_sign == 1 % Sign: positive
        matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPM_filtered_voxel_wise_uncorrected_0_001_positive';
    end
    if effect_sign == -1 % Sign: negative
        matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPM_filtered_voxel_wise_uncorrected_0_001_negative';
    end
end
if analysis_type == 2 % Voxel-wise: p < 0.025, FDR-corrected
    if effect_sign == 1 % Sign: positive
        matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPM_filtered_voxel_wise_FDR_corrected_0_025_positive';
    end
    if effect_sign == -1 % Sign: negative
        matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPM_filtered_voxel_wise_FDR_corrected_0_025_negative';
    end
end
if analysis_type == 3 % Voxel-wise: p < 0.025, FWE-corrected
    if effect_sign == 1 % Sign: positive
        matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPM_filtered_voxel_wise_FWE_corrected_0_025_positive';
    end
    if effect_sign == -1 % Sign: negative
        matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPM_filtered_voxel_wise_FWE_corrected_0_025_negative';
    end
end
if analysis_type == 4 % Cluster-wise: p < 0.025, FWE-corrected
    if effect_sign == 1 % Sign: positive
        matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPM_filtered_cluster_wise_0_001_FWE_0_025_positive';
    end
    if effect_sign == -1 % Sign: negative
        matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name = 'SnPM_filtered_cluster_wise_0_001_FWE_0_025_negative';
    end
end

% ----- Specify the output table type ----- %
matlabbatch{1}.spm.tools.snpm.inference.Report = 'MIPtable';

% ----- Run inference model ----- %
spm_jobman('run', matlabbatch);
