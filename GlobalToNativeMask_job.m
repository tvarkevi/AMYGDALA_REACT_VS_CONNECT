%-----------------------------------------------------------------------
% Job saved on 14-Feb-2020 15:38:21 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = '<UNDEFINED>'; % inputs{1}: (coregistered) anatomical image
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = '<UNDEFINED>'; % inputs{2}: (coregistered) anatomical image
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = '<UNDEFINED>'; % inputs{4}: location of spm12\tpm\TPM.nii 
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                             78 76 85];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'n';
matlabbatch{2}.spm.spatial.coreg.write.ref(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.coreg.write.source = '<UNDEFINED>'; % inputs{4}: roi mask
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'c';
matlabbatch{3}.spm.util.defs.comp{1}.inv.comp{1}.def(1) = cfg_dep('Normalise: Estimate & Write: Deformation (Subj 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','def'));
matlabbatch{3}.spm.util.defs.comp{1}.inv.space = '<UNDEFINED>'; % inputs{5}: (coregistered) anatomical image
matlabbatch{3}.spm.util.defs.out{1}.savedef.ofname = 'inverse_deformation.nii';
matlabbatch{3}.spm.util.defs.out{1}.savedef.savedir.savepwd = 1;
matlabbatch{4}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Deformations: Deformation', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','def'));
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Reslice: Resliced Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = '<UNDEFINED>'; % inputs{6}: bounding box of (coregistered) anatomical image
matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = '<UNDEFINED>'; % inputs{7}: voxel size of (coregistered) anatomical image
matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{4}.spm.spatial.normalise.write.woptions.prefix = 'i';
matlabbatch{5}.spm.spatial.coreg.write.ref = '<UNDEFINED>'; % inputs{8}: (coregistered) anatomical image
matlabbatch{5}.spm.spatial.coreg.write.source(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{5}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{5}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{5}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{5}.spm.spatial.coreg.write.roptions.prefix = 'c';
