function FrameWiseDisplacement(input_filename, input_reference, output_filename)

% Calculate and store framewise displacement.
%
% Input argument(input_filename): Specify the name of the (.txt) file that
% contains the SPM12 motion parameters, e.g.: 
%       input_filename = 'F:\Classification of Amygdala Reactivity (CAR)\Data\NIFTI_MARS\xm13101101\xm13101101_5_1\rp_xm13101101_5_1-0001.txt';
% Input argument(input_reference): Specify the name of the reference image,
% For this experiment, use the first image (not the mean image), e.g.:
%       input_reference = 'F:\Classification of Amygdala Reactivity (CAR)\Data\NIFTI_MARS\xm13101101\xm13101101_5_1\xm13101101_5_1-0011.nii';
% Input argument(output_filename): Specify the name of the output (.txt)
% file wherein the displacements will be stored, e.g.:
%       output_filename = 'F:\Classification of Amygdala Reactivity (CAR)\Data\Nifti MARS\xm13101101\xm13101101_5_1\FD_Jenkinson_xm13101101_5_1-0001.txt';
% Subfunctions: y_FD_Jenkinson


% ----- Obtain framewise displacements (FD) ----- %
FD_Jenkinson = y_FD_Jenkinson(input_filename, input_reference);

% ----- Store FD's into a text file ----- %
fileID = fopen(output_filename, 'w');
for iLine = 1:length(FD_Jenkinson)
    fprintf(fileID, [num2str(FD_Jenkinson(iLine)) '\n']);
end
fclose(fileID);