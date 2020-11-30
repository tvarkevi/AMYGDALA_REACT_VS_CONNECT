# ----- Import modules ----- #
import os
import pandas as pd
import glob
import numpy as np
import nibabel as nb
import scipy.stats
import time
import warnings
from statsmodels.stats.multitest import multipletests


class Experiment:
    
    def __init__(self):
        """ 
        Define parameters and variables of the the experiment and dataset. 
    
        """
        
        # ----- User define directories ----- #
        self.work_dir = input("Enter working directory: ")
        self.data_dir = input("Enter data directory: ")
        
        # ----- Load (user define) data from excel spreadsheet ----- #
        fn_data = input("Enter full directory of excel file with subjects data: " )
        try:
            data_frame = pd.read_excel(fn_data)
            headers = data_frame.columns
            data = data_frame.values
        except FileNotFoundError as e:
            print("Error: ", e)
        
        # ----- Obtain subject variables ----- #
        try:
            # ----- Obtain list with subjects: exclude NaNs ----- #
            self.subjects_list = list(data[:, headers.get_loc("SubjT0ID")])
            self.subjects_PPN = list(data[:, headers.get_loc("PPN")])
            
            # ----- Obtain study identifier ----- #
            self.study_ID = self.subjects_list[0][0:2]
            if 'MARS' in fn_data.split('\\')[-1]:
                this_study_name = 'MARS'
            elif 'BETER' in fn_data.split('\\')[-1]:
                this_study_name = 'BETER'
            elif 'SIM' in fn_data.split('\\')[-1]:
                this_study_name = 'SIM'
            elif 'GRAND' in fn_data.split('\\')[-1]:
                this_study_name = 'GRAND'
            
            # ----- Define data sub-directories ----- #
            self.scans_dir = self.data_dir + '\\NIFTI_' + this_study_name
            self.logs_dir = self.data_dir + '\\LOGS_' + this_study_name
            
            onsets_dir = self.data_dir + '\\ONSETS_' + this_study_name
            if not os.path.exists(onsets_dir):
                os.makedirs(onsets_dir)
            self.onsets_dir = onsets_dir
            
            # ----- Obtain group data ----- #
            self.subjects_group = list(data[:, headers.get_loc("Group")])
            #   Control group = 1
            #   Experiment group = 2
            
            # ----- Obtain other subjects variables from input file ----- #
            for col in range(data.shape[1]):
                setattr(self, "subjects_" + str(headers[col]), list(data[:, col]))
        except UnboundLocalError as e:
            print("Error: ", e)


class Preprocessing(Experiment):
    
    def __init__(self):
        """
        Specify on what scans the analyses should be performed (resting-state
        or emotion task).
        
        """
        
        # ----- Inherit attributes from super class ----- #
        super().__init__()
        
        # ----- Specify the scan type ----- #
        type_of_scans = input("What type of input data do you want to use? Enter REST for resting-state data or EMO for emotion task data: ")
        if (type_of_scans != 'REST') and (type_of_scans != 'EMO'):
            raise(ValueError('This is not a valid input value.'))
        self.scans_type = type_of_scans
        
    def extract_inclusive_FOV_mask(self):
        """
        Extract the inclusive FOV mask for each subject seperately, and for all
        subjects combined.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting inclusive FOV mask ...")
        
        # ----- Define arguments  ----- #
        if self.scans_type == 'REST':
            t2_dir = self.subjects_REST
        elif self.scans_type == 'EMO':
            t2_dir = self.subjects_EMO
        else:
            raise(ValueError('This is not a valid input value.'))
        scans_dir = self.scans_dir + '_' + self.scans_type
        scans_prefix = input("Optional: Enter a prefix for the functional scans (defaults: enter nra for normalized, (coregistered), realigned, and slice-time corrected, or nr for normalized (coregistered) and realigned): ")
        
        # ----- Loop over subjects: voxel-wise correlation map ----- #
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t2_dir = scans_dir + '\\' + this_subject + '\\' + this_subject + t2_dir[iSubject]
            if not os.path.exists(this_subject_t2_dir):
                continue
            
            # ----- Print progress to console ----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\tCreating FOV-mask for subject: " + this_subject)
            
            # ----- Define filename of output FOV mask ----- #
            save_fname = this_subject_t2_dir + '\\FOV_mask_' + this_subject + t2_dir[iSubject] + '.nii'
            if os.path.exists(save_fname):
                os.remove(save_fname) # Delete pre-existing mask
            
            # ----- Find all scans for this subject ----- #
            all_scans = glob.glob(this_subject_t2_dir + '\\' + scans_prefix + this_subject + '*.nii')
            try:
                this_subject_FOV_data = np.zeros(np.shape(nb.load(all_scans[0]).get_fdata()))
            except:
                print("Process failed for subject: " + this_subject)
                continue
            for iScan, this_scan in enumerate(all_scans):
                this_scan_data = nb.load(this_scan).get_fdata()
                
                # ----- Binarize the scan data ----- #
                this_scan_data[np.where(this_scan_data != 0)] = 1
                
                # ----- Add non-zero voxels to mask ----- #
                this_subject_FOV_data += this_scan_data
            
            # ----- Binarize the subject's FOV mask; write to .nii file ----- #
            this_subject_FOV_data[np.where(this_subject_FOV_data < len(all_scans))] = 0
            this_subject_FOV_data[np.where(this_subject_FOV_data == len(all_scans))] = 1
            this_subject_FOV_scan = nb.Nifti1Image(this_subject_FOV_data, nb.load(this_scan).affine)
            nb.save(this_subject_FOV_scan, save_fname)
            
            # ----- Print progress to console ----- #
            print('Code executed with no errors')
    
    def run_preprocessing(self):
        """
        Preprocessing input data (via Matlab).
        
        Scan prefix manual:
            r = realigned
            a = slice-time corrected
            c = coregistered
            c1-c5 = segmentation output
            n = normalized
            s = smoothed
            (e - eroded, not part of this method, see erode_segmentation_mask)
        
        """
        
        # ----- Print progress to console ----- #
        print("\nRunning preprocssing procedure ...")
        
        # ----- Specify the preprocessing type ----- #
        type_of_processing = int(input("What type of preprocessing do you want to perform? Enter 1 for slice-timing correction, 2 for realignment, 3 for coregistration, 4 for normalization, 5 for segmentation, or 6 for smoothing: "))
        if type_of_processing == 1:
            script_name = "SliceTimingCorrection"
        elif type_of_processing == 2:
            script_name = "Realignment"
        elif type_of_processing == 3:
            script_name = "Coregistration"
        elif type_of_processing == 4:
            script_name = "Normalization"
        elif type_of_processing == 5:
            script_name = "Segmentation"
        elif type_of_processing == 6:
            script_name = "Smoothing"
        else:
            raise(ValueError('This is not a valid input value.'))
        
        # ----- Define argument (passed to Matlab) ----- #
        scans_prefix = input("Optional: Enter a prefix for the (T1 or T2) scans (inputs depend on the preprocessing procedure): ")
        scans_dir = self.scans_dir + '_' + self.scans_type
        study_ID = self.study_ID
        working_dir = self.work_dir
        log_fname = self.study_ID + "_" + script_name + "_" + self.scans_type + "_log.txt"
        if os.path.exists(working_dir + '\\' + log_fname):
            os.remove(working_dir + '\\' + log_fname)
        
        # ----- Loop over subjects: Run preprocessing ----- #        
        for iSubject, this_subject in enumerate(self.subjects_list):
            t1_dir = self.subjects_T1[iSubject]
            if self.scans_type == 'REST':
                t2_dir = self.subjects_REST[iSubject]
            elif self.scans_type == 'EMO':
                t2_dir = self.subjects_EMO[iSubject]
            
            # ----- Print progress to console ----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\t" + script_name + " for subject: " + this_subject)
            
            # ----- Print progress to log file ----- #
            fh = open(working_dir + "\\" + log_fname, 'a')
            fh.write(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\t" + script_name + " for subject: " + this_subject + "\n")
            fh.close()
            
            # ----- Run data preprocessing script in Matlab (via shell) ----- #
            matlab_cmd_p1 = "matlab -nodesktop -nosplash -wait -r \"addpath('" + working_dir + "'); "
            matlab_cmd_p2 = "start_up('" + working_dir + "'); "
            matlab_cmd_p3 = script_name + "(" + str(iSubject) + ",'" + this_subject + "','" + t1_dir + "','" + t2_dir + "','" + scans_dir + "','" + scans_prefix + "','" + study_ID + "','" + working_dir + "','" + log_fname + "'); quit"
            try:
                os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
                
                # ----- Print progress to console ----- #
                print('Code executed with no errors')
            except:
                # ----- Print failure to log file ----- #
                fh = open(working_dir + "\\" + log_fname, 'a')
                fh.write("\t\t" + script_name + " failed for subject: " + this_subject + "\n")
                fh.close()
    
    def extract_FD_jenkinson(self):
        """
        Obtain subject-level framewise displacement using the Jenkinson method,
        and write to subject-level directories
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting framewise displacements (FD Jenkinson) ...")
        
        # ----- Define variables ----- #
        refscan = input("Enter the reference scan for framewise displacement computation. Enter either an integer (n'th volume) or M (mean image) (default: enter 1): ")
        
        # ----- Loop over subjects ----- #        
        for iSubject, this_subject in enumerate(self.subjects_list):
            if self.scans_type == 'REST':
                this_subject_t2_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_REST[iSubject]
            elif self.scans_type == 'EMO':
                this_subject_t2_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_EMO[iSubject]
                        
            # ----- Print progress to console ----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\tCalculating framewise displacement for subject: " + this_subject)
            
            # ----- Specify the input file for FD extraction ----- #
            if os.path.exists(this_subject_t2_dir):
                try:
                    input_filename = glob.glob(this_subject_t2_dir + '\\rp_*.txt')[0]
                except:
                    print("Missing data for subject: " + this_subject)
                    continue
            else:
                print("Missing data for subject: " + this_subject)
                continue
            
            # ----- Specify the reference scan and output filename ----- #
            if refscan.isdigit():
                input_reference = glob.glob(this_subject_t2_dir + '\\' + this_subject + '*' + refscan + '.nii')[0]
                output_filename = this_subject_t2_dir + '\\FD_Jenkinson_' + input_reference.split('\\')[-1].split('.')[0] + '.txt'
            elif refscan == 'M':
                input_reference = glob.glob(this_subject_t2_dir + '\\mean' + '*.nii')[0]
                output_filename = this_subject_t2_dir + '\\FD_Jenkinson_' + input_reference.split('\\')[-1].split('.')[0] + '.txt'
            
            # ----- Run y_FD_Jenkinson script in Matlab (via shell) ----- #
            matlab_cmd_p1 = "matlab -automation -wait -r \"addpath('" + self.work_dir + "');"
            matlab_cmd_p2 = "start_up('" + self.work_dir + "'); "
            matlab_cmd_p3 = "FrameWiseDisplacement('" + input_filename + "', '" + input_reference + "', '" + output_filename + "'); quit"
            try:
                os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
                
                # ----- Print progress to console ----- #
                print('Code executed with no errors')
            except:
                # ----- Print failure to console ----- #
                print('Execution of Matlab code failed for subject ' + this_subject)
    
    def erode_segmentation_mask(self):
        """
        Erosion of the white-matter or CSF segmentations (via Matlab).
        
        Scan prefix manual:
            c2 = white-matter segmentation (in anatomical directory)
            c3 = CSF segmentation (in anatomical directory)
        
        """
        
        # ----- Print progress to console ----- #
        print("\nEroding segmentation mask ...")
        
        # ----- Define argument (passed to Matlab) ----- #
        working_dir = self.work_dir
        t1_scans_prefix = input("On what segmentation mask does the erosion need to be applied? Enter a scan prefix (default: enter c2 or c3): ")
        
        # ----- Loop over subjects: Run preprocessing ----- #        
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t1_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]
            
            # ----- Print progress to console ----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\t" + "Eroding segmented data for subject: " + this_subject)
            
            # ----- Find the segmentation file ----- #
            this_subject_file = glob.glob(this_subject_t1_dir + '\\' + t1_scans_prefix + '*.nii')
            if this_subject_file:
                this_subject_file = this_subject_file[0]
            else:
                print('The segmentation file could not be found')
                continue
            
            # ----- Run data preprocessing script in Matlab (via shell) ----- #
            matlab_cmd_p1 = "matlab -nodesktop -nosplash -wait -r \"addpath('" + working_dir + "');"
            matlab_cmd_p2 = "start_up('" + working_dir + "'); "
            matlab_cmd_p3 = "Erosion('" + this_subject_file + "'); quit"
            try:
                os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
                
                # ----- Print progress to console ----- #
                print('Code executed with no errors')
            except:
                # ----- Print failure to console ----- #
                print('Execution of Matlab code failed for subject ' + this_subject)
        
    def run_filtering(self, signal, Fs, cutFreq):
        """
        Filtering according to the ezfilt algorithm, as originally written in
        Matlab by Thomas Gladwin: https://www.tegladwin.com/code.php
        
        """
        
        signal = np.reshape(np.array(signal), [len(signal), 1])
        fcrit = cutFreq
        
        # Create X
        X = np.empty(0)
        t = np.arange(1, signal.shape[0] + 1)
        t = np.divide((t - 0), Fs)
        T = t[-1]
        f0 = 0
        while 1 == 1:
            cos0 = np.cos(np.multiply((2 * np.pi * f0), t))
            if not X.size:
                X = cos0
            else:
                X = np.c_[X, cos0]
            f0 = f0 + 0.5 / T
            if f0 > fcrit:
                break # Break out of loop
        X = np.matrix(X)
        
        # Fit
        beta = (X.T @ X).I @ X.T @ signal
        slow = np.array(X @ beta)
        fast = np.array(signal - slow)
        
        return slow, fast


class EmotionTask(Preprocessing):
    
    def extract_behavioral_data(self):
        """ 
        Extract onsets and durations from task logfile, and write to CSV file. 
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting emotion task behavioral data ...\n")
        
        # ----- Set file attributes ----- #
        n_blocks = 96
        start_trial = 3
        end_trial = 218
        picture_codes = ["777", "888", "999"]
        neu_val = 1
        pos_val = 2
        neg_val = 3
        
        # ----- Loop over all subjects in object ----- #
        # check_subjects = []
        for iSubject, (this_subjID, this_subjPPN) in enumerate(zip(self.subjects_list, self.subjects_PPN)):
                        
            # ----- Obtain logfile for this subject ----- #
            this_logfile = glob.glob(self.logs_dir + "\\" + this_subjID + "*aug2010.log")
            if not this_logfile:
                this_logfile = glob.glob(self.logs_dir + "\\" + this_subjPPN[-3:] + "*aug2010.log")
            if not this_logfile:
                print("This logfile does not exist")
                continue
            
            # ----- Print progress to console ----- #
            print(str(iSubject + 1) + "\tExtracting behavioral log data for subject " + this_subjID + "(" + this_subjPPN + ")")
            
            # ----- Pre-allocate iteration variables ----- #
            n_empty = 0
            miss_val = 0
            
            all_onsets = []
            all_pictures = []
            all_responses = []
            all_RTs = []
            
            # ----- Extract logdata from this logfile ----- #
            fh = open(this_logfile[0], "r")
            for iLine, this_line in enumerate(fh):
                this_line = this_line.rstrip() # Remove new-line ('\n') characters
                this_line = this_line.split('\t') # Read line tab-delimited ('\t')
                
                # ----- Contiue to next line for non-data entries ----- #
                if len(this_line) < 2:
                    if not this_line[0]:
                        n_empty += 1
                    continue # non-data entry: continue to next line
                
                # ----- Extract column header indices ----- #
                if n_empty < 2 and len(this_line) > 2:
                    try:
                        trial_index = this_line.index("Trial")
                        event_type_index = this_line.index("Event Type")
                        picture_index = this_line.index("Code")
                        response_index = this_line.index("Code")
                        onset_index = this_line.index("Time")
                        RT_index = this_line.index("TTime")
                        continue # continue to next line
                    except:
                        continue # continue to next line                            
                
                # ----- Continue to next line for pre=task trials ----- #
                this_trial = int(this_line[trial_index])
                if this_trial < start_trial:
                    continue # continue to next line
                
                # ----- Extract data from valid trials ----- #
                if this_trial%2:
                    if any(c in this_line[picture_index] for c in picture_codes):
                        all_pictures.append(this_line[picture_index])
                        all_onsets.append(float(this_line[onset_index]))
                        consecutive_trial = this_trial + 1
                    if "bedankt" in this_line[picture_index]:
                        break # end of task reached: break out of loop
                if not this_trial%2 and this_trial == consecutive_trial:
                    if "evaluatie" in this_line[picture_index]:
                        all_responses.append(miss_val) # missing value as placeholder
                        all_RTs.append(miss_val) # missing value as placeholder
                    if "Response" in this_line[event_type_index]:
                        all_responses[-1] = this_line[response_index][-1]
                        all_RTs[-1] = this_line[RT_index]
                    if "fill" in this_line[picture_index]:
                        continue # continue to next line
                
                # ----- Continue until end trial reached ----- #
                if this_trial > end_trial:
                    break # end trial reached: break out of loop
                
            fh.close()
            
            # ----- Adjust all onsets for onset of first trial ----- #
            all_onsets = np.asarray(all_onsets)
            all_onsets = list(all_onsets - all_onsets[0])
            
            # ----- Match picture codes to conditions ----- #
            all_conditions = []
            for iPicture, this_picture in enumerate(all_pictures):
                if picture_codes[neu_val - 1] in this_picture:
                    all_conditions.append(neu_val)
                elif picture_codes[pos_val - 1] in this_picture:
                    all_conditions.append(pos_val)
                elif picture_codes[neg_val - 1] in this_picture:
                    all_conditions.append(neg_val)
                else:
                    all_conditions.append(miss_val)
            
            # ----- Extract congruent and incongruent trials ----- #
            all_congruencies = []
            all_incongruencies = []
            for iElement, this_element in enumerate(zip(all_conditions, all_responses)):
                if (this_element[0] == int(this_element[-1])):
                    all_congruencies.append(this_element[-1])
                    all_incongruencies.append(0)
                else: 
                    all_congruencies.append(0)
                    all_incongruencies.append(this_element[-1])
            
            # ----- Create output matrix with logdata ----- #
            column_headers = ['onsets', 'pictures', 'conditions', 'responses', 'congruent', 'incongruent', 'RTs']
            data_matrix = np.array([all_onsets, all_pictures, all_conditions, all_responses, all_congruencies, all_incongruencies, all_RTs]).transpose()
            output_logdata = pd.DataFrame(data_matrix, columns=column_headers)
            output_fn = self.onsets_dir + "\\" + this_subjID + "-" + this_logfile[0].split('\\')[-1][0:-4] + ".csv"
            output_logdata.to_csv(output_fn, index=False)
            
            # ----- Print progress to console ----- #
            if output_logdata.shape[0] < n_blocks:
                # check_subjects.append([iSubject, this_subjID])
                print("Check output for subject " + this_subjID + "(" + this_subjPPN + ")")
            else:
                print("Code executed with no errors")
    
    def extract_confound_regressors(self):
        """ 
        Concatanate the 6 realignment parameters, white-matter, and
        cerebrospinal fluid (9P) into a single confound matrix, to be used in 
        first-level GLM analysis.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting confound regressors ...")
        
        # ----- Loop over subjects: write confound models ----- #
        wm_csf_scans_prefix = input("Optional: Enter a prefix for the white-matter and csf masks (default: enter e for eroded): ")
        t2_scans_prefix = input("Optional: Enter a prefix for the resting-state scans (default: enter nra for normalized, (coregistered), realigned, and slice-time corrected): ")
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_emo_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_EMO[iSubject]
            this_subject_t1_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]
            
            # ----- Print progress to console----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + ' Computing confound model for subject ' + this_subject)
            if not os.path.exists(this_subject_emo_dir):
                print("Missing data for subject " + this_subject)
                continue
            if not os.path.exists(this_subject_t1_dir):
                print("Missing data for subject " + this_subject)
                continue
            
            try:
                # ----- Obtain white-matter and CSF segmentations ----- #
                this_wm_scan_ID = glob.glob(this_subject_t1_dir + '\\' + wm_csf_scans_prefix + 'c2*.nii')
                this_wm_mask = nb.load(this_wm_scan_ID[0]).get_fdata()
                wm_rows, wm_cols, wm_stks = np.where(this_wm_mask > 0) # voxel subscripts
                
                this_csf_scan_ID = glob.glob(this_subject_t1_dir + '\\' + wm_csf_scans_prefix + 'c3*.nii')
                this_csf_mask = nb.load(this_csf_scan_ID[0]).get_fdata()
                csf_rows, csf_cols, csf_stks = np.where(this_csf_mask > 0) # voxel subscripts
                
                # ----- Pre-allocate confound predictors ----- #
                wm_signal = np.array([])
                csf_signal = np.array([])
                
                # ----- Loop over all scans ----- #
                all_rs_scan_IDs = glob.glob(this_subject_emo_dir + '\\' + t2_scans_prefix + this_subject + '*.nii')
                for iScan, this_scan_ID in enumerate(all_rs_scan_IDs):
                    this_img = nb.load(this_scan_ID)
                    this_img_array = this_img.get_fdata()
                    
                    # ----- Extract global mean, white matter, CSF signals ----- #
                    wm_signal = np.append(wm_signal, np.nanmean(this_img_array[wm_rows, wm_cols, wm_stks]))
                    csf_signal = np.append(csf_signal, np.nanmean(this_img_array[csf_rows, csf_cols, csf_stks]))
                
                # ----- Extract motion parameters ----- #
                this_subject_mp_dir = glob.glob(this_subject_emo_dir + '\\rp_*txt')[0]
                motion_pars = np.array([])
                fh = open(this_subject_mp_dir, 'r')
                for iRow, this_line in enumerate(fh):
                    this_line_values = this_line.strip().split()
                    motion_pars = np.r_[motion_pars, np.array([float(v) for v in this_line_values])]
                fh.close()                
                motion_pars = np.reshape(motion_pars, [iRow + 1, len(this_line_values)])
                
                # ----- Compute confound matrix ----- #
                confound_matrix = np.c_[motion_pars, wm_signal, csf_signal]
                
                # ----- Write the confound matrix to (.csv) file ----- #
                savename = this_subject_emo_dir + '\\' + this_subject + self.subjects_EMO[iSubject] + '_confound_regressors.csv'
                pd.DataFrame(confound_matrix).to_csv(savename, index=False, header=False)
                
                # ----- Print progress to console ----- #
                print("Code executed with no errors")
            except:
                # ----- Print progress to console ----- #
                print("Process failed for subject: " + this_subject)
    
    def extract_spike_regressors(self):
        """
        Calculate the number of FD outliers for all subjects based on 
        predefined threshold (recommended: 0.2), and obtain spike regressors 
        for first-level GLM analysis
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting spike regressors ...")
        
        # ----- Define variables ----- #
        threshold = float(input("Enter framewise displament threshold for number-of-outliers calculation (default: 0.5): " ))
        
        # ----- Generate output number of outliers logfile ----- #
        output_fname = self.study_ID + '_Outliers_FD_Jenkinson_EMO_log.txt'
        if os.path.exists(self.work_dir + '\\' + output_fname):
            os.remove(self.work_dir + '\\' + output_fname)
        fh = open(self.work_dir + '\\' + output_fname, 'a')
        fh.write('%s\t%s\t%s\t%s\n' % ('SubjID', 'Group', 'MFD', 'nOutliers'))
        fh.close()
        
        # ----- Loop over subjects ----- #        
        for iSubject, this_subject in enumerate(self.subjects_list):
            
            # ----- Print progress to console ----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\tCalculating the number of outliers, and extracting spike regressors for subject: " + this_subject)
            
            # ----- Find FD Jenkinson file ----- #
            this_subject_directory = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_EMO[iSubject]
            FD_jenkinson = glob.glob(this_subject_directory + '\\FD_Jenkinson*.txt')
            if not FD_jenkinson:
                print("Missing data for subject: " + this_subject)
                continue
            
            # ----- Calculate number of outliers for this subject ----- #
            this_subject_FD_jenkinson = np.loadtxt(FD_jenkinson[0])
            this_subject_MFD = np.nanmean(this_subject_FD_jenkinson)
            this_subject_all_outliers = np.array(np.where(this_subject_FD_jenkinson > threshold))
            if this_subject_all_outliers[0].size == 0:
                this_subject_n_outliers = 0
                
                # ----- Compute empty spike regressors for this subject ----- #
                all_spike_regressors = np.array([])
                
                # ----- Save spike regressor matrix to (.csv) file ----- #
                savename = this_subject_directory + '\\' + this_subject + self.subjects_EMO[iSubject] + '_spike_regressors.csv'
                if os.path.exists(savename):
                    os.remove(savename)
                pd.DataFrame(all_spike_regressors).to_csv(savename, index=False, header=False)
                
                # ----- Write NaN as number of outliers for this subject ----- #
                fh = open(self.work_dir + '\\' + output_fname, 'a')
                fh.write('%s\t%d\t%f\t%d\n' % (this_subject, self.subjects_group[iSubject], this_subject_MFD, this_subject_n_outliers))
                fh.close()
                
                # ----- Print progress to console ----- #
                print('No outliers detected for subject ' + this_subject + ', but code executed with no errors')
            else:
                this_subject_n_outliers = this_subject_all_outliers.shape[-1]
                
                # ----- Compute spike regressors for this subject ----- #
                all_spike_regressors = np.zeros((np.shape(this_subject_FD_jenkinson)[0], this_subject_n_outliers))
                for iSpike, this_spike_index in enumerate(this_subject_all_outliers[-1]):
                    this_spike_regressor = np.zeros(np.shape(this_subject_FD_jenkinson))
                    this_spike_regressor[this_spike_index] = 1
                    all_spike_regressors[:, iSpike] =+ this_spike_regressor
                
                # ----- Save spike regressor matrix to (.csv) file ----- #
                savename = this_subject_directory + '\\' + this_subject + self.subjects_EMO[iSubject] + '_spike_regressors.csv'
                if os.path.exists(savename):
                    os.remove(savename)
                pd.DataFrame(all_spike_regressors).to_csv(savename, index=False, header=False)
                
                # ----- Write subject's number of outliers to file ----- #
                fh = open(self.work_dir + '\\' + output_fname, 'a')
                fh.write('%s\t%d\t%f\t%d\n' % (this_subject, self.subjects_group[iSubject], this_subject_MFD, this_subject_n_outliers))
                fh.close()
                
                # ----- Print progress to console ----- #
                print("Code executed with no errors")
    
    def run_1st_level_analysis(self):
        """ 
        First level statistical analysis of the emotion task data (via Matlab).
        
        """
        
        # ----- Print progress to console ----- #
        print("\nRunning first-level analysis ...")
        
        # ----- Define arguments (passed to Matlab) ----- #
        working_dir = self.work_dir
        onsets_dir = self.onsets_dir
        scans_dir = self.scans_dir + '_' + self.scans_type
        log_fname = self.study_ID + "_1st_level_" + self.scans_type + "_log.txt"
        if os.path.exists(working_dir + "\\" + log_fname):
            os.remove(working_dir + "\\" + log_fname)
        explicit_mask = input("Optional: Enter an explicit mask (filename as listed in the subject's directory, default: press enter to continue): ")
        
        # ----- Loop over subjects: Run preprocessing ----- #        
        t2_scans_prefix = input("Optional: Enter a prefix for the emotion task scans (default: enter nra for normalized, (coregistered), realigned and slice-time corrected): ")
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t1_dir = self.subjects_T1[iSubject]
            this_subject_t2_dir = self.subjects_EMO[iSubject]
            
            # ----- Define the explicit mask [optional] ----- #
            if explicit_mask:
                try:
                    mask_fn = glob.glob(self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + this_subject_t1_dir + '\\' + explicit_mask)[0].split('\\')[-1]
                except:
                    mask_fn = explicit_mask
            else:
                mask_fn = explicit_mask
            
            # ----- Print progress to console ----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\tFirst level analysis for subject: " + this_subject)
            
            # ----- Print progress to log file ----- #
            fh = open(working_dir + "\\" + log_fname, 'a')
            fh.write(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\tFirst level analysis for subject: " + this_subject + "\n")
            fh.close()
            
            # ----- Run data preprocessing script in Matlab (via shell) ----- #
            matlab_cmd_p1 = "matlab -nodesktop -nosplash -wait -r \"addpath('" + working_dir + "');"
            matlab_cmd_p2 = "start_up('" + working_dir + "'); "
            matlab_cmd_p3 = "FirstLevelAnalysis(" + str(iSubject) + ",'" + this_subject + "','" + working_dir + "','" + this_subject_t1_dir + "','" + this_subject_t2_dir + "','" + t2_scans_prefix + "','" + mask_fn + "','" + scans_dir + "','" + onsets_dir + "','" + log_fname + "'); quit"
            try:
                os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
                
                # ----- Print progress to console ----- #
                print("Code executed with no errors")
            except:
                # ----- Print failure to log file ----- #
                fh = open(working_dir + "\\" + log_fname, 'a')
                fh.write("\t\tFirst level analysis failed for subject: " + this_subject + "\n")
                fh.close()
                
                # ----- Print failure to console ----- #
                print("First level analysis failed for subject: " + this_subject)
    
    def create_native_roi_mask(self):
        """
        Transform an ROI mask to native/subject space using inverse
        normalization and coregistration procedures in SPM12.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting native ROI mask ...")
        
        # ----- Specify the ROI (global mask) ----- #
        roi_fn = input("Enter the filename of the ROI mask (as listed in the working directory): ")
        
        # ----- Loop over subjects ----- #
        t1_scans_prefix = input("Optional: Enter a prefix for the T1 native reference image (default: enter c for coregistered): ")
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t1_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]

            # ----- Print progress to console----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + " Computing native roi mask for subject " + this_subject)
            if not os.path.exists(this_subject_t1_dir):
                print("Missing data for subject " + this_subject)
                continue
            
            # ----- Reslice ROI mask from global to native space ----- #
            try:
                template_fn = glob.glob(this_subject_t1_dir + '\\' + t1_scans_prefix + this_subject + '*.nii')[0]
                matlab_cmd_p1 = "matlab -nodesktop -nosplash -wait -r \"addpath('" + self.work_dir + "');"
                matlab_cmd_p2 = "start_up('" + self.work_dir + "'); "
                matlab_cmd_p3 = "GlobalToNativeMask('" + self.work_dir + "','" + template_fn + "','" + roi_fn + "'); quit"
            
                os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
                
                # ----- Print progress to console ----- #
                print("Code executed with no errors")
            except:
                # ----- Print progress to console ----- #
                print("Global to native transformation of ROI mask failed for subject: " + this_subject)
    
    def obtain_global_roi_mask(self):
        """
        Copy and reslice an ROI mask to normalized subject reference space 
        using SPM12.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting global ROI mask ...")
        
        # ----- Specify the ROI (global mask) ----- #
        roi_fn = input("Enter the filename of the ROI mask (as listed in the working directory): ")
        
        # ----- Loop over subjects ----- #
        t1_scans_prefix = input("Optional: Enter a prefix for the T1 reference image (default: enter n for normalized): ")
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t1_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]

            # ----- Print progress to console----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + " Computing native roi mask for subject " + this_subject)
            if not os.path.exists(this_subject_t1_dir):
                print("Missing data for subject " + this_subject)
                continue
            
            # ----- Reslice ROI mask to reference space ----- #
            try:
                template_fn = glob.glob(this_subject_t1_dir + '\\' + t1_scans_prefix + this_subject + '*.nii')[0]
                matlab_cmd_p1 = "matlab -nodesktop -nosplash -wait -r \"addpath('" + self.work_dir + "');"
                matlab_cmd_p2 = "start_up('" + self.work_dir + "'); "
                matlab_cmd_p3 = "CoregisterGlobalMask('" + self.work_dir + "','" + template_fn + "','" + roi_fn + "'); quit"
            
                os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
                
                # ----- Print progress to console ----- #
                print("Code executed with no errors")
            except:
                # ----- Print progress to console ----- #
                print("Reslicing of global ROI mask failed for subject: " + this_subject)
    
    def extract_roi_masked_spm_data(self):
        """
        Extract statistical output from parametric maps using an ROI mask.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting ROI signal ...")
        
        # ----- Specify the ROI (global mask) ----- #
        roi_fn = input("Enter the filename of the ROI mask (as listed in the subject's directory): ")
        mask_type = int(input("Is the mask binary or probabilistic? Enter 1 (binary) or 2 (probabilistic): "))
        if (mask_type != 1) and (mask_type != 2):
            raise ValueError('This is not a valid input value: Enter 1 (for a binary mask) or 2 (for a probabilistic mask).')
        
        # ----- Specify the input scans name ----- #
        input_scans = input("Enter the name of the first level output to apply the ROI mask to (as listed in the subject directory): ")
        
        # ----- Loop over subjects: extract ROI masked data ----- #
        all_lh_roi_means = []
        all_rh_roi_means = []
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t1_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]
            this_subject_emo_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_EMO[iSubject]
            
            # ----- Print progress to console----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + ' Extracting roi regressor for subject ' + this_subject)
            if not os.path.exists(this_subject_emo_dir):
                print("Missing data for subject " + this_subject)
                # ----- Insert NaN for missing subjects ----- #
                all_lh_roi_means.append(np.nan)
                all_rh_roi_means.append(np.nan)
                continue
            if not os.path.exists(this_subject_t1_dir):
                print("Missing data for subject " + this_subject)
                # ----- Insert NaN for missing subjects ----- #
                all_lh_roi_means.append(np.nan)
                all_rh_roi_means.append(np.nan)
                continue
            
            try:
                # ----- Load native ROI mask ----- #
                this_roi_fn = glob.glob(this_subject_t1_dir + '\\' + roi_fn)[0]
                roi_mask = nb.load(this_roi_fn)
                roi_mask_data = roi_mask.get_fdata()
                if mask_type == 1: # binary mask
                    roi_mask_data[np.isnan(roi_mask_data)] = 0 # Set NaN values to zero
                    roi_mask_data[np.where(roi_mask_data < 0.9)] = 0 # Set non-Bolean values to zero
                    roi_mask_data[np.where(roi_mask_data >= 0.9)] = 1 # Set non-Bolean values to zero
                elif mask_type == 2: # probabilistic mask
                    roi_mask_data[np.isnan(roi_mask_data)] = 0 # Set NaN values to zero
                    roi_mask_data[np.where(roi_mask_data < 0.05)] = 0 # Set low probability values to zero
                               
                # ----- Load input scan data ----- #
                this_scan_ID = glob.glob(this_subject_emo_dir + '\\' + 'first_level_analysis' + '\\' + input_scans)[0]
                this_scan = nb.load(this_scan_ID)
                this_scan_data = this_scan.get_fdata()
                this_scan_data[np.isnan(this_scan_data)] = 0 # Set NaN values to zero
                
                # ----- Extract ROI-weighted data  ----- #
                roi_weighted_data = np.multiply(this_scan_data, roi_mask_data)
                
                # ----- Obtain the standard (mni) coordinates of all voxels ----- #        
                all_mni_XYZ = np.array(np.meshgrid(*(range(i) for i in this_scan.shape), indexing='ij'))
                all_mni_XYZ = np.rollaxis(all_mni_XYZ, 0, len(this_scan.shape) + 1)
                all_mni_XYZ = nb.affines.apply_affine(this_scan.affine, all_mni_XYZ)
                
                # ----- Obtain subscripts of roi weighted data ----- #
                roi_XYZ = np.nonzero(roi_weighted_data)
                roi_X = roi_XYZ[0][:] # All row subscripts
                roi_Y = roi_XYZ[1][:] # All column subscripts
                roi_Z = roi_XYZ[2][:] # All stack subscripts
                
                # ----- Concatenate each data value and XYZ into new vectors ----- #
                roi_voxel_values = np.array([])
                roi_voxel_mni_XYZ = np.empty([3,0])
                roi_voxel_weights = np.array([])
                for i, (x, y, z) in enumerate(zip(roi_X, roi_Y, roi_Z)):
                    roi_voxel_values = np.append(roi_voxel_values, roi_weighted_data[x, y, z])
                    roi_voxel_mni_XYZ = np.c_[roi_voxel_mni_XYZ, all_mni_XYZ[x, y, z]]
                    roi_voxel_weights = np.append(roi_voxel_weights, roi_mask_data[x, y, z])
                
                # ----- Split data according to hemisphere ----- #
                lh_roi_voxel_inds = np.where(roi_voxel_mni_XYZ[0,:] < 0)
                lh_roi_voxel_values = roi_voxel_values[lh_roi_voxel_inds[0]]
                # lh_roi_voxel_mni_XYZ = roi_voxel_mni_XYZ[:, lh_roi_voxel_inds[0]]
                lh_roi_voxel_weights = roi_voxel_weights[lh_roi_voxel_inds[0]]
                
                rh_roi_voxel_inds = np.where(roi_voxel_mni_XYZ[0,:] > 0)
                rh_roi_voxel_values = roi_voxel_values[rh_roi_voxel_inds[0]]
                # rh_roi_voxel_mni_XYZ = roi_voxel_mni_XYZ[:, rh_roi_voxel_inds[0]]
                rh_roi_voxel_weights = roi_voxel_weights[rh_roi_voxel_inds[0]]
                    
                # ----- Compute (probability-weighted) mean of each hemisphere's roi weighted data ----- #
                if mask_type == 1:
                    all_lh_roi_means.append(np.nanmean(lh_roi_voxel_values))
                    all_rh_roi_means.append(np.nanmean(rh_roi_voxel_values))
                if mask_type == 2:
                    all_lh_roi_means.append(np.nansum(lh_roi_voxel_values) / np.nansum(lh_roi_voxel_weights)) # weighted mean
                    # all_lh_roi_means.append(np.nansum(lh_roi_voxel_values)) # weighted sum
                    all_rh_roi_means.append(np.nansum(rh_roi_voxel_values) / np.nansum(rh_roi_voxel_weights)) # weighted mean
                    # all_rh_roi_means.append(np.nansum(rh_roi_voxel_values)) # weighted sum
                
                # ----- Print progress to console ----- #
                print("Code executed with no errors")
            except:
                # ----- Insert NaN for failed subjects ----- #
                all_lh_roi_means.append(np.nan)
                all_rh_roi_means.append(np.nan)
                
                # ----- Print progress to console ----- #
                print("Extraction of ROI regressors (left/right) failed for subject: " + this_subject)
        
        # ----- Create output matrix ----- #
        all_subjIDs = np.asarray(self.subjects_list)
        all_lh_roi_means = np.asarray(all_lh_roi_means)
        all_rh_roi_means = np.asarray(all_rh_roi_means)
        output_matrix = np.array([all_subjIDs, all_lh_roi_means, all_rh_roi_means]).transpose()
        
        # ----- Save output to (.csv) file(s) ----- #
        column_headers = ['subjID', 'HemiL', 'HemiR']
        output_matrix = pd.DataFrame(output_matrix, columns=column_headers)
        output_fn = self.work_dir + "\\" + self.study_ID + "_" + self.scans_type + "_Mean_" + input_scans.split('.')[0].capitalize() + ".csv"
        # output_fn = self.work_dir + "\\" + self.study_ID + "_" + self.scans_type + "_Sum_" + input_scans.split('.')[0].capitalize() + ".csv"
        output_matrix.to_csv(output_fn)
        

class RestingState(Preprocessing):
    
    def create_native_roi_mask(self):
        """
        Transform an ROI mask to native/subject space using inverse
        normalization and coregistration procedures in SPM12.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting native ROI mask ...")
        
        # ----- Specify the ROI (global mask) ----- #
        roi_fn = input("Enter the filename of the ROI mask (as listed in the working directory): ")
        
        # ----- Loop over subjects: write confound models ----- #
        t1_scans_prefix = input("Optional: Enter a prefix for the T1 native reference image (default: enter c for coregistered): ")
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t1_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]
            
            # ----- Print progress to console----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + " Computing native roi mask for subject " + this_subject)
            if not os.path.exists(this_subject_t1_dir):
                print("Missing data for subject " + this_subject)
                continue
            
            # ----- Reslice ROI mask from global to native space ----- #
            try:
                template_fn = glob.glob(this_subject_t1_dir + '\\' + t1_scans_prefix + this_subject + '*.nii')[0]
                matlab_cmd_p1 = "matlab -nodesktop -nosplash -wait -r \"addpath('" + self.work_dir + "');"
                matlab_cmd_p2 = "start_up('" + self.work_dir + "'); "
                matlab_cmd_p3 = "GlobalToNativeMask('" + self.work_dir + "','" + template_fn + "','" + roi_fn + "'); quit"
                
                os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
                
                # ----- Print progress to console ----- #
                print("Code executed with no errors")
            except:
                # ----- Print progress to console ----- #
                print("Global to native transformation of ROI mask failed for subject: " + this_subject)
    
    def obtain_global_roi_mask(self):
        """
        Copy and reslice an ROI mask to normalized subject reference space 
        using SPM12.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting global ROI mask ...")
        
        # ----- Specify the ROI (global mask) ----- #
        roi_fn = input("Enter the filename of the ROI mask (as listed in the working directory): ")
        
        # ----- Loop over subjects ----- #
        t1_scans_prefix = input("Optional: Enter a prefix for the T1 reference image (default: enter n for normalized): ")
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t1_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]

            # ----- Print progress to console----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + " Computing native roi mask for subject " + this_subject)
            if not os.path.exists(this_subject_t1_dir):
                print("Missing data for subject " + this_subject)
                continue
            
            # ----- Reslice ROI mask to reference space ----- #
            try:
                template_fn = glob.glob(this_subject_t1_dir + '\\' + t1_scans_prefix + this_subject + '*.nii')[0]
                matlab_cmd_p1 = "matlab -nodesktop -nosplash -wait -r \"addpath('" + self.work_dir + "');"
                matlab_cmd_p2 = "start_up('" + self.work_dir + "'); "
                matlab_cmd_p3 = "CoregisterGlobalMask('" + self.work_dir + "','" + template_fn + "','" + roi_fn + "'); quit"
            
                os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
                
                # ----- Print progress to console ----- #
                print("Code executed with no errors")
            except:
                # ----- Print progress to console ----- #
                print("Reslicing of global ROI mask failed for subject: " + this_subject)
    
    def extract_roi_regressors(self):
        """
        Extract the ROI regressor for each subject, for use in the voxel-wise 
        connectivity analysis.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting ROI signal ...")
        
        # ----- Set filter properties ----- #
        TR = 1.6
        low_freq = 0.01
        high_freq = 0.08
        
        # ----- Specify the ROI (global mask) ----- #
        roi_fn = input("Enter the filename of the ROI mask (as listed in the subject's directory): ")
        mask_type = int(input("Is the mask binary or probabilistic? Enter 1 (binary) or 2 (probabilistic): "))
        if (mask_type != 1) and (mask_type != 2):
            raise ValueError('This is not a valid input value: Enter 1 (for a binary mask) or 2 (for a probabilistic mask).')
        
        # ----- Loop over subjects: write confound models ----- #
        rs_scans_prefix = input("Optional: Enter a prefix for the resting-state scans (default: enter nr for normalized (coregistered) and realigned): ")
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t1_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]
            this_subject_rs_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_REST[iSubject]
            
            # ----- Print progress to console----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + ' Extracting roi regressor for subject ' + this_subject)
            if not os.path.exists(this_subject_rs_dir):
                print("Missing data for subject " + this_subject)
                continue
            if not os.path.exists(this_subject_t1_dir):
                print("Missing data for subject " + this_subject)
                continue
            
            # ----- Load native ROI mask ----- #
            this_roi_fn = glob.glob(this_subject_t1_dir + '\\' + roi_fn)[0]
            roi_mask = nb.load(this_roi_fn)
            roi_mask_data = roi_mask.get_fdata()
            if mask_type == 1: # binary mask
                roi_mask_data[np.isnan(roi_mask_data)] = 0 # Set NaN values to zero
                roi_mask_data[np.where(roi_mask_data < 0.9)] = 0 # Set non-Bolean values to zero
                roi_mask_data[np.where(roi_mask_data >= 0.9)] = 1 # Set non-Bolean values to zero
            elif mask_type == 2: # probabilistic mask
                roi_mask_data[np.isnan(roi_mask_data)] = 0 # Set NaN values to zero
                roi_mask_data[np.where(roi_mask_data < 0.05)] = 0 # Set low probability values to zero
            
            try:
                # ----- Pre-allocate predictor variable ----- #
                lh_roi_signal = np.array([])
                rh_roi_signal = np.array([])
                
                # ----- Loop over all scans ----- #
                all_rs_scan_IDs = glob.glob(this_subject_rs_dir + '\\' + rs_scans_prefix + this_subject + '*.nii')
                for iScan, this_scan_ID in enumerate(all_rs_scan_IDs):
                    this_scan = nb.load(this_scan_ID)
                    this_scan_data = this_scan.get_fdata()
                    this_scan_data[np.isnan(this_scan_data)] = 0 # Set NaN values to zero
                    
                    # ----- Extract ROI-weighted data  ----- #
                    roi_weighted_data = np.multiply(this_scan_data, roi_mask_data)
                    
                    # ----- Obtain the standard (mni) coordinates of all voxels ----- #        
                    all_mni_XYZ = np.array(np.meshgrid(*(range(i) for i in this_scan.shape), indexing='ij'))
                    all_mni_XYZ = np.rollaxis(all_mni_XYZ, 0, len(this_scan.shape) + 1)
                    all_mni_XYZ = nb.affines.apply_affine(this_scan.affine, all_mni_XYZ)
                    
                    # ----- Obtain subscripts of roi weighted data ----- #
                    roi_XYZ = np.nonzero(roi_weighted_data)
                    roi_X = roi_XYZ[0][:] # All row subscripts
                    roi_Y = roi_XYZ[1][:] # All column subscripts
                    roi_Z = roi_XYZ[2][:] # All stack subscripts
                    
                    # ----- Concatenate each data value and XYZ into new vectors ----- #
                    roi_voxel_values = np.array([])
                    roi_voxel_mni_XYZ = np.empty([3,0])
                    roi_voxel_weights = np.array([])
                    for i, (x, y, z) in enumerate(zip(roi_X, roi_Y, roi_Z)):
                        roi_voxel_values = np.append(roi_voxel_values, roi_weighted_data[x, y, z])
                        roi_voxel_mni_XYZ = np.c_[roi_voxel_mni_XYZ, all_mni_XYZ[x, y, z]]
                        roi_voxel_weights = np.append(roi_voxel_weights, roi_mask_data[x, y, z])
                    
                    # ----- Split data according to hemisphere ----- #
                    lh_roi_voxel_inds = np.where(roi_voxel_mni_XYZ[0,:] < 0)
                    lh_roi_voxel_values = roi_voxel_values[lh_roi_voxel_inds[0]]
                    # lh_roi_voxel_mni_XYZ = roi_voxel_mni_XYZ[:, lh_roi_voxel_inds[0]]
                    lh_roi_voxel_weights = roi_voxel_weights[lh_roi_voxel_inds[0]]
                    
                    rh_roi_voxel_inds = np.where(roi_voxel_mni_XYZ[0,:] > 0)
                    rh_roi_voxel_values = roi_voxel_values[rh_roi_voxel_inds[0]]
                    # rh_roi_voxel_mni_XYZ = roi_voxel_mni_XYZ[:, rh_roi_voxel_inds[0]]
                    rh_roi_voxel_weights = roi_voxel_weights[rh_roi_voxel_inds[0]]
                    
                    # ----- Compute (probability-weighted) mean of each hemisphere's roi weighted data ----- #
                    if mask_type == 1:
                        lh_roi_signal = np.append(lh_roi_signal, np.nanmean(lh_roi_voxel_values))
                        rh_roi_signal = np.append(rh_roi_signal, np.nanmean(rh_roi_voxel_values))
                    if mask_type == 2:
                        lh_roi_signal = np.append(lh_roi_signal, (np.nansum(lh_roi_voxel_values) / np.nansum(lh_roi_voxel_weights))) # weighted mean
                        # lh_roi_signal = np.append(lh_roi_signal, np.nansum(lh_roi_voxel_values)) # weighted sum
                        rh_roi_signal = np.append(rh_roi_signal, (np.nansum(rh_roi_voxel_values) / np.nansum(rh_roi_voxel_weights))) # weighted mean
                        # rh_roi_signal = np.append(rh_roi_signal, np.nansum(rh_roi_voxel_values)) # weighted sum
                
                #----- Filter the roi regressors ----- #
                unused, lh_roi_signal = self.run_filtering(lh_roi_signal, 1/TR, low_freq) # high pass
                lh_roi_signal, unused = self.run_filtering(lh_roi_signal, 1/TR, high_freq) # low pass
                
                unused, rh_roi_signal = self.run_filtering(rh_roi_signal, 1/TR, low_freq) # high pass
                rh_roi_signal, unused = self.run_filtering(rh_roi_signal, 1/TR, high_freq) # low pass
                
                # ----- Write the roi regressors to (.csv) file ----- #
                lh_savename = this_subject_rs_dir + '\\' + this_subject + self.subjects_REST[iSubject] + '_lh_roi regressor.csv'
                pd.DataFrame(lh_roi_signal).to_csv(lh_savename, index=False, header=False)
                
                rh_savename = this_subject_rs_dir + '\\' + this_subject + self.subjects_REST[iSubject] + '_rh_roi regressor.csv'
                pd.DataFrame(rh_roi_signal).to_csv(rh_savename, index=False, header=False)
                
                # ----- Print progress to console ----- #
                print("Code executed with no errors")
            except:
                # ----- Print progress to console ----- #
                print("Extraction of ROI regressors (left/right) failed for subject: " + this_subject)
        
    def extract_confound_regressors(self):
        """ 
        Concatanate the 6 realignment parameters, global signal, white-matter, 
        cerebrospinal fluid (9P), the temporal derivatives (18P), quadratic 
        terms (27P), and quadratic terms of the temporal derivatives (36P), 
        into a single confound matrix, to be used in voxel-wise GLM analysis.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting confound regressors ...")
        
        # ----- Set filter properties ----- #
        TR = 1.6
        low_freq = 0.01
        high_freq = 0.08
        
        # ----- Loop over subjects: write confound models ----- #
        gm_scans_prefix = input("Optional: Enter a prefix for the gray matter mask (default: press enter to continue): ")
        wm_csf_scans_prefix = input("Optional: Enter a prefix for the white-matter and csf masks (default: enter e for eroded): ")
        t2_scans_prefix = input("Optional: Enter a prefix for the resting-state scans (default: enter nr for normalized (coregistered) and realigned): ")
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_rs_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_REST[iSubject]
            this_subject_t1_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]
            
            # ----- Print progress to console----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + ' Computing confound model for subject ' + this_subject)
            if not os.path.exists(this_subject_rs_dir):
                print("Missing data for subject " + this_subject)
                continue
            if not os.path.exists(this_subject_t1_dir):
                print("Missing data for subject " + this_subject)
                continue
            
            try:
                # ----- Obtain grey segmentation ----- #
                this_gm_scan_ID = glob.glob(this_subject_t1_dir + '\\' + gm_scans_prefix + 'c1*.nii')
                this_gm_mask = nb.load(this_gm_scan_ID[0]).get_fdata()
                gm_rows, gm_cols, gm_stks = np.where(this_gm_mask > 0) # voxel subscripts
                
                # ----- Obtain global signal ----- #
                this_raw_wm_scan_ID = glob.glob(this_subject_t1_dir + '\\c2*.nii')
                this_raw_wm_mask = nb.load(this_raw_wm_scan_ID[0]).get_fdata()
                raw_wm_rows, raw_wm_cols, raw_wm_stks = np.where(this_raw_wm_mask > 0) # voxel subscripts
                
                this_raw_csf_scan_ID = glob.glob(this_subject_t1_dir + '\\c3*.nii')
                this_raw_csf_mask = nb.load(this_raw_csf_scan_ID[0]).get_fdata()
                raw_csf_rows, raw_csf_cols, raw_csf_stks = np.where(this_raw_csf_mask > 0) # voxel subscripts            
                
                glob_rows = np.concatenate((gm_rows, raw_wm_rows, raw_csf_rows)) # row subscripts
                glob_cols = np.concatenate((gm_cols, raw_wm_cols, raw_csf_cols)) # colum subscripts
                glob_stks = np.concatenate((gm_stks, raw_wm_stks, raw_csf_stks)) # slice subscripts
                
                # ----- Obtain white-matter and CSF segmentations ----- #
                this_wm_scan_ID = glob.glob(this_subject_t1_dir + '\\' + wm_csf_scans_prefix + 'c2*.nii')
                this_wm_mask = nb.load(this_wm_scan_ID[0]).get_fdata()
                wm_rows, wm_cols, wm_stks = np.where(this_wm_mask > 0) # voxel subscripts
                
                this_csf_scan_ID = glob.glob(this_subject_t1_dir + '\\' + wm_csf_scans_prefix + 'c3*.nii')
                this_csf_mask = nb.load(this_csf_scan_ID[0]).get_fdata()
                csf_rows, csf_cols, csf_stks = np.where(this_csf_mask > 0) # voxel subscripts 
                
                # ----- Pre-allocate confound predictors ----- #
                wm_signal = np.array([])
                csf_signal = np.array([])
                glob_signal = np.array([])
                
                # ----- Loop over all scans ----- #
                all_rs_scan_IDs = glob.glob(this_subject_rs_dir + '\\' + t2_scans_prefix + this_subject + '*.nii')
                for iScan, this_scan_ID in enumerate(all_rs_scan_IDs):
                    this_img = nb.load(this_scan_ID)
                    this_img_array = this_img.get_fdata()
                    
                    # ----- Extract global mean, white matter, CSF signals ----- #
                    wm_signal = np.append(wm_signal, np.nanmean(this_img_array[wm_rows, wm_cols, wm_stks]))
                    csf_signal = np.append(csf_signal, np.nanmean(this_img_array[csf_rows, csf_cols, csf_stks]))
                    glob_signal = np.append(glob_signal, np.nanmean(this_img_array[glob_rows, glob_cols, glob_stks]))
                
                # ----- Extract motion parameters ----- #
                this_subject_mp_dir = glob.glob(this_subject_rs_dir + '\\rp_*txt')[0]
                motion_pars = np.array([])
                fh = open(this_subject_mp_dir, 'r')
                for iRow, this_line in enumerate(fh):
                    this_line_values = this_line.strip().split()
                    motion_pars = np.r_[motion_pars, np.array([float(v) for v in this_line_values])]
                fh.close()                
                motion_pars = np.reshape(motion_pars, [iRow + 1, len(this_line_values)])
                
                # ----- Calculate confound expansions ----- #
                motion_pars_dfdx = np.r_[np.zeros((1, np.shape(motion_pars)[-1])), motion_pars[0:-1]] # Friston et al. (1996): Movement-Related effects in fMRI time-series
                # motion_pars_dfdx = np.r_[np.zeros((1, np.shape(motion_pars)[-1])), np.subtract(motion_pars[1:], motion_pars[0:-1])] # Satterthwaite et al. (2013): An improved framework for confound regression and filtering for control of motion artifact in the preprocessing of resting-state functional connectivity data
                motion_pars2 = np.square(motion_pars)
                motion_pars_dfdx2 = np.square(motion_pars_dfdx)
                
                # glob_signal_dfdx = np.r_[0, glob_signal[0:-1]] # Friston et al. (1996): Movement-Related effects in fMRI time-series
                glob_signal_dfdx = np.r_[0, np.subtract(glob_signal[1:], glob_signal[0:-1])] # Satterthwaite et al. (2013): An improved framework for confound regression and filtering for control of motion artifact in the preprocessing of resting-state functional connectivity data
                glob_signal2 = np.square(glob_signal)
                glob_signal_dfdx2 = np.square(glob_signal_dfdx)
                
                # wm_signal_dfdx = np.r_[0, wm_signal[0:-1]] # Friston et al. (1996): Movement-Related effects in fMRI time-series
                wm_signal_dfdx = np.r_[0, np.subtract(wm_signal[1:], wm_signal[0:-1])] # Satterthwaite et al. (2013): An improved framework for confound regression and filtering for control of motion artifact in the preprocessing of resting-state functional connectivity data
                wm_signal2 = np.square(wm_signal)
                wm_signal_dfdx2 = np.square(wm_signal_dfdx)
                
                # csf_signal_dfdx = np.r_[0, csf_signal[0:-1]] # Friston et al. (1996): Movement-Related effects in fMRI time-series
                csf_signal_dfdx = np.r_[0, np.subtract(csf_signal[1:], csf_signal[0:-1])] # Satterthwaite et al. (2013): An improved framework for confound regression and filtering for control of motion artifact in the preprocessing of resting-state functional connectivity data
                csf_signal2 = np.square(csf_signal)
                csf_signal_dfdx2 = np.square(csf_signal_dfdx)
                
                # ----- Specify the (unfiltered) confound model (36P) ----- #
                confound_matrix_raw = np.c_[motion_pars, glob_signal, wm_signal, csf_signal,
                                        motion_pars_dfdx, glob_signal_dfdx, wm_signal_dfdx, csf_signal_dfdx,
                                        motion_pars2, glob_signal2, wm_signal2, csf_signal2,
                                        motion_pars_dfdx2, glob_signal_dfdx2, wm_signal_dfdx2, csf_signal_dfdx2]
                
                # ----- Compute filtered confound matrix ----- #
                confound_matrix = np.zeros(np.shape(confound_matrix_raw))
                for iCol in range(np.shape(confound_matrix_raw)[-1]):
                    this_confound = np.reshape(confound_matrix_raw[:,iCol], [np.shape(confound_matrix_raw)[0], 1])
                    
                    #----- Filter the confound regressor ----- #
                    unused, this_confound = self.run_filtering(this_confound, 1/TR, low_freq) # high pass
                    this_confound, unused = self.run_filtering(this_confound, 1/TR, high_freq) # low pass
                    
                    # ----- Concatenate filtered confound variables ----- #
                    confound_matrix[:,iCol] = confound_matrix[:,iCol] + this_confound.T
                    
                # ----- Write the confound matrix to (.csv) file ----- #
                savename = this_subject_rs_dir + '\\' + this_subject + self.subjects_REST[iSubject] + '_confound_regressors.csv'
                pd.DataFrame(confound_matrix).to_csv(savename, index=False, header=False)
                
                # ----- Print progress to console ----- #
                print("Code executed with no errors")
            except:
                # ----- Print progress to console ----- #
                print("Process failed for subject: " + this_subject)
    
    def extract_spike_regressors(self):
        """
        Calculate the number of FD outliers for all subjects based on 
        predefined threshold (recommended: 0.2), and obtain spike regressors 
        for voxel-wise GLM analysis
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting spike regressors ...")
        
        # ----- Define variables ----- #
        threshold = float(input("Enter framewise displament threshold for number-of-outliers calculation (default: 0.2): " ))
        
        # ----- Generate output number of outliers logfile ----- #
        output_fname = self.study_ID + '_Outliers_FD_Jenkinson_REST_log.txt'
        if os.path.exists(self.work_dir + '\\' + output_fname):
            os.remove(self.work_dir + '\\' + output_fname)
        fh = open(self.work_dir + '\\' + output_fname, 'a')
        fh.write('%s\t%s\t%s\t%s\n' % ('SubjID', 'Group', 'MFD', 'nOutliers'))
        fh.close()
        
        # ----- Loop over subjects ----- #        
        for iSubject, this_subject in enumerate(self.subjects_list):
            
            # ----- Print progress to console ----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\tCalculating the number of outliers, and extracting spike regressors for subject: " + this_subject)
            
            # ----- Find FD Jenkinson file ----- #
            this_subject_directory = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_REST[iSubject]
            FD_jenkinson = glob.glob(this_subject_directory + '\\FD_Jenkinson*.txt')
            if not FD_jenkinson:
                print("Missing data for subject: " + this_subject)
                continue
            
            # ----- Calculate number of outliers for this subject ----- #
            this_subject_FD_jenkinson = np.loadtxt(FD_jenkinson[0])
            this_subject_MFD = np.nanmean(this_subject_FD_jenkinson)
            this_subject_all_outliers = np.array(np.where(this_subject_FD_jenkinson > threshold))
            if this_subject_all_outliers[0].size == 0:
                this_subject_n_outliers = 0
                
                # ----- Compute empty spike regressors for this subject ----- #
                all_spike_regressors = np.array([])
                
                # ----- Save spike regressor matrix to (.csv) file ----- #
                savename = this_subject_directory + '\\' + this_subject + self.subjects_REST[iSubject] + '_spike_regressors.csv'
                if os.path.exists(savename):
                    os.remove(savename)
                pd.DataFrame(all_spike_regressors).to_csv(savename, index=False, header=False)
                
                # ----- Write NaN as number of outliers for this subject ----- #
                fh = open(self.work_dir + '\\' + output_fname, 'a')
                fh.write('%s\t%d\t%f\t%d\n' % (this_subject, self.subjects_group[iSubject], this_subject_MFD, this_subject_n_outliers))
                fh.close()
                
                # ----- Print progress to console ----- #
                print('No outliers detected for subject ' + this_subject + ', but code executed with no errors')
            else:
                this_subject_n_outliers = this_subject_all_outliers.shape[-1]
                
                # ----- Compute spike regressors for this subject ----- #
                all_spike_regressors = np.zeros((np.shape(this_subject_FD_jenkinson)[0], this_subject_n_outliers))
                for iSpike, this_spike_index in enumerate(this_subject_all_outliers[-1]):
                    this_spike_regressor = np.zeros(np.shape(this_subject_FD_jenkinson))
                    this_spike_regressor[this_spike_index] = 1
                    all_spike_regressors[:, iSpike] =+ this_spike_regressor
                
                # ----- Save spike regressor matrix to (.csv) file ----- #
                savename = this_subject_directory + '\\' + this_subject + self.subjects_REST[iSubject] + '_spike_regressors.csv'
                if os.path.exists(savename):
                    os.remove(savename)
                pd.DataFrame(all_spike_regressors).to_csv(savename, index=False, header=False)
                
                # ----- Write subject's number of outliers to file ----- #
                fh = open(self.work_dir + '\\' + output_fname, 'a')
                fh.write('%s\t%d\t%f\t%d\n' % (this_subject, self.subjects_group[iSubject], this_subject_MFD, this_subject_n_outliers))
                fh.close()
                
                # ----- Print progress to console ----- #
                print("Code executed with no errors")
    
    def voxel_wise_connectivity_analysis(self):
        """ 
        Seed-based voxel-wise regrrssion analysis of the resting-state data.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nConducting voxel-wise connectivity analysis ...")
        
        # ----- Pre-allocate output logfile ----- #
        log_fname = self.study_ID + '_Connectivity_' + self.scans_type + '_log.txt'
        if os.path.exists(self.work_dir + '\\' + log_fname):
            os.remove(self.work_dir + '\\' + log_fname)
        
        # ----- Set filter properties ----- #
        TR = 1.6
        low_freq = 0.01
        high_freq = 0.08
        
        # ----- Pre-set to ignore divide-by-zero errors ----- #
        np.seterr(divide='ignore', invalid='ignore')
        
        # ----- Loop over subjects ----- #
        t2_scans_prefix = input("Optional: Enter a prefix for the resting-state scans on (default: enter nr for normalized (coregistered) and realigned): ")
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_rs_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_REST[iSubject]
            this_subject_t1_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]
            
            # ----- Print progress to console----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\tSeed-based connectivity (1st-level) analysis for subject " + this_subject)
            
            # ----- Print progress to log file----- #
            fh = open(self.work_dir + '\\' + log_fname, 'a')
            fh.write(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\tConnectivity analysis for subject: " + this_subject + "\n")
            fh.close()
            
            # ----- Skip subject if exclusion index is zero ----- #
            if self.subjects_Include[iSubject] == 0:
                print('This subject has an include value of zero')
                fh = open(self.work_dir + "\\" + log_fname, 'a')
                fh.write("\t\tThis subject has an include value of zero\n")
                fh.close()
                continue
            
            # ----- Clock start-time ----- #
            start_time = time.time()
            
            try:
                # ----- Obtain grey matter segmentation ----- #
                this_gm_scan_ID = glob.glob(this_subject_t1_dir + '\\' + 'c1*.nii')
                this_gm_mask = nb.load(this_gm_scan_ID[0]).get_fdata()
                gm_rows, gm_cols, gm_stks = np.where(this_gm_mask > 0) # voxel subscripts
                
                gm_coordinates = np.c_[gm_rows, gm_cols, gm_stks]
                
                # ----- Partial volume correction: white matter ----- #
                this_wm_scan_ID = glob.glob(this_subject_t1_dir + '\\' + 'ec2*.nii')
                this_wm_mask = nb.load(this_wm_scan_ID[0]).get_fdata()
                wm_rows, wm_cols, wm_stks = np.where(this_wm_mask > 0) # voxel subscripts
                
                wm_coordinates = np.c_[wm_rows, wm_cols, wm_stks]
                
                wm_overlap_idx = []
                for this_wm_coordinate in wm_coordinates:
                    this_overlap_idx = np.where((gm_coordinates == this_wm_coordinate).all(axis=1))[0]
                    if this_overlap_idx.size == 1:
                        wm_overlap_idx.append(this_overlap_idx[0])
                gm_coordinates = np.delete(gm_coordinates, wm_overlap_idx, axis=0)
                
                # ----- Partial volume correction: cerebrospinal fluid ----- #
                this_csf_scan_ID = glob.glob(this_subject_t1_dir + '\\' + 'ec3*.nii')
                this_csf_mask = nb.load(this_csf_scan_ID[0]).get_fdata()
                csf_rows, csf_cols, csf_stks = np.where(this_csf_mask > 0) # voxel subscripts
                
                csf_coordinates = np.c_[csf_rows, csf_cols, csf_stks]
                
                csf_overlap_idx = []
                for this_csf_coordinate in csf_coordinates:
                    this_overlap_idx = np.where((gm_coordinates == this_csf_coordinate).all(axis=1))[0]
                    if this_overlap_idx.size == 1:
                        csf_overlap_idx.append(this_overlap_idx[0])
                gm_coordinates = np.delete(gm_coordinates, csf_overlap_idx, axis=0)
                
                # ----- Obtain (pre-filtered) roi regressors (left/right) ----- #
                lh_roi_regressor = np.loadtxt(this_subject_rs_dir + '\\' + this_subject + self.subjects_REST[iSubject] + '_lh_roi regressor.csv', delimiter=',')
                rh_roi_regressor = np.loadtxt(this_subject_rs_dir + '\\' + this_subject + self.subjects_REST[iSubject] + '_rh_roi regressor.csv', delimiter=',')
                
                # ----- Obtain (pre-filtered) nuisance variables ----- #
                confound_regressors = np.loadtxt(this_subject_rs_dir + '\\' + this_subject + self.subjects_REST[iSubject] + '_confound_regressors.csv', delimiter=',')
                
                # ----- Obtain spike regressors ----- #
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore") # Catch empty file warning
                    spike_regressors = np.loadtxt(this_subject_rs_dir + '\\' + this_subject + self.subjects_REST[iSubject] + '_spike_regressors.csv', delimiter=',')
                    
                # ----- Load all rs scans for this subject ----- #
                all_rs_scans = []
                
                all_rs_scan_IDs = glob.glob(this_subject_rs_dir + '\\' + t2_scans_prefix + this_subject + '*.nii')
                for iScan, this_rs_scan_ID in enumerate(all_rs_scan_IDs):
                    this_rs_scan = nb.load(this_rs_scan_ID)
                    this_rs_scan_data = this_rs_scan.get_fdata()
                    
                    all_rs_scans.append(this_rs_scan_data)
                
                # ----- Pre-allocate beta map variables  ----- #
                this_subject_b_map_lh = np.zeros(this_rs_scan_data.shape)
                this_subject_t_map_lh = np.zeros(this_rs_scan_data.shape)
                
                this_subject_b_map_rh = np.zeros(this_rs_scan_data.shape)
                this_subject_t_map_rh = np.zeros(this_rs_scan_data.shape)
                
                # ----- Loop over all brain voxels ----- #
                for iVoxel, this_voxel_coordinate in enumerate(gm_coordinates):
                    this_voxel_coordinate = tuple(this_voxel_coordinate)
                    
                    # ----- Obtain the voxel signal (outcome variable) ----- #
                    this_voxel_signal = np.empty(0)
                    for iScan, this_rs_scan in enumerate(all_rs_scans):
                        this_voxel_signal = np.append(this_voxel_signal, all_rs_scans[iScan][this_voxel_coordinate])
                    
                    # ----- Filter the voxel signal (outcome variable) ----- #
                    unused, this_voxel_signal = self.run_filtering(this_voxel_signal, 1/TR, low_freq) # high pass
                    this_voxel_signal, unused = self.run_filtering(this_voxel_signal, 1/TR, high_freq) # low pass
                    
                    y = scipy.stats.zscore(this_voxel_signal, ddof=1)
                    
                    # ----- Run the GLM for this voxel: left hemisphere ----- #
                    X_lh = np.c_[lh_roi_regressor, confound_regressors]
                    if spike_regressors.size == 0:
                        X_lh = np.c_[np.ones(lh_roi_regressor.shape), scipy.stats.zscore(X_lh, ddof=1)]
                    else:
                        X_lh = np.c_[np.ones(lh_roi_regressor.shape), scipy.stats.zscore(X_lh, ddof=1), spike_regressors]
                    B_lh = np.array(np.matrix(X_lh.T @ X_lh).I @ X_lh.T @ y).flatten()
                    
                    MSE_lh = np.sum(np.square(y.flatten() - (X_lh @ B_lh))) / (X_lh.shape[0] - X_lh.shape[1]) # Mean square error
                    VAR_lh = MSE_lh * np.matrix(X_lh.T @ X_lh).I # Variance-covariance matrix
                    
                    this_subject_b_map_lh[this_voxel_coordinate] = B_lh[1]
                    this_subject_t_map_lh[this_voxel_coordinate] = B_lh[1] / np.sqrt(VAR_lh[1, 1])
                    
                    # ----- Run the GLM for this voxel: right hemisphere ----- #
                    X_rh = np.c_[rh_roi_regressor, confound_regressors]
                    if spike_regressors.size == 0:
                        X_rh = np.c_[np.ones(rh_roi_regressor.shape), scipy.stats.zscore(X_rh, ddof=1)]
                    else:
                        X_rh = np.c_[np.ones(rh_roi_regressor.shape), scipy.stats.zscore(X_rh, ddof=1), spike_regressors]
                    B_rh = np.array(np.matrix(X_rh.T @ X_rh).I @ X_rh.T @ y).flatten()
                    
                    MSE_rh = np.sum(np.square(y.flatten() - (X_rh @ B_rh))) / (X_rh.shape[0] - X_rh.shape[1]) # Mean square error
                    VAR_rh = MSE_rh * np.matrix(X_rh.T @ X_rh).I # Variance-covariance matrix
                    
                    this_subject_b_map_rh[this_voxel_coordinate] = B_rh[1]
                    this_subject_t_map_rh[this_voxel_coordinate] = B_rh[1] / np.sqrt(VAR_rh[1, 1])
                    
                # ----- Write subject's b-maps to (.nii) files ----- #
                this_subject_b_map_lh[np.isnan(this_subject_b_map_lh)] = 0 # Set NaN values to zero
                lh_b_nifti_img = nb.Nifti1Image(this_subject_b_map_lh, nb.load(all_rs_scan_IDs[0]).affine)
                lh_b_savename = this_subject_rs_dir + "\\connectivity_b_map_lh_" + this_subject + '.nii'
                nb.save(lh_b_nifti_img, lh_b_savename)
                
                this_subject_b_map_rh[np.isnan(this_subject_b_map_rh)] = 0 # Set NaN values to zero
                rh_b_nifti_img = nb.Nifti1Image(this_subject_b_map_rh, nb.load(all_rs_scan_IDs[0]).affine)
                rh_b_savename = this_subject_rs_dir + "\\connectivity_b_map_rh_" + this_subject + '.nii'
                nb.save(rh_b_nifti_img, rh_b_savename)
                
                # ----- Write subject's t-maps to (.nii) files ----- #
                this_subject_t_map_lh[np.isnan(this_subject_t_map_lh)] = 0 # Set NaN values to zero
                lh_t_nifti_img = nb.Nifti1Image(this_subject_t_map_lh, nb.load(all_rs_scan_IDs[0]).affine)
                lh_t_savename = this_subject_rs_dir + "\\connectivity_t_map_lh_" + this_subject + '.nii'
                nb.save(lh_t_nifti_img, lh_t_savename)
                
                this_subject_t_map_rh[np.isnan(this_subject_t_map_rh)] = 0 # Set NaN values to zero
                rh_t_nifti_img = nb.Nifti1Image(this_subject_t_map_rh, nb.load(all_rs_scan_IDs[0]).affine)
                rh_t_savename = this_subject_rs_dir + "\\connectivity_t_map_rh_" + this_subject + '.nii'
                nb.save(rh_t_nifti_img, rh_t_savename)
                
                # ----- Print progress to logfile ----- #
                fh = open(self.work_dir + '\\' + log_fname, 'a')
                fh.write("\t\tCode executed with no errors" + "\n")
                fh.close()
                
                # ----- Print progress to console----- #       
                print("Code executed with no errors")
            except:
                # ----- Print failure to log file ----- #
                fh = open(self.work_dir + '\\' + log_fname, 'a')
                fh.write("\t\tConnectivity map failed for subject: " + this_subject + "\n")
                fh.close()
                
                # ----- Print failure to console ----- #
                print("Process failed for subject: " + this_subject)
                
            # ----- Clock end-time ----- #
            stop_time = time.time()
            print('Duration: ' + str(stop_time - start_time) + ' seconds')
            
            fh = open(self.work_dir + '\\' + log_fname, 'a')
            fh.write("\t\tDuration: " + str(stop_time - start_time) + " seconds" + "\n")
            fh.close()


class Postprocessing(Preprocessing):
    
    def create_study_FOV_mask(self):
        """
        Extract inclusive FOV mask across subjects.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting study inclusive FOV mask ...")
        
        # ----- Define variables  ----- #
        scans_prefix = input("Optional: Enter a prefix for the functional scans (default: enter FOV_mask_ for the FOV masks): ")
        if self.scans_type == 'REST':
            t2_dir = self.subjects_REST
        elif self.scans_type == 'EMO':
            t2_dir = self.subjects_EMO
        else:
            raise(ValueError('This is not a valid input value.'))
        scans_dir = self.scans_dir + '_' + self.scans_type
        working_dir = self.work_dir
        study_ID = self.study_ID
        working_dir = self.work_dir
        
        # ----- Pre-set to ignore deprecation warning ----- #
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        
        # ----- Obtain template scan ----- #
        template_scan = glob.glob(scans_dir + '\\' + self.subjects_list[0] + '\\' + self.subjects_list[0] + t2_dir[0] + '\\' + scans_prefix + self.subjects_list[0] + t2_dir[0] + '*.nii')[0]
        
        # ----- Loop over subjects: voxel-wise correlation map ----- #
        n_valid_subjects = 0
        container_brain = np.zeros(nb.load(template_scan).shape)
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t2_dir = scans_dir + '\\' + this_subject + '\\' + this_subject + t2_dir[iSubject]
            
            # ----- Print progress to console ----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\tExtracting FOV-mask for subject: " + this_subject)            
            if not os.path.exists(this_subject_t2_dir):
                print('Error: missing data for subject: ' + this_subject)
                continue
            
            # ----- Skip subject if exclusion index is zero ----- #
            if self.subjects_Include[iSubject] == 0:
                print('Warning: This subject has an include value of zero')
                continue
            
            # ----- Find FOV mask for this subject ----- #
            this_subject_FOV_mask = glob.glob(this_subject_t2_dir + '\\' + scans_prefix  + this_subject + t2_dir[iSubject] + '*.nii')
            if not this_subject_FOV_mask:
                print('Error: FOV file missing for subject: ' + this_subject)
                continue
            this_subject_FOV_data = nb.load(this_subject_FOV_mask[0]).get_fdata()
            
            # ----- Binarize the subject's FOV mask ----- #
            this_subject_FOV_data[np.where(this_subject_FOV_data > 0)] = 1
            
            # ----- Add subject's FOV mask to all subjects' FOV mask ----- #
            container_brain += this_subject_FOV_data
            n_valid_subjects += 1
            
            # ----- Print progress to console ----- #
            print('Code executed with no errors')
        
        # ----- Binarize across subjects FOV mask, write to .nii file ----- #
        all_subjects_FOV_data = np.divide(container_brain, n_valid_subjects)
        all_subjects_FOV_data[np.where(all_subjects_FOV_data < 0.9)] = 0
        all_subjects_FOV_scan = nb.Nifti1Image(all_subjects_FOV_data, nb.load(template_scan).affine)
        save_fname = working_dir + '\\' + study_ID + '_Inclusive_FOV_Mask_'+ self.scans_type + '.nii'
        nb.save(all_subjects_FOV_scan, save_fname)
    
    def create_study_grey_matter_mask(self):
        """
        Extract inclusive grey matter mask across subjects.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nExtracting study grey matter mask ...")
        
        # ----- Define variables ----- #
        scans_prefix = input("Optional: Enter a prefix for the T1 scans (default: enter c1n for the grey matter segmentations of the normalized T1 scans): ")
        scans_dir = self.scans_dir + '_' + self.scans_type
        study_ID = self.study_ID
        working_dir = self.work_dir
        
        # ----- Pre-set to ignore deprecation warning ----- #
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        
        # ----- Obtain template scan ----- #
        template_scan = glob.glob(scans_dir + '\\' + self.subjects_list[0] + '\\' + self.subjects_list[0] + self.subjects_T1[0] + '\\' + scans_prefix + self.subjects_list[0] + self.subjects_T1[0] + '*.nii')[0]
        
        # ----- Compute sum of all segmentations (container variable) ----- #
        n_valid_subjects = 0
        container_brain = np.zeros(nb.load(template_scan).shape)
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t1_dir = scans_dir + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]
            
            # ----- Print progress to console ----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\tExtracting grey matter mask for subject: " + this_subject)
            if not os.path.exists(this_subject_t1_dir):
                print('Error: missing data for subject: ' + this_subject)
                continue
            
            # ----- Skip subject if exclusion index is zero ----- #
            if self.subjects_Include[iSubject] == 0:
                print('This subject has an include value of zero')
                continue
            
            # ----- Find grey matter segmentation for this subject ----- #
            this_subject_normalized_gm_mask = glob.glob(this_subject_t1_dir + '\\' + scans_prefix + this_subject + self.subjects_T1[iSubject] + '*.nii')
            if not this_subject_normalized_gm_mask:
                print('Error: grey matter segmentation missing for subject: ' + this_subject)
                continue
            this_subject_normalized_gm_mask = nb.load(this_subject_normalized_gm_mask[0]).get_fdata()
            
            # ----- Binarize the grey matter segmentation mask ----- #
            this_subject_normalized_gm_mask[np.where(this_subject_normalized_gm_mask > 0)] = 1
            
            # ----- Add subject's segmentation to container  ----- #
            container_brain += this_subject_normalized_gm_mask
            n_valid_subjects += 1
            
            # ----- Print progress to console ----- #
            print('Code executed with no errors')
        
        # ----- Create average brain, write to .nii file ----- #
        all_subjects_GM_data = np.divide(container_brain, n_valid_subjects)
        all_subjects_GM_data[np.where(all_subjects_GM_data <= 0.5)] = 0
        all_subjects_GM_data[np.where(all_subjects_GM_data > 0.5)] = 1
        all_subjects_GM_scan = nb.Nifti1Image(all_subjects_GM_data, nb.load(template_scan).affine)
        save_fname = working_dir + '\\' + study_ID + '_Inclusive_GM_Mask_'+ self.scans_type + '.nii'
        nb.save(all_subjects_GM_scan, save_fname)
            
    def create_average_brain(self):
        """
        Compute average brain from normalized T1 scans.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nComputing average brain ...")
        
        # ----- Define variables ----- #
        scans_prefix = input("Optional: Enter a prefix for the T1 scans (default: enter high_res_n for the high-resolution normalized T1 scans): ")
        scans_dir = self.scans_dir + '_' + self.scans_type
        study_ID = self.study_ID
        working_dir = self.work_dir
        
        # ----- Pre-set to ignore deprecation warning ----- #
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        
        # ----- Obtain template scan ----- #
        template_scan = glob.glob(scans_dir + '\\' + self.subjects_list[0] + '\\' + self.subjects_list[0] + self.subjects_T1[0] + '\\' + scans_prefix + self.subjects_list[0] + self.subjects_T1[0] + '*.nii')[0]
        
        # ----- Compute sum of all brains (container variable) ----- #
        n_valid_subjects = 0
        container_brain = np.zeros(nb.load(template_scan).shape)
        for iSubject, this_subject in enumerate(self.subjects_list):
            this_subject_t1_dir = scans_dir + '\\' + this_subject + '\\' + this_subject + self.subjects_T1[iSubject]
            
            # ----- Print progress to console ----- #
            print(str(iSubject + 1) + " (of " + str(len(self.subjects_list)) + ")" + "\tExtracting normalized T1 scan for subject: " + this_subject)
            if not os.path.exists(this_subject_t1_dir):
                print('Error: missing data for subject: ' + this_subject)
                continue
            
            # ----- Skip subject if exclusion index is zero ----- #
            if self.subjects_Include[iSubject] == 0:
                print('This subject has an include value of zero')
                continue
            
            # ----- Find T1 scan for this subject ----- #
            this_subject_normalized_t1 = glob.glob(this_subject_t1_dir + '\\' + scans_prefix + this_subject + self.subjects_T1[iSubject] + '*.nii')
            if not this_subject_normalized_t1:
                print('Error: normalized T1 scan missing for subject: ' + this_subject)
                continue
            this_subject_normalized_t1 = nb.load(this_subject_normalized_t1[0]).get_fdata()
            
            # ----- Add subject's T1 scan to container  ----- #
            container_brain += this_subject_normalized_t1
            n_valid_subjects += 1
            
            # ----- Print progress to console ----- #
            print('Code executed with no errors')
        
        # ----- Create average brain, write to .nii file ----- #
        average_brain_data = np.divide(container_brain, n_valid_subjects)
        average_brain_scan = nb.Nifti1Image(average_brain_data, nb.load(template_scan).affine)
        save_fname = working_dir + '\\' + study_ID + '_Average_Brain_'+ self.scans_type + '.nii'
        nb.save(average_brain_scan, save_fname)
    
    def motion_correction_benchmark_qc_fc(self):
        """
        Calculate the voxel-wise correlations betweeen the mean framewise 
        displacement (MFD) and functional connectivity (FC) values, as a motion 
        correction benchmark, across subjects per dataset.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nComputing QC-FC benchmarking ...")
        
        # ----- Define variables ----- #
        mfd_fname = input("Enter the filename of the study MFD data (as listed in the working directory): ")
        gm_mask = input("Enter the filename of the study grey matter mask (as listed in the working directory): ")
        scans_prefix = input("Optional: Enter a prefix for the functional connectivity maps (default: enter connectivity_b_map_lh_ or connectivity_b_map_rh_ for the connectivity maps): ")
        if self.scans_type == 'REST':
            t2_dir = self.subjects_REST
        elif self.scans_type == 'EMO':
            t2_dir = self.subjects_EMO
        else:
            raise(ValueError('This is not a valid input value.'))
        working_dir = self.work_dir
        study_ID = self.study_ID
        
        # ----- Clock start-time ----- #
        start_time = time.time()
        print('\nInitiating computations ...')
        
        # ----- Read subjects data from spreadsheet ----- #
        included_subjects = np.where(np.array(self.subjects_Include) == 1)[0]
        subjects_list = np.array(self.subjects_list)[included_subjects]
        t2_dir = np.array(t2_dir)[included_subjects]
        
        # ----- Extract subjects MFD data ----- #
        mfd_fcontents = pd.read_csv(mfd_fname, sep='\t', header=0)
        mfd_subjects = np.array(mfd_fcontents['SubjID'])
        mfd_data = np.zeros(included_subjects.shape)
        for iSubject, this_subject in enumerate(mfd_subjects):
            mfd_data[np.where(subjects_list == this_subject)[0]] = mfd_fcontents['MFD'][iSubject]
        
        # ----- Obtain study inclusive grey matter mask ----- #
        this_gm_scan_ID = glob.glob(working_dir + '\\' + gm_mask)
        this_gm_mask = nb.load(this_gm_scan_ID[0]).get_fdata()
        gm_rows, gm_cols, gm_stks = np.where(this_gm_mask > 0) # voxel subscripts
        
        # ----- Pre-allocate outcome variables ----- #
        this_study_corr_matrix = np.zeros(nb.load(glob.glob(working_dir + '\\' + gm_mask)[0]).shape)
        this_study_p_matrix = np.zeros(nb.load(glob.glob(working_dir + '\\' + gm_mask)[0]).shape)
        
        # ----- Loop over all brain voxels ----- #
        perc = 0
        for iVoxel, unused in enumerate(gm_rows):
            this_voxel_coordinate = (gm_rows[iVoxel], gm_cols[iVoxel], gm_stks[iVoxel])
            
            # ----- Loop over subjects: obtain voxel connectivity vector ----- #
            this_voxel_connectivity_vector = np.zeros(len(subjects_list))
            for iSubject, this_subject in enumerate(subjects_list):
                this_subject_t2_dir = self.scans_dir + '_' + self.scans_type + '\\' + this_subject + '\\' + this_subject + t2_dir[iSubject]
                
                # ----- Find connectivity map for this subject ----- #
                this_subject_connectivity_map = glob.glob(this_subject_t2_dir + '\\' + scans_prefix + this_subject + '*.nii')
                if not this_subject_connectivity_map:
                    print('Error: Connectivity map missing for subject: ' + this_subject)
                    continue
                this_subject_connectivity_map = nb.load(this_subject_connectivity_map[0]).get_fdata()
                
                # ----- Obtain voxel connectivity value for this subject  ----- #
                this_voxel_connectivity_vector[iSubject] = this_subject_connectivity_map[this_voxel_coordinate]
                
            # ----- Correlate beta value and MFD vectors for this voxel ----- #
            this_corr_value, this_p_value = scipy.stats.spearmanr(this_voxel_connectivity_vector, mfd_data, nan_policy='omit')
            this_study_corr_matrix[this_voxel_coordinate] = this_corr_value
            this_study_p_matrix[this_voxel_coordinate] = this_p_value
            
            # ----- Print progress to console ----- #
            if not iVoxel % round(len(gm_rows) / 100):
                perc += 1
                if perc > 100:
                    continue
                print("Progress report: " + str(perc) + "%")
        
        # ----- Save QC-FC correlation and p-maps to (.nii) file ----- #
        this_study_corr_matrix_scan = nb.Nifti1Image(this_study_corr_matrix, nb.load(this_gm_scan_ID[0]).affine)
        save_fname = working_dir + '\\' + study_ID + '_QC-FC_Corr_' + scans_prefix + self.scans_type + '.nii'
        nb.save(this_study_corr_matrix_scan, save_fname)
        
        this_study_p_matrix_scan = nb.Nifti1Image(this_study_p_matrix, nb.load(this_gm_scan_ID[0]).affine)
        save_fname = working_dir + '\\' + study_ID + '_QC-FC_Prob_' + scans_prefix + self.scans_type + '.nii'
        nb.save(this_study_p_matrix_scan, save_fname)
        
        this_study_sig_matrix = np.zeros(nb.load(glob.glob(working_dir + '\\' + gm_mask)[0]).shape)
        this_study_sig_matrix[np.where(np.logical_and(this_study_p_matrix > 0, this_study_p_matrix < 0.05))] = 1
        this_study_sig_matrix_scan = nb.Nifti1Image(this_study_sig_matrix, nb.load(this_gm_scan_ID[0]).affine)
        save_fname = working_dir + '\\' + study_ID + '_QC-FC_Sig_' + scans_prefix + self.scans_type + '.nii'
        nb.save(this_study_sig_matrix_scan, save_fname)
        
        # ----- Correct for multiple comparisons via FDR ----- #
        reject, pvalscorr = multipletests(this_study_p_matrix[gm_rows, gm_cols, gm_stks], method='fdr_bh', is_sorted=False)[:2] # Benjamini & Hochberg, 1995
        
        # ----- Compute number of significant voxels and save to (.txt) file ----- #
        log_fname = self.work_dir + '\\' + study_ID + '_QC-FC_NSig_' + scans_prefix + self.scans_type + '.txt'
        if os.path.exists(log_fname):
            os.remove(log_fname)
        fh = open(log_fname, 'a')
        fh.write('Number of significant voxels (raw): ' + str(len(np.where(this_study_sig_matrix == 1)[0])) + ' of ' + str(len(gm_rows)) + "\n")
        fh.write('Number of significant voxels (FDR corrected): ' + str(len(np.where(pvalscorr < 0.05)[0])) + ' of ' + str(len(pvalscorr)) + "\n")
        fh.close()
        
        # ----- Clock end-time ----- #
        stop_time = time.time()
        print('Number of significant voxels (raw): ' + str(len(np.where(this_study_sig_matrix == 1)[0])) + ' of ' + str(len(gm_rows)))
        print('Number of significant voxels (FDR-corrected): ' + str(len(np.where(pvalscorr < 0.05)[0])) + ' of ' + str(len(pvalscorr)))
        print('Duration: ' + str(stop_time - start_time) + ' seconds')
        print('Code executed with no errors')
    
    def motion_correction_benchmark_discriminability(self):
        """
        Conduct discriminability analysis as a motion benchmark. Compute the 
        difference in functional connectivity between the high-motion and low
        motion (non-psychiatric) subjects.
        
        """
        
        # ----- Print progress to console ----- #
        print("\nConducting discriminability analysis ...")
        
        # ----- Define variables ----- #
        mfd_fname = input("Enter the filename of the study MFD data (as listed in the working directory): ")
        gm_mask = input("Enter the filename of the study grey matter mask (as listed in the working directory): ")
        scans_prefix = input("Optional: Enter a prefix for the functional connectivity maps (default: enter connectivity_b_map_lh_ or connectivity_b_map_rh_ for the connectivity maps): ")
        scans_dir = self.scans_dir + '_' + self.scans_type
        if self.scans_type == 'REST':
            t2_dir = self.subjects_REST
        elif self.scans_type == 'EMO':
            t2_dir = self.subjects_EMO
        else:
            raise(ValueError('This is not a valid input value.'))
        study_ID = self.study_ID
        working_dir = self.work_dir
        
        # ----- Clock start-time ----- #
        start_time = time.time()
        
        # ----- Extract data from MFD file ----- #
        mfd_fcontents = pd.read_csv(mfd_fname, sep='\t', header=0)
        mfd_subjects = np.array(mfd_fcontents['SubjID'])
        
        # ----- Extract subjects with an include value of one ----- #
        included_subjects = np.where(np.array(self.subjects_Include) == 1)[0]
        all_included_subjects_list = np.array(self.subjects_list)[included_subjects]
        all_included_subjects_group = np.array(self.subjects_group)[included_subjects]
        all_included_subjects_t2_dir = np.array(t2_dir)[included_subjects]
        all_included_subjects_mfd = np.zeros(included_subjects.shape)
        all_included_subjects_n_out = np.zeros(included_subjects.shape)
        for iSubject, this_subject in enumerate(mfd_subjects):
            all_included_subjects_mfd[np.where(all_included_subjects_list == this_subject)[0]] = mfd_fcontents['MFD'][iSubject]
            all_included_subjects_n_out[np.where(all_included_subjects_list == this_subject)[0]] = mfd_fcontents['nOutliers'][iSubject]
        
        # ----- Extract control participants ----- #
        control_subjects = np.where(np.array(all_included_subjects_group) == 1)[0]
        all_control_subjects_list = all_included_subjects_list[control_subjects]
        all_control_subjects_t2_dir = all_included_subjects_t2_dir[control_subjects]
        all_control_subjects_mfd = all_included_subjects_mfd[control_subjects]
        all_control_subjects_n_out = all_included_subjects_n_out[control_subjects]
        
        # ----- Split control participants according to median of MFD ----- #
        median_mfd = np.median(all_control_subjects_mfd)
        
        low_motion_subjects_list = all_control_subjects_list[np.where(all_control_subjects_mfd < median_mfd)]
        low_motion_subjects_t2_dir = all_control_subjects_t2_dir[np.where(all_control_subjects_mfd < median_mfd)]
        low_motion_subjects_mfd = all_control_subjects_mfd[np.where(all_control_subjects_mfd < median_mfd)]
        low_motion_subjects_n_out = all_control_subjects_n_out[np.where(all_control_subjects_mfd < median_mfd)]
        
        high_motion_subjects_list = all_control_subjects_list[np.where(all_control_subjects_mfd > median_mfd)]
        high_motion_subjects_t2_dir = all_control_subjects_t2_dir[np.where(all_control_subjects_mfd > median_mfd)]
        high_motion_subjects_mfd = all_control_subjects_mfd[np.where(all_control_subjects_mfd > median_mfd)]
        high_motion_subjects_n_out = all_control_subjects_n_out[np.where(all_control_subjects_mfd > median_mfd)]
        
        # ----- Print group motion characteristics to console ----- #
        print('\nDescriptives: ')
        print('Median MFD = ' + str(median_mfd))
        print('\nLow motion group (< median): ')
        print('N = ' + str(len(low_motion_subjects_list)))
        print('Mean MFD = ' + str(np.mean(low_motion_subjects_mfd)))
        print('Mean number of outliers = ' + str(np.mean(low_motion_subjects_n_out)))
        print('\nHigh motion group (> median): ')
        print('N = ' + str(len(high_motion_subjects_list)))
        print('Mean MFD = ' + str(np.mean(high_motion_subjects_mfd)))
        print('Mean number of outliers: ' + str(np.mean(high_motion_subjects_n_out)))
        
        # ----- Obtain template scan ----- #
        template_scan = glob.glob(scans_dir + '\\' + high_motion_subjects_list[0] + '\\' + high_motion_subjects_list[0] + high_motion_subjects_t2_dir[0] + '\\' + scans_prefix + high_motion_subjects_list[0] + '*.nii')[0]
        
        # ----- Pre-set to ignore deprecation warning ----- #
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        
        # ----- Pre-allocate statistical variables: sample size, sum ----- #
        n_low = len(low_motion_subjects_list) # sample size
        n_high = len(high_motion_subjects_list) # sample size
        
        S_low = np.zeros(nb.load(template_scan).shape) # sum of betas
        S_high = np.zeros(nb.load(template_scan).shape) # sum of betas
        
        SS_low = np.zeros(nb.load(template_scan).shape) # Sum of squares
        SS_high = np.zeros(nb.load(template_scan).shape) # Sum of squares
        
        # ----- Loop over subjects: low motion group ----- #
        all_low_motion_subjects_connectivity_maps = [] # Container variable
        for iSubject, this_subject in enumerate(low_motion_subjects_list):
            this_subject_t2_dir = scans_dir + '\\' + this_subject + '\\' + this_subject + low_motion_subjects_t2_dir[iSubject]
            
            # ----- Find connectivity map for this subject ----- #
            this_subject_connectivity_map = glob.glob(this_subject_t2_dir + '\\' + scans_prefix + this_subject + '*.nii')
            this_subject_connectivity_map = nb.load(this_subject_connectivity_map[0]).get_fdata()
            
            # ----- Set NaN values to zero ----- #
            this_subject_connectivity_map[np.isnan(this_subject_connectivity_map)] = 0 # Set NaN values to zero
            
            # ----- Add subject's connectivity map to container variable ----- #
            all_low_motion_subjects_connectivity_maps.append(this_subject_connectivity_map)
            
            # ----- concatenate statistical matrices ----- #
            S_low += this_subject_connectivity_map # Sum of connectivity values
        
        # ----- Loop over subjects: high motion group ----- #
        all_high_motion_subjects_connectivity_maps = [] # Container variable
        for iSubject, this_subject in enumerate(high_motion_subjects_list):
            this_subject_t2_dir = scans_dir + '\\' + this_subject + '\\' + this_subject + high_motion_subjects_t2_dir[iSubject]
            
            # ----- Find connectivity map for this subject ----- #
            this_subject_connectivity_map = glob.glob(this_subject_t2_dir + '\\' + scans_prefix + this_subject + '*.nii')
            this_subject_connectivity_map = nb.load(this_subject_connectivity_map[0]).get_fdata()
            
            # ----- Set NaN values to zero ----- #
            this_subject_connectivity_map[np.isnan(this_subject_connectivity_map)] = 0 # Set NaN values to zero
            
            # ----- Add subject beta map to container variable ----- #
            all_high_motion_subjects_connectivity_maps.append(this_subject_connectivity_map)
            
            # ----- concatenate statistical matrices ----- #
            S_high += this_subject_connectivity_map # Sum of connectivity values
        
        # ----- Loop over subjects: low motion group ----- #
        M_low = np.divide(S_low, n_low) # mean of connectivity values
        for iSubject, this_subject_connectivity_map in enumerate(all_low_motion_subjects_connectivity_maps):
            this_subject_d2_map = np.square(np.subtract(this_subject_connectivity_map, M_low)) # Deviance from mean beta
            
            SS_low += this_subject_d2_map # Sum of squared deviances
            
        S2_low = np.divide(SS_low, (n_low - 1)) # Variance of connectivity values
        
        # ----- Loop over subjects: high motion group ----- #
        M_high = np.divide(S_high, n_high) # mean of connectivity values
        for iSubject, this_subject_connectivity_map in enumerate(all_high_motion_subjects_connectivity_maps):
            this_subject_d2_map = np.square(np.subtract(this_subject_connectivity_map, M_high)) # Deviance from mean beta
            
            SS_high += this_subject_d2_map # Sum of squared deviances
        
        S2_high = np.divide(SS_high, (n_high - 1)) # Variance of connectivity values
        
        # ----- Calculate between-group t-map ----- #
        d_map = np.subtract(M_high, M_low)
        t_map = np.divide(np.subtract(M_high, M_low), np.sqrt(np.add(np.divide(S2_high, n_high), np.divide(S2_low, n_low))))
        
        # ----- Apply study grey matter mask to t-map ----- #
        this_gm_scan_ID = glob.glob(working_dir + '\\' + gm_mask)
        this_gm_mask = nb.load(this_gm_scan_ID[0]).get_fdata()
        gm_rows, gm_cols, gm_stks = np.where(this_gm_mask > 0) # voxel subscripts
        
        d_map = np.multiply(d_map, this_gm_mask)
        t_map = np.multiply(t_map, this_gm_mask)
        
        # ----- Obtain p-values ----- #
        df = len(low_motion_subjects_list) + len(high_motion_subjects_list) - 2 # n1 + n2 -2; independent samples t-test
        pvals = scipy.stats.t.sf(np.abs(t_map[gm_rows, gm_cols, gm_stks]), df) * 2 # Multiplied by two for one-sided p-values (i.e., high motion group > low motion group)
        print('\nNumber of significant voxels (raw): ' + str(len(np.where(pvals < 0.05)[0])) + ' of ' + str(len(pvals)))
        
        # ----- Correct for multiple comparisons via FDR ----- #
        reject, pvalscorr = multipletests(pvals, method='fdr_bh', is_sorted=False)[:2] # Benjamini & Hochberg, 1995
        print('Number of significant voxels (FDR-corrected): ' + str(len(np.where(pvalscorr < 0.05)[0])) + ' of ' + str(len(pvalscorr)))
         
        # ----- Write d- and t-map to (.nii) file ----- #
        d_scan = nb.Nifti1Image(d_map, nb.load(template_scan).affine)
        save_d_fname = working_dir + '\\' + study_ID + '_Discriminability_D_Map_'+ scans_prefix + self.scans_type + '.nii'
        nb.save(d_scan, save_d_fname)
        
        t_scan = nb.Nifti1Image(t_map, nb.load(template_scan).affine)
        save_t_fname = working_dir + '\\' + study_ID + '_Discriminability_T_Map_'+ scans_prefix + self.scans_type + '.nii'
        nb.save(t_scan, save_t_fname)
        
        # ----- Write results to (.txt) file ----- #
        log_fname = self.work_dir + '\\' + study_ID + '_Discriminability_NSig_' + scans_prefix + self.scans_type + '.txt'
        if os.path.exists(log_fname):
            os.remove(log_fname)
        fh = open(log_fname, 'a')
        fh.write('Descriptives: ' + '\n\n')
        fh.write('Median MFD = ' + str(median_mfd) + '\n\n')
        fh.write('Low motion group (< median): \n')
        fh.write('N = ' + str(len(low_motion_subjects_list)) + '\n')
        fh.write('Mean MFD = ' + str(np.mean(low_motion_subjects_mfd)) + '\n')
        fh.write('Mean number of outliers = ' + str(np.mean(low_motion_subjects_n_out)) + '\n\n')
        fh.write('High motion group (> median): \n')
        fh.write('N = ' + str(len(high_motion_subjects_list)) + '\n')
        fh.write('Mean MFD = ' + str(np.mean(high_motion_subjects_mfd)) + '\n')
        fh.write('Mean number of outliers: ' + str(np.mean(high_motion_subjects_n_out)) + '\n\n')
        fh.write('Number of significant voxels (raw): ' + str(len(np.where(pvals < 0.05)[0])) + ' of ' + str(len(pvals)) + '\n')
        fh.write('Number of significant voxels (FDR-corrected): ' + str(len(np.where(pvalscorr < 0.05)[0])) + ' of ' + str(len(pvalscorr)))
        fh.close()
        
        # ----- Clock end-time ----- #
        stop_time = time.time()
        print('\nDuration: ' + str(stop_time - start_time) + ' seconds')
        print('Code executed with no errors')
        
        
class GroupAnalysis(Postprocessing):
    
    def run_2nd_level_analysis(self):
        """ 
        Parametric second-level statistical analysis of the emotion task data
        (via Matlab).
        
        """
        
        # ----- Print progress to console ----- #
        print("\nRunning parametric second-level analysis ...")
        
        # ----- Define arguments (passed to Matlab) ----- #
        working_dir = self.work_dir
        study_ID = self.study_ID
        scans_dir = self.scans_dir + '_' + self.scans_type
        if self.scans_type == 'REST':
            t2_dir = self.subjects_REST
        elif self.scans_type == 'EMO':
            t2_dir = self.subjects_EMO
        else:
            raise(ValueError('This is not a valid input value.'))
        first_level_map = input("Enter the name (without file-extension) of the first-level statistical map that should be used as inputs for the analysis (e.g. con_0001 or connectivity_b_map_lh_): ")
        explicit_mask = input("Optional: Enter an explicit mask (filename as listed in the working directory, default: press enter to continue): ")
        if explicit_mask:
            try:
                mask_fn = glob.glob(working_dir + '\\' + explicit_mask)[0].split('\\')[-1]
            except:
                mask_fn = explicit_mask
        else:
            mask_fn = explicit_mask
        
        # ----- Read subjects data from spreadsheet ----- #
        included_subjects = np.where(np.array(self.subjects_Include) == 1)[0]
        subjects_list = np.array(self.subjects_list)[included_subjects]
        t2_dir = np.array(t2_dir)[included_subjects]
        
        # ----- Loop over subjects: obtain scan filenames ----- #
        scan_fnames = []
        for iSubject, this_subject in enumerate(subjects_list):
            this_subject_t2_dir = scans_dir + '\\' + this_subject + '\\' + this_subject + t2_dir[iSubject]
            
            # ----- Extract scan filename for this subject ----- #
            if self.scans_type == 'EMO':
                this_subject_contrast_map = glob.glob(this_subject_t2_dir + '\\first_level_analysis\\' + first_level_map + '*.nii')
            if self.scans_type == 'REST':
                this_subject_contrast_map = glob.glob(this_subject_t2_dir + '\\' + first_level_map + '*.nii')
            scan_fnames.append(this_subject_contrast_map[0])
        
        # ----- Write scan filenames to (.txt) file ----- #
        container_file = working_dir + '\\scan_filenames.txt'
        with open(container_file, 'w') as fh:
            for item in scan_fnames:
                fh.write("%s\n" % item)
        fh.close()
        
        # ----- Run data preprocessing script in Matlab (via shell) ----- #
        matlab_cmd_p1 = "matlab -nodesktop -nosplash -wait -r \"addpath('" + working_dir + "');"
        matlab_cmd_p2 = "start_up('" + working_dir + "'); "
        matlab_cmd_p3 = "SecondLevelAnalysis('" + working_dir + "','" + study_ID + "','" + first_level_map + "','" + container_file + "','" + mask_fn + "'); quit"
        try:
            os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
        except:
            # ----- Print failure to console ----- #
            print("Parametric second-level analysis failed for this first-level map:\t" + first_level_map)
        
        # ----- Delete (.txt) file with scan filenames ----- #
        os.remove(container_file)
    
    def voxel_matched_regression_analysis(self):
        """ 
        Voxel-matched linear regression analysis, following the procedures 
        described by Mennes, M., Kelly, C., Zuo, X. N., Di Martino, A., 
        Biswal, B., Castellanos, F. X., & Milham, M. P. (2010). Inter-individual
        differences in resting state functional connectivity predict task-induced 
        BOLD activity. Neuroimage, 50(4), 1690–1701.
        
        Landscape-based cluster permutation analysis is conducted on the data 
        using the procedures described by Gladwin, T. E., Vink, M., & Mars, R. 
        B. (2016). A landscape-based cluster analysis using recursive search 
        instead of a threshold parameter. MethodsX, 3, 477-482.
        
        """
        
        # ----- Define arguments (passed to Matlab) ----- #
        workdir = self.work_dir
        subjects_fn = input("Enter the name of the text-file that contains the subjects list (e.g. voxel_matched_regression_analysis.txt): ")
        rs_dir = self.data_dir + '\\NIFTI_GRAND_REST'
        emo_dir = self.data_dir + '\\NIFTI_GRAND_EMO'
        connectivity_prefix = input("Enter a prefix for the functional connectivity maps that should be used as inputs for the analysis (e.g. connectivity_b_map_lh_ or connectivity_b_map_rh_): ")
        reactivity_fn = input("Enter the name of the contrast statistical maps that should be used as inputs for the analysis (e.g. con_0001.nii, con_0002.nii, con_0003.nii, con_0004.nii, or con_0005.nii): ")
        mask_fn = input("Enter the filename of a mask to extract the voxel indices from (as listed in the working directory, default: GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii): ")
        fwhm = int(input("Enter the full width half maximum (FWHM) at which you want to smooth the data (default: enter 6): "))
        n_perms = int(input("Enter the number of permutations you want to conduct (default: 2000): "))
        n_tests = int(input("Enter the number of tests that are conducted in total; i.e., the FWER-correction factor (e.g., 30): "))
        
        # ----- Clock start-time ----- #
        start_time = time.time()
        
        # ----- Print progress to console ----- #
        print("\nConducting voxel-matched regression with landscape-based permutation cluster analysis ...\n")
        
        # ----- Run data preprocessing script in Matlab (via shell) ----- #
        matlab_cmd_p1 = "matlab -nosplash -wait -r \"addpath('" + workdir + "');"
        matlab_cmd_p2 = "start_up('" + workdir + "'); "
        matlab_cmd_p3 = "VoxelMatchedRegressionAnalysis('" + workdir + "','" + subjects_fn + "','" + rs_dir + "','" + emo_dir + "','" + connectivity_prefix + "','" + reactivity_fn + "','" + mask_fn + "'," + str(fwhm) + "," + str(n_perms) + "," + str(n_tests) + "); quit"
        try:
            os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
        except:
            # ----- Print failure to console ----- #
            print("voxel-matched regression with landscape-based cluster permutation analysis failed for this contrast:\t" + reactivity_fn)
        
        # ----- Clock end-time ----- #
        stop_time = time.time()
        print('Duration: ' + str(stop_time - start_time) + ' seconds')
        print('Done!')
    
    def roi_based_regression_analysis(self):
        """ 
        Voxel-by-voxel ROI-based linear regression analysis, loosely
        following the procedures described by Mennes, M., Kelly, C., 
        Zuo, X. N., Di Martino, A., Biswal, B., Castellanos, F. X., & 
        Milham, M. P. (2010). Inter-individual differences in resting state
        functional connectivity predict task-induced BOLD activity. 
        Neuroimage, 50(4), 1690–1701.
        
        Landscape-based cluster permutation analysis is conducted on the data 
        using the procedures described by Gladwin, T. E., Vink, M., & Mars, R. 
        B. (2016). A landscape-based cluster analysis using recursive search 
        instead of a threshold parameter. MethodsX, 3, 477-482.
        
        """
        
        # ----- Define arguments (passed to Matlab) ----- #
        workdir = self.work_dir
        rs_dir = self.data_dir + '\\NIFTI_GRAND_REST'
        connectivity_prefix = input("Enter a prefix for the functional connectivity maps that should be used as inputs for the analysis (e.g. connectivity_b_map_lh_ or connectivity_b_map_rh_): ")
        reactivity_fn = input("Enter the name of the contrast statistical maps that should be used as inputs for the analysis (e.g. con_0001.nii, con_0002.nii, con_0003.nii, con_0004.nii, or con_0005.nii): ")
        reactivity_hemi = int(input("Specify the hemisphere to extract the ROI reactivity data from. Enter 1 for the left hemisphere, or 2 for the right hemisphere: "))
        if (reactivity_hemi != 1) and (reactivity_hemi != 2):
            raise ValueError('This is not a valid input value for reactivity_hemi. Enter 1 for the left hemisphere, or 2 for the right hemisphere.')
        mask_fn = input("Enter the filename of a mask to extract the voxel indices from (as listed in the working directory, default: GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii): ")
        fwhm = int(input("Enter the full width half maximum (FWHM) at which you want to smooth the data (default: enter 6): "))
        n_perms = int(input("Enter the number of permutations you want to conduct (default: 2000): "))
        n_tests = int(input("Enter the number of tests that are conducted in total; i.e., the FWER-correction factor (e.g., 30): "))
        
        # ----- Clock start-time ----- #
        start_time = time.time()
        
        # ----- Print progress to console ----- #
        print("\nConducting ROI-based regression with landscape-based permutation cluster analysis ...\n")
        
        # ----- Run data preprocessing script in Matlab (via shell) ----- #
        matlab_cmd_p1 = "matlab -nosplash -wait -r \"addpath('" + workdir + "');"
        matlab_cmd_p2 = "start_up('" + workdir + "'); "
        matlab_cmd_p3 = "ROIBasedRegressionAnalysis('" + workdir + "','" + rs_dir + "','" + connectivity_prefix + "','" + reactivity_fn + "','" + str(reactivity_hemi) + "','" + mask_fn + "'," + str(fwhm) + "," + str(n_perms) + "," + str(n_tests) + "); quit"
        try:
            os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
        except:
            # ----- Print failure to console ----- #
            print("ROI-based regression with landscape-based cluster permutation analysis failed for this contrast:\t" + reactivity_fn)
        
        # ----- Clock end-time ----- #
        stop_time = time.time()
        print('Duration: ' + str(stop_time - start_time) + ' seconds')
        print('Done!')
    
    def run_snpm_specification_and_computation(self):
        """ 
        Second-level statistical analysis of resting-state vs. emotion task 
        data: specification and computation steps (via Matlab).
        
        """
        
        # ----- Print progress to console ----- #
        print("\nRunning non-parametric second-level analysis ...")
        
        # ----- Define arguments (passed to Matlab) ----- #
        working_dir = self.work_dir
        output_dir = working_dir + '\\' + input("Enter a name of the directory where the output files will be written to (in the working directory): ")
        hemi = input("What hemisphere do you want to run the analysis on? Enter l (left) or r (right): ")
        data_file = input("Enter the name of a text file that contains the predictor-of-interest (and covariate) data, as well as the connectivity map filenames (note that this file needs to be manually prepared beforehand) (e.g. SnPM_Input_SpmT_0001_HemiL.txt): ")
        if (hemi != 'l') and (hemi != 'r'):
            raise(ValueError('This is not a valid input value, enter l (left) or r (right).'))
        if hemi == 'l':
            scans_dir = self.data_dir + '\\SnPM_HemiL'
        elif hemi == 'r':
            scans_dir = self.data_dir + '\\SnPM_HemiR'
        if not os.path.exists(scans_dir):
            raise(FileNotFoundError("This directory does not exist."))
        include_covariate = input("Optional: Enter the name of the covariate you want to include (default: press enter to continue): ")
        n_perms = input("Enter the number of permutations that you want to conduct (e.g. 10000): ")
        explicit_mask = input("Optional: Enter an explicit mask (filename as listed in the working directory, default: GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii): ")
        
        # ----- Run data 2nd analusis script in Matlab (via shell) ----- #
        matlab_cmd_p1 = "matlab -nosplash -wait -r \"addpath('" + working_dir + "');"
        matlab_cmd_p2 = "start_up('" + working_dir + "'); "
        matlab_cmd_p3 = "SnPMSpecificationAndComputation('" + working_dir + "','" + output_dir + "','" + scans_dir + "','" + data_file + "','" + include_covariate + "'," + n_perms + ",'" + explicit_mask + "'); quit"
        try:
            os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
        except:
            # ----- Print failure to console ----- #
            print("Second level analysis failed")
        
        # ----- Return to working directory ----- #
        os.chdir(working_dir)
    
    def run_snpm_inference(self):
        """ 
        Second-level statistical analysis of resting-state vs. emotion task 
        data: inference step (via Matlab).
        
        """
        
        # ----- Print progress to console ----- #
        print("\nRunning non-parametric second-level inference ...")
        
        # ----- Define arguments (passed to Matlab) ----- #
        working_dir = self.work_dir
        analysis_dir = working_dir + '\\' + input("Enter the output directory of the SnPM analysis (as listed in the working directory): ")
        if not os.path.exists(analysis_dir):
            raise FileNotFoundError("This directory does not exist.")
        analysis_type = int(input("What type of analysis do you want to perform? Enter 1 for voxel-wise uncorrected, 2 for voxel-wise FDR-corrected, 3 for voxel-wise FWE-corrected, or 4 for cluster-wise FWE-corrected: "))
        if (analysis_type != 1) and (analysis_type != 2) and (analysis_type != 3) and (analysis_type != 4):
            raise ValueError('This is not a valid input value: Enter 1 (voxel-wise, uncorrected), 2 (voxel-wise, FDR-corrected), 3 (voxel-wise, FWE-corrected), or 4 (cluster-wise, FWE-corrected).')
        effect_sign = int(input("What type of effects do you want to analyze? Enter 1 for positive effects, or -1 for negative effects: "))
        if (effect_sign != 1) and (effect_sign != -1):
            raise ValueError('This is not a valid input value: Enter 1 (positive effects) or -1 (negative effects).')
        
        # ----- Run data 2nd analusis script in Matlab (via shell) ----- #
        matlab_cmd_p1 = "matlab -nodesktop -nosplash -wait -r \"addpath('" + working_dir + "');"
        matlab_cmd_p2 = "start_up('" + working_dir + "'); "
        matlab_cmd_p3 = "SnPMInference('" + analysis_dir + "'," + str(analysis_type) + "," + str(effect_sign) + "); quit"
        try:
            os.system(matlab_cmd_p1 + matlab_cmd_p2 + matlab_cmd_p3)
        except:
            # ----- Print failure to console ----- #
            print("Second level analysis failed")
        
        # ----- Return to working directory ----- #
        os.chdir(working_dir)


# ----- Run code from main method ----- #
if __name__ == "__main__":
    # ======================================================================= #
    # ANALYSIS PIPELINE OF THE RESTING-STATE DATA: REPEAT FOR MARS AND BETER  #
    # ======================================================================= #
    # ----- Preperation of the resting-state data ----- #
    # my_experiment = RestingState()
    # ----- Preprocessing of the resting-state data ----- #
    # my_experiment.run_preprocessing() # Realignment (inputs: none)
    # my_experiment.extract_FD_jenkinson() # FD Jenkinson (inputs: 1)
    # my_experiment.run_preprocessing() # Coregistration (inputs: r)
    # my_experiment.run_preprocessing() # Normalization (inputs: r)
    # my_experiment.extract_inclusive_FOV_mask() # Compute FOV mask (inputs: nr)
    # my_experiment.run_preprocessing() # Segmentation (inputs: n)
    # my_experiment.erode_segmentation_mask() # Erosion of white-matter segmentation (inputs: c2)
    # my_experiment.erode_segmentation_mask() # Erosion of CSF segmentation (inputs: c3)
    # ----- 1st-level (connectivity) analysis of the resting state data ----- #
    # my_experiment.obtain_global_roi_mask() # Copy and reslice ROI mask to subject reference space (inputs: n)
    # my_experiment.extract_roi_regressors() # Extract ROI regressors (inputs: rAmygdala_total_binary_mask_thr_0_5.nii, 1, nr)
    # my_experiment.extract_confound_regressors() # Extract confound regressors (inputs: none, e, nr)
    # my_experiment.extract_spike_regressors() # Extract spike regressors (inputs: 0.2)
    # my_experiment.voxel_wise_connectivity_analysis() # Run voxel-wise connectivity analysis (inputs: nr)
    # ======================================================================= #
    # ANALYSIS PIPELINE OF THE EMOTION TASK DATA: REPEAT FOR MARS AND BETER   #
    # ======================================================================= #
    # ----- Preperation of the emotion task data ----- #
    # my_experiment = EmotionTask()
    # ----- Preprocessing of the emotion task data ----- #
    # my_experiment.run_preprocessing() # Slice-timing correction (inputs: none)
    # my_experiment.run_preprocessing() # Realignment (inputs: r)
    # my_experiment.extract_FD_jenkinson() # FD Jenkinson (inputs: 1)
    # my_experiment.run_preprocessing() # Coregistration (inputs: ra)
    # my_experiment.run_preprocessing() # Normalization (inputs: ra)
    # my_experiment.extract_inclusive_FOV_mask() # Compute FOV mask (inputs: nra)
    # my_experiment.run_preprocessing() # Segmentation (inputs: n)
    # my_experiment.erode_segmentation_mask() # Erosion of white-matter segmentation (inputs: c2)
    # my_experiment.erode_segmentation_mask() # Erosion of CSF segmentation (inputs: c3)
    # ----- First-level analysis of the emotion task data ----- #
    # my_experiment.extract_behavioral_data() # Extraction of behavioral task data
    # my_experiment.obtain_global_roi_mask() # Copy and reslice ROI mask to subject reference space (inputs: n)
    # my_experiment.extract_confound_regressors() # Extract confound regressors (inputs: none, e, nra)
    # my_experiment.extract_spike_regressors() # Extract spike regressors (inputs: 0.5)
    # my_experiment.run_1st_level_analysis() # First-level analysis (inputs: none, nra)
    # my_experiment.extract_roi_masked_spm_data() # Extract ROI-masked SPM data (inputs: rAmygdala_total_binary_mask_thr_0_5.nii, 1, spmT_000X.nii/con_000X.nii)
    # ======================================================================= #
    #  POST-PROCESSING PIPELINE: REPEAT FOR MARS AND BETER                    #
    # ======================================================================= #
    # my_experiment = Postprocessing()
    # ----- Compute study-specific masks ----- #
    # my_experiment.create_study_FOV_mask() # Compute study FOV mask (inputs: FOV_mask_)
    # my_experiment.create_study_grey_matter_mask() # Compute study grey matter mask (inputs: c1n)
    # my_experiment.create_average_brain() # Create an average brain (inputs: n)
    # ----- Motion correction benchmarks ----- #
    # my_experiment.motion_correction_benchmark_qc_fc() # Compute QC-FC motion correction benchmarks (inputs: connectivity_b_map_lh_)
    # my_experiment.motion_correction_benchmark_qc_fc() # Compute QC-FC motion correction benchmarks (inputs: connectivity_b_map_rh_)
    # my_experiment.motion_correction_benchmark_discriminability() # Compute discriminability motion correction benchmark (inputs: connectivity_b_map_lh_)
    # my_experiment.motion_correction_benchmark_discriminability() # Compute discriminability motion correction benchmark (inputs: connectivity_b_map_rh_)
    # ======================================================================= #
    # SECOND-LEVEL ANALYSIS PIPELINE: CONDUCT EACH STEP ONLY ONCE             #
    # ======================================================================= #
    my_experiment = GroupAnalysis()
    # ----- Second-level analysis: resting-state data ----- #
    my_experiment.run_2nd_level_analysis() # Second-level analysis (inputs: connectivity_b_map_lh)
    # my_experiment.run_2nd_level_analysis() # Second-level analysis (inputs: connectivity_b_map_rh)
    # ----- Second-level analysis: emotion task data ----- #
    # my_experiment.run_2nd_level_analysis() # Second-level analysis (inputs: con_0001)
    # my_experiment.run_2nd_level_analysis() # Second-level analysis (inputs: con_0002)
    # my_experiment.run_2nd_level_analysis() # Second-level analysis (inputs: con_0003)
    # my_experiment.run_2nd_level_analysis() # Second-level analysis (inputs: con_0004)
    # my_experiment.run_2nd_level_analysis() # Second-level analysis (inputs: con_0005)
    # ----- ROI-based regression analysis: resting-state vs. emotion task data ----- #
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0001.nii, 1, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0001.nii, 2, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0001.nii, 1, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0001.nii, 2, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0002.nii, 1, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0002.nii, 2, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0002.nii, 1, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0002.nii, 2, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0003.nii, 1, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0003.nii, 2, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0003.nii, 1, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0003.nii, 2, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0004.nii, 1, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0004.nii, 2, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0004.nii, 1, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0004.nii, 2, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0005.nii, 1, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0005.nii, 2, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0005.nii, 1, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.roi_based_regression_analysis() # ROI-based connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0005.nii, 2, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # ----- Voxel-matched regression analysis: resting-state vs. emotion task data ----- #
    # my_experiment.voxel_matched_regression_analysis() # Voxel-matched connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0001.nii, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.voxel_matched_regression_analysis() # Voxel-matched connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0002.nii, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.voxel_matched_regression_analysis() # Voxel-matched connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0003.nii, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.voxel_matched_regression_analysis() # Voxel-matched connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0004.nii, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.voxel_matched_regression_analysis() # Voxel-matched connectivity vs. reactivity analysis (inputs: connectivity_t_map_lh_, spmT_0005.nii, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.voxel_matched_regression_analysis() # Voxel-matched connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0001.nii, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.voxel_matched_regression_analysis() # Voxel-matched connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0002.nii, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.voxel_matched_regression_analysis() # Voxel-matched connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0003.nii, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.voxel_matched_regression_analysis() # Voxel-matched connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0004.nii, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    # my_experiment.voxel_matched_regression_analysis() # Voxel-matched connectivity vs. reactivity analysis (inputs: connectivity_t_map_rh_, spmT_0005.nii, GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii, 6, 2000)
    
