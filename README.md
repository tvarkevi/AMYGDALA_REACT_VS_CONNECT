# AMYGDALA_REACT_VS_CONNECT
Pipeline to predict emotion task reactivity of the amygdala using resting-state connectivity of the amygdala as seed-region.

> Preliminary note: The Python module [amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py) used here assumes the scripts and analysis files are stored in one and the same working directory. The emotion task (EMO) and resting-state (REST) fMRI data needs to have been converted to NIFTI format, and stored in a data directory (data_dir) organized such that each combination of dataset (MARS, BETER) and scan type (EMO, REST) has a different subdirectory (i.e., NIFTI_MARS_REST, NIFTI_MARS_EMO, NIFTI_BETER_REST, and NIFTI_BETER_EMO) within the data directory (e.g. data_dir > NIFTI_MARS_REST). Within each of these subdirectories, each subject has to have its own folder, with the name of that folder corresponding to the subject identifier (e.g. xm13101101). The different types of scans (T1, REST, EMO) for each subject need to be stored as subdirectories within that subject directory, and the identifiers of these subdirectories need to be indicated by suffixes that are added to the main subject identifier (e.g. data_dir > NIFTI_BETER_EMO > be1a031010 > **be1a031010_5_1**). The first two letters of each subject identifier specify the dataset whereto the subject belongs; 'xm' indicates the MARS dataset (e.g. xm13101101), and 'be' indicates the BETER dataset (e.g. be1a031010). The base directory of the data also needs to contain, for each dataset (MARS, BETER) seperately, a directory (LOGS_MARS, LOGS_BETER) that contains the log files of the emotion task (e.g. data_dir > LOGS_BETER > be1a031010-3amygdala_hippo_aug2010.log). These log files always need to end with **aug2010.log**, otherwise the program will not find them. Finally, an Excel (.xlsx) file needs to be made that specifies, for each subject, the subject identifier (e.g. xm13101101), the subject number (e.g. MARS011), group membership (e.g. 1 or 0) and the specific suffixes added to the main subject ID in order to differentiate the different scan subdirectories for that subject (e.g. \_4_1, \_5_1, etc.).

## 1. Setting up the experiment

The Python module [amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py) that is used to run the analyses described here consists of multiple classes. The parent class that is inherited by all the other subclasses is called Experiment. This class has a (constructor) method called \__init__ that constructs all the attributes needed by (and specific to) the experiment. Three user inputs need to be specified in the console:
1. The full working directory that contains all the scripts and auxilliary files needed for the analysis.
2. The full data directory containing each subject identifier as a seperate folder, in which each of the relevant scans of the subjects are specified by subdirectories in the following manner: data_dir > NIFTI_BETER_REST > be1a031010 > be1a031010_5_1 (see also the preliminary note at the beginning of this user manual).
3. The full directory of the Excel file that contains the subject information. The first row in this file needs to contain at least the following column headers: SubjID (or SubjT0ID, containing the subject identifiers), PPN (containing the subject number), Group (containing the group memberships), T1 (containing the anatomical scan directory suffixes for each subject), REST (containing the resting-state directory suffixes for each subject), and EMO (containing the task fMRI directory suffixes for each subject). Other columns containing e.g. questionnaire data can also be specified. Each row in the Excel file specifies a single subject, with the first two letters specifying the dataset that subject belongs to; 'xm' for the MARS dataset, and 'be' for the BETER dataset (see also the preliminary note at the beginning of this user manual)

To initialize the Experiment class along with all its attributes, define the working directory and run the following code in the console:

```
working_dir = r'O:\Project directory\Analysis'

import os
os.chdir(working_dir)

import amygdala_project as Amy
my_experiment = Amy.Experiment()
```
This code generates an output argument called *my_experiment* which contains attributes that define the parameters of the experiment (e.g. *my_experiment.subjects_list* contains a list of all the subject identifiers).

> Note: Since the Experiment class is inherited by all the other subclasses, it does not need to be defined seperately for the analysis to be conducted; it is called automatically when initializing the subclasses later on. Hence, initializing Experiment before initializing the other classes is not strictly necessary.

## 2. Data preprocessing

Many of the data preprocessing methods in the [amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py) module are conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command.

The following preprocessing steps are supported by the pipeline:
1. Realignment of the functional scans (section 2.1)
2. Extraction of an inclusive mask of the functional images (section 2.2)
3. Calculation of framewise displacements (section 2.3)
4. Slice-timing correction of the functional scans (section 2.4)
5. Coregistration of the anatomical scan to the mean functional image (section 2.5)
6. Segmentation of the anatomical image (section 2.6)
7. Erosion of the segmentation volumes (section 2.7)
8. Normalization of the beta (connectivity) maps (or other functional images; section 2.8)
9. Smoothing of the normalized functional images (section 2.9)
10. Filtering of the timeseries data (section 2.10)

### 2.1 Realignment

Realignment of the functional data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the realignment process are [Realignment.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Realignment.m) and [Realignment_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Realignment_job.m)

To conduct the realignment process, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since this class inherits from the main Experiment class, the same three user inputs as described earlier (see section 1) need to be entered. Furthermore, the program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data .
2. The specific preprocessing step that needs to be conducted. Enter 1 for realignment.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. This should be skipped at this stage. Simply press the enter key.

> Note, the above block of code assumes that the [amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py) module has already been imported. If this is not the case, run the following code beforehand:

```
working_dir = r'O:\Project directory\Analysis'

import os
os.chdir(working_dir)
```

### 2.2 Inclusive mask

Extraction of a subject-specific inclusive voxel mask is supported by module via the *extract_inclusive_FOV_mask* method of the Preprocessing class. If inclusive masks should be conducted, enter the following code in the console:

```
my_experiment = Amy.Preprocessing()
my_experiment.extract_inclusive_FOV_mask()
```

Since this class inherits from the main Experiment class, the same three user inputs as described earlier (see section 1) need to be entered. Furthermore, the program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data .
2. An optional prefix to indicate the exact scan identifiers to base the extraction on. Enter r for the realigned scans.

### 2.3 Slice-timing correction

Slice-timing correction of the functional data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the slice-timing correction are [SliceTimingCorrection.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/SliceTimingCorrection.m) and [SliceTimingCorrection_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/SliceTimingCorrection_job.m)

To conduct the slice-timing correction, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since this class inherits from the main Experiment class, the same three user inputs as described earlier (see section 1) need to be entered. Furthermore, the program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data .
2. The specific preprocessing step that needs to be conducted. Enter 2 for slice-timing correction.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. Slice-timing correction needs to be performed on the realigned functional scans. Enter r for the realigned scans.


