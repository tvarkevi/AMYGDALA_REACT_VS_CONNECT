# AMYGDALA_REACT_VS_CONNECT
Pipeline to predict emotion task reactivity of the amygdala using resting-state connectivity of the amygdala as seed-region.

> Preliminary note: The Python module [amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py) used here assumes the following steps having been taken before running the analyses:
> 1. Many of the procedures described in this manual make use of [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. Framewise displacement in this pipeline is calculated using the [Jenkinson](https://pdfs.semanticscholar.org/291c/9d66cf886ed27d89dc03d42a4612b66ab007.pdf) method using code of the [DPARSF](http://rfmri.org/DPARSF) software package. The bandpass filter used in many of the methods described here was originally developed by Thomas Gladwin; the ezfilt.m script can be found at: https://www.tegladwin.com/code.php. Each of these programs have to be installed in order for this pipeline to work properly.
> 2. The necessary analysis scripts and auxilliary files need to stored in one and the same working directory. The emotion task (EMO) and resting-state (REST) fMRI data needs to have been converted to NIFTI format, and stored in a (seperate) data directory (data_dir) organized such that each combination of dataset (MARS, BETER) and scan type (EMO, REST) has a different subdirectory (i.e., NIFTI_MARS_REST, NIFTI_MARS_EMO, NIFTI_BETER_REST, and NIFTI_BETER_EMO) within that data directory (e.g. data_dir > NIFTI_MARS_REST). 
> 3. Within each of the NIFTI data subdirectories, each subject has to have its own subfolder, with the name of that subfolder corresponding to the identifier of that subject (e.g. xm13101101). The different types of scans (T1, REST, EMO) for each subject need to be stored as subdirectories within that subject directory, and the identifiers of these subdirectories need to be indicated by suffixes that are added to the main subject identifier (e.g. data_dir > NIFTI_BETER_EMO > be1a031010 > **be1a031010_5_1**). 
> 4. The first two letters of each subject identifier specify the dataset where to the subject belongs; 'xm' indicates the MARS dataset (e.g. xm13101101), and 'be' indicates the BETER dataset (e.g. be1a031010). 
> 5. The base directory of the data needs to contain, for each dataset (MARS, BETER) seperately, a directory (LOGS_MARS, LOGS_BETER) that contains the log files of the emotion task (e.g. data_dir > LOGS_BETER > be1a031010-3amygdala_hippo_aug2010.log). These log files need to end with **aug2010.log**, otherwise the program will not be able to identify them. 
> 6. An Excel (.xlsx) file needs to be present in the working directory which specifies, for each subject, the subject identifier (e.g. xm13101101), the subject number (e.g. MARS011), group membership (e.g. 1 or 0) and the specific suffixes added to the main subject ID in order to differentiate the different scan subdirectories for that subject (e.g. \_4_1, \_5_1, etc.).

Table of contents
* 1. [Setting up the experiment](#1. Setting up the experiment)

## 1. Setting up the experiment <a name="Setting up the experiment"></a>

The Python module [amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py) that is used to run the analyses described here consists of multiple classes. The parent class that is inherited by all the other subclasses is called Experiment. This class has a (constructor) method called \_\_init__ that constructs all the attributes needed by (and specific for) the experiment in question. Three user inputs need to be specified in the console when initializing the Experiment class:
1. The full working directory that contains all the scripts and auxilliary files needed for the analysis.
2. The full base directory of the data, containing each subject identifier as a seperate folder, in which each of the relevant scans of the subjects are specified by subdirectories in the following manner: data_dir > NIFTI_BETER_REST > be1a031010 > be1a031010_5_1 (see also point 3 of the preliminary note at the beginning of this manual).
3. The full directory of the Excel file that contains the subject information. The first row in this file needs to contain at least the following column headers: SubjID (or SubjT0ID, containing the subject identifiers), PPN (containing the subject number), Group (containing the group memberships), T1 (containing the anatomical scan directory suffixes for each subject), REST (containing the resting-state directory suffixes for each subject), and EMO (containing the task fMRI directory suffixes for each subject). Other columns containing e.g. questionnaire data can also be specified. Each row in the Excel file specifies a single subject, with the first two letters specifying the dataset that subject belongs to; 'xm' for the MARS dataset, and 'be' for the BETER dataset (see also point 6 of the preliminary note at the beginning of this user manual)

To initialize the Experiment class and all its attributes, define the working directory and run the following code in the console:

```
working_dir = r'O:\Project directory\Analysis'

import os
os.chdir(working_dir)

import amygdala_project as Amy
my_experiment = Amy.Experiment()
```
This code generates an output argument called *my_experiment* which contains attributes that define the parameters of the experiment (e.g. *my_experiment.subjects_list* contains a list of all the subject identifiers).

> Note: Since the Experiment class is inherited by all the other subclasses, it is not strictly necessary to initialize this method on its own, before conducting any of the other procedures described here. That is, it is called automatically when initializing the other subclasses of the pipeline.

## 2. Data preprocessing

Many of the data preprocessing methods in the [amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py) module are conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command.

The following preprocessing steps are supported by the pipeline:
1. Realignment of the functional data (section 2.1)
2. Extraction of an inclusive mask, based on the functional images (section 2.2)
3. Calculation of framewise displacements (section 2.3)
4. Slice-timing correction of the functional data (section 2.4)
5. Coregistration of the anatomical scan to the mean functional image (section 2.5)
6. Segmentation of the anatomical image (section 2.6)
7. Erosion of the segmentation data (section 2.7)
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

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The specific preprocessing step that needs to be conducted. Enter 1 for realignment.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. This option should be skipped at this stage; simply press the enter key to continue.

> Note, the above block of code assumes that the [amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py) module has already been imported. If this is not the case, run the following code beforehand in the console:

```
working_dir = r'O:\Project directory\Analysis'

import os
os.chdir(working_dir)

import amygdala_project as Amy
```

### 2.2 Inclusive mask extraction

Extraction of a subject-specific inclusive voxel mask is supported by the [amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py) module via the *extract_inclusive_FOV_mask* method of the Preprocessing class. Enter the following code in the console to execute this process:

```
my_experiment = Amy.Preprocessing()
my_experiment.extract_inclusive_FOV_mask()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. An optional prefix to indicate the exact scan identifiers to base the extraction of the inclusive mask on. Enter r for the realigned functional images.
3. The last input asks whether an across-subjects inclusive mask needs to be constructed. Enter 1 for yes or 0 for no.

### 2.3 Framewise displacement

The extraction of framewise displacement (FD) data is conducted using [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB function that is responsible for the extraction process is called [FrameWiseDisplacement.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/FrameWiseDisplacement.m), and makes use of the y_FD_Jenkinson.m script of the [DPARSF](http://rfmri.org/DPARSF) software package. The DPARSF function calculates the framewise displacement data using the [Jenkinson](https://pdfs.semanticscholar.org/291c/9d66cf886ed27d89dc03d42a4612b66ab007.pdf) method, based on the motion parameter estimates generated by the realignment procedure in [SPM12]((https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)) (see section 2.1). 

Enter the following code in the console to extract the FD Jenkinson data:

```
my_experiment = Amy.Preprocessing()
my_experiment.extract_FD_jenkinson()
```

> Note: It is essential that both [SPM12]((https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and [DPARSF](http://rfmri.org/DPARSF) are added to the MATLAB path in order for the *extract_FD_jenkinson* method to work properly.

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The reference scan that is used to calculate the framewise displacement data. Enter 1 for the first image, or M for the mean image.

### 2.4 Slice-timing correction

Slice-timing correction of the functional data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are used for the slice-timing correction are [SliceTimingCorrection.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/SliceTimingCorrection.m) and [SliceTimingCorrection_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/SliceTimingCorrection_job.m)

To conduct the slice-timing correction, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The specific preprocessing step that needs to be conducted. Enter 2 for slice-timing correction.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. Slice-timing correction needs to be performed on the realigned functional scans. Enter r for the realigned scans.

### 2.5 Coregistration

Coregistration of the raw anatomical image to the mean functional image is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the coregistration are [Coregistration.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Coregistration.m) and [Coregistration_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Coregistration_job.m)

To conduct the coregistration procedure, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The specific preprocessing step that needs to be conducted. Enter 3 for coregistration.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. The raw anatomical volume needs to be used as source image at this stage. Simply press the enter key to continue.

### 2.6 Segmentation 

Segmentation of the anatomical data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the segmentation are [Segmentation.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Segmentation.m) and [Segmentation_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Segmentation_job.m)

To conduct the segmentation procedure, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The specific preprocessing step that needs to be conducted. Enter 4 for segmentation.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. The segmentation can be conducted on either the raw or coregistered anatomical data. Simply press the enter key to conduct the segmentation on the raw T1 image, or enter c for the coregistered anatomical image.

### 2.7 Erosion 

Erosion of the segmentation data is supported by the [amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py) module, via the *erode_segmentation_mask* method of the Preprocessing class. This method uses the MATLAB function [Erosion.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Erosion.m) via a shell call command. 

Enter the following code in the console to execute this process:

```
my_experiment = Amy.Preprocessing()
my_experiment.erode_segmentation_mask()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. An optional prefix to indicate the exact scan identifiers on which the erosion needs to be applied. In theory, the erosion procedure can be conducted on all outputs of the segmentation procedure (c1, c2, c3, c4, c5). To perform the erosion on the white-matter or CSF segmentations, enter c2 or c3, respectively.

### 2.8 Normalization

Normalization of the functional data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the normalization procedure are [Normalization.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Normalization.m) and [Normalization_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Normalization_job.m)

To conduct the normalization procedure, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The specific preprocessing step that needs to be conducted. Enter 5 for normalization.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. In theory, the normalization can be performed on any series of functional images. Enter b_map_ for the raw beta (connectivity) maps.

### 2.9 Smoothing

Smoothing of the functional data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the smoothing procedure are [Smoothing.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Smoothing.m) and [Smoothing_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Smoothing_job.m)

To conduct the smoothing procedure, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The specific preprocessing step that needs to be conducted. Enter 6 for smoothing.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. In theory, the smoothihng procedure can be applied on any series of functional images. Enter nb_map_ for the normalized beta (connectivity) maps.

### 2.10 Filtering

Filtering of the timeseries data is conducted via the *run_filtering* method of the Preprocessing class. This method utilizes a filtering procedure *ezfilt.m* developed by Thomas Gladwin (https://www.tegladwin.com/), rewritten into this pipeline. Note that this filtering method is integrated in many of the more advanced procedures described here, and does not need to be initialized seperately. Nevertheless, it *is* possible to call the method manually from the console. For instance, to bandpass filter a given signal (your_signal) at 0.01-0.08 Hz, at a given repetition time (your_TR), enter the following code:

```
TR = your_TR

my_experiment = Amy.Preprocessing()

unused, your_signal = my_experiment.run_filtering(your_signal, 1/TR, 0.01)
your_signal, unused = my_experiment.run_filtering(your_signal, 1/TR, 0.08)
```

## 3. Emotion task

The emotion task described here is based on the paradigm of [Van Buuren et al. (2011)](https://www.sciencedirect.com/science/article/pii/S000632231100254X) (see also [Heesink et al., 2018](https://www.cambridge.org/core/journals/european-psychiatry/article/neural-activity-during-the-viewing-of-emotional-pictures-in-veterans-with-pathological-anger-and-aggression/C3179D960B2B1DA0C0C33D34B5FCA2D6)).

The analysis of the emotion task data is conducted via the EmotionTask subclass. This class inherits all attributes and methods of the Preprocessing class, which itself (in turn) inherits from the main Experiment class. The first-level analysis of the emotion task data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via shell call commands.

The following preprocessing pipeline is recommended for the emotion task data:
1. Realignment of the raw functional images (see section 2.1)
2. Extraction of the FD Jenkinson data (see section 2.3)
3. Slice-timing correction of the realigned functional images (see section 2.4)
4. Coregistration of the anatomical image to the mean functional image (see section 2.5)
5. Extraction of a (subject-level) inclusve field-of-view (FOV) voxel mask from the realigned functional images (see section 2.2)
6. Segmentation of the coregistered anatomical image (see section 2.6)
7. Erosion of the white-matter and CSF segmentations (see section 2.7)

As part of the first-level analysis pipeline of the emotion task data, the following steps need to be performed after the preprocessing steps:
1. Extraction of the behavioral task data from the logfiles (see section 3.1)
2. Extraction of the confound regressors (see section 3.2)
3. Extraction of the spike regressors (see section 3.3)
4. First-level analysis in [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (see section 3.4)
5. Construction of a (subject-level) region-of-interest (ROI) mask in native space (see section 3.5)
6. Extraction (across subjects) of the ROI-masked statistical paramatric data (see section 3.6)

### 3.1 Behavioral data

Before the first-level analysis of the emotion task data can be conducted, the task data needs to be extracted from the raw logfiles. As described in the preliminary note of this manual, the raw logfiles need to be copied to the data directory (data_dir), in a subfolder named after the dataset in question; i.e., data_dir > LOGS_MARS, or data_dir > LOGS_BETER. Furthermore, all individual logfiles need to end with **aug2010.log**, otherwise the program will not be able to find them (see point 5 of the preliminary note of this manual). 

The extraction of the behavioral data from the raw logfiles is performed by the *extract_behavioral_data* method of the EmotionTask class. This method creates a comma seperated (.csv) file that contains the onsets, picture identifiers, conditions, responses, and reaction times (RT) of all emotional picture trials. This CSV file is written to a new subdirectory in the data directory (data_dir), named after the dataset in question: i.e., data_dir > ONSETS_MARS, or data_dir > ONSETS_BETER.

In order to extract the task data from the raw logfiles, enter the following code in a console:

```
my_experiment = Amy.EmotionTask()
my_experiment.extract_behavioral_data()
```

> Note, the above block of code assumes that the [amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py) module has already been imported. If this is not the case, run the following code beforehand in the console:

```
working_dir = r'O:\Project directory\Analysis'

import os
os.chdir(working_dir)

import amygdala_project as Amy
```

### 3.2 Confound regressors

Before the first-level analysis of the emotion task data can be conducted, the confound regressor model needs to be extracted from the data. The confound regressors are extracted using the *extract_confound_regressors* method of the EmotionTask class. This method computes a total of 36 nuisance parameters, following recent recommendations by [Satterthwaite et al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3811142/), [Ciric et al. (2017)](https://www.sciencedirect.com/science/article/pii/S1053811917302288), [Parkes et al. (2018)](https://www.sciencedirect.com/science/article/pii/S1053811917310972), [Ciric et al. (2018)](https://www.nature.com/articles/s41596-018-0065-y), and [Satterthwaite et al. (2019)](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.23665): i.e., the six realignment parameters (6P). their derivatives (12P), quadratic terms (18P), and quadratic terms of their derivatives (24P), as well as the white-matter, CSF, and global mean signal (27P), their derivatives (30P), quadratic terms (33P), and quadratic terms of their derivatives (36P).

To execute the *extract_confound_regressors* method, enter the following code in the console:

```
my_experiment = Amy.EmotionTask()
my_experiment.extract_confound_regressors()
```

Since the EmotionTask subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter EMO for the emotion task data.
2. An optional prefix to indicate the gray matter scan to base the confound model on. It is recommended that this option is skipped at this stage. Simply press the enter key to continue.
3. An optional prefix to indicate the white matter and CSF scans to base the confound model on. It is recommended that the eroded white-matter and CSF segmentations are used at this stage. Enter e for the eroded white-matter and CSF segmentations.
4. An optional prefix to indicate the exact functinoal scan identifiers on which the analysis needs to be performed. It is recommended that the slice-time corrected and realigned functional images are used at this stage. Enter ar for the slice-time corrected realigned functional scans.

The confound regressor process creates an output CSV file in the emotion task scan directory called (e.g.) data_dir > NIFTI_MARS_EMO > xm13101101 > xm13101101_3_1 > **xm13101101_3_1_confound_regressors.csv**. This file can be used to define the nuisance regressors of the GLM defined in the first-level analysis described in section 3.4.

### 3.3 Spike regressors

Before the first-level analysis of the emotion task data can be conducted, the spike regressors need to be extracted from the FD Jenkinson data. The spike regressors are extracted using the *extract_spike_regressors* method of the EmotionTask class. The number and identity of the spike regressors computed by this method are defined by how many and what frames exceed a given pre-specified threshold (e.g. 0.2 mm).

To execute the *extract_spike_regressors* method, enter the following code in the console:

```
my_experiment = Amy.EmotionTask()
my_experiment.extract_spike_regressors()
```

Since the EmotionTask subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter EMO for the emotion task data.
2. The FD Jenkinson threshold above which individual frames will be flagged as outliers. A threshold of 0.2 mm is recommended at this stage.

The spike regressor process creates an output CSV file in the emotion task scan directory called (e.g.) data_dir > NIFTI_MARS_EMO > xm13101101 > xm13101101_3_1 > **xm13101101_3_1_spike_regressors.csv**. This file can be used to define the GLM in the first-level analysis described in section 3.4. The process also creates an output text files in the working directory, containing, for all subjects of the dataset in question (MARS, BETER), the number of outliers and mean framewise displacement (MFD). This file is stored in the working directory as either **MARS_Outliers_FD_Jenkinson_EMO.txt** or **BETER_Outliers_FD_Jenkinson_EMO.txt**.

### 3.4 First-level analysis

The first-level analysis of the emotion task data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts responsible for the first-level procedure are [FirstLevelAnalysis.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/FirstLevelAnalysis.m) and [FirstLevelAnalysis_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/FirstLevelAnalysis_job.m). The first-level analysis step of the pipeline conducts the following three procedures:
1. Model specification, in which Neutral, Positive, and Negative (picture) conditions are defined, respectively, with the task onsets and durations derived from the output files generated in section 3.1. The confound regressors and spike regressors generated in sections 3.2 and 3.3, respectively, are entered as nuisance variables at this stage.
2. Model estimation, using the SPM.mat file (e.g., data_dir > NIFTI_MARS_EMO > xm13101101 > xm13101101_3_1 > LEV1 > **SPM.mat**) generated in the model specification step.
3. Contrast specification, using the following contrasts (respectively): Negative vs. Neutral, Positive vs. Neutral, and Negative vs. Positive.

To perform the first-level analysis, enter the following code in a console:

```
my_experiment = Amy.EmotionTask()
my_experiment.run_1st_level_analysis()
```

Since the EmotionTask subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter EMO for the emotion task data.
2. An optional explicit mask to base the first-level analysis on. It is recommended that this option be skipped at this stage. Simply press the enter key to continue.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. It is recommended that the slice-time corrected and realigned functional data is used at this stage. Enter ar for the slice-time corrected and realigned functional (emotion task) images.

The output beta, contrast, and t-maps of the first-level procedure are stored in a subdirectory called LEV1, in the emotion task folder of each subject; e.g., data_dir > NIFTI_MARS_EMO > xm13101101 > xm13101101_3_1 > **LEV1**. In this folder can be found the beta-maps of each of the predictors specified in the model, as well as the specified contrast maps and corresponding t-maps.

### 3.5 Native ROI mask

AFter the first-level analysis has been performed, the statistical parametric data of the region-of-interest (ROI), in this case the amygdala, needs to be extracted. The first step towards accomplishing this aim is generating an ROI mask in native space for each subject in the dataset.

The *create_native_roi_mask* module utilizes [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the procedure itself are called [GlobalToNativeMask.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/GlobalToNativeMask.m) and [GlobalToNativeMask_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/GlobalToNativeMask_job.m). The method inverse normalizes the specified ROI mask, via the following 5-step procedure:
1. The coregistered (to the mean functional image) anatomical image is normalized to NMI space.
2. The ROI mask is coregistered to the normalized anatomical image generated in step 1.
3. The deformation matrix yielded by step 1 is used to create an inverse (MNI to native) deformation field.
4. The inverse deformation field yielded by step 3 is used to reslice the T1-coregistered ROI mask yielded by step 2 into native space.
5. The native ROI mask generated in step 4 is coregistered to the anatomical image that was used as input in step 1.

To perform the global-to-native ROI mask transformation, enter the following code in a console:

```
my_experiment = Amy.EmotionTask()
my_experiment.create_native_roi_mask()
```

Since the EmotionTask subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter EMO for the emotion task data.
2. The name of the ROI mask used for the analyses, as listed in the working directory. It is recommended that a probability map of the amygdala is used. 
3. An optional prefix to indicate the exact scan identifier of the anatomical image to base the inverse normalization on. It is recommended that the anatomical image coregistered to the mean functional image is used for this procedure. Enter c for the coregistered anatomical image.

The result of this procedure is a set of files and images that are written to the T1 directory of each subject. The native ROI mask used for further analysis steps is indicated by the prefix **cic** followed by the original filename of the ROI mask; e.g. cicAmygdala_total_probability_map.nii.

### 3.6 ROI-masked SPM data

AFter the first-level analysis has been performed, the statistical parametric data of the region-of-interest (ROI), in this case the amygdala, needs to be extracted. The second step towards accomplishing this goal is masking the statistical parametric output generated by the first-level analysis (see section 3.4), by using the native ROI mask generated by the previous stage (section 3.5). 

The *extract_roi_masked_spm_data* method computes the probability weighted mean (or sum) of the statistical parametric maps of each subject, seperately for each hemisphere. To execute the *extract_roi_masked_spm_data* method, enter the following code in the console:

```
my_experiment = Amy.EmotionTask()
my_experiment.extract_roi_masked_spm_data()
```

Since the EmotionTask subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter EMO for the emotion task data.
2. The name of the ROI mask used for the analyses, as listed in the T1 directory of each subject. The ROI mask needs to be resliced to native space (see section 3.5). It is recommended that the output file of section 3.5 (e.g. cicAmygdala_total_probability_map.nii) is used
3. The nature of the ROI mask specified in step 2. Enter 1 for a binary ROI mask, or 2 for a probabilistic one. It is recommended that a probability map of the amygdala is used
4. The name of the statistical parametric map generated by the first-level analysis (see section 3.4), on which the ROI mask should be applied (for each subject). The con_0001.nii file corresponds to the contrast Negative vs. Neutral, con_0002.nii to the contrast Positive vs. Neutral, and con_0003.nii to the contrast Negative vs. Positive.

The probability weighted mean (or summed) SPM values are written to an output CSV file in the working directory. The file is named after the specific dataset entered into the analyses (MARS, BETER), the summary statistic used to mask the parametric maps (Mean, Sum), and the specific contrast map on which the ROI mask was applied (e.g., con_0001); for example: **MARS_EMO_Mean_Con_0001.csv**

## 4. Resting-state

The analysis of the resting-state data is conducted using the RestingState subclass. This class inherits all attributes and methods of the Preprocessing class, which itself (in turn) inherits from the main Experiment class. The seed-based connectivity (first-level) analysis of the resting-state data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via shell call commands.

The following preprocessing pipeline is recommended for the resting-state data:
1. Realignment of the raw functional images (see section 2.1)
2. Extraction of the FD Jenkinson data (see section 2.3)
3. Coregistration of the anatomical image to the mean functional image (see section 2.5)
4. Extraction of a (subject-level) inclusve field-of-view (FOV) voxel mask from the realigned functional images (see section 2.2)
5. Segmentation of the coregistered anatomical image (see section 2.6)
6. Erosion of the white-matter and CSF segmentations (see section 2.7)

As part of the first-level seed-based connectivity pipeline of the resting-state data, the following steps need to be performed after the preprocessing steps:
1. Construction of a (subject-level) native ROI mask (see section 4.1)
2. Extraction of the region-of-interest (ROI) regressors (see section 4.2)
3. Extraction of the confound regressors (see section 4.3)
4. Extraction of the spike regressors (see section 4.4)
5. Voxel-wise seed-based functional connectivity analysis (see section 4.5)
6. Normalization of the beta (connectivity) maps (see section 2.8)
7. Optional: smoothing of the normalized beta (connectivity) maps (see section 2.9)

### 4.1 Native ROI mask

Before the seed-based connectivity analysis can be performed, the region-of-interest (ROI) mask needs to be transformed to native space. The *create_native_roi_mask* module used for this procedure ultilizes [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the procedure itself are [GlobalToNativeMask.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/GlobalToNativeMask.m) and [GlobalToNativeMask_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/GlobalToNativeMask_job.m). This method inverse normalizes the ROI mask specified, via the following 5-step procedure:
1. The coregistered (to the mean functional image) anatomical image is normalized to NMI space.
2. The ROI mask is coregistered to the normalized anatomical image generated in step 1.
3. The deformation matrix yielded by step 1 is used to create an inverse (MNI to native) deformation field.
4. The inverse deformation field yielded by step 3 is used to reslice the T1-coregistered ROI mask yielded by step 2 into native space.
5. The native ROI mask generated in step 4 is coregistered to the anatomical image that was used as input in step 1.

To perform the global-to-native ROI mask transformation, enter the following code in a console:

```
my_experiment = Amy.RestingState()
my_experiment.create_native_roi_mask()
```

Since the RestingState subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter REST for the resting-state data.
2. The name of the ROI mask used for the analysis, as listed in the working directory. It is recommended that a probability map of the amygdala is used. 
3. An optional prefix to indicate the exact scan identifier of the anatomical image on which to base the inverse normalization. It is recommended that the coregistered (to the mean functional image) anatomical image is selected for this procedure. Enter c for the coregistered anatomical image.

The result of this procedure is a set of files and images that are written to the T1 directory of each subject in the data directory. The native ROI mask used for further analysis steps is indicated by the prefix **cic**, followed by the original filename of the ROI mask; e.g. cicAmygdala_total_probability_map.nii.


### 4.2 ROI regressors

Before the connectivity analysis of the resting-state data can be conducted, the ROI regressors (left/right hemisphere) need to be extracted from the data. The ROI regressors are extracted using the *extract_roi_regressors* method of the RestingState class. This method computes one ROI regressor per hemisphere, based on a pre-specified ROI mask, which can be either probabilistic or binary.

To execute the *extract_roi_regressors* method, enter the following code in the console:

```
my_experiment = Amy.RestingState()
my_experiment.extract_roi_regressors()
```

Since the RestingState subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter REST for the resting-state data.
2. The name of the ROI mask used for the analysis, as listed in the working directory. It is recommended that a probability map of the amygdala is used.
3. The nature of the ROI mask specified in step 2. Enter 1 for a binary ROI mask, or 2 for a probabilistic one.
4. An optional prefix to indicate what version of the ROI mask should be used for the analyses. It is recommended that the mask is transformed to subject/native space using the procedure described in section 4.1, which generates an output native ROI mask in the T1 directory of each subject. Enter cic to utilize this output (of section 4.1) native ROI mask.
5. An optional prefix to indicate the exact functinoal scan identifiers on which the analysis needs to be performed. It is recommended that the realigned functional images are used at this stage. Enter r for the realigned functional scans.

The ROI regressor process creates two output CSV files in the resting-state scan directory called (e.g.) data_dir > NIFTI_MARS_REST > xm13101101 > xm13101101_5_1 > **xm13101101_5_1_lh_roi regressor.csv** and **xm13101101_5_1_rh_roi regressor.csv**. These files can be used to define the predictor-of-interest signal used in the GLM specified in the voxel-wise connectivity analysis, as described in section 4.5.

### 4.3 Confound regressors

Before the connectivity analysis of the resting-state data can be conducted, the confound regressor model needs to be extracted from the data. The confound regressors are extracted using the *extract_confound_regressors* method of the RestingState class. This method computes a total of 36 nuisance parameters, following recent recommendations by [Satterthwaite et al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3811142/), [Ciric et al. (2017)](https://www.sciencedirect.com/science/article/pii/S1053811917302288), [Parkes et al. (2018)](https://www.sciencedirect.com/science/article/pii/S1053811917310972), [Ciric et al. (2018)](https://www.nature.com/articles/s41596-018-0065-y), and [Satterthwaite et al. (2019)](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.23665): i.e., the six realignment parameters (6P). their derivatives (12P), quadratic terms (18P), and quadratic terms of their derivatives (24P), as well as the white-matter, CSF, and global mean signal (27P), their derivatives (30P), quadratic terms (33P), and quadratic terms of their derivatives (36P).

To execute the *extract_confound_regressors* method, enter the following code in the console:

```
my_experiment = Amy.RestingState()
my_experiment.extract_confound_regressors()
```

Since the RestingState subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter REST for the resting-state data.
2. An optional prefix to indicate the gray matter scan to base the confound model on. It is recommended that this option is skipped at this stage. Simply press the enter key to continue.
3. An optional prefix to indicate the white matter an CSF masks to base the confound model on. It is recommended that the eroded white-matter and CSF segmentations are used at this stage. Enter e for the eroded white-matter and CSF segmentations.
4. An optional prefix to indicate the exact functinoal scan identifiers on which the analysis will to be performed. It is recommended that the realigned functional images are used at this stage. Enter r for the realigned functional scans.

The confound regressor process creates an output CSV file in the resting-state scan directory called (e.g.) data_dir > NIFTI_MARS_REST > xm13101101 > xm13101101_5_1 > **xm13101101_5_1_confound_regressors.csv**. This file can be used to define the nuisance regressors of the GLM defined in voxel-wise connectivity analysis model described in section 4.5.

### 4.4 Spike regressors

Before the connectivity analysis of the resting-state data can be conducted, the spike regressors need to be extracted from the FD Jenkinson data. The spike regressors are extracted using the *extract_spike_regressors* method of the RestingState class. The number and identity of the spike regressors computed by this method are defined by how many and what frames exceed a given pre-specified threshold (e.g. 0.2 mm).

To execute the *extract_spike_regressors* method, enter the following code in the console:

```
my_experiment = Amy.RestingState()
my_experiment.extract_spike_regressors()
```

Since the RestingState subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter REST for the resting-state data.
2. The FD Jenkinson threshold above which individual frames will be flagged as outliers. A threshold of 0.2 mm is recommended at this stage.

The spike regressor process creates an output CSV file in the resting-state scan directory called (e.g.) data_dir > NIFTI_MARS_REST > xm13101101 > xm13101101_5_1 > **xm13101101_5_1_spike_regressors.csv**. This file can be used to define the GLM in the first-level analysis described in section 4.5. The process also creates an output text file in the working directory, containing, for all subjects of the dataset in question (MARS, BETER), the number of outliers and mean framewise displacement (MFD). This file is stored in the working directory as either **MARS_Outliers_FD_Jenkinson_REST.txt** or **BETER_Outliers_FD_Jenkinson_REST.txt**.

### 4.5 Voxel-wise connectivity analysis

The voxel-wise seed-based functional connectivity analysis of the resting-state data is conducted via the *voxel_wise_connectivity_analysis* method. This method performs two GLMs for each voxel within the grey matter mask of the brain.. In the first GLM, the signal of the left ROI regressor is entered as predictor-of-interest, along with the 36 confound model (see section 4.3) and spike regressors (see section 4.4) as covariates-of-no-interest, and the signal of each voxel (individually) as outcome variable. In the second GLM, the signal of the right ROI regressor is entered as predictor-of-interest, again with the 36 confound model (see section 4.3) and spike regressors (see section 4.4) as covariates-of-no-interest, and the signal of each voxel (individually) as outcome variable. 

The GLM (beta) coefficients of the predictors-of-interest (left/right ROI) are extracted by the method, and stored as NIFTI files that can be found in the resting-state directory of each subject (e.g. **b_map_lh_xm13101101.nii**, **b_map_rh_xm13101101.nii**). The beta-values contained within each of these files represent a quantification of the relationship of the left/right ROI signals, in this case the left and right amygdalae, with the signal of each individual gray matter voxel, corrected for the nuisance effects of motion, white-matter, CSF, and global signal (36 confound model), with censoring of the frames flagged as framewise displacement outliers (spike regressors).

To perform the voxel-wise seed-based connectivity analysis, enter the following code in a console:

```
my_experiment = Amy.RestingState()
my_experiment.voxel_wise_connectivity_analysis()
```

Since the RestingState subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter REST for the resting-state data.
2. An optional prefix to indicate the exact scan identifiers on which the conncetivity analysis needs to be conducted. It is recommended that the realigned functional data is used at this stage. Enter r for the realigned functional (resting-state) images.

