# AMYGDALA_REACT_VS_CONNECT
Pipeline to predict emotion task reactivity of the amygdala using resting-state connectivity of the amygdala as seed-region.

> Preliminary note: The Python module [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) used here assumes the following steps having been undertaken before running the analyses:
> 1. Many of the procedures described in this manual make use of [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. Framewise displacement in this pipeline is calculated via the [Jenkinson](https://pdfs.semanticscholar.org/291c/9d66cf886ed27d89dc03d42a4612b66ab007.pdf) method using code of the [DPARSF](http://rfmri.org/DPARSF) software package. The bandpass filter used in many of the methods described here was originally developed by Thomas Gladwin; the ezfilt.m script can be found at: https://www.tegladwin.com/code.php. **Each of these programs have to be installed in a subfolder called Programs that is placed within the working directory in order for this pipeline to work properly.** The reason this last step is important is because many of the methods described here employ a MATLAB script called [start_up](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/start_up.m), which adds the above programs to the MATLAB path, but in order to do so, requires each of the programs to have been placed in the working_dir > Programs subfolder. 
> 2. The necessary analysis scripts and auxilliary files need to stored in one and the same working directory. The emotion task (EMO) and resting-state (REST) fMRI data needs to have been converted to NIFTI format, and stored in a (seperate) data directory (data_dir) organized such that each combination of dataset (MARS, BETER) and scan type (EMO, REST) has a different subdirectory (i.e., NIFTI_MARS_REST, NIFTI_MARS_EMO, NIFTI_BETER_REST, and NIFTI_BETER_EMO) within that data directory (e.g. data_dir > NIFTI_MARS_REST). 
> 3. Within each of the NIFTI data subdirectories, each subject has to have its own subfolder, with the name of that subfolder corresponding to the identifier of that subject (e.g. xm13101101). The different types of scans (T1, REST, EMO) for each subject need to be stored as subdirectories within that subject directory, and the identifiers of these subdirectories need to be indicated by suffixes that are added to the main subject identifier (e.g. data_dir > NIFTI_BETER_EMO > be1a031010 > **be1a031010_5_1**). 
> 4. The first two letters of each subject identifier specify the dataset where to the subject belongs; 'xm' indicates the MARS dataset (e.g. xm13101101), and 'be' indicates the BETER dataset (e.g. be1a031010). 
> 5. The base directory of the data needs to contain, for each dataset (MARS, BETER) seperately, a directory (LOGS_MARS, LOGS_BETER) that contains the log files of the emotion task (e.g. data_dir > LOGS_BETER > be1a031010-3amygdala_hippo_aug2010.log). These log files need to end with **aug2010.log**, otherwise the program will not be able to identify them. 
> 6. An Excel (.xlsx) file needs to be present in the working directory which specifies, for each subject, the subject identifier (e.g. xm13101101), the subject number (e.g. MARS011), group membership (e.g. 1 or 0) and the specific suffixes added to the main subject ID in order to differentiate the different scan subdirectories for that subject (e.g. \_4_1, \_5_1, etc.).

Table of contents:
1. [Setting up the experiment](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#1-setting-up-the-experiment)
2. [Data preprocessing](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#2-data-preprocessing)
    1. [Slice-timing correction](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#21-slice-timing-correction)
    2. [Realignment](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#22-realignment)
    3. [Framewise displacement](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#23-framewise-displacement)
    4. [Coregistration](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#24-coregistration)
    5. [Normalization](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#25-normalization)  
    6. [Inclusive mask extraction](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#26-inclusive-mask-extraction)
    7. [Segmentation](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#27-segmentation)
    8. [Erosion](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#28-erosion)
3. [Emotion task](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#3-emotion-task)
    1. [Behavioral data](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#31-behavioral-data)
    2. [Confound regressors](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#32-confound-regressors)
    3. [Spike regressors](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#33-spike-regressors)
    4. [First-level analysis](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#34-first-level-analysis)
    5. [Reslice ROI mask](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#35-reslice-roi-mask)
    6. [ROI-masked SPM data](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#36-roi-masked-spm-data)
4. [Resting-state](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#4-resting-state)
    1. [Reslice ROI mask](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#41-reslice-roi-mask)
    2. [ROI regressors](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#42-roi-regressors)
    3. [Confound regressors](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#43-confound-regressors)
    4. [Spike regressors](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#44-spike-regressors)
    5. [Voxel-wise connectivity analysis](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#45-voxel-wise-connectivity-analysis)
5. [Post-processing](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#5-post-processing)
    1. [Inclusive FOV mask](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#51-inclusive-fov-mask)
    3. [Inclusive grey matter mask](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#52-inclusive-grey-matter-mask)
    2. [Average brain](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#53-average-brain)
    4. [QC-FC motion correction benchmark](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#54-qc-fc-motion-correction-benchmark)
    5. [Discriminability motion correction benchmark](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#55-discriminability-motion-correction-benchmark)
6. [Group-level analysis](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#6-group-level-analysis)
    1. [Second-level analysis](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#61-second-level-analysis)
    2. [Voxel-matched regression analysis](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#62-voxel-matched-regression-analysis)
    3. [ROI-based regression analysis](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#63-roi-based-regression-analysis)
    4. [Statistical non-parametric analysis: specification and computation](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#64-statistical-non-parametric-analysis-specification-and-computation)
    5. [Statistical non-parametric analysis: inference](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#65-statistical-non-parametric-analysis-inference)

## 1. Setting up the experiment

The Python module [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) that is used to run the analyses described here consists of multiple classes. The parent class that is inherited by all the other subclasses is called Experiment. This class has a (constructor) method called \_\_init__ that constructs all the attributes needed by (and specific for) the experiment in question. Three user inputs need to be specified in the console when initializing the Experiment class:
1. The full working directory that contains all the scripts and auxilliary files needed for the analysis.
2. The full base directory of the data, containing each subject identifier as a seperate folder, in which each of the relevant scans of the subjects are specified by subdirectories in the following manner: data_dir > NIFTI_BETER_REST > be1a031010 > be1a031010_5_1 (see also point 3 of the preliminary note at the beginning of this manual).
3. The full directory of the Excel file that contains the subject information. The first row in this file needs to contain at least the following column headers: SubjID (or SubjT0ID, containing the subject identifiers), PPN (containing the subject number), Group (containing the group memberships), T1 (containing the anatomical scan directory suffixes for each subject), REST (containing the resting-state directory suffixes for each subject), and EMO (containing the task fMRI directory suffixes for each subject). Other columns containing e.g. questionnaire data can also be specified. Each row in the Excel file specifies a single subject, with the first two letters specifying the dataset that subject belongs to; 'xm' for the MARS dataset, and 'be' for the BETER dataset (see also point 6 of the preliminary note at the beginning of this user manual)

To initialize the Experiment class and all its attributes, define the working directory and run the following code in the console:

```
working_dir = r'O:\Project directory\Analysis'

import os
os.chdir(working_dir)

import amygdala_recon as Amy
my_experiment = Amy.Experiment()
```
This code generates an output argument called *my_experiment* which contains attributes that define the parameters of the experiment (e.g. *my_experiment.subjects_list* contains a list of all the subject identifiers).

> Note: Since the Experiment class is inherited by all the other subclasses, it is not strictly necessary to initialize this method on its own, before conducting any of the other procedures described here (although it is possible). It is called automatically when initializing the other subclasses of the pipeline.

> Note: Later on, for the post-processing and group analysis steps, the above-mentioned Excel file that contains the subject information also needs to contain an additional column that detail whether or not a subject will be included in the second-level analysis (column header: Include).

## 2. Data preprocessing

Many of the data preprocessing methods in the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module are conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command.

The following preprocessing steps are supported by the pipeline:
1. Slice-timing correction of the functional data ([section 2.1](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#21-slice-timing-correction))
2. Realignment of the functional data ([section 2.2](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#22-realignment))
3. Calculation of framewise displacements ([section 2.3](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#23-framewise-displacement))
4. Coregistration of the functional images to the anatomical image ([section 2.4](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#24-coregistration))
5. Normalization of the functional data ([section 2.5](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#25-normalization))
6. Extraction of an inclusive mask, based on the functional images ([section 2.6](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#26-inclusive-mask-extraction))
7. Segmentation of the anatomical image ([section 2.7](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#27-segmentation))
8. Erosion of the segmentation data ([section 2.8](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#28-erosion))

### 2.1 Slice-timing correction

Slice-timing correction of the functional data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are used for the slice-timing correction are [SliceTimingCorrection.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/SliceTimingCorrection.m) and [SliceTimingCorrection_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/SliceTimingCorrection_job.m)

To conduct the slice-timing correction, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The specific preprocessing step that needs to be conducted. Enter 1 for slice-timing correction.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. This option should be skipped at this stage; simply press the enter key to continue.

### 2.2 Realignment

Realignment of the functional data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the realignment process are [Realignment.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Realignment.m) and [Realignment_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Realignment_job.m)

To conduct the realignment process, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The specific preprocessing step that needs to be conducted. Enter 2 for realignment.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. It is recommended that the motion parameters are estimated from the raw unprocessed functional images, while the realignment itself is performed on the slice-time corrected images. Enter a for the slice-time corrected scans, or simply press enter to continue.

> Note, the above block of code assumes that the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module has already been imported. If this is not the case, run the following code beforehand in the console:

```
working_dir = r'O:\Project directory\Analysis'

import os
os.chdir(working_dir)

import amygdala_recon as Amy
```

### 2.3 Framewise displacement

The extraction of framewise displacement (FD) data is conducted using [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB function that is responsible for the extraction process is called [FrameWiseDisplacement.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/FrameWiseDisplacement.m), and makes use of the y_FD_Jenkinson.m script of the [DPARSF](http://rfmri.org/DPARSF) software package. The DPARSF function calculates the framewise displacement data using the [Jenkinson](https://pdfs.semanticscholar.org/291c/9d66cf886ed27d89dc03d42a4612b66ab007.pdf) method, based on the motion parameter estimates generated by the realignment procedure in [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (see section 2.1). 

Enter the following code in the console to extract the FD Jenkinson data:

```
my_experiment = Amy.Preprocessing()
my_experiment.extract_FD_jenkinson()
```

> Note: It is essential that both [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and [DPARSF](http://rfmri.org/DPARSF) are added to the MATLAB path in order for the *extract_FD_jenkinson* method to work properly.

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The reference scan that is used to calculate the framewise displacement data. Enter 1 for the first image (recommended), or M for the mean image.

### 2.4 Coregistration

Coregistration of the functional images to the raw anatomical image is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the coregistration are [Coregistration.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Coregistration.m) and [Coregistration_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Coregistration_job.m)

To conduct the coregistration procedure, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The specific preprocessing step that needs to be conducted. Enter 3 for coregistration.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. The raw anatomical volume needs to be used as reference image at this stage, with the mean EPI as a source image; the other functional scans need remain in alignment with the mean image. The appropriate prefix of these other scans should be entered at this stage: Enter ra for the realigned slice-time corrected functional images, or r for the realigned (non-slice-timing corrected) images. 

### 2.5 Normalization

Normalization of the functional data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the normalization procedure are [Normalization.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Normalization.m) and [Normalization_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Normalization_job.m)

To conduct the normalization procedure, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The specific preprocessing step that needs to be conducted. Enter 4 for normalization.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. In theory, the normalization can be performed on any series of functional images. Enter ra for the realigned slice-time corrected (coregistered) functional images, or r for the realigned (non-slice-timing corrected) (coregistered) images.

### 2.6 Inclusive mask extraction

Extraction of a subject-specific inclusive voxel mask is supported by the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module via the *extract_inclusive_FOV_mask* method of the Preprocessing class. Enter the following code in the console to execute this process:

```
my_experiment = Amy.Preprocessing()
my_experiment.extract_inclusive_FOV_mask()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. An optional prefix to indicate the exact scan identifiers to base the extraction of the inclusive mask on. Enter nra for the realigned, slice-time corrected, and (coregistered) normalized functional images, or nr for the realigned (non-slice-timing corrected, but coregistered) normalized images.


### 2.7 Segmentation 

Segmentation of the anatomical data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the segmentation are [Segmentation.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Segmentation.m) and [Segmentation_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Segmentation_job.m)

To conduct the segmentation procedure, enter the following code in a console:

```
my_experiment = Amy.Preprocessing()
my_experiment.run_preprocessing()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The specific preprocessing step that needs to be conducted. Enter 5 for segmentation.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. The segmentation can be conducted on either the raw or normalized anatomical data. Enter n for the normalized anatomical image, or simply press enter to use the raw anatomical scan.

### 2.8 Erosion 

Erosion of the segmentation data is supported by the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module, via the *erode_segmentation_mask* method of the Preprocessing class. This method uses the MATLAB function [Erosion.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/Erosion.m) via a shell call command. 

Enter the following code in the console to execute this process:

```
my_experiment = Amy.Preprocessing()
my_experiment.erode_segmentation_mask()
```

Since the Preprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the preprocessing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. An optional prefix to indicate the exact scan identifiers on which the erosion needs to be applied. In theory, the erosion procedure can be conducted on all outputs of the segmentation procedure (c1, c2, c3, c4, c5). To perform the erosion on the white-matter or CSF segmentations, enter c2 or c3, respectively.

## 3. Emotion task

The emotion task described here is based on the paradigm of [Van Buuren et al. (2011)](https://www.sciencedirect.com/science/article/pii/S000632231100254X) (see also [Heesink et al., 2018](https://www.cambridge.org/core/journals/european-psychiatry/article/neural-activity-during-the-viewing-of-emotional-pictures-in-veterans-with-pathological-anger-and-aggression/C3179D960B2B1DA0C0C33D34B5FCA2D6)).

The analysis of the emotion task data is conducted via the EmotionTask subclass. This class inherits all attributes and methods of the Preprocessing class, which itself (in turn) inherits from the main Experiment class. The first-level analysis of the emotion task data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via shell call commands.

The following preprocessing pipeline is recommended for the emotion task data:
1. Slice-timing correction of the raw functional images (see [section 2.1](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#21-slice-timing-correction))
2. Realignment of the (slice-time corrected) functional images (see [section 2.2](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#22-realignment))
3. Extraction of the FD Jenkinson data (see [section 2.3](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#23-framewise-displacement))
4. Coregistration of the functional images to the anatomical image (see [section 2.4](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#24-coregistration))
5. Normalization of the realigned (slice-time corrected) and coregistered functional images (see [section 2.5](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#25-normalization))
6. Extraction of a (subject-level) inclusve field-of-view (FOV) voxel mask from the preprocessed functional images (see [section 2.6](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#26-inclusive-mask-extraction))
7. Segmentation of the anatomical image (see [section 2.7](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#27-segmentation))
8. Erosion of the white-matter and CSF segmentations (see [section 2.8](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#28-erosion))

As part of the first-level analysis pipeline of the emotion task data, the following steps need to be performed after the preprocessing steps:
1. Extraction of the behavioral task data from the logfiles (see [section 3.1](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#31-behavioral-data))
2. Extraction of the confound regressors (see [section 3.2](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#32-confound-regressors))
3. Extraction of the spike regressors (see [section 3.3](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#33-spike-regressors))
4. First-level analysis in [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (see [section 3.4](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#34-first-level-analysis))
5. Construction of a region-of-interest (ROI) mask in subject reference space (see [section 3.5](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#35-reslice-roi-mask))
6. Extraction (across subjects) of the ROI-masked statistical paramatric data (see [section 3.6](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#36-roi-masked-spm-data))

### 3.1 Behavioral data

Before the first-level analysis of the emotion task data can be conducted, the task data needs to be extracted from the raw logfiles. As described in the preliminary note of this manual, the raw logfiles need to be copied to the data directory (data_dir), in a subfolder named after the dataset in question; i.e., data_dir > LOGS_MARS, or data_dir > LOGS_BETER. Furthermore, all individual logfiles need to end with **aug2010.log**, otherwise the program will not be able to find them (see point 5 of the preliminary note of this manual). 

The extraction of the behavioral data from the raw logfiles is performed by the *extract_behavioral_data* method of the EmotionTask class. This method creates a comma seperated (.csv) file that contains the onsets, picture identifiers, conditions, responses, and reaction times (RT) of all emotional picture trials, as well as two columns that are masked by congruency (condition == response, otherwise zero) and incongruency (condition /= response, otherwise zero). This CSV file is written to a new subdirectory in the data directory (data_dir), named after the dataset in question: i.e., data_dir > ONSETS_MARS, or data_dir > ONSETS_BETER.

In order to extract the task data from the raw logfiles, enter the following code in a console:

```
my_experiment = Amy.EmotionTask()
my_experiment.extract_behavioral_data()
```

> Note, the above block of code assumes that the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module has already been imported. If this is not the case, run the following code beforehand in the console:

```
working_dir = r'O:\Project directory\Analysis'

import os
os.chdir(working_dir)

import amygdala_recon as Amy
```

### 3.2 Confound regressors

Before the first-level analysis of the emotion task data can be conducted, the confound regressor model needs to be extracted from the data. The confound regressors are extracted using the *extract_confound_regressors* method of the EmotionTask class. This method computes a total of 8 nuisance parameters: i.e., the six realignment parameters (6P). and the white-matter and CSF signals (8P).

To execute the *extract_confound_regressors* method, enter the following code in the console:

```
my_experiment = Amy.EmotionTask()
my_experiment.extract_confound_regressors()
```

Since the EmotionTask subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter EMO for the emotion task data.
2. An optional prefix to indicate the white matter and CSF scans to base the confound model on. It is recommended that the eroded white-matter and CSF segmentations are used at this stage. Enter e for the eroded white-matter and CSF segmentations.
3. An optional prefix to indicate the exact functional scan identifiers on which the analysis needs to be performed. It is recommended that the realigned, slice-time corrected, and (coregistered) normalized functional images are used at this stage. Enter nra for the realigned, slice-time corrected, and (coregistered) normalized functional scans.

The confound regressor process creates an output CSV file in the emotion task scan directory called (e.g.) data_dir > NIFTI_MARS_EMO > xm13101101 > xm13101101_3_1 > **xm13101101_3_1_confound_regressors.csv**. This file can be used to define the nuisance regressors of the GLM defined in the first-level analysis described in section 3.4.

### 3.3 Spike regressors

Before the first-level analysis of the emotion task data is conducted, the spike regressors can, but do not have to be extracted from the FD Jenkinson data (this step is optional). The spike regressors are extracted using the *extract_spike_regressors* method of the EmotionTask class. The number and identity of the spike regressors computed by this method are defined by how many and what frames exceed a given pre-specified threshold (e.g. 0.5 mm).

To execute the *extract_spike_regressors* method, enter the following code in the console:

```
my_experiment = Amy.EmotionTask()
my_experiment.extract_spike_regressors()
```

Since the EmotionTask subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter EMO for the emotion task data.
2. The FD Jenkinson threshold above which individual frames will be flagged as outliers. A threshold of 0.5 mm is recommended at this stage (see [Siegel et al.,2013](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.22307)).

The spike regressor process creates an output CSV file in the emotion task scan directory called (e.g.) data_dir > NIFTI_MARS_EMO > xm13101101 > xm13101101_3_1 > **xm13101101_3_1_spike_regressors.csv**. This file can be used to define the GLM in the first-level analysis described in section 3.4. The process also creates an output text files in the working directory, containing, for all subjects of the dataset in question (MARS, BETER), the number of outliers and mean framewise displacement (MFD). This file is stored in the working directory as either **MARS_Outliers_FD_Jenkinson_EMO.txt** or **BETER_Outliers_FD_Jenkinson_EMO.txt**.

### 3.4 First-level analysis

The first-level analysis of the emotion task data is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts responsible for the first-level procedure are [FirstLevelAnalysis.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/FirstLevelAnalysis.m) and [FirstLevelAnalysis_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/FirstLevelAnalysis_job.m). The first-level analysis step of the pipeline conducts the following three procedures:
1. Model specification, in which Neutral, Positive, and Negative (picture) conditions are defined, respectively, with the task onsets and durations derived from the output files generated in section 3.1. The confound regressors (and spike regressors, if so desired) generated in section 3.2 (and 3.3) are entered as nuisance variables at this stage.
2. Model estimation, using the SPM.mat file (e.g., data_dir > NIFTI_MARS_EMO > xm13101101 > xm13101101_3_1 > first_level_analysis > **SPM.mat**) generated in the model specification step.
3. Contrast specification, using the following contrasts (respectively): Pictures vs. Baseline, Negative vs. Neutral, Positive vs. Neutral, and Negative vs. Positive, Negative + Positive vs. Neutral.

To perform the first-level analysis, enter the following code in a console:

```
my_experiment = Amy.EmotionTask()
my_experiment.run_1st_level_analysis()
```

Since the EmotionTask subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter EMO for the emotion task data.
2. An optional explicit mask to base the first-level analysis on. It is recommended that this option be skipped at this stage. Simply press the enter key to continue.
3. An optional prefix to indicate the exact scan identifiers on which the preprocessing needs to be conducted. It is recommended that the realigned, slice-time corrected, and (coregistered) normalized functional data is used at this stage. Enter nra for the realigned, slice-time corrected, and (coregistered) normalized functional (emotion task) images.

The output beta, contrast, and t-maps of the first-level procedure are stored in a subdirectory called first_level_analysis, in the emotion task folder of each subject; e.g., data_dir > NIFTI_MARS_EMO > xm13101101 > xm13101101_3_1 > **first_level_analysis**. In this folder can be found the beta-maps of each of the predictors specified in the model, as well as the specified contrast maps and corresponding t-maps.

### 3.5 Reslice ROI mask

AFter the first-level analysis has been performed, the statistical parametric data of the region-of-interest (ROI), in this case the amygdala, needs to be extracted. The first step towards accomplishing this aim is generating an ROI mask in the reference space of each subject in the dataset. This can be achieved in either of two ways: 
- Option 1: The analyses are conducted mainly in native space, and therefore the ROI mask needs to be resliced to the native reference space of each subject; this can be accomplished via the *create_native_roi_mask* method. 
- Option 2: The analyses are conducted mainly in standard MNI space, and therefore the ROI mask needs to be resliced to the standard reference space of each subject; this can be achieved via the *obtain_global_roi_mask* method.

The *create_native_roi_mask* method utilizes [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the procedure itself are called [GlobalToNativeMask.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/GlobalToNativeMask.m) and [GlobalToNativeMask_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/GlobalToNativeMask_job.m). The method inverse normalizes the specified ROI mask, via the following 5-step procedure:
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
2. The name of the ROI mask used for the analyses, as listed in the working directory. A binary mask or probability map of the amygdala should be used. 
3. An optional prefix to indicate the exact scan identifier of the anatomical image to base the inverse normalization on. It is recommended that the anatomical image coregistered to the mean functional image is used for this procedure. Enter c for the coregistered anatomical image.

The result of this procedure is a set of files and images that are written to the T1 directory of each subject. The native ROI mask used for further analysis steps is indicated by the prefix **cic** followed by the original filename of the ROI mask; e.g. cicAmygdala_total_binary_mask_thr_0_5.nii.

The *obtain_global_roi_mask* method also utilizes [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the procedure itself are called [CoregisterGlobalMask.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/CoregisterGlobalMask.m) and [CoregisterGlobalMask_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/CoregisterGlobalMask_job.m). The method reslices the standard ROI mask to the standard reference space of each subject (i.e., with the same voxel size).

To reslice the ROI mask to standard reference space, enter the following code in a console:

```
my_experiment = Amy.EmotionTask()
my_experiment.obtain_global_roi_mask()
```

Since the EmotionTask subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter EMO for the emotion task data.
2. The name of the ROI mask used for the analyses, as listed in the working directory. A binary mask or probability map of the amygdala should be used. 
3. An optional prefix to indicate the exact scan identifier of the anatomical image to base the inverse normalization on. It is recommended that the normalized anatomical image is used for this procedure. Enter n for the normalized anatomical image.

The output of this procedure is a resliced version of the ROI mask in standard reference space for each subject. The resliced ROI mask used for further analysis steps is indicated by the prefix **r** followed by the original filename of the ROI mask; e.g. rAmygdala_total_binary_mask_thr_0_5.nii.

### 3.6 ROI-masked SPM data

AFter the first-level analysis has been performed, the statistical parametric data of the region-of-interest (ROI), in this case the amygdala, needs to be extracted. The second step towards accomplishing this goal is masking the statistical parametric output generated by the first-level analysis (see section 3.4), by using the resliced ROI mask generated by the previous stage (section 3.5). 

The *extract_roi_masked_spm_data* method computes the ROI mean of the statistical parametric maps of each subject, separately for each hemisphere. To execute the *extract_roi_masked_spm_data* method, enter the following code in the console:

```
my_experiment = Amy.EmotionTask()
my_experiment.extract_roi_masked_spm_data()
```

Since the EmotionTask subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter EMO for the emotion task data.
2. The name of the ROI mask used for the analyses, as listed in the T1 directory of each subject. The ROI mask needs to be resliced to native space (see section 3.5). It is recommended that the output file of section 3.5 (e.g. rAmygdala_total_binary_mask_thr_0_5.nii, or cicAmygdala_total_binary_mask_thr_0_5.nii) is used
3. The nature of the ROI mask specified in step 2. Enter 1 for a binary ROI mask, or 2 for a probabilistic one.
4. The name of the statistical parametric map generated by the first-level analysis (see section 3.4), on which the ROI mask should be applied (for each subject). The con_0001.nii/spmT_0001.nii files correspond to the contrast Pictures (Neutral + Positive + Negative) vs. Baseline, con_0002.nii/spmT_0002.nii to the contrast Negative vs. Neutral, con_0003.nii/spmT_0003.nii to the contrast Positive vs. Neutral, and con_0004.nii/spmT_0004.nii to the contrast Negative vs. Positive, and con_0005.nii/spmT_0005.nii to the contrast Negative + Positive vs. Neutral.

The ROI-weighted subjects means of SPM values are written to an output CSV file in the working directory. The file is named after the specific dataset entered into the analyses (MARS, BETER), the summary statistic used to mask the parametric maps (Mean, Sum), and the specific contrast map on which the ROI mask was applied (e.g., con_0001); for example: **MARS_EMO_Mean_Con_0001.csv** or **MARS_EMO_Mean_spmT_0001.csv**

## 4. Resting-state

The analysis of the resting-state data is conducted using the RestingState subclass. This class inherits all attributes and methods of the Preprocessing class, which itself (in turn) inherits from the main Experiment class.

The following preprocessing pipeline is recommended for the resting-state data:
1. Realignment of the raw functional images (see [section 2.2](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#22-realignment))
2. Extraction of the FD Jenkinson data (see [section 2.3](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#23-framewise-displacement))
3. Coregistration of the functional images to the anatomical image (see [section 2.4](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#24-coregistration))
4. Normalization of the realigned and coregistered functional images (see [section 2.5](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#25-normalization))
5. Extraction of a (subject-level) inclusve field-of-view (FOV) voxel mask from the preprocessed functional images (see [section 2.6](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#26-inclusive-mask-extraction))
6. Segmentation of the anatomical image (see [section 2.7](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#27-segmentation))
7. Erosion of the white-matter and CSF segmentations (see [section 2.8](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#28-erosion))

As part of the first-level seed-based connectivity pipeline of the resting-state data, the following steps need to be performed after the preprocessing steps:
1. Construction of a region-of-interest (ROI) mask in subject reference space (see [section 3.5](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#41-reslice-roi-mask))
2. Extraction of the region-of-interest (ROI) regressors (see [section 4.2](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#42-roi-regressors))
3. Extraction of the confound regressors (see [section 4.3](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#43-confound-regressors))
4. Extraction of the spike regressors (see [section 4.4](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#44-spike-regressors))
5. Voxel-wise seed-based functional connectivity analysis (see [section 4.5](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#45-voxel-wise-connectivity-analysis))

### 4.1 Reslice ROI mask

Before the seed-based connectivity analysis can be performed, the region-of-interest (ROI) mask, in this case the amygdala, needs to be resliced to the reference space of each subject in the dataset. This can be achieved in either of two ways: 
- Option 1: The analyses are conducted mainly in native space, and therefore the ROI mask needs to be resliced to the native reference space of each subject; this can be accomplished via the *create_native_roi_mask* method. 
- Option 2: The analyses are conducted mainly in standard MNI space, and therefore the ROI mask needs to be resliced to the standard reference space of each subject; this can be achieved via the *obtain_global_roi_mask* method.

The *create_native_roi_mask* method ultilizes [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the procedure itself are [GlobalToNativeMask.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/GlobalToNativeMask.m) and [GlobalToNativeMask_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/GlobalToNativeMask_job.m). This method inverse normalizes the ROI mask specified, via the following 5-step procedure:
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
2. The name of the ROI mask used for the analysis, as listed in the working directory. A binary mask or probability map of the amygdala should be used. 
3. An optional prefix to indicate the exact scan identifier of the anatomical image on which to base the inverse normalization. It is recommended that the coregistered (to the mean functional image) anatomical image is selected for this procedure. Enter c for the coregistered anatomical image.

The result of this procedure is a set of files and images that are written to the T1 directory of each subject in the data directory. The native ROI mask used for further analysis steps is indicated by the prefix **cic**, followed by the original filename of the ROI mask; e.g. cicAmygdala_total_binary_mask_thr_0_5.nii.

The *obtain_global_roi_mask* method also utilizes [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts that are responsible for the procedure itself are called [CoregisterGlobalMask.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/CoregisterGlobalMask.m) and [CoregisterGlobalMask_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/CoregisterGlobalMask_job.m). The method reslices the standard ROI mask to the standard reference space of each subject (i.e., with the same voxel size).

To reslice the ROI mask to standard reference space, enter the following code in a console:

```
my_experiment = Amy.EmotionTask()
my_experiment.obtain_global_roi_mask()
```

Since the EmotionTask subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter REST for the resting-state data.
2. The name of the ROI mask used for the analyses, as listed in the working directory. A binary mask or probability map of the amygdala should be used. 
3. An optional prefix to indicate the exact scan identifier of the anatomical image to base the inverse normalization on. It is recommended that the normalized anatomical image is used for this procedure. Enter n for the normalized anatomical image.

The output of this procedure is a resliced version of the ROI mask in standard reference space for each subject. The resliced ROI mask used for further analysis steps is indicated by the prefix **r** followed by the original filename of the ROI mask; e.g. rAmygdala_total_binary_mask_thr_0_5.nii.

### 4.2 ROI regressors

Before the connectivity analysis of the resting-state data can be conducted, the ROI regressors (left/right hemisphere) need to be extracted from the data. The ROI regressors are extracted using the *extract_roi_regressors* method of the RestingState class. This method computes one ROI regressor per hemisphere, based on a pre-specified ROI mask, which can be either probabilistic or binary.

To execute the *extract_roi_regressors* method, enter the following code in the console:

```
my_experiment = Amy.RestingState()
my_experiment.extract_roi_regressors()
```

Since the RestingState subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter REST for the resting-state data.
2. The name of the ROI mask used for the analysis, as listed in the working directory. A binary mask or probability map of the amygdala should be used. 
3. The nature of the ROI mask specified in step 2. Enter 1 for a binary ROI mask, or 2 for a probabilistic one.
4. An optional prefix to indicate what version of the ROI mask should be used for the analyses. It is recommended that the output file of section 4.1 (e.g. rAmygdala_total_binary_mask_thr_0_5.nii, or cicAmygdala_total_binary_mask_thr_0_5.nii) is used. Enter r for the mask in standard MNI (subject) space, or cic for the mask in native (subject) space
5. An optional prefix to indicate the exact functional scan identifiers on which the analysis needs to be performed. It is recommended that the normalized (coregistered) and realigned functional images are used at this stage. Enter nr for the normalized (coregistered) and realigned functional scans.

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
4. An optional prefix to indicate the exact functional scan identifiers on which the analysis will to be performed. It is recommended that the normalized (coregistered) and realigned functional images are used at this stage. Enter nr for the normalized (coregistered) and realigned functional scans.

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

The voxel-wise seed-based functional connectivity analysis of the resting-state data is conducted via the *voxel_wise_connectivity_analysis* method. This method performs two GLMs for each voxel within the grey matter mask of the brain. In the first GLM, the signal of the left ROI regressor is entered as predictor-of-interest, along with the 36 confound model (see section 4.3) and spike regressors (see section 4.4) as covariates-of-no-interest, and the signal of each voxel (individually) as outcome variable. In the second GLM, the signal of the right ROI regressor is entered as predictor-of-interest, again with the 36 confound model (see section 4.3) and spike regressors (see section 4.4) as covariates-of-no-interest, and the signal of each voxel (individually) as outcome variable. 

The GLM (beta) coefficients of the predictors-of-interest (left/right ROI) are extracted by the method and converted to t-maps, both of which are stored as NIFTI files that can be found in the resting-state directory of each subject (e.g. **connectivity_b_map_lh_xm13101101.nii** and **connectivity_b_map_rh_xm13101101.nii**,  alongside **connectivity_t_map_rh_xm13101101.nii** and **connectivity_t_map_rh_xm13101101.nii**). The b-values contained within the b-maps represent a quantification of the relationship of the left/right ROI signals, in this case the left and right amygdalae, with the signal of each individual gray matter voxel, corrected for the nuisance effects of motion, white-matter, CSF, and global signal (36 confound model), with censoring of the frames flagged as framewise displacement outliers (spike regressors). The t-maps contain the b-values corrected for the uncertainty of the estimate by having divided by its own standard error (effectivity converting them into t-values).

To perform the voxel-wise seed-based connectivity analysis, enter the following code in a console:

```
my_experiment = Amy.RestingState()
my_experiment.voxel_wise_connectivity_analysis()
```

Since the RestingState subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter REST for the resting-state data.
2. An optional prefix to indicate the exact scan identifiers on which the connectivity analysis needs to be conducted. It is recommended that the normalized (coregistered) and realigned functional data is used at this stage. Enter nr for the normalized (coregistered) and realigned functional (resting-state) images.

## 5. Post-processing

After the first-level analysis of both the emotion task and resting-state data, a number of post-processing steps need to be performed before the group-level analysis can be conducted. The post-processing procedures are conducted via the Postprocessing subclass. This class inherits all attributes and methods of the Preprocessing class, which itself (in turn) inherits from the main Experiment class.
As before, many of these post-processing methods are conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command.

The following post-processing steps are supported by the pipeline:
1. Creation of an across-subjects study-specific FOV mask (see [section 5.1](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#51-inclusive-fov-mask))
2. Creation of an across-subjects study-specific grey matter mask (see [see section 5.2](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#52-inclusive-grey-matter-mask))
3. Creation of an across-subjects study-specific average (T1) brain (see [see section 5.3](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#53-average-brain))
4. Extraction of QC-FC motion correction bencmark (see [section 5.4](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#54-qc-fc-motion-correction-benchmark))
5. Extraction of discriminability motion correction benchmark (see [section 5.5](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#55-discriminability-motion-correction-benchmark))

### 5.1 Inclusive FOV mask

Extraction of a study-specific inclusive voxel (i.e., FOV) mask is supported by the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module via the *create_study_FOV_mask* method of the Postprocessing class. Enter the following code in the console to execute this process:

```
my_experiment = Amy.Postprocessing()
my_experiment.create_study_FOV_mask()
```

Since the Postprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the post-processing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. A prefix to indicate the exact scan identifiers to base the extraction of the inclusive mask on. Enter FOV_mask_ for the FOV mask images.

The *create_study_FOV_mask* method creates an output NIFTI image in the working directory called (e.g.) working_dir > **MARS_Inclusive_FOV_Mask_REST.nii** or **BETER_Inclusive_FOV_Mask_REST.nii**. This file can be used to mask the input images of the second-level analysis, the procedures of which are described in section 6.

### 5.2 Inclusive grey matter mask

Extraction of a study-specific (inclusive) grey matter mask is supported by the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module via the *create_study_grey_matter_mask* method of the Postprocessing class. Enter the following code in the console to execute this process:

```
my_experiment = Amy.Postprocessing()
my_experiment.create_study_grey_matter_mask()
```

Since the Postprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the post-processing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. An optional prefix to indicate the exact scan identifiers to base the extraction of the inclusive mask on. Enter c1n for grey matter segmentations of the normalized T1 images.

The *create_study_grey_matter_mask* method creates an output NIFTI image in the working directory called (e.g.) working_dir > **MARS_Inclusive_GM_Mask_REST.nii** or **BETER_Inclusive_GM_Mask_REST.nii**. This file can be used to mask the input images of the second-level analysis (e.g. alongside the study-specific FOV mask), the procedures of which are described in section 6.

### 5.3 Average brain

Extraction of a study-specific average T1 image is supported by the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module via the *create_average_brain* method of the Postprocessing class. Enter the following code in the console to execute this process:

```
my_experiment = Amy.Postprocessing()
my_experiment.create_average_brain()
```

Since the Postprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the post-processing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. An optional prefix to indicate the exact scan identifiers to base the average T1 scan on. Enter n for the normalized T1 images.

The *create_average_brain* method creates an output NIFTI image in the working directory called (e.g.) working_dir > **MARS_Average_Brain_REST.nii** or **BETER_Average_Brain_REST.nii**.

> Note: It is recommended that the *create_average_brain* method is applied to T1 scans that have been resliced (i.e., normalized) using a high resolution, e.g., 1-by-1-by-1 voxel size, rather than the lower resolution based on the functional data, generated by the normalization process described in [section 2.5](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#25-normalization)

### 5.4 QC-FC motion correction benchmark

In order to assess the success of the motion correction procedures described above, i.e., the realignment method described in section 2.2, and the confound and spike regression methods described in sections 3.2, 3.3, 4.3, and 4.4, a number of motion correction benchmarks can be extracted from the output of the first-level (connectivity) analysis, as recommended by [Ciric et al. (2017)](https://www.sciencedirect.com/science/article/pii/S1053811917302288) and [Parkes et al. (2018)](https://www.sciencedirect.com/science/article/pii/S1053811917310972). One such motion correction benchmark is the QC-FC correlation map, which is created by a voxel-wise across-subjects correlation analysis that involves the mean framewise displacement data vs. the (voxel-wise) functional connectivity data.

The extraction of this QC-FC motion correction bencmark is supported by the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module via the *motion_correction_benchmark_qc_fc* method of the Postprocessing class. Enter the following code in the console to execute this process:

```
my_experiment = Amy.Postprocessing()
my_experiment.motion_correction_benchmark_qc_fc()
```

Since the Postprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the post-processing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The filename as listed in the working directory of the output textfile, created via the procedures described in section 4.4, that contains the mean framewise displacement (MFD) data and number of flagged outlier frames (nOutliers) across subjects, for the resting-state data (e.g. BETER_Outliers_FD_Jenkinson_REST_log.txt or MARS_Outliers_FD_Jenkinson_REST_log.txt).
3. The filename as listed in the working directory of the study-specific grey matter mask created in section 5.2. Enter the filename of the study-specific grey matter mask you want to use (e.g. MARS_Inclusive_GM_Mask_REST.nii or BETER_Inclusive_GM_Mask_REST.nii).
4. A prefix to indicate the exact scan identifiers to base the motion correction benchmark on. Enter either connectivity_b_map_lh_ or connectivity_b_map_rh_ for the left or right normalized beta maps (respectively).

The *motion_correction_benchmark_qc_fc* method creates a number of output NIFTI images in the working directory that detail, voxel-wise and across subjects, the correlation coefficients between the mean framewise displacement data and functional connectivity values (e.g. working_dir > **BETER_QC-FC_Corr_connectivity_b_map_lh_REST.nii**), associated p-values (e.g. working_dir > **BETER_QC-FC_Prob_connectivity_b_map_lh_REST.nii**), and a binary map that indicates whether or not these p-values exceed an alpha of < 0.05 uncorrected (e.g. working_dir > **BETER_QC-FC_Sig_connectivity_b_map_lh_REST.nii**). The number of voxels with a value of 1 in this latter file, i.e., those with a p-value < 0.05 uncorrected, is also calculated by the method, as well as the number of voxels significant _after_ FDR-correction ([Benjamini & Hochberg, 1995](https://www.jstor.org/stable/pdf/2346101.pdf?casa_token=6sYKUuXNtwIAAAAA:0-4q3v3-IU9zPV4hGp5IOLVdTccoKf_dudfBuT-VZQhT8JJ8-d77z7XwEo9cKL21utTmvp8w8FEytbHdTCkm3wmWoAT6LrRwJj09U7FHu2OUmCAg4WDg)), all of which are written to both the console, and an output text file stored in the working directory (e.g. working_dir > **BETER_QC-FC_NSig_connectivity_b_map_lh_REST.txt**).

### 5.5 Discriminability motion correction benchmark

The second motion correction benchmark that is supported by the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module is a group analysis of the control group, median split according to the mean framewise displacement data. That is to say, the non-psychiatric control subjects are median-split based on the subject-level mean framewise displacement values; a group analysis is then conducted on the connectivity maps. If the motion correction pipeline was succesful in mitigating the effects of head motion, group differences in functional connectivity would be expected to be minimal at this stage.

This so-called discriminability analysis is supported by the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module via the *motion_correction_benchmark_discriminability* method of the Postprocessing class. Enter the following code in the console to execute this process:

```
my_experiment = Amy.Postprocessing()
my_experiment.motion_correction_benchmark_discriminability()
```

Since the Postprocessing subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans on which the post-processing needs to be conducted. Enter REST for resting-state data or EMO for emotion task data.
2. The filename as listed in the working directory of the output textfile, created via the procedures described in section 4.4, that contains the mean framewise displacement (MFD) data and number of flagged outlier frames (nOutliers) across subjects, for the resting-state data (e.g. BETER_Outliers_FD_Jenkinson_REST_log.txt or MARS_Outliers_FD_Jenkinson_REST_log.txt).
3. The filename as listed in the working directory of the study-specific grey matter mask created in section 5.2. Enter the filename of the study-specific grey matter mask you want to use (e.g. MARS_Inclusive_GM_Mask_REST.nii or BETER_Inclusive_GM_Mask_REST.nii).
4. A prefix to indicate the exact scan identifiers to base the motion correction benchmark on. Enter either connectivity_b_map_lh_ or connectivity_b_map_rh_ for the left or right normalized beta maps (respectively).

The *motion_correction_benchmark_discriminability* method creates an output NIFTI image in the working directory that contains the voxel-wise the independent samples t-test scores that assess the difference in mean connectivity values between the high and low motion groups (e.g. working_dir > **BETER_Discriminability_T_Map_connectivity_b_map_lh_REST.nii**). If the number of significant voxels are low at this stage, as indicated by the number of voxels that exceed a given t-threshold at n1 + n2 - 2 degrees of freedom, the motion correction procedures are likely to have been succesful. The number of voxels with a p-value < 0.05 uncorrected, is also calculated by the method, as well as the number of voxels significant _after_ FDR-correction ([Benjamini & Hochberg, 1995](https://www.jstor.org/stable/pdf/2346101.pdf?casa_token=6sYKUuXNtwIAAAAA:0-4q3v3-IU9zPV4hGp5IOLVdTccoKf_dudfBuT-VZQhT8JJ8-d77z7XwEo9cKL21utTmvp8w8FEytbHdTCkm3wmWoAT6LrRwJj09U7FHu2OUmCAg4WDg)), both of which are again stored to both the console, and an output text file written to the working directory (e.g. working_dir > **BETER_Discriminability_NSig_connectivity_b_map_lh_REST.txt**).
> Note that a _one-sided_ t-test is performed at a voxel-level since head motion is expected to inflate (rather than deflate) measures of functional connectivity. This means that an increase in functional connectivity would be expected in the high vs. low motion group if one assumes that motion artefacts were influential, and thus a one-sided independent samples t-test can be used.

## 6. Group-level analysis

The prediction of the emotion task reactivity of the amygdala using the seed-based resting-state connectivity of the amygdala is conducted via the GroupAnalysis subclass of the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module. This class inherits all attributes and methods of the Postprocessing class, which itself (in turn) inherits from the Preprocessing subclass, and by extension, the main Experiment class.

The following group-level analysis steps are supported by the pipeline:
1. Standard group-level analysis of the task reactivity or resting-state connectivity data (see [section 6.1](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#61-second-level-parametric-analysis))
2. Voxel-matched regression analysis of the connectivity vs. reactivity data (see [section 6.2](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#62-voxel-matched-regression-analysis))
2. ROI-based regression analysis of the connectivity vs. reactivity data (see [section 6.3](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#63-roi-based-regression-analysis))

### 6.1 Second-level analysis

Before the resting-state connectivity vs. task reactivity analysis can be conducted, the group-level effects of both the resting-state and task reactivity data can be extracted via standard statistical (second-level) analysis. The second-level analysis method supported by the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module is conducted using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html), via a shell call command. The MATLAB scripts responsible for the second-level procedure are [SecondLevelAnalysis.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/FirstLevelAnalysis.m) and [SecondLevelAnalysis_job.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/FirstLevelAnalysis_job.m). 

To perform the standard parametric second-level analysis, enter the following code in a console:

```
my_experiment = Amy.GroupAnalysis()
my_experiment.run_2nd_level_analysis()
```

Since the GroupAnalysis subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter REST for the resting-state data or EMO for the emotion task data.
2. The filename, as listed in the subjects' first-level output directories, of the statistical map that is to be used for the analysis (without file extension; e.g. con_0001 or con_0002, or connectivity_b_map_lh_ or connectivity_b_map_rh_).
3. An optional explicit mask to base the second-level analysis on. Simply press the enter key to continue or enter the filename of the study-specific voxel mask you want to use (e.g.GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii).

The output beta, contrast, and t-maps of the second-level procedure are stored in a subfolder of the working directory that is named after the study ID (MARS, BETER), and contrast image used for the analysis (con_0001, con_0002, con_0003, con_0004, con_0005); e.g., working_dir > **BETER_second_level_analysis_con_0001**. In this folder can be found the beta-maps of each of the predictors specified in the model, as well as the specified contrast maps and corresponding t-maps.

### 6.2 Voxel-matched regression analysis

Voxel-matched connectivity vs. reactivity regression analysis using landscape-based clustering with permutation testing is supported by the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module via the *voxel_matched_regression_analysis* method, and follows the approach described by [Mennes et al. (2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2839004/). This method conducts a simple linear regression model across participants for each locus in standard MNI space, wherein the connectivity value of each voxel coordinate across subjects is used as predictor-of-interest, and the reactivity value spatially matched to that voxel is used as outcome variable. The beta-map yielded by this approach is then converted into a t-map and subjected to landscape-based cluster analysis with permutation/randomisation testing, following the procedures described by [Gladwin, Vink, and Mars (2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4950168/).

The voxel-matched regression analysis method of the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module is conducted via a shell call command in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html). The MATLAB script responsible for the voxel-matched regression procedures is [VoxelMatchedRegressionAnalysis.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/FirstLevelAnalysis.m).

Enter the following code in the console to execute this process:

```
my_experiment = Amy.GroupAnalysis()
my_experiment.voxel_matched_regression_analysis()
```

Since the GroupAnalysis subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter REST for the resting-state data or EMO for the emotion task data.
2. The name of a text file listed in the working directory that contains the list of subjects that are to be included in the analysis (e.g. voxel_matched_regression_analysis.txt). Note that this file needs to be manually prepared beforehand.
3. A prefix to indicate the exact scan identifiers to base the voxel-matched regression analyses on; i.e., on the predictor side of the equation. Enter either connectivity_t_map_lh_ or connectivity_map_t_rh_ for the left or right t-standardized connectivity maps (respectively).
4. The filename, as listed in the subjects' first-level output directories, of the contrast images used for the voxel-matched regression analyses; i.e., on the outcome side of the equation (e.g. spmT_0001.nii, spmT_0002.nii, spmT_0003.nii, spmT_0004.nii, spmT_0005.nii).
5. A **non-optional** explicit mask file to use for the voxel-matched regression analysis, as listed in the working directory (e.g. GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii).
6. The full-width half-maximum (FWHM) at which the input connectivity and reactivity maps are to be smoothed, prior to the actual voxel-matched regression analysis. It is recommended that a FWHM of 6 cubic mm is used at this stage.
7. The number of permutation/randomisations that are to be conducted during the landscape-based cluster analysis. It is recommended that 2000 permutations are entered at this stage.

The *voxel_matched_regression_analysis* method creates a new subdirectory within the main working directory, which is named after both the connectivity prefix (entered in step 3) and reactivity contrast image (entered in step 4) used in the analysis, and which contains the output files and images of the landscape-based cluster analysis (e.g. working_dir > **voxel_matched_regression_connectivity_t_map_lh_vs_spmT_0001**). The *R.mat* file written to this directory is the essential output yielded by the analysis; the other files and images contained by the directory are derived from this main output file.

### 6.3 ROI-based regression analysis

The ROI-based regression analysis ia a variation on the main voxel-matched regression method described in section 6.2, which in itself is based on the approach originally documented by [Mennes et al. (2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2839004/). Rather than matching the connectivity value of each locus to the reactivity value at that corresponding coordinate in standard MNI space, the ROI-based regression analysis conducts a simple linear regression model for each voxel, wherein the connectivity value of that voxel across subjects is still used as predictor-of-interest (as in section 6.2), but rather than matching the reactivity value of the corresponding locus in MNI space (as in section 6.2), the same across-subjects vector of mean reactivity values of a given ROI (the output of [section 3.6](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/README.md#36-roi-masked-spm-data)) is used as outcome variable for each voxel, and thus for each regression model. The beta-map yielded by this approach is again converted into a t-map and subjected to landscape-based cluster analysis with permutation/randomisation testing, following the procedures described by [Gladwin, Vink, and Mars (2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4950168/).

The ROI-based connectivity vs. reactivity regression analysis with landscape-based clustering using permutation testing is supported by the [amygdala_recon.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_recon.py) module via the *roi_based_regression_analysis* method. This method uses a shell call command in [MATLAB R2016b](https://nl.mathworks.com/products/matlab.html); the MATLAB script responsible for the ROI-matched regression procedure is called [ROIMatchedRegressionAnalysis.m](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/FirstLevelAnalysis.m). Note that the method assumes that a text file is prepared manually beforehand, in which the subject identifiers to be included in the analysis, are listed alongside the ROI reactivity values for the left hemisphere and right hemisphere, respectively (in that order; with no colunn headers).

Enter the following code in the console to execute this process:

```
my_experiment = Amy.GroupAnalysis()
my_experiment.roi_matched_regression_analysis()
```

Since the GroupAnalysis subclass inherits from the main Experiment class, the same three user inputs as described above (see section 1) need to be entered. The program will ask for the following additional inputs to be specified in the console:
1. The type of scans to be used for the analysis (this input-dependent attribute is inherited from the \_\_init__ method of the Preprocessing class). Enter REST for the resting-state data or EMO for the emotion task data.
2. A prefix to indicate the exact scan identifiers to base the voxel-matched regression analyses on; i.e., on the predictor side of the equation. Enter either connectivity_t_map_lh_ or connectivity_t_map_rh_ for the left or right normalized connectivity maps (respectively).
3. The filename, as listed in the subjects' first-level output directories, of the contrast images used for the ROI-matched regression analyses; i.e., on the outcome side of the equation (e.g. spmT_0001.nii, spmT_0002.nii, spmT_0003.nii, spmT_0004.nii, spmT_0005.nii). 
> Note that a text file needs to be manually prepared before the method is initialized, in which the subject identifiers are listed alongside the ROI reactivity values of the left and right hemisphere, respectively (in that order; with no colunn headers), and this file needs to be named such that the name of the input contrast image corresponds to the filename of the text file that is used an an implicit input by the method, in the following manner: **ROI_SpmT_0001_based_regression_analysis.txt** (using spmT_0001.nii as an example). The name of the input contrast image thus needs to correspond to the filename of the text file that is used an an implicit input by the method.
4. The hemisphere to base the reactivity values on, as extracted from the text file noted in the text above. Enter 1 for the left hemisphere or 2 for the right hemisphere.
> Note that the method assumes that the text file that is used as an implicit input for the analysis contains one column of subject identifiers, followed by one column of ROI reactivity values for the left hemisphere for that contrast, and one column of ROI reactivity values for the right hemisphere for that same contrast. These columns of data need to be entered into the text file in that specific order: subject identifiers, followed by reactivity values of the left hemisphere, followed by reactivity values of the right hemisphere, with no column headers.
5. A **non-optional** explicit mask file to use for the voxel-matched regression analysis, as listed in the working directory (e.g. GRAND_Inclusive_FOV_GM_Mask_REST_EMO.nii).
6. The full-width half-maximum (FWHM) at which the input connectivity and reactivity maps are to be smoothed, prior to the actual voxel-matched regression analysis. It is recommended that a FWHM of 6 cubic mm is used at this stage.
7. The number of permutation/randomisations that are to be conducted during the landscape-based cluster analysis. It is recommended that 2000 permutations are entered at this stage.

The *roi_based_regression_analysis* method creates a new subdirectory within the main working directory, which is named after the connectivity prefix (entered in step 2), reactivity contrast image (entered in step 3), and hemisphere of reactivity values (entered in step 4) used for the analysis, and which contains the output files and images of the landscape-based cluster analysis (e.g. working_dir > **roi_based_regression_connectivity_t_map_lh_vs_spmT_0001_HemiL**). The *R.mat* file written to this directory is the essential output yielded by the analysis; the other files and images contained by the directory are derived from this main output file.


*Tim Varkevisser*
