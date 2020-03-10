# AMYGDALA_REACT_VS_CONNECT
Pipeline to predict emotion task reactivity of the amygdala using resting-state connectivity of the amygdala as seed-region.

> Preliminary note: The Python module ([amygdala_project.py](https://github.com/tvarkevi/AMYGDALA_REACT_VS_CONNECT/blob/master/amygdala_project.py)) used here assumes the scripts and analysis files are stored in one and the same working directory. The task and resting-state fMRI data of the experiment needs to have been converted to NIFTI format, and stored in a data directory where the different types of data (T1, resting-state, task data) are stored as subdirectories specified by suffixes that added to subject identifier (e.g. subject001 > subject001_1_1). Finally, an Excel (.xlsx) file needs to be made that specifies, for each subject, the subject identifier, and specific suffixes added to the subject ID to differentiate the subdirectories.

## 1. Setting up the experiment

