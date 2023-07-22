# Reproduction of Data in [Cascarina (2023)](add link when published)

Important disclaimer: all code and resulting data have not yet been peer reviewed and are subject to change.

This directory contains all necessary information and code to reproduce the data in [Cascarina (2023)](add link when published). Below is a list of dependencies and associated version numbers used for testing and data analysis:

| Package | Version |
| ----------- | ----------- |
| Python | 3.7.15 | 
| Numpy | 1.21.5 |
| Scipy | 1.7.3 |
| Pandas | 1.3.5 |
| Matplotlib | 3.5.3 |
| Seaborn | 0.12.0 |
| Biopython | 1.79 |
| goatools | 1.0.2 |
| tqdm | 4.64.1 |
| upsetplot | 0.8.0 |

To reproduce figures and tables in the paper:
1. Download the required data from Zenodo (NOTE: this is currently only available as a privately shared link provided in the cover letter for peer review purposes. If you do not have this link, the code shared here will fail).
2. In a single, main directory, extract all files from zipped archives in separate directories with names identical to those of the zipped archives. For example, the files contained in "Archaea.zip" should be unzipped into a folder called "Archaea". At the end of this process, the following directories should exist and should be kept in a single location:
    - Archaea
    - Archaea_SCRAMBLED
    - Bacteria
    - Bacteria_SCRAMBLED
    - Eukaryota
    - Eukaryota_SCRAMBLED
    - Viruses
    - Viruses_SCRAMBLED
    - GOAfiles
    - SecondaryLCDs_by_LCDcategory
    - Pfam_Data
    - Observed_vs_Scrambled_LCDfrequency_Statistics
    - GOterm_Results_C-rich_LCDs_ModelOrganisms
3. Download all files from this Github repository. Extract all files into the main directory from Step 2: unlike the extracted files above, these files should not exist in a separate directory.
4. Run the "Reproducibility_BATCH.bat" file from the command line to run all analyses in succession.
    - NOTE: Some of the commands in the batch file may take up to 1 day to run. The first two commands in particular may take ~15-17hrs on a basic desktop computer. However, these files can be run concurrently, so you can run them in parallel and delete them in the batch file prior to running the batch file if preferred. These two commands must finish execution before running the remaining commands in the batch file.

## License info
All code in this directory is subject to the terms of the GPLv3 license (see LICENSE in this directory). If you use data, code, or information contained in this directory or the associated manuscript, please cite:

INCLUDE REFERENCE WHEN PUBLISHED
