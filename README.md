# group_23_project

# Analysis of Breast Cancer Proteomes
Final project for the R for Bio-Data Science course (22160) - Group 23

## Data

The original, raw data used for this project comes from the following study:

Mertins, P., Mani, D., Ruggles, K. et al. Proteogenomics connects somatic mutations to signalling in breast cancer. Nature 534, 55–62 (2016)
https://doi.org/10.1038/nature18003

In the study, the influence of DNA mutations on protein expression in breast cancer was examined.

For this project the following files were used:

-   **77_cancer_proteomes_CPTAC_itraq.csv**: includes iTRAQ proteome profiling of 77 breast cancer samples generated by the Clinical Proteomic Tumor Analysis Consortium (NCI/NIH). The file contains a collection of proteins, with their respective expression levels, found across different breast cancer patients 
    
-   **PAM50_proteins.csv**: PAM50 classification system (relevant for further analysis)

-   **clinical_data_breast_cancer.csv**: clinical information related to the samples
    
### Data download 

Data is available in [Kaggle](https://www.kaggle.com/datasets/piotrgrabo/breastcancerproteomes/data?select=77_cancer_proteomes_CPTAC_itraq.csv), or can be downloaded directly from [GitHub](https://github.com/BCPP/BreastCancerProteomes/tree/master).

The data files were not included in this GitHub repository, as specified in the project requirements. As the data could not be retrieved automatically from any web server, if you wish to execute this project end-to-end, you must download the corresponding files and add them to the `/data/_raw/` directory.

## Executing the analysis 

The R folder contains all the .qmd scripts created for this project, which break down the analysis into easy-to-follow steps. To execute the entire analysis in one go, please run `00_all.qmd`.

## Contributors

Ana Pastor Mediavilla, @AnaPasMed
Amanda Jimenez Lobato, @AmandaJimLob99
Carlos de Santiago León, @CarlosSanti00
Laura Figueiredo Tor, @LFT18
Monika Karolina Borkowska, @mondaymon
