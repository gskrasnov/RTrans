# RTrans
## a versatile RNA-Seq gene expression data analysis pipeline

### Main features:
   - Differential expression: intra-group comparison, ANOVA/GLM, non-paramentric tests
   - Gene set enrichment analysis (Gene Ontology, KEGG, Reactome, WikiPathways, Disease Ontology, DisGeNET, Network of Cancer Genes)
   - KEGG pathways visualization
   - Pathway-centric gene expression profiles
   - Detailed analysis of pre-defined functional groups of genes
   - Excel reports with embedded heatmaps and sparklines; heatmaps, MDS/PCA, etc.

## Quick start
### Requirements

Before running RTrans, ensure you have installed:
   - python3 (xlsxwriter, numpy packages are needed)
   - R >= 3.6 (all needed packages will be installed at first run)

### Installing RTrans
Next, clone RTrans GitHub repository:

`git clone https://github.com/gskrasnov/RTrans.git`

move to the RTrans folder:

`cd RTrans`


### Preparing counts data
Before launching RTrans, you need to prepare gene expression data. Typically, these are files with RNA-Seq read counts per gene derived with:
   - HTSeq-counts, featureCounts (Subread package), STAR (during reads maping)
   - eXpress, Salmon, kallisto
The read counts files should containg two columns: 1) Ensembl gene IDs, 2) read counts per gene. File naming should correspond to sample naming (e.g. `sample1.txt`, `sample2.txt`). By default, counts files should be placed to 'counts_data' folder.

Optionally, RTrans accepts read counts matrix (see the Tutorials)

### Specifying RTrans startup parameters

RTrans parameters should be specified in `RTrans.parameters.xlsx` Excel workbook.
Main options can be set in the "Basic parameters", "Sample setup", "Heatmaps annotation" sheets.

##### "Basic parameters" sheet:
Pay a special attention to the parameters:
- *Species* 
- *counts dir* (counts_data by default)
- *counts suffix* (.txt by default)

##### "Sample setup" sheet:
Here you can specify groups to be compared. One column - one comparison. One row - one sample.
**"0"** means control group (e.g. healthy individuals)
**"1"** means target group (e.g. patients)

Sample names should be present at the 1st column

Alternatively, you can specify here ANOVA covariates and describe ANOVA/GLM model ("Models to test" parameter at "Basic parameters" sheet) like this: Age + Gender + Disease

##### "Heatmaps annotation" sheet:
Paste here any sample info (like disease status, age, gender, etc.) to be included in heatmaps. One row - one sample. Sample names should be present at the 1^st^ column and should be the same as in the  "Sample setup" sheet

That's almost all. Now you are ready to laucn RTrans.

###  run RTrans
After filling up parameters \*.xlsx file, you can launch Rtrans. Go to the directory with `RTrans.parameters.xlsx` and `counts_data` folder and type:

` Rscript <path to RTrans.R>`

Rtrans will scan parameteres \*.xlsx file, counts data, determine CPU count, available RAM and adjust the number of threads for parallel execution.

All output results will be placed in the current folder. Typically, one analysis required 1-2 hours (including differential expression, gene set enrichment tests, pathway visualization, plots, etc.)
