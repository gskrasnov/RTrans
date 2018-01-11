
This scrpit is aimed at the analysis of RNA-Seq data:
   - Differential expression: GLM multivariate testing (multiple groups and variates) or t-test (two groups)
   - Two types of Gene Ontology, KEGG, Reactome enrichment analyses
   - Visualization of KEGG pathways (MAPK, PI3K/mTOR, p53, etc.)
   - Pathway-centric differential expression profiles
   - Detailed analysis of pre-defined functional groups of genes
   - Creating heatmaps, PCA plots, Excel reports etc.

###  To run RTrans you need:
   - derive read counts per gene files, which are generated with HTSeq-count, RSEM or featureCount. By default, *.counts files should be placed to 'counts' folder.
   - set parameters. See the example in the 'RTrans parameters example.xlsx'. Provide this file to "Startup.Data(...)" method
   - install python3 with xlsxwriter module to create Excel reports


