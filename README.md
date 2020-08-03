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
   - R >= 3.6 (needed packages will be installed at first run)

Next, clone RTrans GitHub repository:
`
git clone https://github.com/gskrasnov/RTrans.git
`


###  To run RTrans you need:
   - derive read counts per gene files, which are generated with HTSeq-count, RSEM or featureCount. By default, *.counts files should be placed to 'counts' folder.
   - set parameters. See the example in the 'RTrans parameters example.xlsx'. Provide this file to "Startup.Data(...)" method
   - install python3 with xlsxwriter module to create Excel reports


