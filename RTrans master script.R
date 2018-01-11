####################        RTrans        #######################################################
##
##  This scrpit is aimed at the analysis of RNA-Seq data:
##   - Differential expression: GLM multivariate testing (multiple groups and variates) or t-test (two groups)
##   - Two types of Gene Ontology, KEGG, Reactome enrichment analyses
##   - Visualization of KEGG pathways (MAPK, PI3K/mTOR, p53, etc.)
##   - Pathway-centric differential expression profiles
##   - Detailed analysis of pre-defined functional groups of genes
##   - Creating heatmaps, PCA plots, Excel reports etc.
##
##  To run RTrans you need:
##   - derive read counts per gene files, which are generated with HTSeq-count or featureCount. By default, *.counts files should be placed to 'counts' folder.
##   - set parameters. See the example in the 'RTrans parameters example.xlsx'. Provide this file to "Startup.Data(...)" method
##   - install python3 with xlsxwriter module to create Excel reports
##
#################################################################################################

# loading RTrans functions and objects
if (!('rstudioapi' %in% rownames(installed.packages()))) install.packages("rstudioapi")
source(sort(list.files(path = sprintf("%s/RTrans.base",dirname(rstudioapi::getActiveDocumentContext()$path)),pattern = "RTrans.methods.*R",full.names=T),decreasing = T)[1])

# Loading parameters from Excel workbook and preprocessing data
Startup.Data = Prepare(Parameters.xlsx.file.name = '{script_dir}/RTrans.parameters for Prassolov.xlsx')
GLM.model.list = Startup.Data$Pars$Models.to.Test
GLM.model = GLM.model.list[1]  # for debug only

for (GLM.model in GLM.model.list){
  Analysis.Data = Analyze.GLM(Startup.Data, GLM.model)
  Create.Heatmaps(Startup.Data, GLM.model = GLM.model)
  inDetails(Startup.Data, GLM.model = GLM.model)

  KEGG.classic.Enrichment(Startup.Data, GLM.model = GLM.model)
  KEGG.trends.Enrichment(Startup.Data, GLM.model = GLM.model)
  KEGG.Pathways.Visualization(Startup.Data, GLM.model = GLM.model)
  KEGG.Expression.Profiles.Custom(Startup.Data, GLM.model = GLM.model)
  
  Reactome.classic.Enrichment(Startup.Data, GLM.model = GLM.model)
  Reactome.trends.Enrichment(Startup.Data, GLM.model = GLM.model)
  Reactome.Expression.Profiles.Custom(Startup.Data, GLM.model = GLM.model)

  GO.classic.Enrichment(Startup.Data, GLM.model = GLM.model, GO.types = c('BP','CC','MF'))
  GO.trends.Enrichment(Startup.Data, GLM.model = GLM.model, GO.types = c('BP','CC','MF'))
  #### GO.Expression.Profiles.Custom(Startup.Data, GLM.model = GLM.model)
  
  topGO.Enrichment(Startup.Data, GLM.model = GLM.model, GO.types = c('BP'))
  topGO.Expression.Profiles.Custom(Startup.Data, GLM.model = GLM.model, GO.type = 'BP')
  topGO.Expression.Profiles.Enriched(Startup.Data, GLM.model = GLM.model, GO.type = 'BP')
}


# Summarize.GLM.results(Startup.Data, GLM.models = GLM.model.list, out.dir = )
# Summarize.inDetails.GLM.results(Startup.Data, GLM.models = GLM.model.list)
# Summarize.KEGG.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.model.list)
# Summarize.KEGG.DE.info(Startup.Data, GLM.models = GLM.model.list)
# Summarize.Reactome.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.model.list)
# Summarize.GO.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.model.list)
# Summarize.topGO.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.model.list)

Summarize.results.by.Analyses.groups(Startup.Data)


################################
################################

# if you want to analyze custom gene set with LogFC and p-value entered from external file, use the following steps:
if (!('rstudioapi' %in% rownames(installed.packages()))) install.packages("rstudioapi")
source(sort(list.files(path = sprintf("%s/RTrans.base",dirname(rstudioapi::getActiveDocumentContext()$path)),pattern = "RTrans.methods.*R",full.names=T),decreasing = T)[1])

Startup.Data = Prepare(Parameters.xlsx.file.name = '{script_dir}/RTrans parameters example.xlsx',ext.DE.data.files = c('ext.DE.data.example.txt'), ext.DE.data.species = 'hsa')

for (f in Startup.Data$ext.DE.data.files){
  Analysis.Data = Process.ext.DE.data(Startup.Data,f)
  # ....... and then run GO, Reactome, KEGG analyses as usual
}
Summarize.GLM.results(Startup.Data)
# ....... run summarizing as usual

