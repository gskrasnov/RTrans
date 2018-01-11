setClass("RTransParameters",
         representation("list")
)

setClass("RTransStartupData",
         representation("list")
)

RTransAnalysisData.prototype = vector(mode="list")
RTransAnalysisData.prototype$GLM.analysis.performed = F
RTransAnalysisData.prototype$KEGG.classic.Enrichment.performed = F
RTransAnalysisData.prototype$GO.Enrichment.performed = F
RTransAnalysisData.prototype$Heatmaps.created = F
RTransAnalysisData.prototype$Pathway.Visualization.performed = F
RTransAnalysisData.prototype$GO.Expression.Profiles.Created = F
RTransAnalysisData.prototype$Enriched.Pathways.Classic = c()
RTransAnalysisData.prototype$Enriched.Pathways.Trended = c()
RTransAnalysisData.prototype$Enrich.Results.Summary.by.DB = vector(mode="list")

setClass("RTransAnalysisData",
         representation("list"),
         prototype = prototype(RTransAnalysisData.prototype)
)

sum.mod = function(x){
  x=x[!is.na(x)]; return(length(x[x]))
}


Verify.path = function(path, max.chars = 254){
  path.intact = path
  path.quote = NULL
  if(startsWith(path,"'"))  path.quote = "'"
  if(startsWith(path,'"'))  path.quote = '"'
  if(!is.null(path.quote)){
    path %<>% gsub(path.quote,'', ., fixed = TRUE)
    max.chars %<>% -2
  }
  
  path %<>% gsub('\\','/', ., fixed = TRUE)
  if(startsWith(path,'/')) {
    full.path = path
    full.path.is.provided = TRUE
  } else if(substr(path, 2, 2) == ':') {
    full.path = path
    full.path.is.provided = TRUE
  } else {
    sprintf('%s/%s',getwd(),path) %>% gsub('\\','/', ., fixed = TRUE) -> full.path
    full.path.is.provided = FALSE
  }
  
  if(nchar(full.path) <= max.chars) return(path.intact)
  
  strsplit(full.path, '/', fixed = TRUE)[[1]] -> full.path.split
  tail(full.path.split,1) -> file.name
  head(full.path.split,-1) %>% paste(., collapse = '/') -> dir.name
  
  if(nchar(dir.name) >= max.chars + 4){
    warning(sprintf('Cannot truncate file name "%s" since dir length (%d) exceeds %d symbols', full.path, nchar(dir.name), max.chars + 4))
    return(path.intact)
  }
  
  ext.pos = tail(gregexpr('.', file.name, fixed = TRUE)[[1]], 1)
  if(ext.pos < 0) ext.pos = nchar(file.name)
  file.name.base = substr(file.name, 1, ext.pos - 1)
  file.name.ext = substr(file.name, ext.pos, nchar(file.name))
  max.chars.for.file.name.base = max.chars - nchar(dir.name) - nchar(file.name.ext) - 1
  if(full.path.is.provided){
    path.corrected = sprintf("%s/%s[...]%s", dir.name, substr(file.name.base, 1, max.chars.for.file.name.base - 5), file.name.ext)
  } else {
    path.corrected = sprintf("%s[...]%s", substr(file.name.base, 1, max.chars.for.file.name.base - 5), file.name.ext)
  }
  if(!is.null(path.quote)) path.corrected = sprintf("%s%s%s", path.quote, path.corrected, path.quote)
  return(path.corrected)
}

png.mod <- function(filename, width = 480, height = 480, pointsize = 12, res = NA, ...){
  png(filename = Verify.path(filename), width = width, height = height, pointsize = pointsize, res = res, ...)
}

write.table.mod <- function(x, file = "", append = FALSE, quote = TRUE, sep = " ", ...){
  write.table(x = x, file = Verify.path(file), append = append, quote = quote, sep = sep, ...)
}

Calc.min.high.CPM.samples.from.total.samples.count = function(total.samples,min.group.size=NULL){
  if(is.null(min.group.size)) return (as.integer((sqrt(total.samples) + sqrt(max(0,total.samples-5)))/2*1.5 + total.samples/8))
  if(min.group.size == 8) return (5)
  if(min.group.size == 7) return (5)
  if(min.group.size == 6) return (4)
  if(min.group.size == 5) return (3)
  if(min.group.size == 4) return (3)
  if(min.group.size == 3) return (2)
  if(min.group.size == 2) return (2)
  if(min.group.size == 1) return (1)
  return(as.integer(min(
    ((sqrt(total.samples) + sqrt(max(0,total.samples-5)))/2*1.5 + total.samples/7)*1.2,
    min.group.size/1.45)))
  # return(as.integer(min(
  #   ((sqrt(total.samples)/1.4 + sqrt(max(0,total.samples-5)))/2*1.4 + total.samples/7)*1.2,
  #   min.group.size/1.8)))
}

# load.or.install.packages.std <- function(pkg.list, force.re.install=FALSE){
#   if(force.re.install){
#     install.packages(pkgs = pkg.list, dependencies = TRUE)
#     return()
#   }
#   
#   for (pkg in pkg.list) {
#     if(force.re.install){
#       eval(parse(text = sprintf("install.packages(\"%s\", dep = TRUE)",pkg)))
#     } else {
#       if (!(pkg %in% rownames(installed.packages())))  eval(parse(text = sprintf("install.packages(\"%s\")",pkg)))
#     }
#   }
#   for (pkg in pkg.list) eval(parse(text = sprintf("suppressPackageStartupMessages(library(\"%s\"))",pkg)))
# }
# 
# load.or.install.packages.bio <- function(pkg.list, force.re.install=FALSE){
#   if(force.re.install){
#     source("https://bioconductor.org/biocLite.R")
#     biocLite(pkgs = pkg.list, dependencies = TRUE)
#     return()
#   }
#   
#   if (any(!(pkg.list %in% rownames(installed.packages())))){
#     source("https://bioconductor.org/biocLite.R")
#     for (pkg in pkg.list){
#       if(force.re.install){
#         eval(parse(text = sprintf("biocLite(\"%s\", dep = TRUE)", pkg)))
#       } else {
#         if (!(pkg %in% rownames(installed.packages()))) eval(parse(text = sprintf("biocLite(\"%s\")", pkg)))
#       }
#     }
#   }
#   for (pkg in pkg.list) eval(parse(text = sprintf("suppressPackageStartupMessages(library(\"%s\"))",pkg)))
# }

load.or.install.packages.std <- function(pkg.list, force.re.install=FALSE){
  if(force.re.install){
    install.packages(pkgs = pkg.list, dependencies = TRUE)
    return()
  }
  
  
  pkg.to.install = pkg.list[!(pkg.list %in% rownames(installed.packages()))]
  if(length(pkg.to.install) > 0){
    install.packages(pkg.to.install)
  }
  
  for (pkg in pkg.list) {
    if(force.re.install){
      eval(parse(text = sprintf("install.packages(\"%s\", dep = TRUE)",pkg)))
    } else {
      if (!(pkg %in% rownames(installed.packages())))  eval(parse(text = sprintf("install.packages(\"%s\")",pkg)))
    }
  }
  for (pkg in pkg.list) eval(parse(text = sprintf("suppressPackageStartupMessages(library(\"%s\"))",pkg)))
}

load.or.install.packages.bio <- function(pkg.list, force.re.install=FALSE){
  if(force.re.install){
    source("https://bioconductor.org/biocLite.R")
    biocLite(pkgs = pkg.list, dependencies = TRUE)
    return()
  }
  
  if (any(!(pkg.list %in% rownames(installed.packages())))){
    source("https://bioconductor.org/biocLite.R")
    for (pkg in pkg.list){
      if(force.re.install){
        eval(parse(text = sprintf("biocLite(\"%s\", dep = TRUE)", pkg)))
      } else {
        if (!(pkg %in% rownames(installed.packages()))) eval(parse(text = sprintf("biocLite(\"%s\")", pkg)))
      }
    }
  }
  for (pkg in pkg.list) eval(parse(text = sprintf("suppressPackageStartupMessages(library(\"%s\"))",pkg)))
}


Extract.number.from.text <- function (x,Species) {
  comment.start = regexpr("\\#",x)[1]
  if (comment.start > 0) x =substring(x,1,comment.start-1)
  x = gsub("\\;","",gsub("\\,","",gsub("\t","",gsub(" ","",gsub(Species,"",x)))))
  return(x)
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

Cleanup.Model.Name = function(x){
  x %<>% gsub("\\:","(d)",.) %>% gsub("\\*","(m)",.) %>% gsub("\\","-",.,fixed = TRUE) %>% gsub("/","-",.,fixed = TRUE) %>% gsub("|","-",.,fixed = TRUE)
  return(x)
}


Weighted.Canberra.dist = function(x,y,min.value.to.calc.dist = 30,Canberra.weight.power = 0.4){
  if (min.value.to.calc.dist > 0){
    keep = apply(cbind(x,y),1,max) > min.value.to.calc.dist
    x.filtered = x[keep]
    y.filtered = y[keep]
    sums.counts = x.filtered + y.filtered
    weights = sums.counts^(Canberra.weight.power - 1)
    return (sum(abs(x.filtered - y.filtered)*weights)/sum(weights))
  } else {
    sums.counts = x + y
    weights = sums.counts^(Canberra.weight.power - 1)
    return (sum(abs(x - y)*weights)/sum(weights))
  }
}

Get.KEGG.pathway.list.from.file <- function(file.name){
  conn <- file(file.name,open="r")
  KEGG.lines = readLines(conn)
  close(conn)
  
  KEGG.lines = KEGG.lines[!(KEGG.lines %in% "")]
  for (n in 1:length(KEGG.lines)){
    KEGG.lines[n] = Extract.number.from.text(KEGG.lines[n])
  }
  KEGG.lines = KEGG.lines[!(KEGG.lines %in% "")]
  KEGG.lines = as.numeric(KEGG.lines)
  KEGG.lines = KEGG.lines[!(is.na(KEGG.lines))]
  return(KEGG.lines)
}

Eliminate.Redundant.Predictors.Info = function(predictor.values){
  
  #predictor.values = Analysis.Data$predictor.values
  if(is.vector(predictor.values)){
    return(as.data.frame(predictor.values))
  }
  
  na.preds = apply(predictor.values, MARGIN = 2, FUN = function(x){ all(is.na(x))  })
  if(sum(!na.preds) == 1){
    tmp = as.data.frame(predictor.values[,!na.preds])
    colnames(tmp) = colnames(predictor.values)[!na.preds]
    return(tmp)
  }
  predictor.values = predictor.values[,!na.preds]
  
  na.counts = apply(predictor.values, MARGIN = 2, FUN = function(x){ sum(is.na(x))  })
  predictor.values = predictor.values[,order(na.counts)]
  
  if(dim(as.data.frame(predictor.values))[2] == 1){
    return(as.data.frame(predictor.values))
  }
  
  pred.count = dim(predictor.values)[2]
  res.array = data.frame(array(dim=c(0,pred.count)))
  col_nx = 1
  col_ny = 2
  for (col_nx in 1:pred.count){
    col.data_x = predictor.values[,col_nx]
    res = apply(predictor.values, 2, function(col.data_y){
      #print(col.data_y)
      all(apply(cbind(col.data_x, col.data_y),1, function(pred.pair){
        x = pred.pair[1]
        y = pred.pair[2]
        if(is.na(x) && is.na(y)) return(TRUE)
        if(is.na(x) && !is.na(y)) return(FALSE)
        if(!is.na(x) && is.na(y)) return(TRUE)
        return(x==y)
      }))
    })
    res[col_nx] = FALSE
    res.array = rbind(res.array,res)
    #print(res)
  }
  
  colnames(res.array) = colnames(predictor.values)
  rownames(res.array) = colnames(predictor.values)
  
  repeat{
    omitted.yes = FALSE
    ny = 1
    for(ny in 1:ncol(res.array)){
      keep = !res.array[ny,]
      if(any(!keep)){
        if(sum(keep) == 1){
          return(predictor.values[,colnames(res.array)[keep]])
        }
        res.array = res.array[keep, keep]
        #res.array = res.array[, keep]
        omitted.yes = TRUE
        break
      }
    }
    if (!omitted.yes) break
  }
  
  return(predictor.values[,colnames(res.array)])
}


mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }

  dict = pattern
  names(dict) = replacement
  dict = rev(dict[order(nchar(dict), dict)])
  pattern = dict
  names(pattern) = NULL
  replacement = names(dict)

  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, fixed = TRUE)
  }
  return(result)
}


Replace.Preds.with.Abbrebiations <- function(text, Pars){
  dict = rev(Pars$abbreviations_to_available.preds[order(nchar(Pars$abbreviations_to_available.preds), Pars$abbreviations_to_available.preds)])
  return(mgsub(dict,names(dict),text))
}
Revert.Preds.from.Abbrebiations <- function(text, Pars){
  dict = rev(Pars$available.preds_to_abbreviations[order(nchar(Pars$available.preds_to_abbreviations), Pars$available.preds_to_abbreviations)])
  return(mgsub(dict, names(dict),text))
}


Get.Default.Parameters = function(){
  Default.Parameters = new('RTransParameters')
  
  # General parameters
  Default.Parameters[['DE.package']] = 'edgeR'
  Default.Parameters[['Models.to.Test']] = NA
  Default.Parameters[['Species']] = 'hsa'
  Default.Parameters[['add.MDS.dims.as.predictors']] = FALSE
  Default.Parameters[['Normalize.predictors']] = TRUE
  Default.Parameters[['results.dir']] = '{script_dir}' 
  Default.Parameters[['counts.dir']] = '{script_dir}/counts'
  Default.Parameters[['suppl.data.dir']] = '{script_dir}/RTrans.base'
  Default.Parameters[['Create.Excel.results']] = TRUE
  Default.Parameters[['Bypass.Completed.steps']] = TRUE
  Default.Parameters[['RNA.Seq.norm.method']] = 'TMM'
  Default.Parameters[['min.samples.with.sufficient.CPM']] = '{auto}'
  Default.Parameters[['sufficient.CPM.to.analyze']] = 2
  Default.Parameters[['use.exact.test.in.binary.predictors']] = TRUE
  Default.Parameters[['use.QLfit']] = TRUE
  
  
  # Calculating distance matrices
  Default.Parameters[['Create.Distance.matrices']] = '{auto}'
  Default.Parameters[['Distance.method']] = 'canberra.weighted'
  Default.Parameters[['Minkowski.power']] = 0.5
  Default.Parameters[['Canberra.weight.power']] = 0.4
  Default.Parameters[['Min.CPM.to.include.in.dist']] = 25
  
  #Heatmaps
  Default.Parameters[['Create.Heatmaps']] =  TRUE
  Default.Parameters[['Top.genes.to.include.in.heatmaps.list']] = c(50,500,5000,25000)
  
  # Gene ontology gene set enrichments analysis (GSEA)
  Default.Parameters[['Perform.GO.Enrich__with.topGO']] = TRUE
  Default.Parameters[['GO.Enrich__with.topGO___mode.list']] = c('upreg', 'downreg')
  Default.Parameters[['GO.Enrich__with.topGO___max.DE.genes.list']] = c(40, 80, 250, 500, 1000, 2000)
  Default.Parameters[['GO.Enrich__with.topGO___gene.min.Score.threshold']] = 0
  Default.Parameters[['GO.Enrich__with.topGO___gene.max.PValue.threshold']] = 0.05
  Default.Parameters[['GO.Enrich__with.topGO___max.Genes.in.term.to.show']] = 100
  Default.Parameters[['GO.Enrich__with.topGO___gene.min.abs.LogFC.threshold']] = 0.3
  
  # Parameters of creating gene expression profiles for Enriched GO terms + Custom GO terms (see 'GO terms (DE profiles)' sheet)
  Default.Parameters[['Create.topGO.Expression.Profiles']] = TRUE
  Default.Parameters[['topGO.Expression.Profiles___Max.gene.PValue.list']] = c(0.01, 0.05, 1.00)
  Default.Parameters[['topGO.Expression.Profiles___Min.gene.logCPM.list']] = c(0, 3, 5, 7)
  Default.Parameters[['topGO.Expression.Profiles___Max.gene.PValue.for.summary.across.models']] = 1
  Default.Parameters[['topGO.Expression.Profiles___Min.gene.logCPM.for.summary.across.models']] = 3 
  Default.Parameters[['topGO.Expression.Profiles___Sort.LogFCs']] = T
  Default.Parameters[['topGO.Expression.Profiles___Maximal.genes.in.terms.to.visualize']] = 120
  Default.Parameters[['topGO.Expression.Profiles___remove.redundant.GO.terms']] = TRUE
  
  # Additionally, parameters for the selection of enriched GO terms
  Default.Parameters[['topGO.Enriched.Expression.Profiles___Max.terms.to.visualize']] = 400
  Default.Parameters[['topGO.Enriched.Expression.Profiles___Max.Pvalue.threshold']] = 0.01
  Default.Parameters[['topGO.Enriched.Expression.Profiles___Max.FDR.threshold']] = 1
  Default.Parameters[['topGO.Enriched.Expression.Profiles___scoring.adjustment.to.gene.list.size']] = T
  Default.Parameters[['topGO.Enriched.Expression.Profiles___scoring.adj.list.range.optimal.size']] = 100
  Default.Parameters[['topGO.Enriched.Expression.Profiles___scoring.adj.list.range.start']] = 40
  Default.Parameters[['topGO.Enriched.Expression.Profiles___scoring.adj.list.range.end']] = 350
  Default.Parameters[['topGO.Enriched.Expression.Profiles___limit.Pvalues.to']] = 1e-30
  Default.Parameters[['topGO.Enriched.Expression.Profiles___sum.scores.with.power']] = 2.8
  Default.Parameters[['topGO.Enriched.Expression.Profiles___minimal.score']] = 0
  Default.Parameters[['topGO.Enriched.Expression.Profiles___desired.GO.terms.count__to.discard.score']] = 40
  Default.Parameters[['topGO.Enriched.Expression.Profiles___minimal.term.size']] = 5
  Default.Parameters[['topGO.Enriched.Expression.Profiles___sort.terms.by']] = 'name'
  
  
  # GSEA with clusterProfiler
  Default.Parameters[['Perform.classic.Enrich__with.clusterProfiler']] = TRUE
  Default.Parameters[['Perform.GO.classic.Enrich__with.clusterProfiler']] = TRUE
  Default.Parameters[['Perform.KEGG.classic.Enrich__with.clusterProfiler']] = TRUE
  Default.Parameters[['Perform.Reactome.classic.Enrich__with.clusterProfiler']] = TRUE
  
  Default.Parameters[['classic.Enrich__with.clusterProfiler___mode.list']] = c('upreg', 'downreg')
  Default.Parameters[['classic.Enrich__with.clusterProfiler___max.DE.genes.list']] = c(40, 80, 250, 500)
  Default.Parameters[['classic.Enrich__with.clusterProfiler___gene.min.Score.threshold']] = 0
  Default.Parameters[['classic.Enrich__with.clusterProfiler___gene.max.PValue.threshold']] = 0.05
  Default.Parameters[['classic.Enrich__with.clusterProfiler___gene.min.abs.LogFC.threshold']] = 0.3
  Default.Parameters[['classic.Enrich__with.clusterProfiler___term.max.Pvalue.threshold']] = 0.05
  Default.Parameters[['classic.Enrich__with.clusterProfiler___term.max.FDR.threshold']] = 1.0
  Default.Parameters[['classic.Enrich__with.clusterProfiler___term.max.Qvalue.threshold']] = 1.0
  Default.Parameters[['classic.Enrich__with.clusterProfiler___generate.Summary.Plot']] = TRUE
  
  Default.Parameters[['trends.Enrich__with.clusterProfiler___term.max.Pvalue.threshold']] = 0.05
  Default.Parameters[['trends.Enrich__with.clusterProfiler___term.max.FDR.threshold']] = 1.0
  Default.Parameters[['trends.Enrich__with.clusterProfiler___term.max.Qvalue.threshold']] = 1.0
  
  # Parameters of creating gene expression profiles for Enriched GO terms + Custom GO terms (see 'GO terms (DE profiles)' sheet)
  Default.Parameters[['clusterProfiler.Expression.Profiles___Max.gene.PValue.list']] = c(0.01, 0.05, 1.00)
  Default.Parameters[['clusterProfiler.Expression.Profiles___Min.gene.logCPM.list']] = c(0, 3, 5, 7)
  Default.Parameters[['clusterProfiler.Expression.Profiles___Max.gene.PValue.for.summary.across.models']] = 1
  Default.Parameters[['clusterProfiler.Expression.Profiles___Min.gene.logCPM.for.summary.across.models']] = 3 
  Default.Parameters[['clusterProfiler.Expression.Profiles___Sort.LogFCs']] = T
  Default.Parameters[['clusterProfiler.Expression.Profiles___Maximal.genes.in.terms.to.visualize']] = 120
  Default.Parameters[['clusterProfiler.Expression.Profiles___remove.redundant.entries']] = TRUE
  
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.terms.to.visualize']] = 300
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___balance.between.Classic.and.Trends.Enrich']] = 0.5
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.Pvalue.threshold__classic.test']] = 0.01
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.Pvalue.threshold__trends.test']] = 0.05
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.FDR.threshold__classic.test']] = 1
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.FDR.threshold__trends.test']] = 1
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___scoring.adjustment.to.gene.list.size']] = TRUE
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___scoring.adj.list.range.optimal.size']] = 100
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___scoring.adj.list.range.start']] = 40
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___scoring.adj.list.range.end']] = 350
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___limit.Pvalues.to']] = 1e-30
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___sum.scores.with.power']] = 2.7
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___minimal.score']] = 40
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___desired.DB.entries.count__to.discard.score']] = 40
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___minimal.DB.entry.size']] = 5
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___sort.terms.by']] = 'name'
  
  
  for (db in c('GO','KEGG','Reactome')){
    for (item in names(Default.Parameters)[startsWith(names(Default.Parameters), 'classic.Enrich__with.clusterProfiler___')]){
      Default.Parameters[[sprintf('%s.%s',db,item)]] = Default.Parameters[[item]]
    }
  }
  
  for (db in c('GO','KEGG','Reactome')){
    for (item in names(Default.Parameters)[startsWith(names(Default.Parameters), 'trends.Enrich__with.clusterProfiler___')]){
      Default.Parameters[[sprintf('%s.%s',db,item)]] = Default.Parameters[[item]]
    }
  }
  
  for (db in c('GO','KEGG','Reactome')){
    for (item in names(Default.Parameters)[startsWith(names(Default.Parameters), 'clusterProfiler.Enriched.Expression.Profiles___')]){
      Default.Parameters[[sprintf('%s.%s',db,item)]] = Default.Parameters[[item]]
    }
  }
  
  for (db in c('GO','KEGG','Reactome')){
    for (item in names(Default.Parameters)[startsWith(names(Default.Parameters), 'clusterProfiler.Expression.Profiles___')]){
      Default.Parameters[[sprintf('%s.%s',db,item)]] = Default.Parameters[[item]]
    }
  }
  
  
  Default.Parameters[['GO.Enrich__with.clusterProfiler___Ontology.list']] = c('BP','MF','CC')
  
  Default.Parameters[['Create.clusterProfiler.Expression.Profiles']] = TRUE
  Default.Parameters[['Create.clusterProfiler.Custom.GO.Expression.Profiles']] = TRUE
  Default.Parameters[['Create.clusterProfiler.Custom.KEGG.Expression.Profiles']] = TRUE
  Default.Parameters[['Create.clusterProfiler.Custom.Reactome.Expression.Profiles']] = TRUE
  Default.Parameters[['Create.clusterProfiler.Enriched.GO.Expression.Profiles']] = TRUE
  Default.Parameters[['Create.clusterProfiler.Enriched.KEGG.Expression.Profiles']] = TRUE
  Default.Parameters[['Create.clusterProfiler.Enriched.Reactome.Expression.Profiles']] = TRUE
  
  
  
  # KEGG Pathways visualization
  Default.Parameters[['Perform.Pathway.Visualization']] = TRUE
  Default.Parameters[['Pathway.Visualization___CPM.aware']] = TRUE
  Default.Parameters[['Pathway.Visualization___max.gene.PValue']] = 0.05
  Default.Parameters[['Pathway.Visualization___LogFC.limits']] = 3
  
  #misc
  Default.Parameters[['Use.Info.Table']] = T
  
  return(Default.Parameters)
}


Eliminate.non.informative.lines = function(table){
  table = table[!(apply(table,MARGIN = 1, FUN = function(x) { any(x %in% NA) })),]
  commented.lines = lapply(X = table[,1],FUN = function (x) substring(x,first=1,last=1)) %in% '#'
  return(table[!commented.lines,])
}

endswith = function(x,ending){
  return(substring(x,nchar(x)-nchar(ending)+1,nchar(x)) == ending)
}

Delete.symbols.in.text <- function(x,symbols) {
  for (symbol in symbols){
    while (TRUE){
      if (grepl(pattern = symbol,x = x,fixed=T)){
        x = gsub(symbol,"",x,fixed=T)
      } else {
        break
      }
    }
  }
  return (x)
}

Delete.duplicated.symbols.in.text <- function(x,symbols) {
  for (symbol in symbols){
    while (TRUE){
      if (grepl(pattern = sprintf('%s%s',symbol,symbol),x = x,fixed=T)){
        x = gsub(pattern = sprintf('%s%s',symbol,symbol),symbol,x,fixed=T)
      } else {
        break
      }
    }
  }
  return (x)
}


Text.to.str.vector = function(x,sep = c('\\','/',' ',';')){
  for (symbol in sep) x = gsub(symbol,',',x,fixed=T)
  x = Delete.duplicated.symbols.in.text(x,',')
  return(strsplit(x,',')[[1]])
}

Text.to.int.vector = function(x){
  for (symbol in c('\\','/',' ',';')) x = gsub(symbol,',',x,fixed=T)
  x = Delete.duplicated.symbols.in.text(x,',')
  return(as.numeric(as.character(strsplit(x,',')[[1]])))
}

Read.Parameters <- function(Parameters.xlsx.file.name = NULL, prepare.for.ext.DE.data = FALSE){
  #!!!! as.numeric(as.character(testdata$x)), not as.numeric(testdata$x) 
  
  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")
  if (!('rstudioapi' %in% rownames(installed.packages()))) install.packages("rstudioapi")
  
  tryCatch(expr = {
    script.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
  }, error = function (err){
    script.dir = NULL
    warning('rstudioapi failed to locate script directory')
  })
  if (is.null(Parameters.xlsx.file.name)){
    Parameters.xlsx.file.name = sprintf("%s/RTrans.parameters.xlsx",dirname(rstudioapi::getActiveDocumentContext()$path))
    if (file.exists(Parameters.xlsx.file.name)){
      cat(sprintf('\nFile with parameters found: %s\n',Parameters.xlsx.file.name))
    } else stop(sprintf('Parameters.xlsx.file.name argument is not specified. No parameters *.xlsx file has been found in the default location '))
  }
  
  if(!is.null(script.dir))  Parameters.xlsx.file.name = gsub('\\{script_dir\\}',script.dir,Parameters.xlsx.file.name)
  
  
  if (!file.exists(Parameters.xlsx.file.name)) stop(sprintf('file %s does not exist',Parameters.xlsx.file.name))
  if (substring(Parameters.xlsx.file.name,nchar(Parameters.xlsx.file.name)-4,nchar(Parameters.xlsx.file.name)) != '.xlsx') stop('Please fill up input parameters in Excel workbook (see an example in ./RTrans.base/Parameters.example.xlsx)')
  
  
  # reading main worksheet 'Parameters'
  
  tryCatch(expr = { par.table = suppressWarnings(readxl::read_excel(Parameters.xlsx.file.name,sheet='Parameters',col_names = FALSE))
  }, error = function (err){
    stop(sprintf('Sheet "Parameters" is not found in workbook %s. Use an example in ./RTrans.base/Parameters.example.xlsx',Parameters.xlsx.file.name))
  })
  
  Pars = Get.Default.Parameters()
  
  tryCatch(expr = { schema = suppressWarnings(readxl::read_excel(Parameters.xlsx.file.name,sheet='Sample setup',col_names = TRUE))
  }, error = function (err){
    stop(sprintf('Sheet "Sample setup" is not found in workbook %s. Use an example in ./RTrans.base/Parameters.example.xlsx',Parameters.xlsx.file.name))
  })
  
  available.preds = colnames(schema)[-1]
  available.preds_to_abbreviations = sprintf('pred_%.4d',1:length(available.preds))
  names(available.preds_to_abbreviations) = available.preds
  abbreviations_to_available.preds = available.preds
  names(abbreviations_to_available.preds) = sprintf('pred_%.4d',1:length(available.preds))
  
  par.table = data.frame(par.table)
  par.table = par.table[,1:2]
  par.table = Eliminate.non.informative.lines(par.table)
  par.table[,1] = sapply(par.table[,1],function(x) gsub('"',"",gsub('"','',x)))
  par.table[,2] = sapply(par.table[,2],function(x) gsub('"',"",gsub('"','',x)))
  par.table[,1] = sapply(par.table[,1],function(x) gsub("'","",gsub('"','',x)))
  par.table[,2] = sapply(par.table[,2],function(x) gsub("'","",gsub('"','',x)))
  colnames(par.table) = c('variable','value')
  i = 1
  for (i in 1:dim(par.table)[1]){
    var = par.table[i,1]
    val = par.table[i,2]
    if (!(var %in% names(Pars))){
      message(sprintf('Warning. Unknown parameter "%s". Omitting',var))
      next
    }
    if (var == 'DE.package' | var == 'Method'){
      val = casefold(val)
      if (val %in% c('edger','deseq')){
        Pars[[var]] = val
      } else if ('spearman' %in% val){
        Pars[[var]] = val
      } else if ('pearson' %in% val){
        Pars[[var]] = 'pearson'
      } else if ('cor' %in% val){
        Pars[[var]] = 'combined_cor'
      } else {
        stop('Incorrect parameter "Method/DE.package". Only DESeq/edgeR/spearman/pearson/combined_cor are allowed')
      }
    } else if (var == 'Models.to.Test'){
      Pars[[var]] = val
    } else if (var == 'Species'){
      val = casefold(val)
      if (val %in% c('hsa','has','h.sapiens','homo sapiens','h. sapiens','sapiens','human')){
        val = 'hsa'
      } else if (val %in% c('mmu','m.musculus','mus musculus','m. musculus','musculus','mouse')){
        val = 'mmu'
      } else if (val %in% c('dme','d.melanogaster','drosophila melanogaster','d. melanogaster','melanogaster','drosophila','fly','fruitfly')){
        val = 'dme'
      }
      Pars[[var]]  = val
    } else if (var %in% c('counts.dir','results.dir','suppl.data.dir')){
      Pars[[var]] = gsub('\\','/',val,fixed=T)
    }  else if (val == '{auto}'){
      next
    } else if (var == 'Distance.method'){
      val = casefold(val)
      if (!(val %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski",'1-cor','canberra.weighted'))){
        stop('Incorrect parameter "Distance.method". These values are allowed: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski","1-cor", "canberra.weighted" (see an example in ./RTrans.base/Parameters.example.xlsx)')
      }
      Pars[[var]] = val
    } else if (var == 'RNA.Seq.norm.method'){
      if (!(val %in% c("TMM","upperquartile","none","RLE"))){
        stop('Incorrect parameter "RNA.Seq.norm.method". These values are allowed: "TMM","upperquartile","none","RLE" (see an example in ./RTrans.base/Parameters.example.xlsx)')
      }
      Pars[[var]] = val
    } else if(val == '"' || val == ''){
      #var = 'GO.classic.Enrich__with.clusterProfiler___mode.list'
      template = paste(unlist(strsplit(x = var, split = '.', fixed = TRUE))[-1], collapse = '.')
      if (!(template %in% names(Pars))){
        message(sprintf('Warning. Unknown template %s for parameter "%s". Omitting', template, var))
        next
      }
      Pars[[var]] = Pars[[template]]
    } else if (var == 'topGO.Enriched.Expression.Profiles___sort.terms.by'){
      val = casefold(val)
      if (!(val %in% c("name","id","score","none"))){
        stop('Incorrect parameter "topGO.Enriched.Expression.Profiles___sort.terms.by". These values are allowed: "name","id","score","none" (see an example in ./RTrans.base/Parameters.example.xlsx)')
      }
      # if (val == 'p') val = 'pvalue'
      Pars[[var]] = val
    } else if (var == 'min.samples.with.sufficient.CPM'){
      if (val == '{auto}') { Pars[[var]] = val
      } else {
        val = as.numeric(as.character(val))
        if (is.na(var)){
          stop(sprintf('Incorrect parameter "%s". Only numeric values or "{auto}" are allowed (see an example in ./RTrans.base/Parameters.example.xlsx)',var))
        }
        Pars[[var]] = val
      }
      
    } else if (endswith(var,'ntology.list')) {
      val = toupper(val)
      val = Text.to.str.vector(val)
      if (length(val) == 0 || any(!(val %in% c('BP', 'MF', 'CC')))){
        stop(sprintf('Incorrect parameter "%s". Only "BP", "MF" or "CC" ontologies are allowed (see an example in ./RTrans.base/Parameters.example.xlsx)',var))
      }  
      Pars[[var]] = val
    } else if (endswith(var,'_mode.list')) {
      val = casefold(val)
      val = Text.to.str.vector(val)
      if (length(val) == 0 || any(!(val %in% c('upreg', 'downreg', 'up', 'down')))){
        stop(sprintf('Incorrect parameter "%s". Only "upreg" or "downreg" values are allowed (see an example in ./RTrans.base/Parameters.example.xlsx)',var))
      }  
      Pars[[var]] = val
    } else if (endswith(var,'.list')) {
      val = Text.to.int.vector(val)
      if(any(is.na(val))){
        stop(sprintf('Incorrect parameter "%s". Only numeric values are allowed, e.g. "20,40,80,500" (see an example in ./RTrans.base/Parameters.example.xlsx)',var))
      }
      Pars[[var]] = val
    } else if (typeof(Pars[[var]]) == "logical"){
      val = casefold(val)
      if (!(val %in% c('true','t','f','false','yes','y','n','no','on','off'))){
        stop(sprintf('Incorrect parameter "%s". Only logical values are allowed, e.g. True/False (see an example in ./RTrans.base/Parameters.example.xlsx)',var))
      }
      val = (substring(val,1,1) == 't' || substring(val,1,1) == 'y' || val == 'on')
      Pars[[var]] = val
    } else if (typeof(Pars[[var]]) %in% c("numeric","double","float")){
      val = as.numeric(as.character(val))
      if (is.na(var)){
        stop(sprintf('Incorrect parameter "%s". Only numeric values are allowed (see an example in ./RTrans.base/Parameters.example.xlsx)',var))
      }
      Pars[[var]] = val
    } else if (typeof(Pars[[var]]) %in% c("character")){
      Pars[[var]] = val
    } else {
      print(var)
      print(val)
      print(typeof(Pars[[var]]))
      print(Pars[[var]])
      stop('Something strange.... needs debug')
    }
  }
  
  if(!is.null(script.dir)){
    Pars$suppl.data.dir = gsub('\\{script_dir\\}',script.dir,Pars$suppl.data.dir)
    Pars$results.dir = gsub('\\{script_dir\\}',script.dir,Pars$results.dir)
    Pars$counts.dir = gsub('\\{script_dir\\}',script.dir,Pars$counts.dir)
  }
  #Pars$Use.Info.Table = T
  
  # testing "Models.to.test" parameter
  if(!prepare.for.ext.DE.data){
    if (!('Models.to.Test' %in% names(Pars))){
      stop('Mandatory parameter "Models.to.Test" is missing (see an example in ./RTrans.base/Parameters.example.xlsx)')
    }
  }
  # looking for unfilled parameters
  for (n in names(Pars)[!(names(Pars) %in% par.table[,1])]){
    message(sprintf('Warning. Parameter %s is missing. Setting to default\n\n',n))
  }
  
  Pars$available.preds_to_abbreviations = available.preds_to_abbreviations
  Pars$abbreviations_to_available.preds = abbreviations_to_available.preds
  Pars$Models.to.Test %<>% Replace.Preds.with.Abbrebiations(., Pars) %>% gsub('"', '', ., fixed=TRUE) %>%
    Text.to.str.vector(. ,sep = c(',',';')) %>% sapply(., Delete.prefix.in.model) %>% sapply(., trim) %>%
    Revert.Preds.from.Abbrebiations(., Pars)
  
  return (Pars)
}


Get.BiomaRt.table <- function(complete.gene.list, Pars, DB.data, forced.maRt.table = NULL, use.official.gene.symbol = FALSE){
  hashmd5 = substr(digest(complete.gene.list),1,8)
  if (is.null(forced.maRt.table)) { General.maRt.table_file.name = sprintf("%s/maRt.table.%s.%s.tsv",Pars$suppl.data.dir,Pars$Species,hashmd5)
  } else General.maRt.table_file.name = forced.maRt.table
  
  if (file.exists(General.maRt.table_file.name)) {
    General.maRt.table = read.table(file=General.maRt.table_file.name, header = T, sep = '\t', stringsAsFactors = F)
  } else {
    
    if(Pars$Species %in% names(DB.data$Taxons.by.KEGG.codes)){
      taxon = DB.data$Taxons.by.KEGG.codes[Pars$Species]
    } else {
      taxon = Pars$Species
    }
    taxon = tolower(taxon)
    # taxon = 'Homo sapiens hhhh'
    tmp = unlist(strsplit(x = taxon, split = ' '))
    dataset.name = tolower(sprintf('%s%s_gene_ensembl', substr(tmp[1],1,1), tmp[2]))

    cat('\nQuerying biomaRt for gene info and saving this data to disk...\n')
    mart <- useMart("ensembl", dataset=dataset.name) #, host="www.ensembl.org"
    # if (Pars$Species == 'mmu') mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl") #, host="www.ensembl.org"
    # if (Pars$Species == 'dme') mart <- useMart("ensembl", dataset="dmelanogaster_gene_ensembl") #, host="www.ensembl.org"
    #curl.handle = getCurlHandle()
    
    needed.attributes = c("ensembl_gene_id","external_gene_name", "description","gene_biotype","entrezgene")
    if (Pars$Species == 'dme') needed.attributes = c(needed.attributes,'flybasename_gene','flybase_gene_id')
    if(use.official.gene.symbol){
      General.maRt.table = getBM(attributes=needed.attributes,filters="external_gene_name",values=complete.gene.list, mart=mart)
      General.maRt.table = General.maRt.table[!(duplicated(General.maRt.table[,"external_gene_name"])),]
      rownames(General.maRt.table) = General.maRt.table[,"external_gene_name"]
    } else {
      General.maRt.table = getBM(attributes=needed.attributes,filters="ensembl_gene_id",values=complete.gene.list, mart=mart)
      ##write.table(listAttributes(mart = mart),"ttt.tsv",sep = '\t')  ## listing all the attributes
      General.maRt.table = General.maRt.table[!(duplicated(General.maRt.table[,"ensembl_gene_id"])),]
      rownames(General.maRt.table) = General.maRt.table[,"ensembl_gene_id"]
    }
    General.maRt.table[,"description"] = gsub(pattern = "\\[Source:.*", replacement = "", x = General.maRt.table[,"description"], ignore.case = T,perl = FALSE)
    write.table.mod(x= General.maRt.table,file = General.maRt.table_file.name,sep = '\t')
  }
  # org.Dm.eg.db::
  return(General.maRt.table)
}


# Organisms.info.file.name = 'RTrans.base/Organisms.info.xlsx'
Read.Organisms.Info <- function(Organisms.info.file.name = NULL){

  DB.data = vector(mode = "list")
  
  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")
  if (!('rstudioapi' %in% rownames(installed.packages()))) install.packages("rstudioapi")
  
  if (is.null(Organisms.info.file.name)){
    Organisms.info.file.name = sprintf("%s/RTrans.base/Organisms.info.xlsx",dirname(rstudioapi::getActiveDocumentContext()$path))
    if (file.exists(Organisms.info.file.name)){
      cat(sprintf('\nFile with organisms info found: %s\n',Organisms.info.file.name))
    } else stop(sprintf('Organisms.info.file.name argument is not specified. No organism info *.xlsx file has been found in the default location '))
  }
  
  tryCatch(expr = {
    script.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
  }, error = function (err){
    script.dir = NULL
    warning('rstudioapi failed to locate script directory')
  })
  if(!is.null(script.dir))  Organisms.info.file.name = gsub('\\{script_dir\\}',script.dir,Organisms.info.file.name)
  
  
  if (!file.exists(Organisms.info.file.name)) stop(sprintf('file %s does not exist',Organisms.info.file.name))
  if (substring(Organisms.info.file.name,nchar(Organisms.info.file.name)-4,nchar(Organisms.info.file.name)) != '.xlsx') stop('Please fill up organisms info in Excel workbook')
  
  
  # reading main worksheet 'KEGG organism codes'
  
  tryCatch(expr = { KEGG.code.table = suppressWarnings(readxl::read_excel(Organisms.info.file.name,sheet='KEGG organism codes',col_names = TRUE))
  }, error = function (err){
    stop(sprintf('Sheet "KEGG organism codes" is not found in workbook %s',Organisms.info.file.name))
  })
  
  KEGG.code.table = as.data.frame(KEGG.code.table)
  KEGG.code.table = KEGG.code.table[!(KEGG.code.table[,1] %in% NA),]
  KEGG.code.table = KEGG.code.table[!startsWith(KEGG.code.table[,1],"#"),]
  
  KEGG.codes.by.alias.tmp.part1 = KEGG.code.table[,c('Aliases','KEGG organism code')]
  KEGG.codes.by.alias.tmp.part2 = KEGG.code.table[,c('Taxon','KEGG organism code')]
  colnames(KEGG.codes.by.alias.tmp.part1) = c('Aliases','KEGG organism code')
  colnames(KEGG.codes.by.alias.tmp.part2) = c('Aliases','KEGG organism code')
  KEGG.codes.by.alias.tmp = rbind(KEGG.codes.by.alias.tmp.part1, KEGG.codes.by.alias.tmp.part2)
  colnames(KEGG.codes.by.alias.tmp) = c('Aliases','KEGG organism code')
  
  KEGG.codes.by.alias.tmp = KEGG.codes.by.alias.tmp[!(KEGG.codes.by.alias.tmp[,'Aliases'] %in% NA),]
  KEGG.codes.by.alias.tmp = KEGG.codes.by.alias.tmp[!(duplicated(KEGG.codes.by.alias.tmp[,'Aliases'])),]
  KEGG.codes.by.alias = KEGG.codes.by.alias.tmp[,'KEGG organism code']
  names(KEGG.codes.by.alias) = tolower(KEGG.codes.by.alias.tmp[,'Aliases'])
  
  KEGG.code.table = KEGG.code.table[!(duplicated(KEGG.code.table[,'KEGG organism code'])),]
  rownames(KEGG.code.table) = KEGG.code.table[,'KEGG organism code']
  
  Taxons.by.KEGG.codes = KEGG.code.table[,'Taxon']
  names(Taxons.by.KEGG.codes) = KEGG.code.table[,'KEGG organism code']
  
  DB.data[['KEGG.code.table']] = KEGG.code.table
  DB.data[['KEGG.codes.by.alias']] = KEGG.codes.by.alias
  DB.data[['Taxons.by.KEGG.codes']] = Taxons.by.KEGG.codes
  
  tryCatch(expr = { Bioconductor.code.table = suppressWarnings(readxl::read_excel(Organisms.info.file.name,sheet='Bioconductor DBs',col_names = TRUE))
  }, error = function (err){
    stop(sprintf('Sheet "Bioconductor DBs" is not found in workbook %s',Organisms.info.file.name))
  })
  
  Bioconductor.code.table = as.data.frame(Bioconductor.code.table)
  Bioconductor.code.table = Bioconductor.code.table[!(duplicated(Bioconductor.code.table[,'KEGG organism code'])),]
  Bioconductor.Org.DBs.by.KEGG.codes = Bioconductor.code.table[,'Bioconductor DB']
  names(Bioconductor.Org.DBs.by.KEGG.codes) = Bioconductor.code.table[,'KEGG organism code']
  
  DB.data[['Bioconductor.Org.DBs.by.KEGG.codes']] = Bioconductor.Org.DBs.by.KEGG.codes
  
  tryCatch(expr = { Reactome.code.table = suppressWarnings(readxl::read_excel(Organisms.info.file.name,sheet='Reactome names',col_names = TRUE))
  }, error = function (err){
    stop(sprintf('Sheet "Reactome names" is not found in workbook %s',Organisms.info.file.name))
  })
  
  Reactome.code.table = as.data.frame(Reactome.code.table)
  Reactome.code.table = Reactome.code.table[!(duplicated(Reactome.code.table[,'KEGG organism code'])),]
  Reactome.names.by.KEGG.codes = Reactome.code.table[,'Reactome name']
  names(Reactome.names.by.KEGG.codes) = Reactome.code.table[,'KEGG organism code']
  
  DB.data[['Reactome.names.by.KEGG.codes']] = Reactome.names.by.KEGG.codes
  
  return(DB.data)
}

Read.DB.entries.from.Excel = function(file.name, sheet.name){
  #file.name = 'C:/Users/gskra/OneDrive/Moskalev flies-2016/RTrans.parameters.xlsx'
  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")
  if (substring(file.name,nchar(file.name)-4,nchar(file.name)) != '.xlsx') stop('Please fill up input parameters in Excel workbook (see an example in ./RTrans.base/Parameters.example.xlsx)')
  if (!file.exists(file.name)) stop(sprintf('file %s does not exist',file.name))
  
  # reading worksheet
  tryCatch(expr = { par.table = suppressWarnings(readxl::read_excel(file.name,sheet=sheet.name, col_names = FALSE))
  }, error = function (err){
    stop(sprintf('Sheet "%s" is not found in workbook %s. Use an example in ./RTrans.base/Parameters.example.xlsx', sheet.name, file.name))
  })
  
  par.table = data.frame(par.table)
  par.table = par.table[,1]
  par.table = par.table[!is.na(par.table)]
  commented.lines = sapply(X = par.table,FUN = function (x) substring(x,first=1,last=1)) %in% '#'
  par.table = par.table[!commented.lines]
  return(par.table)
}

Read.Custom.GO.terms = function(file.name){
  Read.DB.entries.from.Excel(file.name, 'GO terms (DE profiles)')
}

Read.Custom.Reactome.pathways = function(file.name){
  Read.DB.entries.from.Excel(file.name, 'Reactome (DE-profiles)')
}

Read.Custom.KEGG.pathways <- function(file.name, Species, sheet.name = 'KEGG (visualization)'){
  #file.name = 'C:/Users/gskra/OneDrive/Moskalev flies-2016/RTrans.parameters.xlsx'
  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")
  if (substring(file.name,nchar(file.name)-4,nchar(file.name)) != '.xlsx') stop('Please fill up input parameters in Excel workbook (see an example in ./RTrans.base/Parameters.example.xlsx)')
  if (!file.exists(file.name)) stop(sprintf('file %s does not exist',file.name))
  
  # reading worksheet
  tryCatch(expr = { par.table = suppressWarnings(readxl::read_excel(file.name,sheet=sheet.name,col_names = F))
  }, error = function (err){
    stop(sprintf('Sheet "%s" is not found in workbook %s. Use an example in ./RTrans.base/Parameters.example.xlsx', sheet.name, file.name))
  })
  
  par.table = data.frame(par.table)
  par.table = par.table[,1]
  par.table = par.table[!is.na(par.table)]
  commented.lines = sapply(X = par.table,FUN = function (x) substring(x,first=1,last=1)) %in% '#'
  par.table = par.table[!commented.lines]
  final.list = as.numeric(as.character(par.table))
  if (any(is.na(final.list)))  message(sprintf('Please provide only numerical KEGG id (e.g. 5133 instead of hsa05133) in "%s" sheet. Such values:  %s', sheet.name, toString(par.table[is.na(final.list)])))
  final.list = final.list[!is.na(final.list)]
  
  return(sprintf("%s%.5d",Species,final.list))
}

Read.KEGG.pathways.to.omit <- function(file.name,Species){
  #file.name = 'C:/Users/gskra/OneDrive/Moskalev flies-2016/RTrans.parameters.xlsx'
  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")
  if (substring(file.name,nchar(file.name)-4,nchar(file.name)) != '.xlsx') stop('Please fill up input parameters in Excel workbook (see an example in ./RTrans.base/Parameters.example.xlsx)')
  if (!file.exists(file.name)) stop(sprintf('file %s does not exist',file.name))
  
  # reading worksheet 'Omit KEGG pathways'
  tryCatch(expr = { par.table = readxl::read_excel(file.name,sheet='Omit KEGG pathways',col_names = F)
  }, error = function (err){
    stop(sprintf('Sheet "Omit KEGG pathways" is not found in workbook %s. Use an example in ./RTrans.base/Parameters.example.xlsx',file.name))
  })
  
  par.table = data.frame(par.table)
  par.table = par.table[,1]
  par.table = par.table[!is.na(par.table)]
  commented.lines = sapply(X = par.table,FUN = function (x) substring(x,first=1,last=1)) %in% '#'
  par.table = par.table[!commented.lines]
  final.list = as.numeric(as.character(par.table))
  if (any(is.na(final.list))) message(sprintf('Please provide only numerical KEGG id (e.g. 5133 instead of hsa05133) in "Omit KEGG pathways" sheet. Such values:  %s',toString(par.table[is.na(final.list)])))
  final.list = final.list[!is.na(final.list)]
  
  return(sprintf("%s%.5d",Species,final.list))
}

Convert.GO.term.names_to_GO.term.IDs = function(terms, GO.term.ID.by.name, remove.redundant = TRUE, prefix = ''){
  not_found_terms_count = 0
  terms__converted = c()
  for (x in terms){
    x = toupper(x)
    if (!(substring(x,1,3) %in% 'GO:')){
      if (x %in% names(GO.term.ID.by.name)){
        x = GO.term.ID.by.name[x]
      } else {
        message(sprintf('\nWarning. Custom GO term "%s" is not found', x))
        not_found_terms_count = not_found_terms_count + 1
      }
    }
    if (!(remove.redundant) | (!(x %in% terms__converted))){
      terms__converted = c(terms__converted, x)
    }
  }
  if (not_found_terms_count > 0)   message(sprintf('\n%sTotal custom GO %d terms were not found', prefix, not_found_terms_count))
  return(terms__converted)
}

# Parameters.file.name = 'C:/Users/gskra/Documents/NA_5FU-2/RTrans.parameters for NA_5FU-3.xlsx'
# sheet.name = 'GO terms (inDetails)'

Read.DB.entries.from.Excel__MultiSet = function(file.name, sheet.name){
  #file.name = 'C:/Users/gskra/OneDrive/Moskalev flies-2016/RTrans.parameters.xlsx'
  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")
  if (substring(file.name,nchar(file.name)-4,nchar(file.name)) != '.xlsx') stop('Please fill up input parameters in Excel workbook (see an example in ./RTrans.base/Parameters.example.xlsx)')
  if (!file.exists(file.name)) stop(sprintf('file %s does not exist',file.name))
  
  # reading worksheet
  tryCatch(expr = { DB.entry.table = suppressWarnings(readxl::read_excel(file.name,sheet=sheet.name, col_names = TRUE))
  }, error = function (err){
    stop(sprintf('Sheet "%s" is not found in workbook %s. Use an example in ./RTrans.base/Parameters.example.xlsx', sheet.name, file.name))
  })
  
  DB.entry.table = as.data.frame(DB.entry.table)
  DB.entry.vector = 
  lapply(as.list(DB.entry.table),function(x){ 
    x = x[!is.na(x)]
    x = x[!startsWith(x, '#')]
    return(x)
  })
  
  DB.entry.vector.clean = vector(mode = 'list')
  #names(DB.entry.vector.clean)
  for(entry in names(DB.entry.vector)){
    if(startsWith(entry, '#')) next
    DB.entry.vector.clean[[entry]] = DB.entry.vector[[entry]]
  }
  
  return(DB.entry.vector.clean)
  
  # par.table = par.table[,1]
  # par.table = par.table[!is.na(par.table)]
  # return(par.table)
}



Read.Completed.steps.status <- function(step.name,Completed.steps.file){
  if (file.exists(Completed.steps.file)){
    Completened.steps = as.vector(t(read.table(file = Completed.steps.file)))
  } else{
    Completened.steps = c()
  }
  if(step.name %in% Completened.steps){
    #cat(sprintf('\nStep %s had been alrealy accomplished. Bypassing.\n',step.name))
    return(TRUE)
  }
  return(FALSE)
}

Add.Completed.step <- function(step.name,Completed.steps.file){
  if (file.exists(Completed.steps.file)){
    Completened.steps = as.vector(t(read.table(file = Completed.steps.file)))
  } else {
    Completened.steps = c()
  }
  Completened.steps = c(Completened.steps, step.name)
  write.table.mod(x=Completened.steps, file = Completed.steps.file,sep='\t')
}


Delete.spaces.in.model <- function(model) {
  repeat{
    if(startsWith(model," ")){
      model = substr(model,2,nchar(model))
    } else {
      break
    }
  }
  
  repeat{
    if(endsWith(model," ")){
      model = substr(model, 1, nchar(model)-1)
    } else {
      break
    }
  }
  
  for (symbol in c('+',':','~','*')){
    while (TRUE){
      if (grepl(pattern = sprintf(" %s",symbol),x = model,fixed=T) | grepl(pattern = sprintf("%s ",symbol),x = model,fixed=T)){
        model = gsub(sprintf(" %s",symbol),symbol,gsub(sprintf("%s ",symbol),symbol,model,fixed=T),fixed=T)
      } else {
        break
      }
    }
  }
  return (model)
}

Delete.prefix.in.model <- function(model) {
  for (symbol in c('~')){
    while (TRUE){
      if (grepl(pattern = sprintf(" %s",symbol),x = model,fixed=T) | grepl(pattern = sprintf("%s ",symbol),x = model,fixed=T)){
        model = gsub(sprintf(" %s",symbol),symbol,gsub(sprintf("%s ",symbol),symbol,model,fixed=T),fixed=T)
      } else {
        break
      }
    }
  }
  model = gsub('~','',model,fixed=T)
  return (model)
}

Extract.predictors <- function (model) {
  model = Delete.spaces.in.model(model)
  model = gsub('\\~','',model)
  for (symbol in c('\\+','\\:','\\~','\\*')){
    model = gsub(symbol,"\\+",model,fixed=T)
  }
  tryCatch(expr = {   model = data.frame(strsplit(model,split = '\\+'),stringsAsFactors = F)[,1]
  }, error = function (err){
    model = strsplit(as.character(model),split = '[',fixed = T)[[1]]
  })
  return(model[!duplicated(model)])
  #return(levels(factor(model)))
}

replace.spec.symbols <- function (predictor.name) {
  predictor.name = gsub("+","&plus&",predictor.name,fixed=T)
  predictor.name = gsub(":","&colon&",predictor.name,fixed=T)
  predictor.name = gsub("~","&tilde&",predictor.name,fixed=T)
  predictor.name = gsub("*","&asterisk&",predictor.name,fixed=T)
  return(predictor.name)
}

revert.spec.symbols <- function (predictor.name) {
  predictor.name = gsub("&plus&","+",predictor.name,fixed=T)
  predictor.name = gsub("&colon&",":",predictor.name,fixed=T)
  predictor.name = gsub("&tilde&","~",predictor.name,fixed=T)
  predictor.name = gsub("&asterisk&","*",predictor.name,fixed=T)
  return(predictor.name)
}



Prepare = function(Parameters.xlsx.file.name = NULL, forced.maRt.table = NULL, force.re.install = FALSE, disable.excel = FALSE, ext.DE.data.files = NULL,
                   ext.DE.data.species = NULL, ext.DE.data.PValue.if.absent = 0.01, ext.DE.data.logCPM.if.absent = 5.0, ext.DE.data.Score.if.absent = 'auto'){
  #Parameters.xlsx.file.name  ='C:/Users/gskra/OneDrive/Evgeniev/RTrans.parameters.xlsx'
  if (!('rstudioapi' %in% rownames(installed.packages()))) install.packages("rstudioapi")
  if (!('magrittr' %in% rownames(installed.packages()))) install.packages("magrittr")
  library(magrittr)
  tryCatch(expr = {
    script.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
  }, error = function (err){
    script.dir = NULL
    warning('rstudioapi failed to locate script directory')
  })
  
  #ext.DE.data.files=c('ext.DE.data.example.txt')
  if (is.null(Parameters.xlsx.file.name)){
    Parameters.xlsx.file.name = sprintf("%s/RTrans.parameters.xlsx",dirname(rstudioapi::getActiveDocumentContext()$path))
    if (file.exists(Parameters.xlsx.file.name)){
      cat(sprintf('\nFile with parameters found: %s\n',Parameters.xlsx.file.name))
    } else stop(sprintf('Parameters.xlsx.file.name argument is not specified. No parameters *.xlsx file has been found in the default location '))
  }
  if (!is.null(script.dir))  Parameters.xlsx.file.name = gsub('\\{script_dir\\}',script.dir,Parameters.xlsx.file.name)
  if (!is.null(script.dir) & !is.null(forced.maRt.table)) forced.maRt.table = gsub('\\{script_dir\\}',script.dir,forced.maRt.table)
  
  #Startup.Data = new('RTransData')
  Pars = Read.Parameters(Parameters.xlsx.file.name, !is.null(ext.DE.data.files))
  if(!is.null(ext.DE.data.species)){
    Pars$Species = ext.DE.data.species
  }
  
  DB.data = Read.Organisms.Info(Organisms.info.file.name = sprintf('%s/Organisms.info.xlsx',Pars$suppl.data.dir))
  
  Startup.Data = new('RTransStartupData')
  Startup.Data$DB.data = DB.data
  
  Pars$Species = tolower(Pars$Species)
  if (!(Pars$Species %in% DB.data$KEGG.codes.by.alias)){
    if(Pars$Species %in% names(DB.data$KEGG.codes.by.alias)){
      Pars$Species = DB.data$KEGG.codes.by.alias[Pars$Species]
    } else {
      stop(sprintf('Organism "%s" is not found among available taxons and KEGG abbreviations. Please verify RTrans parameters *.xlsx-file',Pars$Species))
    }
  }
  
  Startup.Data$org.DB.name = NULL
  if(Pars$Species %in% names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
    Startup.Data$org.DB.name = Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes[[Pars$Species]]
  }
  
  # DB.data$Reactome.names.by.KEGG.codes
  # DB.data$Bioconductor.Org.DBs.by.KEGG.codes
  # DB.data$KEGG.codes.by.alias
  
  # installing/loading packages
  cat('\nLoading or installing packages...\n')

  std.pkg.list = c('gplots','chemometrics','ggdendro','digest',"lme4","readxl","dendextend",'ggplot2', 'dplyr','stringr','forcats','pheatmap', 'magrittr', 'utils')
  bio.pkg.list = c()
  if(!is.null(Startup.Data$org.DB.name))  bio.pkg.list = c(Startup.Data$org.DB.name)
  bio.pkg.list = c(bio.pkg.list, "Rgraphviz","edgeR","DESeq2","FactoMineR","biomaRt","topGO","clusterProfiler", "DOSE","ReactomePA","reactome.db",'AnnotationDbi','GOSemSim')

  if(force.re.install){
    install.packages(std.pkg.list, dep = TRUE)
    source("https://bioconductor.org/biocLite.R")
    biocLite(bio.pkg.list, dep = TRUE)
    update.packages(ask = FALSE)
  } else {
    std.pkg.list.to.install = std.pkg.list[!(std.pkg.list %in% installed.packages())]
    if(length(std.pkg.list.to.install > 0)) install.packages(std.pkg.list.to.install)
    
    bio.pkg.list.to.install = bio.pkg.list[!(bio.pkg.list %in% installed.packages())]
    if(length(bio.pkg.list.to.install > 0)){
      source("https://bioconductor.org/biocLite.R")
      biocLite(bio.pkg.list.to.install)
    }
  }
  
  for(lib in c(std.pkg.list, bio.pkg.list))  eval(parse(text = sprintf('suppressPackageStartupMessages(library(%s))',lib)))
  
  if (Pars$Pathway.Visualization___CPM.aware){
    pathviewmod.Installed.successfully = TRUE
    if (!('pathviewmod' %in% rownames(installed.packages()))){
      pathviewmod.Installed.successfully = FALSE
      tryCatch(expr = {
        install.packages(sprintf('%s/pathview_mod',Pars$suppl.data.dir), repos = NULL, type="source")
        pathviewmod.Installed.successfully = TRUE
      }, error = function (err){
        pathviewmod.Installed.successfully = FALSE
      })
    }
    
    tryCatch(expr = { suppressPackageStartupMessages(library('pathviewmod'))
    }, error = function (err){
      pathviewmod.Installed.successfully = FALSE
    })
    
    if (!pathviewmod.Installed.successfully | !('pathviewmod' %in% rownames(installed.packages()))){
      message('Warning. RTrans is provided with a modified, CPM-aware version of "pathview". It gives much better results for the analysis of massive gene expression profiles changes, e.g. age-associated. RTrans was unable to install the modified package, "pathviewmod" to your system. CPM-aware pathway visualization will be disabled.')
      Pars$Pathway.Visualization___CPM.aware = F
    }
  }
  
  if (!Pars$Pathway.Visualization___CPM.aware){
    load.or.install.packages.bio(c("pathview"))
    suppressPackageStartupMessages(library('pathview'))
  }
  
  
  cat('\nChecking input *.counts files and pre-processing expression data...\n')
  
  Pars$RefSeq.Info.File = sprintf('%s/%s.Ensembl.Genes.Info.txt',Pars$suppl.data.dir,Pars$Species)
  Pars$Pathway.names.table.File = sprintf('%s/KEGG.pathway.names.tsv',Pars$suppl.data.dir)
  Pars$GO.full.descriptions.File = sprintf('%s/Gene Ontology terms names and descriptions.tsv',Pars$suppl.data.dir)
  
  
  if (Pars$Use.Info.Table & file.exists(Pars$RefSeq.Info.File)) {
    Startup.Data$Info.table = read.table(Pars$RefSeq.Info.File,stringsAsFactors = FALSE, blank.lines.skip = FALSE,sep = '\t')
  } else {
    Pars$Use.Info.Table = FALSE
    Startup.Data$Info.table = NA
  }
  
  # Loading lists of custom GO terms and KEGG pathways
  topGO.Expression.Profiles___Custom.GO.terms = Read.Custom.GO.terms(Parameters.xlsx.file.name)
  ## GO.clusterProfiler.Expression.Profiles___Custom.DB.entries is exactly topGO.Expression.Profiles___Custom.GO.terms (assigned further)
  
  Custom.Pathways.to.Visualize = Read.Custom.KEGG.pathways(Parameters.xlsx.file.name, Pars$Species)
  Omit.pathways = Read.KEGG.pathways.to.omit(Parameters.xlsx.file.name,Pars$Species)
  
  KEGG.clusterProfiler.Expression.Profiles___Custom.DB.entries = Read.Custom.KEGG.pathways(Parameters.xlsx.file.name, Pars$Species, sheet.name = 'KEGG (DE profiles)')
  
  R.DBe = Read.DB.entries.from.Excel(Parameters.xlsx.file.name, 'Reactome (DE profiles)')
  #x =  DB.entries[1]
  R.Sp = sapply(R.DBe, function(x){
    found.obj = gregexpr("-",x)[[1]]
    if(length(found.obj) == 2){
      return(tolower(substr(x, found.obj[1] + 1, found.obj[2] - 1)))
    } else {
      return ('-')
    }
    })
  
  
  R.Sp.incorrect = sum(Pars$Species != R.Sp)
  R.Sp.correct = sum(Pars$Species == R.Sp)
  if(R.Sp.incorrect > 0){
    cat(sprintf('%d Reactome pathways with unexpected species ID have been detected: %s.  Species for these IDs will be replaced with "%s"', R.Sp.incorrect, paste(levels(as.factor(R.Sp[Pars$Species != R.Sp])), collapse = ', '), Pars$Species))
    #x = R.DBe[1]
    R.DBe = sapply(R.DBe, function(x){
      found.obj = gregexpr("-",x)[[1]]
      if(length(found.obj) == 2){
        return(
          sprintf('%s%s%s', substr(x, 1, found.obj[1]),
                toupper(Pars$Species),
                substr(x, found.obj[2], nchar(x)))
        )
      } else {
        return (x)
      }
    })
    names(R.DBe) = NULL
  }
  Reactome.clusterProfiler.Expression.Profiles___Custom.DB.entries = R.DBe
  
  ## Loading pathway names
  # if (Pars$Perform.Pathway.GSEA | Pars$Perform.Pathway.Visualization){
  Pathway.names.table = read.table(Pars$Pathway.names.table.File,sep='\t',stringsAsFactors = F,quote = "")
  colnames(Pathway.names.table) = c("KEGG ID","Pathway name")
  Pathway.names.table[,"KEGG ID"] = sprintf("%s%.5d",Pars$Species,Pathway.names.table[,"KEGG ID"])
  rownames(Pathway.names.table) = Pathway.names.table[,"KEGG ID"]
  Pathway.names.table[,"Pathway name"] = gsub(pattern = "'",replacement = "",x = Pathway.names.table[,"Pathway name"])
  Pathway.names.table[,"Pathway name"] = gsub(pattern = "/",replacement = "-",x = Pathway.names.table[,"Pathway name"])
  Startup.Data$Pathway.names.table = Pathway.names.table
  # }
  
  ## Downloading KEGG info using clusterProfiler functions - code is copied from clusterProfiler sources
  source(sprintf('%s/kegg.clusterProfiler.R',Pars$suppl.data.dir))
  if(Pars$Species == 'dme'){
    Startup.Data$KEGG_DATA <- prepare_KEGG(species = Pars$Species, keyType = 'ncbi-geneid')
  } else {
    Startup.Data$KEGG_DATA <- prepare_KEGG(species = Pars$Species, keyType = 'kegg')
  }
  
  
  #
  
  ## Retrieving names of all GO terms
  #Startup.Data$all.GO.terms <- Term(GOTERM)
  #
  
  # Loading GO names and descriptions
  # if (Pars$Perform.GO.Enrich__with.topGO | length(topGO.Expression.Profiles___Custom.GO.terms) > 0){
    
  GO.full.descriptions = read.table(Pars$GO.full.descriptions.File,sep='\t',check.names = FALSE,quote = "",stringsAsFactors = F,header = F,comment.char = '#')
  colnames(GO.full.descriptions) = c('GO term','name','full description')
  rownames(GO.full.descriptions) = GO.full.descriptions[,'GO term']
  Startup.Data$GO.full.descriptions = GO.full.descriptions
  
  tmp = GO.full.descriptions[!(duplicated(GO.full.descriptions[,'name'])),]
  GO.term.ID.by.name = tmp[,'GO term']
  names(GO.term.ID.by.name) = toupper(tmp[,'name'])
  Startup.Data$GO.term.ID.by.name = GO.term.ID.by.name
  
  topGO.Expression.Profiles___Custom.GO.terms = 
    Convert.GO.term.names_to_GO.term.IDs(terms = topGO.Expression.Profiles___Custom.GO.terms,
                                       GO.term.ID.by.name = GO.term.ID.by.name,
                                       remove.redundant = Pars$topGO.Expression.Profiles___remove.redundant.GO.terms,
                                       prefix = '[GO terms (DE profiles)]: ')
  
  
  GO.terms.inDetails = Read.DB.entries.from.Excel__MultiSet(Parameters.xlsx.file.name, 'GO terms (inDetails)')
  GO.terms.inDetails = lapply(GO.terms.inDetails, function(x){
    Convert.GO.term.names_to_GO.term.IDs(terms = x, GO.term.ID.by.name = GO.term.ID.by.name,
                                         remove.redundant = Pars$topGO.Expression.Profiles___remove.redundant.GO.terms,
                                         prefix = '[GO terms (inDetails)]: ')
  })
  Startup.Data$GO.terms.inDetails = GO.terms.inDetails
  rm(GO.term.ID.by.name)
  
  Analyses.groups = Read.DB.entries.from.Excel__MultiSet(Parameters.xlsx.file.name, 'Analyses groups')
  if(length(Analyses.groups) == 0){
    Analyses.groups[['main analysis']] = Pars$Models.to.Test
  }

  Completed.steps.file = sprintf("%s/Completed.steps.list",Pars$results.dir)
  
  setwd(Pars$results.dir)
  
  ##### Pre-creating complete gene list and retrieving biomaRt infos for these genes
  Startup.Data$use.ext.DE.data = FALSE
  if(is.null(ext.DE.data.files)){
    tryCatch(expr = { schema = suppressWarnings(readxl::read_excel(Parameters.xlsx.file.name,sheet='Sample setup',col_names = TRUE))
    }, error = function (err){
      stop(sprintf('Sheet "Sample setup" is not found in workbook %s. Use an example in ./RTrans.base/Parameters.example.xlsx',Parameters.xlsx.file.name))
    })
    schema %<>% as.data.frame
    commented.lines = sapply(schema[,1], FUN = function (x) substring(x,first=1,last=1)) %in% '#'
    schema = schema[!commented.lines,]
    schema = schema[,!apply(schema, 2, function(x) all(is.na(x)))]
    if (!("Sample names" %in% colnames(schema))) stop(sprintf('Column "Sample names" is not found in the workbook %s. Use an example in ./RTrans.base/Parameters.example.xlsx',Parameters.xlsx.file.name))
    rownames(schema) = schema$'Sample names'
    
    #colnames(schema) = replace.spec.symbols(colnames(schema))
  
    #if (!("Sample names" %in% colnames(schema))) stop(sprintf('Column "Sample names" is not found in the workbook %s. Use an example in ./RTrans.base/Parameters.example.xlsx',Parameters.xlsx.file.name))
    
    if ("File names" %in% colnames(schema)){
      file.names = sprintf("%s/%s",Pars$counts.dir,schema$'File names')
    } else {
      file.names = sprintf("%s/%s.counts",Pars$counts.dir,schema$'Sample names')
    }
    if (any(!file.exists(file.names))){
      #warning(sprintf('Warning; The following files do not exist: %s',toString(file.names[!file.exists(file.names)])))
      cat(sprintf('\nTotal %d of %d *.counts files do not exist: %s',
                  sum.mod(!file.exists(file.names)), length(file.names), toString(
                    sapply(strsplit(file.names[!file.exists(file.names)], split='/'), function(x) { tail(x,1)})
                    )))
      if (sum.mod(file.exists(file.names)) == 0) stop('No *.counts files found. Fix sample setup')
      schema=schema[file.exists(file.names),]
      rownames(schema) = schema$'Sample names'
      file.names = file.names[file.exists(file.names)]
    }
    samples_counts=suppressMessages(readDGE(file.names))
    
    if (length(Pars$Models.to.Test) == 1){
      if (Pars$Models.to.Test=='{auto}')    Pars$Models.to.Test = colnames(schema)[!(colnames(schema) %in% c('Sample names','File names'))]
    }
    Pars$Models.to.Test = unique(Pars$Models.to.Test[Pars$Models.to.Test != ""])
    
    # ## creating predictors substitution dict, in order to avoiding recognizing predictors
    # ## like "HCT116, p53- vs p53+, 4h" as a GLM model "~ <HCT116, p53- vs p53> + <, 4h>"
    # available.preds = colnames(schema)[-1]
    # available.preds_to_abbreviations = sprintf('pred_%d',1:length(available.preds))
    # names(available.preds_to_abbreviations) = available.preds
    # abbreviations_to_available.preds = available.preds
    # names(abbreviations_to_available.preds) = sprintf('pred_%d',1:length(available.preds))
    # Pars$available.preds_to_abbreviations = available.preds_to_abbreviations
    # Pars$abbreviations_to_available.preds = abbreviations_to_available.preds


    #sapply(Pars$Models.to.Test, function(x){ r = Replace.Preds.with.Abbrebiations(x,Pars)  })
    #sapply(d, function(x){ r = Revert.Preds.from.Abbrebiations(x,Pars)  })
    
    # checking for the presence of all needed predictors
    all.preds = c()
    conv.models = sapply(Pars$Models.to.Test, function(x){ r = Replace.Preds.with.Abbrebiations(x, Pars) })
    for (tmp.Pred.set in conv.models)  all.preds = append(all.preds, Extract.predictors(tmp.Pred.set))
    all.preds = levels(factor(all.preds))
    all.preds = sapply(all.preds, function(x){ r = Revert.Preds.from.Abbrebiations(x,Pars) })
    if(sum.mod(!(all.preds %in% colnames(schema))) > 0){
      stop(sprintf('The following predictors from GLM models ("Parameters" sheet) are not found in "Sample Setup" sheet: %s. Please check "Sample Setup" row names or "Models.to.Test" parameter in "Parameters" sheet',toString(all.preds[!(all.preds %in% colnames(schema))])))
    }
    
    Analyses.groups %>% unlist %>% .[!duplicated(.)] -> models.from.Analyses.groups

    all.preds.from.Analyses.groups = c()
    conv.models = sapply(models.from.Analyses.groups, function(x){ r = Replace.Preds.with.Abbrebiations(x, Pars) })
    for (tmp.Pred.set in conv.models)  all.preds.from.Analyses.groups = append(all.preds.from.Analyses.groups, Extract.predictors(tmp.Pred.set))
    all.preds.from.Analyses.groups = levels(factor(all.preds.from.Analyses.groups))
    all.preds.from.Analyses.groups = sapply(all.preds.from.Analyses.groups, function(x){ r = Revert.Preds.from.Abbrebiations(x,Pars) })
    if(sum.mod(!(all.preds.from.Analyses.groups %in% colnames(schema))) > 0){
      stop(sprintf('The following predictors from GLM models ("Analyses groups" sheet) are not found in "Sample Setup" sheet: %s. Please check "Sample Setup" row names or GLm models listed in "Analyses groups" sheet',toString(all.preds.from.Analyses.groups[!(all.preds.from.Analyses.groups %in% colnames(schema))])))
    }

    noint = rownames(samples_counts$counts) %in% c("__no_feature","__ambiguous",'__too_low_aQual','__not_aligned','__alignment_not_unique')
    complete.gene.list = sort(rownames(samples_counts$counts)[!noint],decreasing = F)
    General.maRt.table = Get.BiomaRt.table(complete.gene.list, Pars, DB.data, forced.maRt.table = NULL)
  } else {
    Startup.Data$use.ext.DE.data = TRUE
    Pars$Models.to.Test = ext.DE.data.files
    Startup.Data$ext.DE.data.files = ext.DE.data.files
    Startup.Data$ext.DE.data = vector(mode = "list")
    
    complete.Ensembl.gene.list = c()
    for(f in ext.DE.data.files){
      DE.data = read.table(f,stringsAsFactors = FALSE,sep = '\t',header = T,check.names = FALSE)
      colnames(DE.data) = tolower(colnames(DE.data))
      if(!('gene' %in%  colnames(DE.data))){
        stop(sprintf('Column "gene" should be present in the external expression data file "%s"',f))
      }
      if(!('logfc' %in%  colnames(DE.data))){
        stop(sprintf('Column "logFC" should be present in the external expression data file "%s"',f))
      }
      if(!('pvalue' %in%  colnames(DE.data))){
        cat(sprintf('Column "PValue" is not present in the external expression data file "%s". Assuming p-value = %g\n',f,ext.DE.data.PValue.if.absent))
        # DE.data = DE.data[,-3]
        new.colnames = c(colnames(DE.data),'pvalue')
        DE.data = cbind(DE.data,data.frame(rep(ext.DE.data.PValue.if.absent,dim(DE.data)[1])))
        colnames(DE.data) = new.colnames
      }
      if(!('fdr' %in%  colnames(DE.data))){
        cat(sprintf('Column "fdr" is not present in the external expression data file "%s". Deriving FDR from p-values\n',f))
        new.colnames = c(colnames(DE.data),'fdr')
        DE.data = cbind(DE.data,data.frame(p.adjust(DE.data[,'pvalue'],method = 'BH')))
        colnames(DE.data) = new.colnames
      }
      if(!('logcpm' %in%  colnames(DE.data))){
        cat(sprintf('Column "logCPM" is not present in the external expression data file "%s". Assuming logCPM = %g\n',f,ext.DE.data.logCPM.if.absent))
        new.colnames = c(colnames(DE.data),'logcpm')
        DE.data = cbind(DE.data,data.frame(rep(ext.DE.data.logCPM.if.absent,dim(DE.data)[1])))
        colnames(DE.data) = new.colnames
      }
      if(!('score' %in%  colnames(DE.data))){
        if (ext.DE.data.Score.if.absent == 'auto'){
          cat(sprintf('Column "Score" is not present in the external expression data file "%s". Deriving Score from P-values and logFC\n',f))
          Scores = ((-1)*log2(DE.data[,"pvalue"] + 1e-300))^1.5* 10 *DE.data[,"logfc"]
          Scores = 10*abs(Scores)**0.5
          new.colnames = c(colnames(DE.data),'score')
          DE.data = cbind(DE.data,data.frame(Scores))
          colnames(DE.data) = new.colnames
        } else {
          cat(sprintf('Column "Score" is not present in the external expression data file "%s". Assuming Score = %g\n',f,ext.DE.data.Score.if.absent))
          new.colnames = c(colnames(DE.data),'score')
          DE.data = cbind(DE.data,data.frame(rep(ext.DE.data.Score.if.absent,dim(DE.data)[1])))
          colnames(DE.data) = new.colnames
        }
      }
      rownames(DE.data) = DE.data[,'gene']
      complete.gene.list = DE.data[,'gene']
      
      ## detecting if there is ENSEMBL gene accessions
      ENS.genes.count = length(grep(pattern = 'ENS',x = complete.gene.list,ignore.case = TRUE))
      if (ENS.genes.count/length(complete.gene.list) < 0.33){
        cat(sprintf('Transforming official gene symbols to Ensembl code (file "%s")\n',f))
      }
      
      maRt.table = Get.BiomaRt.table(complete.gene.list, Pars, DB.data, use.official.gene.symbol = TRUE)
      cat(sprintf('Total %d of %d official gene symbols cannot be mapped to Ensembl (usually most of them are LOC.... c3orf...) \n',sum.mod(!(DE.data[,'gene'] %in% rownames(maRt.table))),dim(DE.data)[1]))
      
      present = DE.data[,'gene'] %in% rownames(maRt.table)
      DE.data = DE.data[present,]

      new.colnames = c(colnames(DE.data),'ensembl_gene_id')
      DE.data = cbind(DE.data,data.frame(maRt.table[DE.data[,'gene'],'ensembl_gene_id']))
      colnames(DE.data) = new.colnames
      
      cat(sprintf('%d of %d mapped Ensembl gene IDs are duplicated\n',sum.mod(duplicated(DE.data[,'ensembl_gene_id'])),dim(DE.data)[1]))
      DE.data = DE.data[!duplicated(DE.data[,'ensembl_gene_id']),]
      Startup.Data$ext.DE.data[[f]] = DE.data
      
      complete.Ensembl.gene.list = c(complete.Ensembl.gene.list,as.character(DE.data[,'ensembl_gene_id']))
    }
    
#    Pars$Models.to.Test
    complete.Ensembl.gene.list = complete.Ensembl.gene.list[!duplicated(complete.Ensembl.gene.list)]
    General.maRt.table = Get.BiomaRt.table(complete.Ensembl.gene.list, Pars, DB.data, forced.maRt.table = NULL)
  }
  
  tmp.maRt.table = General.maRt.table[!(General.maRt.table[,'ensembl_gene_id'] %in% NA),]
  tmp.maRt.table = tmp.maRt.table[!duplicated(tmp.maRt.table[,'ensembl_gene_id']),]
  Ensembl.names__to__gene.symbols = tmp.maRt.table[,'external_gene_name']
  names(Ensembl.names__to__gene.symbols) = sapply(tmp.maRt.table[,'ensembl_gene_id'],function(y) { as.character(y) })
  
  
  
  if(length(Analyses.groups) == 0){
    cat(sprintf('\nmodel group "main analyses" has been created. %d GLM models are added to this group\n', length(Pars$Models.to.Test)))
    Analyses.groups[['main analysis']] = Pars$Models.to.Test
  }
  
  Analyses.groups %>% unlist %>% .[!duplicated(.)] -> models.from.Analyses.groups
  
  absent.models = models.from.Analyses.groups[!(models.from.Analyses.groups %in% Pars$Models.to.Test)]
  if(length(absent.models) > 0){
    cat(sprintf('\n%d GLM models which are pre-defined within the groups of analyses (see "Analyses group" sheet) are absent in Models.to.test parameter ("Parameters" sheet): %s\n', length(absent.models), toString(absent.models)))
    cat(sprintf('\nThese %d models will be added\n', length(absent.models)))
    Pars$Models.to.Test = c(Pars$Models.to.Test, absent.models)
  }
  
  remaining.models = Pars$Models.to.Test[!(Pars$Models.to.Test %in% models.from.Analyses.groups)]
  if(length(remaining.models) > 0){
    cat(sprintf('\nThe following GLM models are absent within the pre-defined analyses groups (see "Analyses groups" sheet): %s\n', toString(remaining.models)))
    cat(sprintf('\nmodel group "remaining analyses" has been created. These %d models are added to this group\n', length(remaining.models)))
    Analyses.groups[['remaining analyses']] = remaining.models
  }
  
  names(Pars$Models.to.Test) = NULL
  
  Pars$Analyses.groups = Analyses.groups

  
  #####
  ##### Adding MDS coordinates of samples as predictors, if this option (add.MDS.dims.as.predictors) is et TRUE
  
  if ((Pars$add.MDS.dims.as.predictors | Pars$DE.package=='edger') & is.null(ext.DE.data.files)){
    counts = samples_counts$counts
    # excluding meta-tags and low expression genes
    noint = rownames(counts) %in% c("__no_feature","__ambiguous",'__too_low_aQual','__not_aligned','__alignment_not_unique')
    cpms = cpm(counts)
    if (Pars$min.samples.with.sufficient.CPM == '{auto}') {
      min.samples.with.sufficient.CPM = Calc.min.high.CPM.samples.from.total.samples.count(dim(schema)[1])
    } else min.samples.with.sufficient.CPM = Pars$min.samples.with.sufficient.CPM
    #keep = rowSums(cpms > Pars$sufficient.CPM.to.analyze) >= min.samples.with.sufficient.CPM & !noint
    keep = rowSums(cpms > min(Pars$sufficient.CPM.to.analyze,1)) >= min(min.samples.with.sufficient.CPM,3) & !noint
    counts = counts[keep,]
    colnames(counts) = schema$'Sample names'
    d = DGEList(counts=counts)
    d = calcNormFactors(d,method=Pars$RNA.Seq.norm.method) #,"TMM","upperquartile","none"
    if(Pars$DE.package=='edger'){
      write.table.mod(cpm(d, normalized.lib.sizes = TRUE),sprintf('CPMs_%s.tsv',Pars$DE.package),sep='\t',quote = F)
    }
  }
  
  if (Pars$add.MDS.dims.as.predictors & is.null(ext.DE.data.files)){
    MDS.info = plotMDS(d,dim.plot = c(2,3))
    two.dim.MDS.info = cbind(MDS.info$x,MDS.info$y)
    MDS.info = plotMDS(d,dim.plot = c(1,2))
    Cm.three.dim.MDS.info = cbind(MDS.info$x,two.dim.MDS.info)
    Cm.three.dim.MDS.info = data.frame(Cm.three.dim.MDS.info)
    colnames(Cm.three.dim.MDS.info) = c('MDS dim1','MDS dim2','MDS dim3')
    Pars$Models.to.Test = c(Pars$Models.to.Test,colnames(Cm.three.dim.MDS.info))
    Startup.Data$Cm.three.dim.MDS.info = Cm.three.dim.MDS.info
    rm(counts)
    suppressWarnings(rm(noint))
    rm(cpms)
    rm(keep)
    rm(d)
  }
  if(is.null(ext.DE.data.files))  rm(samples_counts)
  
  if(Pars$DE.package=='deseq' & is.null(ext.DE.data.files)) {
    cond = rep(c(0,1),as.integer(dim(schema)[1]/2))
    cond = c(cond,rep(0, dim(schema)[1] - length(cond))) ### artificial 'conditions'
    sampleTable <- data.frame(sampleName = schema$'Sample names',
                              fileName = basename(file.names),
                              condition = factor(cond))
    dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = dirname(file.names)[1],
                                      design= ~ condition)
    
    dds <- estimateSizeFactors(dds)
    norm.counts = cpm(counts(dds,normalized=T))
    write.table.mod(norm.counts,sprintf('CPMs_%s.tsv',Pars$DE.package),sep='\t',quote = F)
    
  }
  
  # testing for the presence of Python
  if (disable.excel) Pars$Create.Excel.results = FALSE
  if (Pars$Create.Excel.results & !disable.excel){
    if(suppressWarnings(system('python3 -c "import xlsxwriter"',ignore.stdout = T, ignore.stderr = T)) == 0){
      Pars$python.bin = "python3"
    } else if (suppressWarnings(system('python -c "import xlsxwriter"',ignore.stdout = T, ignore.stderr = T)) == 0){
      Pars$python.bin = "python"
    } else {
      stop('Python3 with xlsxwriter package is required for creating Excel workbooks containing reports. If python is already installed, open command line and type smth. like "python -m pip install xlsxwriter" (or "python3 -m pip install xlsxwriter"). Otherwise, set "Create.Excel.results" as False in startup RTrans parameters xlsx-file or run "Prepare(...)" with additional "disable.excel=TRUE" argument')
      #Pars$Create.Excel.results = FALSE
    }
    
  }
  
  Startup.Data$Pars = Pars
  Startup.Data$General.maRt.table = General.maRt.table
  Startup.Data$Ensembl.names__to__gene.symbols = Ensembl.names__to__gene.symbols
  Startup.Data$topGO.Expression.Profiles___Custom.GO.terms = topGO.Expression.Profiles___Custom.GO.terms
  Startup.Data$GO.clusterProfiler.Expression.Profiles___Custom.DB.entries = topGO.Expression.Profiles___Custom.GO.terms
  Startup.Data$KEGG.clusterProfiler.Expression.Profiles___Custom.DB.entries = KEGG.clusterProfiler.Expression.Profiles___Custom.DB.entries
  Startup.Data$Reactome.clusterProfiler.Expression.Profiles___Custom.DB.entries = Reactome.clusterProfiler.Expression.Profiles___Custom.DB.entries
  Startup.Data$Custom.Pathways.to.Visualize = Custom.Pathways.to.Visualize
  Startup.Data$Omit.pathways = Omit.pathways
  Startup.Data$Completed.steps.file = Completed.steps.file
  Startup.Data$Parameters.xlsx.file.name = Parameters.xlsx.file.name
  
  Startup.Data$Startup.hashmd5 = digest(toString(c(Startup.Data$Pars,Startup.Data$General.maRt.table,Startup.Data$topGO.Expression.Profiles___Custom.GO.terms,Startup.Data$Custom.Pathways.to.Visualize,Startup.Data$Omit.pathways)))
  
  cat('\nPreparing completed. Now you can run GLM/DE analysis\n')
  return(Startup.Data)
  
  #####
}

#Startup.Data = startup.ata
#GLM.model = "Tg_status [cells]"

Analyze.GLM = function(Startup.Data, GLM.model, bypass.if.completed = FALSE, forced.parameters = NULL){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  ####################################################################
  #### Analysing GLM/differential expression
  
  GLM.model = gsub("~","",GLM.model,fixed = TRUE)
  
  cat(sprintf("\nProcessing model  \"~ %s\"...",GLM.model))
  setwd(Pars$results.dir)
  
  Analysis.Data = new("RTransAnalysisData")
  
  
  ## extracting names of all used individual predictors
  conv.GLM.model = Replace.Preds.with.Abbrebiations(GLM.model,Pars)
  Current.Predictors.list = Extract.predictors(conv.GLM.model)
  Current.Predictors.list = sapply(Current.Predictors.list, function(x){ r = Revert.Preds.from.Abbrebiations(x,Pars) })
  
  #Analysis.name = sprintf("Associations with %s",GLM.model)
  Analysis.name = Cleanup.Model.Name(GLM.model)
  Analysis.Data$Analysis.name = Analysis.name
  Analysis.Data$GLM.model = GLM.model
  
  tryCatch(expr = { schema = suppressWarnings(readxl::read_excel(Startup.Data$Parameters.xlsx.file.name ,sheet='Sample setup',col_names = T))
  }, error = function (err){
    stop(sprintf('Sheet "Sample setup" is not found in workbook %s. Use an example in ./RTrans.base/Parameters.example.xlsx',Startup.Data$Parameters.xlsx.file.name))
  })
  schema %<>% as.data.frame
  commented.lines = sapply(X = schema[,1], FUN = function (x) substring(x,first=1,last=1)) %in% '#'
  schema = schema[!commented.lines,]
  schema = schema[,!apply(schema, 2, function(x) all(is.na(x)))]
  if (Pars$add.MDS.dims.as.predictors) schema = cbind(schema, Startup.Data$Cm.three.dim.MDS.info)
  rownames(schema) = schema$'Sample names'
  
  
  #cor.test(data.matrix(schema[,'Gclc']),data.matrix(schema[,'Mst57D']),method = 'spearman')
  

  # all.preds = c()
  # conv.substrate = c(Current.Predictors.list, Pars$Models.to.Test)
  # conv.substrate = sapply(conv.substrate, function(x){ r = Replace.Preds.with.Abbrebiations(x,Pars) })
  # for (tmp.Pred.set in conv.substrate)  all.preds = append(all.preds, Extract.predictors(tmp.Pred.set))
  # all.preds = levels(factor(all.preds))
  # all.preds = sapply(all.preds, function(x){ r = Revert.Preds.from.Abbrebiations(x,Pars) })
  # if(sum.mod(!(all.preds %in% colnames(schema))) > 0){
  #   stop(sprintf('The following predictors are not found in "Sample Setup" sheet: "%s". Please check "Sample Setup" row names, "Models.to.Test" parameter in "Parameters" sheet or currently used GLM.model',toString(all.preds[!(all.preds %in% colnames(schema))])))
  # }
  # 
  
  # names(all.preds) = NULL
  # all.preds
  # colnames(schema)
  
  
  for (Pred in Current.Predictors.list){
    Discard = (data.frame(schema,check.names = FALSE)[,Pred] %in% 'discard') | is.na(data.frame(schema,check.names = FALSE)[,Pred])
    schema = schema[!Discard,]
  }
  rownames(schema) = schema$'Sample names'
  
  tryCatch(expr = {
    Components = data.matrix(schema[,Current.Predictors.list]) 
  }, error = function (err){
    Components = data.matrix(as.numeric(schema[,Current.Predictors.list]))
  })
  

  ## converting 'character' to 'double'
  if (typeof(Components) == 'character'){
    if (dim(Components)[2] == 1){
      Components = data.matrix(apply(Components,MARGIN = 1,as.numeric))
      colnames(Components) = Current.Predictors.list
    } else {
      nm = colnames(Components)
      Components = t(data.matrix(apply(Components,MARGIN = 1,as.numeric)))
      colnames(Components) = nm
    }
  }
  
  if (Pars$Normalize.predictors) {
    predictors.range = apply(X = Components, MARGIN = 2, FUN = max) - apply(X = Components, MARGIN = 2, FUN = min)
    for (i in 1:dim(Components)[2]) {
      Components[,i] = Components[,i] / predictors.range[i]
    }
  }
  colnames(Components) = Current.Predictors.list
  rownames(Components) = rownames(schema)
  Predictors.count = length(Current.Predictors.list)
  Main.Component = as.numeric(Components[,Predictors.count])
  
  
  #### testing if this step has been preiously completed
  step.hashmd5 = digest(cbind(schema$'Sample names',Components)) ## + make test if there are ready files
  Analysis.Data$step.hashmd5 = step.hashmd5
  
  if (bypass.if.completed & Read.Completed.steps.status(sprintf("%s.%s.%s.complete",GLM.model,Startup.Data$Startup.hashmd5,step.hashmd5),Startup.Data$Completed.steps.file)){
    message(sprintf('\nGLM model "%s" processing was previously completed. Bypassing...',GLM.model))
    return(Analysis.Data)
  }
  ####
  
  cat(sprintf("\n%s:   loading and processing RNA-Seq data...",GLM.model))
  
  # setting design matrix
  GLM.formula = sprintf("~ %s", GLM.model)
  #Current.Predictors.list = sprintf('&&&%s&&&',Current.Predictors.list)
  #Current.Predictors.list = mgsub("&&&","",Current.Predictors.list)
  GLM.formula = mgsub(Current.Predictors.list,
                      sprintf('as.numeric(Components[,"%s"])',Pars$available.preds_to_abbreviations[Current.Predictors.list]),
                      GLM.formula,fixed = TRUE)
  
  GLM.formula = Revert.Preds.from.Abbrebiations(GLM.formula, Pars)
  #for (Pred in Current.Predictors.list){
  #  GLM.formula = gsub(Pred,sprintf('as.numeric(Components[,"%s"])',Pred),GLM.formula,fixed = T)
  #}
  
  Single.binary.predictor = FALSE
  Single.binary.predictor__minimal.group.size = NULL
  if (length(Current.Predictors.list) == 1){
    if(length(levels(as.factor(as.numeric(Components[,Pred])))) == 2){
      cat(sprintf('\nRTrans has detected that GLM model has only one binary predictor'))
      Single.binary.predictor__minimal.group.size = min(table(as.numeric(Components[,Pred])))
      Single.binary.predictor = TRUE
    }
  }
  
  # loading files with counts
  
  if ("File names" %in% colnames(schema)){
    file.names = sprintf("%s/%s",Pars$counts.dir,schema$'File names')
  } else {
    file.names = sprintf("%s/%s.counts",Pars$counts.dir,schema$'Sample names')
  }
  if (sum.mod(!file.exists(file.names)) > 0){
    #warning(sprintf('Warning; The following files do not exist: %s',toString(file.names[!file.exists(file.names)])))
    cat(sprintf('\nTotal %d of %d *.counts files do not exist',sum.mod(!file.exists(file.names)),length(file.names)))
    if (sum.mod(file.exists(file.names)) == 0) stop(sprintf('No *.counts files found. Fix "Sample setup" worksheet and check directory "%s" ',Pars$counts.dir))
    schema=schema[file.exists(file.names),]
    rownames(schema) = schema$'Sample names'
    # print(dim(schema))
    tmp = colnames(Components)
    Components = data.matrix(Components[file.exists(file.names),])
    Main.Component = Main.Component[file.exists(file.names)]
    colnames(Components) = tmp
    file.names = file.names[file.exists(file.names)]
  }
  

  eval(parse(text = sprintf("design <- model.matrix(%s)",GLM.formula)))
  
  current.GLM.results.dir = sprintf("%s/~ %s, results",Pars$results.dir,Analysis.name)
  Analysis.Data$current.GLM.results.dir = current.GLM.results.dir
  dir.create(current.GLM.results.dir, showWarnings = FALSE)
  setwd(current.GLM.results.dir)
  
  Main.Predictor.name = colnames(Components)[Predictors.count]
  Analysis.Data$Main.Predictor.name = Main.Predictor.name
  Analysis.Data$Main.Component = Main.Component
  
  #use.exact.test.in.binary.predictors = FALSE
  
  if(!(Pars$DE.package %in% c('deseq'))){
    samples_counts=suppressMessages(readDGE(file.names))
    counts = samples_counts$counts
    # naming 'counts'
    colnames(counts) = schema$'Sample names'
    
    # calculating stats and excluding meta-tags
    meta.tags = c("__no_feature","__ambiguous",'__too_low_aQual','__not_aligned','__alignment_not_unique')
    meta.tags = meta.tags[meta.tags %in% rownames(counts)]
    if (length(meta.tags) > 0){
      mapping.stats.table = array(dim = c(dim(counts)[2],length(meta.tags) + 2))
      rownames(mapping.stats.table) = schema$'Sample names'
      colnames(mapping.stats.table) = c(meta.tags,'total hits','useful reads')
      for (mt in meta.tags){
        mapping.stats.table[schema$'Sample names',mt] = counts[mt,schema$'Sample names']
      }
      mapping.stats.table[schema$'Sample names','total hits'] = apply(counts,2,sum)
      
      noint = rownames(counts) %in% meta.tags
      counts = counts[!noint,]
      mapping.stats.table[schema$'Sample names','useful reads'] = apply(counts,2,sum)
      write.table.mod(mapping.stats.table,'mapping stats.tsv',sep='\t')
    }
    
    # excluding low expression genes
    cpms = cpm(counts)
    if (Pars$min.samples.with.sufficient.CPM == '{auto}') {
      if (Single.binary.predictor){
        min.samples.with.sufficient.CPM = Calc.min.high.CPM.samples.from.total.samples.count(dim(schema)[1],Single.binary.predictor__minimal.group.size)
      } else {
        min.samples.with.sufficient.CPM = Calc.min.high.CPM.samples.from.total.samples.count(dim(schema)[1])
      }
    } else min.samples.with.sufficient.CPM = Pars$min.samples.with.sufficient.CPM
    
    cat(sprintf('\nSetting CPM limit %.1f for at least %d samples',Pars$sufficient.CPM.to.analyze,min.samples.with.sufficient.CPM))
    
    keep = rowSums(cpms > Pars$sufficient.CPM.to.analyze) >= min.samples.with.sufficient.CPM
    cat(sprintf('\n%d of %d genes passed CPM threshold\n',sum.mod(keep),length(keep)))
    counts = counts[keep,]
  
    # plotting histograms - before normalization
      
    plot.subtitle = ''
    if (length(levels(factor(Main.Component))) == 2) plot.subtitle = sprintf('color indicates %s (quantity; blue-orange)',Main.Predictor.name)
    if (length(levels(factor(Main.Component))) == 3) plot.subtitle = sprintf('color indicates %s (quantity; blue-green-orange)',Main.Predictor.name)
    if (length(levels(factor(Main.Component))) >= 4) plot.subtitle = sprintf('color indicates %s (quantity; blue-green-orange gradient)',Main.Predictor.name)
    
    col.variety <- colorRampPalette(c("#0b00c4","#02babc","#25ad00","#cfcd00","#ff8400"))(512)
    Main.Component.mod = Main.Component - min(Main.Component)
    cl = col.variety[1 + as.integer(Main.Component.mod/max(Main.Component.mod)*511)]
    
    cpms = cpm(counts)
    png.mod(filename = sprintf('CPM density, non-normalized.png'),
        units="in",width=8.23,height=6.3,pointsize=12,
        res=600)
    plot(density(log(cpms[,1])),col=cl[1],xlim=c(-4,10),main = 'Log2(CPM) density, non-normalized', sub = plot.subtitle)
    for (x in 2:dim(cpms)[2]){
      lines(density(log(cpms[,x])),col=cl[x])
    }
    dev.off()
    
    
    # creating DGEList object
    if(Single.binary.predictor & Pars$use.exact.test.in.binary.predictors){
      d = DGEList(counts=counts, group = Main.Component)
    } else {
      d = DGEList(counts=counts)
    }
    
    
    # Calculating Norm Factors
    #RNA.Seq.norm.method = 'RLE'  #,"TMM","upperquartile","none","RLE"
    d = calcNormFactors(d,method=Pars$RNA.Seq.norm.method) #,"TMM","upperquartile","none","RLE"
    Analysis.Data$edgeR_d = d
    
    # plotting histograms - after normalization
    cpms = cpm(d)
    png.mod(filename = sprintf('CPM density, normalized (%s).png',Pars$RNA.Seq.norm.method),
        units="in",width=8.23,height=6.3,pointsize=12,
        res=600)
    plot(density(log(cpms[,1])),col=cl[1],xlim=c(-4,10),main = sprintf('Log2(CPM) density, normalized (%s)',Pars$RNA.Seq.norm.method), sub = plot.subtitle)
    for (x in 2:dim(cpms)[2]){
      lines(density(log(cpms[,x])),col=cl[x])
    }
    dev.off()
  }
  
  simple.DE.mode = FALSE
  
  if(Pars$DE.package=='edger'){
    
    # Creating MDS/PCA plots
    
    #col.variety = c("darkgreen","blue","red","darkgoldenrod2","chartreuse2","darkorchid2","black","gray59","lightsteelblue1","orange4","sienna3","bisque3")
    ### gradient: blue-teal-green-yellow-orange
    col.variety<-colorRampPalette(c("#0b00c4","#02babc","#25ad00","#cfcd00","#ff8400"))(length(levels(factor(Main.Component))))
    
    plot.title = sprintf('MDS plot, %s',Main.Predictor.name)
    plot.subtitle = ''
    if (length(levels(factor(Main.Component))) == 2) plot.subtitle = sprintf('color indicates %s (blue-orange)',Main.Predictor.name)
    if (length(levels(factor(Main.Component))) == 3) plot.subtitle = sprintf('color indicates %s (blue-green-orange)',Main.Predictor.name)
    if (length(levels(factor(Main.Component))) >= 4) plot.subtitle = sprintf('color indicates %s (blue-green-orange gradient)',Main.Predictor.name)
    
    
    if(dim(d$counts)[2] > 2){
      tryCatch(expr = {
        plotMDS(d,dim.plot = c(1,2), col = col.variety[factor(Main.Component)])
        title(main = plot.title, sub = plot.subtitle)
        Create.MDS <<- TRUE
      }, error = function (err){
        warning('MDS plots creating failed')
        Create.MDS <<- FALSE
      })
    } else {
      Create.MDS <<- FALSE
      msg = sprintf('MDS plots will not be created bacause only %d samples present',dim(d$counts)[2])
      cat(sprintf('\n%s\n', msg))
      #warning(msg)
    }
    
    
    if(Create.MDS){
      pdf(file = 'MDS.dim.1-2.pdf',
          width=18.23,height=16.3,pointsize=12)
      MDS.info = plotMDS(d,dim.plot = c(1,2), col = col.variety[factor(Main.Component)])
      two.dim.MDS.info = cbind(MDS.info$x,MDS.info$y)
      title(main = plot.title, sub = plot.subtitle)
      dev.off()
      
      
      if(dim(d)[2] > 3){
        pdf(file = 'MDS.dim.2-3.pdf',
            width=18.23,height=16.3,pointsize=12)
        MDS.info = plotMDS(d,dim.plot = c(2,3), col = col.variety[factor(Main.Component)])
        title(main = plot.title, sub = plot.subtitle)
        dev.off()
        three.dim.MDS.info = cbind(two.dim.MDS.info,MDS.info$y)
        colnames(three.dim.MDS.info) = c('MDS dim1','MDS dim2','MDS dim3')
        write.table.mod(x = three.dim.MDS.info,file = 'MDS info.tsv',sep = '\t')
      }
      
      png.mod(filename = 'MDS.dim.1-2.png',units="in",res=600,
          width=8.23,height=6.3,pointsize=12)
      plotMDS(d,dim.plot = c(1,2), col = col.variety[factor(Main.Component)])
      title(main = plot.title, sub = plot.subtitle)
      dev.off()
      
      if(dim(d)[2] > 3){
        png.mod(filename = 'MDS.dim.2-3.png',units="in",res=600,
            width=8.23,height=6.3,pointsize=12)
        a = plotMDS(d,dim.plot = c(2,3), col = col.variety[factor(Main.Component)])
        title(main = plot.title, sub = plot.subtitle)
        dev.off()
      }
    }
    
    ## Calculating distance matrices
    C.Create.Distance.matrices = Pars$Create.Distance.matrices
    if (Pars$Create.Distance.matrices == '{auto}'){
      C.Create.Distance.matrices = (dim(cpms)[2] <= 50)
    }
    
    if (C.Create.Distance.matrices){
      if (Pars$Distance.method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")){
        dist.matrix = data.frame(dist(cpms,Pars$Distance.method))
      } else {
        N = dim(cpms)[2]
        dist.matrix = data.frame(array(dim = c(N,N)))
        colnames(dist.matrix) = colnames(cpms)
        rownames(dist.matrix) = colnames(cpms)
        
        #F = colnames(cpms)[1] # this is for debug only
        ready = 0
        for (F in colnames(cpms)){
          #S = colnames(cpms)[2]
          for (S in colnames(cpms)){
            if (Pars$Distance.method == 'canberra.weighted') dist.matrix[F,S] = Weighted.Canberra.dist(cpms[,F],cpms[,S],Pars$Min.CPM.to.include.in.dist,Pars$Canberra.weight.power)
            if (Pars$Distance.method == '1-cor') dist.matrix[F,S] = 1 - suppressWarnings(cor.test(cpms[,F],cpms[,S],method = 'spearman')$estimate)
            ready = ready+1
            cat(sprintf('\rCalculating dist matrices. Completed %.2f%%',ready/N/N*100))
            
            #Weighted.Canberra.dist(c(1,2,3,4,5,6,7,8,9,10),c(1,20,30,4,5,6,7,8,9,10),min.value.to.calc.dist = 0.99)
          }
        }
      }

      tmp=as.dist(dist.matrix)
      dend=hclust(tmp)
      dend=as.dendrogram(dend)
      
      dend.colors = col.variety[factor(Main.Component)]
      names(dend.colors) = schema$'Sample names'
      labels_colors(dend) = dend.colors[labels(dend)]
      png.mod(filename = sprintf('Clustering, dist as %s.png',Pars$Distance.method),units="in",res=600,
          width=8.23,height=6.3,pointsize=12)
      plot(dend,main = sprintf('%s, clustering dendrogram [%s]',GLM.model,Pars$Distance.method),
           sub = sprintf('distances are calculated with "%s" method',Pars$Distance.method))
      dev.off()
      ##hang=-1
      
      
      #Current.Predictors.list = c('Age','Genotype','Gender')
      info.block = data.frame(as.data.frame(schema)[colnames(cpms),Current.Predictors.list])
      colnames(info.block) = Current.Predictors.list
      info.block.mod = cbind(array(dim = c(dim(info.block)[2],dim(info.block)[2])),t(info.block))
      colnames(info.block.mod) = colnames(cbind(info.block,dist.matrix))
      class(info.block.mod) <- "numeric"
      write.table.mod(x=rbind(info.block.mod,
                          cbind(info.block,dist.matrix)),
                  file='distance matrix.tsv',sep="\t",na = "")
    }
    
    # evaluating dispersion
    # if(TTest.instead.of.GLM){
    #   estimateGLMCommonDisp
    # } else {
    #d=estimateDisp(d,design)
	  # print(dim(design))
	  # print(dim(d$counts))
    
    if(dim(design)[1] < 3){
      simple.DE.mode = TRUE
      msg = sprintf('Sample setup with only %d samples is found. Dispersion, p-value, FDR and LR values cannot be calculated. Score-based P-value "evaluation" will be performed instead.',dim(design)[1])
      #cat(msg)
      warning(msg)
      if(!(all(Main.Component == c(0,1)) || all(Main.Component == c(1,0)))){
        stop('Only comparison of two samples is supported in "simple" mode')
      }
      cpm.add = 0.4
      norm.counts = cpm(d, normalized.lib.sizes = TRUE)
      LogFCs = log2((norm.counts[,2] + cpm.add)/(norm.counts[,1] + cpm.add))
      LogCPMs = log((apply(norm.counts,1,sum)/dim(norm.counts)[2]) + 0.01)
      
      #p.values = 10^((-1)*(((2*abs(LogFCs) + 2 + LogCPMs)/3)**1.4))
      p.values = 10^((-1) * abs(LogFCs)^0.5 * (abs(LogFCs) + 1 + 1.4*LogCPMs)/2.5  )
      #plot(density(p.values))
      p.values = sapply(p.values, function(x) { min(1,x)})
      FDRs = p.adjust(p.values,method = 'BH')

      Scores = ((-1)*log2(p.values + 1e-300))^1.5* 10 * LogFCs
      Scores = 10*abs(Scores)**0.5
      
      LRs = Scores/5
      
      Stats.table = cbind(as.data.frame(LogFCs), as.data.frame(LogCPMs), as.data.frame(LRs), as.data.frame(p.values), as.data.frame(FDRs))
      colnames(Stats.table) = c('logFC','logCPM','LR','PValue','FDR')
    } else {
      simple.DE.mode = FALSE
      
      if(Single.binary.predictor & Pars$use.exact.test.in.binary.predictors){
        d = estimateDisp(d)
      } else {
        d=estimateGLMCommonDisp(d,design)
        d=estimateGLMTrendedDisp(d,design)
        d=estimateGLMTagwiseDisp(d,design)
      }
      # maximum sipersion should be taken in order to preserve reliability
      
      # Plotting Biological coefficient of variation plot
      png.mod(filename = 'MeanVar.png',units="in",res=600,
          width=8.23,height=6.3,pointsize=12)
      plotMeanVar(d, show.tagwise.vars = TRUE, NBline = TRUE)
      dev.off()
      png.mod(filename = 'BCV.png',units="in",res=600,
          width=8.23,height=6.3,pointsize=12)
      plotBCV(d)
      dev.off()
      
      if(Single.binary.predictor & Pars$use.exact.test.in.binary.predictors){
        ## exact test
        et <- exactTest(d)
        tt <- topTags(et,n=nrow(d))
        norm.counts = cpm(d, normalized.lib.sizes = TRUE)
        Stats.table = tt$table
        Stats.table %<>% mutate(LR = -log2(PValue + 1e-250))
        rownames(Stats.table) = rownames(tt$table)
      } else if (Pars$use.QLfit){
        ## GLM approximation
        fit <- glmQLFit (d,design)
        lrt <- glmQLFTest(fit)
        #lrt <- glmLRT(fit,coef=1)
        tt <- topTags(lrt,n=nrow(d))
        norm.counts = cpm(d, normalized.lib.sizes = TRUE)
        Stats.table = tt$table
      } else {
        ## GLM approximation
        fit <- glmFit(d,design)
        lrt <- glmLRT(fit)
        #lrt <- glmLRT(fit,coef=1)
        tt <- topTags(lrt,n=nrow(d))
        norm.counts = cpm(d, normalized.lib.sizes = TRUE)
        Stats.table = tt$table
      }
      
    }    
  }  else if(Pars$DE.package=='deseq') {
    cond = as.numeric(Components[,Current.Predictors.list[1]])
    sampleTable <- data.frame(sampleName = schema$`Sample names`,
                              fileName = basename(file.names),
                              condition = cond)
    dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = dirname(file.names)[1],
                                      design= ~ condition)
    
    dds <- dds[ rowSums(counts(dds)) > 1, ]

    dds <- DESeq(dds)
    res <- results(dds)
    plotMA(res, main="DESeq2", ylim=c(-2,2))
    res = res[order(res$pvalue),]
    
    norm.counts = cpm(counts(dds[rownames(res),],normalized=T))
    
    Stats.table = cbind(
      data.frame(res[,'log2FoldChange']),
      data.frame(log2(apply(norm.counts,1,mean))),
      data.frame(res[,'stat']*10),
      res[,c('pvalue','padj')])
    colnames(Stats.table) = c('logFC','logCPM','LR','PValue','FDR')
    tt = Stats.table
  }  else {
    # cond = as.numeric(Components[,Current.Predictors.list[1]])
    # sampleTable <- data.frame(sampleName = schema$`Sample names`,
    #                           fileName = basename(file.names),
    #                           condition = cond)
    # dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
    #                                   directory = dirname(file.names)[1],
    #                                   design= ~ condition)
    # 
    # dds <- dds[ rowSums(counts(dds)) > 1, ]
    # 
    # 
    # # excluding low expression genes
    # cpms = cpm(counts(dds))
    # if (Pars$min.samples.with.sufficient.CPM == '{auto}') { min.samples.with.sufficient.CPM = min(4,max(1,as.integer(dim(schema)[1]/2.3 + 0.5)))
    # } else min.samples.with.sufficient.CPM = Pars$min.samples.with.sufficient.CPM
    # 
    # keep = rowSums(cpms > Pars$sufficient.CPM.to.analyze) >= min.samples.with.sufficient.CPM
    # cat(sprintf('\n%d of %d genes passed CPM threshold',sum.mod(keep),length(keep)))
    # dds = dds[keep,]
    # dds <- DESeq(dds)
    #norm.counts = cpm(counts(dds,normalized=T))
    
    norm.counts = cpm(d, normalized.lib.sizes = TRUE)
    
    trimming = 0.03
    
    Cor.table = data.matrix(array(dim=c(dim(norm.counts)[1],9)))
    colnames(Cor.table) = c('Spearman r','Spearman P','Pearson r','Pearson P','logFC','logCPM','LR','PValue','FDR')
    #Pars$DE.package = 'combined_cor'
    for (j in 1:dim(norm.counts)[1]){
      if(Pars$DE.package %in% c('combined_cor', 'spearman')){
        ct = suppressWarnings(cor.test(Main.Component,norm.counts[j,],method = 'spearman',alternative = 'two.sided'))
        Cor.table[j,'Spearman r'] = ct$estimate
        Cor.table[j,'Spearman P'] = ct$p.value
      }
      if(Pars$DE.package %in% c('combined_cor', 'pearson')){
        ct = suppressWarnings(cor.test(Main.Component,norm.counts[j,],method = 'pearson',alternative = 'two.sided'))
        Cor.table[j,'Pearson r'] = ct$estimate
        Cor.table[j,'Pearson P'] = ct$p.value
      }
    }
    
    
    rel_sd_trim = function(x) sd_trim(x,trim = trimming,const = F)/mean(x)+0.001
    
    #library(chemometrics)
    if(Pars$DE.package == 'spearman'){
      Cor.table[,'PValue'] = Cor.table[,'Spearman P']
      Cor.table[,'logFC'] = 2.5*Cor.table[,'Spearman r'] * log(1 + 4*apply(norm.counts,1,rel_sd_trim))
      #Cor.table[,'logFC'] = log((((1 + 1.5*abs(Cor.table[,'Spearman r']))**2
        #*apply(norm.counts,1,rel_sd_trim))**0.33)*3,2)*sign(Cor.table[,'Spearman r'])
    } else if(Pars$DE.package == 'pearson'){
      Cor.table[,'PValue'] = Cor.table[,'Pearson P']
      Cor.table[,'logFC'] = 2.5*Cor.table[,'Pearson r'] * log(1 + 4*apply(norm.counts,1,rel_sd_trim))
      #Cor.table[,'logFC'] = log((((1 + 1.5*abs(Cor.table[,'Pearson r']))**2
      #                            *apply(norm.counts,1,rel_sd_trim))**0.33)*3,2)*sign(Cor.table[,'Pearson r'])
    } else {  ## combined_cor
      integrate.corr.rs <- function(x){
        if(length(levels(factor(sign(x)))) > 1) return(100)
        return (((abs(x[1])*abs(x[1])*abs(x[2]))^0.33)*sign(x[1]))
      }
      Rs = apply(cbind(Cor.table[,'Spearman r'],Cor.table[,'Pearson r']),1,integrate.corr.rs)  ## 'joint' Rs
      Cor.table[,'PValue'] = apply(cbind(Cor.table[,'Spearman P'],Cor.table[,'Pearson P']),1,max)
      Cor.table[,'logFC'] = 2.5*Rs * log(1 + 4*apply(norm.counts,1,rel_sd_trim))
      #Cor.table[,'logFC'] = log((((1 + 1.5*abs(Rs))**2
      #                            *apply(norm.counts,1,rel_sd_trim))**0.33)*3,2)*sign(Rs)
    }
    Cor.table[,'LR'] = abs(Cor.table[,'logFC']) * abs(log(Cor.table[,'PValue'],2))*2
    Cor.table[,'logCPM'] = log2(apply(norm.counts,1,mean))
    Cor.table[,'FDR'] = p.adjust(Cor.table[,'PValue'],method = 'BH')

    #res <- results(dds)
    #plotMA(res, main="DESeq2", ylim=c(-2,2))
    #res = res[order(res$pvalue),]

    rownames(Cor.table) = rownames(norm.counts)
    Stats.table = Cor.table[order(Cor.table[,'LR'],decreasing = T),c('logFC','logCPM','LR','PValue','FDR')]
    tt = Stats.table
  }
  
  ##############################################################
  ############ Creating signle output results table
  
  # combining FRD, LR, P-value info and Read Count Info
  #Scores = abs(tt$table['logFC'])**0.7 * log(tt$table['FDR'])*(-1)
  #colnames(Scores) = c('Score')
  
  
  cat(sprintf("\n%s:   calculating correlations, scores and writing expression tables...",Analysis.Data$GLM.model))
  
  #all.preds = c()
  #for (tmp.Pred.set in c(Analysis.Data$GLM.model,Pars$Models.to.Test))  all.preds = append(all.preds,Extract.predictors(tmp.Pred.set))
  
  #all.preds = levels(as.factor(all.preds))
  
  if (Pars$add.MDS.dims.as.predictors) {
    Analysis.Data$predictor.values = suppressWarnings(as.data.frame(schema, row.names = rownames(schema))[colnames(norm.counts),])
    Analysis.Data$predictor.values.nr = Eliminate.Redundant.Predictors.Info(Analysis.Data$predictor.values)
    all.preds.nr = colnames(Analysis.Data$predictor.values.nr)
    
    Counts.Table = suppressWarnings(rbind(t(as.data.frame(schema, row.names = rownames(schema))[colnames(norm.counts), all.preds.nr]),norm.counts[rownames(Stats.table),]))
    # Counts.Table = suppressWarnings(rbind(t(data.matrix(schema)[colnames(norm.counts),all.preds]),norm.counts[rownames(tt),]))
    #} else Counts.Table = rbind(t(as.data.frame(schema,row.names = rownames(schema))[colnames(norm.counts),all.preds]),norm.counts[rownames(tt),])
  } else {
    all.preds = colnames(schema)[!(colnames(schema) %in% 'Sample names')]
    Analysis.Data$predictor.values = suppressWarnings(data.matrix(schema)[colnames(norm.counts), all.preds])
    if (length(all.preds) == 1){
      Analysis.Data$predictor.values = as.data.frame(Analysis.Data$predictor.values)
      colnames(Analysis.Data$predictor.values) = c(all.preds)
    }
    Analysis.Data$predictor.values.nr = Eliminate.Redundant.Predictors.Info(Analysis.Data$predictor.values)
    
    all.preds.nr = colnames(Analysis.Data$predictor.values.nr)
    Counts.Table = suppressWarnings(rbind(t(data.matrix(schema)[colnames(norm.counts), all.preds.nr]), norm.counts[rownames(Stats.table),]))
  }
  
  #colnames(schema)
  ###Stats.table = cbind(tt$table,Scores)
  #head(Counts.Table)
  
  
  # predictor.values = Analysis.Data$predictor.values

  a = array(dim = c(length(all.preds.nr) ,dim(Stats.table)[2]))
  rownames(a) = all.preds.nr
  colnames(a) = colnames(Stats.table)
  
  Stats.table.ext = rbind(a,Stats.table)
  
  ## correlations with predictors
  Cor.table = data.matrix(array(dim=c(dim(norm.counts)[1],4)))
  colnames(Cor.table) = c('Spearman r','Spearman P','Pearson r','Pearson P')
  if(!simple.DE.mode){
    for (j in 1:dim(norm.counts)[1]){
      ct = suppressWarnings(cor.test(Main.Component,norm.counts[rownames(tt)[j],],method = 'spearman',alternative = 'two.sided'))
      Cor.table[j,1] = ct$estimate
      Cor.table[j,2] = ct$p.value
      ct = suppressWarnings(cor.test(Main.Component,norm.counts[rownames(tt)[j],],method = 'pearson',alternative = 'two.sided'))
      Cor.table[j,3] = ct$estimate
      Cor.table[j,4] = ct$p.value
    }
  }
  #length(norm.counts[rownames(tt)[j],])
  a = array(dim = c(length(all.preds.nr) ,4))
  
  Cor.table.ext = rbind(a,Cor.table)
  
  
  # genes.anno <- getBM(attributes=c(attributes,"gene_biotype"), filters="ensembl_gene_id",
  #                     values=rownames(tt), mart=mart, curl = curl.handle)
  # rownames(genes.anno) = genes.anno[,"ensembl_gene_id"]
  
  Mart.table = Startup.Data$General.maRt.table[rownames(Stats.table),c("external_gene_name","gene_biotype","description")]
  mart.absent = !(rownames(Stats.table) %in% rownames(Startup.Data$General.maRt.table))
  for (n in 1:length(mart.absent)){
    if (mart.absent[n]){
      rownames(Mart.table)[n] = rownames(Stats.table)[n]
    }
  }
  
  #sum(grepl("^NA",rownames(Mart.table)))
  
  if (Pars$Use.Info.Table){
    Mart.table = cbind(Mart.table,data.frame(Startup.Data$Info.table[rownames(Stats.table),"RefSeq_Summary"]))
    colnames(Mart.table) = c("Gene name","Biotype","description","RefSeq_Summary")
  }  else{
    colnames(Mart.table) = c("Gene name","Biotype","description")
  }
  a = array(dim = c(length(all.preds.nr) ,dim(Mart.table)[2]))
  rownames(a) = all.preds.nr
  colnames(a) = colnames(Mart.table)
  Mart.table.ext = rbind(a,Mart.table)
  
  
  
  
  #head(Stats.table)
  ## Calculating scores
  if(!simple.DE.mode){
    if (Predictors.count > 1){
      Scores = (-1)*log2(Stats.table[,"PValue"] + 1e-300)*(Stats.table[,"LR"]+10)*((0.2+abs(Cor.table[,"Spearman r"]))^0.4)*Stats.table[,"logFC"]
    } else{
      Scores = (-1)*log2(Stats.table[,"PValue"] + 1e-300)*(Stats.table[,"LR"]+10)*Cor.table[,"Spearman r"]
    }
    #Scores = (-1)*log2(Stats.table[,"PValue"])*Stats.table[,"logFC"]
    Scores = 10*abs(Scores)**0.5
  }
  #med = median(Scores)
  #Scores = array(Scores,dim = c(length(Scores),1))
  #Scores = apply(Scores,1,function(x) max(0,x-med*5))
  Scores = array(Scores,dim = c(length(Scores),1))
  a = array(dim = c(length(all.preds.nr) ,dim(Scores)[2]))
  Scores.table.ext = rbind(a,Scores)
  colnames(Scores.table.ext) = c('Score')
  
  
  #sum(grepl("^NA",rownames(ResTable.gene.part)))
  
  ResTable = cbind(Mart.table.ext,Scores.table.ext,Stats.table.ext,Cor.table.ext,Counts.Table)
  colnames(ResTable) = gsub(pattern = "Pearson.",replacement = "Pearson ",x = colnames(ResTable),fixed = T)
  colnames(ResTable) = gsub(pattern = "Spearman.",replacement = "Spearman ",x = colnames(ResTable),fixed = T)
  colnames(ResTable) = gsub(pattern = "Gene.name",replacement = "Gene name",x = colnames(ResTable),fixed = T)

  ResTable.gene.part = ResTable[(length(all.preds.nr)+1):dim(ResTable)[1],]
  ResTable.gene.part = ResTable.gene.part[order(-ResTable.gene.part[,"Score"], ResTable.gene.part[,"PValue"]), ]
  ResTable.sorted.by.Score = rbind(ResTable[1:length(all.preds.nr),],ResTable.gene.part)
  Analysis.Data$ResTable.pred.part = ResTable[1:length(all.preds.nr),]
  #colnames(ResTable.sorted.by.Score)
  
  output.file.name = Verify.path(sprintf("%s_%s_combined.txt",Analysis.Data$Analysis.name,Pars$DE.package))
  write.table.mod(ResTable.sorted.by.Score, file = output.file.name, sep='\t',na = '', col.names = colnames(ResTable.sorted.by.Score))
  if (Pars$Create.Excel.results){
    CL = sprintf('%s "%s/DE.results.to.Excel.py"', Pars$python.bin, Pars$suppl.data.dir)
    CL = gsub('/','\\',CL,fixed = T)
    CL = paste(c(CL,Verify.path(sprintf('"%s, DE results.xlsx"',Analysis.Data$Analysis.name)),sprintf('"%s"',c(output.file.name))),collapse = ' ')
    system(CL)
    if (!file.exists(Verify.path(sprintf('%s, DE results.xlsx',Analysis.Data$Analysis.name)))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (DE results) for ~%s FAILED',Analysis.Data$Analysis.name)) }
    
  }
  
  ResTable.gene.part = ResTable.gene.part[!is.na(ResTable.gene.part[,"Score"]),]
  sorted.norm.counts = norm.counts[rownames(ResTable.gene.part),]
  write.table.mod(sorted.norm.counts,'sorted.norm.counts.tsv',sep='\t')
  
  #cor.test(counts['FBgn0004242',],counts['FBgn0013972',],method = 'spearman')

  Analysis.Data$ResTable.gene.part = ResTable.gene.part
  saveRDS(ResTable.gene.part, file='DE.results.db')
  Analysis.Data$sorted.norm.counts = sorted.norm.counts
  Analysis.Data$Main.Component = Main.Component
  Analysis.Data$schema = schema
  Analysis.Data$Current.Predictors.list = Current.Predictors.list
  
  suppressWarnings(rm(counts))
  suppressWarnings(rm(noint))
  suppressWarnings(rm(cpms))
  suppressWarnings(rm(keep))
  suppressWarnings(rm(d))
  suppressWarnings(rm(ResTable))
  rm(ResTable.gene.part)
  rm(sorted.norm.counts)
  rm(norm.counts)
  
  Analysis.Data$GLM.analysis.performed = TRUE
  saveRDS(Analysis.Data, file = sprintf('%s, analysis data.db',Analysis.Data$Analysis.name))
  cat('\nDifferential expression/GLM analysis completed.\n')
  return(invisible(Analysis.Data))
}


# Create.Heatmaps.OLD = function(Startup.Data,Analysis.Data,bypass.if.completed = F,forced.parameters = NULL){
#   if (!Analysis.Data$GLM.analysis.performed) stop('Create.Heatmaps can be processed only after GLM/DE analysis. It is not performed.')
#   if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
#   } else Pars = forced.parameters
#   
#   if (bypass.if.completed & Read.Completed.steps.status(sprintf("%s.%s.%s.heatmaps",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
#     message(sprintf('\nCreating heatmaps for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
#     return()
#   }
#   
#   #############################################################
#   ########### Creating Heatmaps
#   
#   setwd(Analysis.Data$current.GLM.results.dir)
#   
#   cat(sprintf("\n%s:   creating heatmaps...",Analysis.Data$GLM.model))
#   #sorted.norm.counts = ResTable.gene.part[,14:dim(ResTable.gene.part)[2]]
#   rn.refseq = Analysis.Data$ResTable.gene.part[rownames(Analysis.Data$sorted.norm.counts),"Gene name"]
#   rn.ensembl = rownames(Analysis.Data$sorted.norm.counts)
#   no.refseq.name = is.na(rn.refseq)
#   
#   for (i in 1:length(no.refseq.name)){
#     if (no.refseq.name[i]){
#       rn.refseq[i] = rn.ensembl[i]
#     }
#   }
#   
#   dup = duplicated(rn.refseq)
#   
#   for (i in 1:length(dup)){
#     if (dup[i]){
#       rn.refseq[i] = sprintf('%s, %s',rn.refseq[i],rn.ensembl[i])
#     }
#   }
#   
#   sorted.norm.counts.for.heatmaps = Analysis.Data$sorted.norm.counts
#   rownames(sorted.norm.counts.for.heatmaps) = rn.refseq
#   
#   
#   
#   #x=50
#   for (x in Pars$Top.genes.to.include.in.heatmaps.list){
#     x = min(x,dim(sorted.norm.counts.for.heatmaps)[1])
#     if (x == dim(sorted.norm.counts.for.heatmaps)[1]){
#       Heatmap_Title = sprintf("%s, all genes. CPM**0.5 ",Analysis.Data$GLM.model)
#       png.mod(filename = sprintf('%s.heatmap.all.genes.png',Analysis.Data$Analysis.name,x),
#           units="in", width=7, height=7, pointsize=12, res=300)
#     } else {
#       Heatmap_Title = sprintf("%s, top %d DE genes. CPM**0.5 ",Analysis.Data$GLM.model,x)
#       png.mod(filename = sprintf('%s.heatmap.%d.top.genes.png',Analysis.Data$Analysis.name,x),
#           units="in", width=7, height=7, pointsize=12, res=300)
#     }
#     
#     
#     hmcols<-colorRampPalette(c("white","yellow","orange","red","blue"))(256)
#     limit = quantile(sorted.norm.counts.for.heatmaps[1:x,],0.983)
#     par(cex.main=0.6)
#     max.chars = max(nchar(colnames(sorted.norm.counts.for.heatmaps)))
#     heatmap.2(main = Heatmap_Title, data.matrix(pmin(sorted.norm.counts.for.heatmaps[1:x,],limit))**0.5,Colv = FALSE,Rowv = T,scale = 'none',
#               col = hmcols, key=TRUE, symkey=FALSE,
#               density.info="none", trace="none",
#               cexRow=0.55/x*50,cexCol=min(1, 0.8/max.chars*10),dendrogram='row')
#     dev.off()
#     
#     
#     if (x == dim(sorted.norm.counts.for.heatmaps)[1]){
#       Heatmap_Title = sprintf("%s, all genes. CPM z-score ",Analysis.Data$GLM.model)
#       png.mod(filename = sprintf('%s.heatmap.norm.all.genes.png',Analysis.Data$Analysis.name,x),
#           units="in", width=7, height=7, pointsize=12, res=300)
#     } else {
#       Heatmap_Title = sprintf("%s, top %d DE genes. CPM z-score",Analysis.Data$GLM.model,x)
#       png.mod(filename = sprintf('%s.heatmap.norm.%d.top.genes.png',Analysis.Data$Analysis.name,x),
#           units="in", width=7, height=7, pointsize=12, res=300)
#     }
#     
#     
#     intra.norm = sorted.norm.counts.for.heatmaps/rowSums(sorted.norm.counts.for.heatmaps)
#     par(cex.main=0.6)
#     hmcols<-colorRampPalette(c("#7372b2","#78abf8",'white',"#ffb31f","#ff6521"))(512)
#     heatmap.2(main = Heatmap_Title, intra.norm[1:x,]**1,Colv = FALSE,Rowv = T,scale = 'row',
#               col = hmcols, key=TRUE, symkey=FALSE,
#               density.info="none", trace="none",
#               cexRow=0.55/x*50,cexCol=min(1, 0.8/max.chars*10),dendrogram='row')
#     dev.off()
#   }
#   Add.Completed.step(sprintf("%s.%s.%s.heatmaps",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
#   Analysis.Data$Heatmaps.created = T
#   cat('\nCreating heatmaps completed.\n')
# }

Create.Heatmaps = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL,
                           Analysis.Data.RDS.file = NULL, bypass.if.completed = FALSE, forced.parameters = NULL,
                           limit.to.genes = NULL, out.dir = NULL, max.Genes.counts.list = NULL, cluster.samples = FALSE,
                           sort.by.predictor = NULL,
                           annotation.logFC.range = 4, annotation.log10P.range = 12){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  if(is.null(GLM.model))  GLM.model = Analysis.Data$GLM.model
  
  if (!Analysis.Data$GLM.analysis.performed) stop('Create.Heatmaps can be processed only after GLM/DE analysis. It is not performed.')

  if (bypass.if.completed & Read.Completed.steps.status(sprintf("%s.%s.%s.heatmaps",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
    message(sprintf('\nCreating heatmaps for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
    return()
  }

  
  if (cluster.samples && !is.null(sort.by.predictor)){
    msg = '\nCreate.Heatmaps(...):   Both sort.by.predictor and cluster.samples parameters are provided. Clustering samples will not be applied\n'
    cat(msg)
    warning(msg)
    cluster.samples = FALSE
  }
  
    
  #############################################################
  ########### Creating Heatmaps
  
  if(is.null(out.dir))  out.dir = Analysis.Data$current.GLM.results.dir
  
  cat(sprintf("\n~ %s:   creating heatmaps...\n", GLM.model))
  #sorted.norm.counts = ResTable.gene.part[,14:dim(ResTable.gene.part)[2]]
  
  if(!is.null(limit.to.genes)){
    keep = (rownames(Analysis.Data$sorted.norm.counts) %in% limit.to.genes)
  } else {
    keep = rep(TRUE, ncol(Analysis.Data$sorted.norm.counts))
  }
  
  if(sum(keep) == 0){
    cat('\nNo genes to render\n')
    return(invisible())
  }
  
  rn.refseq = Analysis.Data$ResTable.gene.part[rownames(Analysis.Data$sorted.norm.counts)[keep],"Gene name"]
  rn.ensembl = rownames(Analysis.Data$sorted.norm.counts)[keep]
  no.refseq.name = is.na(rn.refseq)
  
  for (i in 1:length(no.refseq.name)){
    if (no.refseq.name[i]){
      rn.refseq[i] = rn.ensembl[i]
    }
  }
  
  dup = duplicated(rn.refseq)
  
  for (i in 1:length(dup)){
    if (dup[i]){
      rn.refseq[i] = sprintf('%s, %s',rn.refseq[i],rn.ensembl[i])
    }
  }
  
  if(sum(keep) < 2){
    cat('\nHeatmap will not be created since <2 genes are given\n')
    return(invisible())
  }
  sorted.norm.counts.for.heatmaps = Analysis.Data$sorted.norm.counts[keep,]
  rownames(sorted.norm.counts.for.heatmaps) = rn.refseq
  
  corresponding.Ensembl.IDs = rownames(Analysis.Data$sorted.norm.counts[keep,])
  genes.annotation = as.data.frame(Analysis.Data$ResTable.gene.part[corresponding.Ensembl.IDs, c('logFC','logCPM')])
  
  if(is.null(annotation.logFC.range) | 'auto' %in% annotation.logFC.range){
    annotation.logFC.range = (mean(abs(genes.annotation[1:50,'logFC'])) + quantile(abs(genes.annotation[,'logFC']), 0.99)) / 2
    annotation.logFC.range %<>% min(.,6) %>% max(., 1.5) %>% round(.,1)
  }
  if(is.null(annotation.log10P.range) | 'auto' %in% annotation.log10P.range){
    annotation.log10P.range = 15
    #annotation.log10P.range = (mean(abs(genes.annotation[1:50,'logFC'])) + quantile(abs(genes.annotation[,'logFC']), 0.99)) / 2
    #annotation.log10P.range %<>% min(.,6) %>% max(., 1.5) %>% round(.,1)
  }
  
  genes.annotation[,'logFC'] %<>% pmin(. , annotation.logFC.range) %>% pmax (., -annotation.logFC.range)
  genes.annotation[,'logCPM'] %<>% pmin(. , 10) %>% pmax (., 0)
  genes.annotation = cbind(genes.annotation, (-1)*log10(Analysis.Data$ResTable.gene.part[corresponding.Ensembl.IDs, 'PValue'] + 1e-100))
  colnames(genes.annotation)[3] <- '-log10(p)'
  genes.annotation[,'-log10(p)'] %<>% pmin(. , annotation.log10P.range) %>% pmax (., 0)

  #head(Analysis.Data$ResTable.gene.part$'Spearman r')
  Spearman_r.is.present = FALSE
  if('Spearman r' %in% colnames(Analysis.Data$ResTable.gene.part)){
    if(any(!is.na(Analysis.Data$ResTable.gene.part[corresponding.Ensembl.IDs, 'Spearman r']))){
      genes.annotation = cbind(genes.annotation, Analysis.Data$ResTable.gene.part[corresponding.Ensembl.IDs, 'Spearman r'])
      colnames(genes.annotation)[4] <- 'Spearman_r'
      Spearman_r.is.present = TRUE
    }
  }
  
  rownames(genes.annotation) = rownames(sorted.norm.counts.for.heatmaps)

  #sample.order = 1:ncol(sorted.norm.counts.for.heatmaps)
  
  if(!is.null(sort.by.predictor)){
    sample.annotations = as.data.frame(Analysis.Data$predictor.values)
    if(sort.by.predictor %in% colnames(sample.annotations)){
      sample.order = order(sample.annotations[,sort.by.predictor])
      sorted.norm.counts.for.heatmaps = sorted.norm.counts.for.heatmaps[,sample.order]
    } else if (!msg.sent){
      msg = sprintf('\nPredictor %s is not found among sample annotations. Sorting will not be applied\n', sort.by.predictor)
      cat(msg)
      warning(msg)
    }
  }
  
  if(is.null(max.Genes.counts.list))  max.Genes.counts.list = Pars$Top.genes.to.include.in.heatmaps.list
  max.Genes.counts.list = sapply(max.Genes.counts.list, function(x) { min(x, nrow(sorted.norm.counts.for.heatmaps)) })
  max.Genes.counts.list = max.Genes.counts.list[!duplicated(max.Genes.counts.list)]
  
  
  #cluster.samples = TRUE
  x=max.Genes.counts.list[length(max.Genes.counts.list)]
  x=120
  for (x in max.Genes.counts.list){
    cat(sprintf("\r  processing heatmap for top %d DE genes...",x))
    x = min(x,dim(sorted.norm.counts.for.heatmaps)[1])
    
    #genes.annotation = rownames(sorted.norm.counts.for.heatmaps[1:x,])
    
    hmcols<-colorRampPalette(c("white","yellow","orange","red","blue"))(256)
    limit = quantile(sorted.norm.counts.for.heatmaps[1:x,],0.983)
    par(cex.main=0.4)
    max.chars = max(nchar(colnames(sorted.norm.counts.for.heatmaps)))
    
    sample.annotations = as.data.frame(Analysis.Data$predictor.values.nr)
    #sample.annotations = sample.annotations[,!endsWith(colnames(sample.annotations),' D')]
    
    keep = !(apply(sample.annotations,2,function(x) { length(levels(as.factor(x)))}) %in% 1)
    sample.annotations.clf = sample.annotations[, keep]
    
    if(sum(keep) == 1){
      sample.annotations.clf = as.data.frame(sample.annotations.clf)
      colnames(sample.annotations.clf) = colnames(sample.annotations)[which(keep)]
      rownames(sample.annotations.clf) = rownames(sample.annotations)
    }
    
    
    current.genes.annotation = genes.annotation[1:x,]
    logFC.colors <- colorRampPalette(c("#4875e6","white","#f95136"))(512)
    Spearman_r.colors <- colorRampPalette(c("#0e8bec","white","#f13004"))(512)
    P.colors <- colorRampPalette(c("white","#f1ab04","#69cd0f"))(512)
    logCPM.colors <- colorRampPalette(c("white", "#ecab0e"))(512)
    
    logFC.color.coord.start <- (length(logFC.colors) * (min(current.genes.annotation$logFC) + annotation.logFC.range) / 2 / annotation.logFC.range) %>% round %>% min(., length(logFC.colors) - 1) %>% max(., 0)
    logFC.color.coord.end <- (length(logFC.colors) * (max(current.genes.annotation$logFC) + annotation.logFC.range) / 2 / annotation.logFC.range) %>% round %>% min(., length(logFC.colors)) %>% max(., 0)
    P.color.coord.start <- (length(P.colors) * (min(current.genes.annotation$`-log10(p)`) + 0) / annotation.log10P.range) %>% round %>% min(., length(P.colors) - 1) %>% max(., 0)
    P.color.coord.end <- (length(P.colors) * (max(current.genes.annotation$`-log10(p)`) + 0) / annotation.log10P.range) %>% round %>% min(., length(P.colors)) %>% max(., 0)
    logCPM.color.coord.start <- (length(logCPM.colors) * (min(current.genes.annotation$logCPM) + 0) / 10) %>% round %>% min(., length(logCPM.colors) - 1) %>% max(., 0)
    logCPM.color.coord.end <- (length(logCPM.colors) * (max(current.genes.annotation$logCPM) + 0) / 10) %>% round %>% min(., length(logCPM.colors)) %>% max(., 0)
    if(Spearman_r.is.present){
      Speaman_r.color.coord.start <- (length(Spearman_r.colors) * (min(current.genes.annotation$Spearman_r) + 1) / 2) %>% round %>% min(., length(Spearman_r.colors) - 1) %>% max(., 0)
      Speaman_r.color.coord.end <- (length(Spearman_r.colors) * (max(current.genes.annotation$Spearman_r) + 1) / 2) %>% round %>% min(., length(Spearman_r.colors)) %>% max(., 0)
    }
    
    ann.colors = list(
      logFC = logFC.colors [logFC.color.coord.start : logFC.color.coord.end],
      logCPM = logCPM.colors[logCPM.color.coord.start : logCPM.color.coord.end]
    )
    ann.colors[['-log10(p)']] = P.colors[P.color.coord.start : P.color.coord.end]
    if(Spearman_r.is.present){
      ann.colors[['Spearman_r']] = Spearman_r.colors[Speaman_r.color.coord.start : Speaman_r.color.coord.end]
    }
    
    width.min = 8
    width.max = 19
    hm.width = max(width.min, min(width.max, dim(sorted.norm.counts.for.heatmaps)[2]*17/29))
    hm.height = min(16, max(8, 8*(ncol(sample.annotations) + ncol(genes.annotation))/9))
    
    if (x == dim(sorted.norm.counts.for.heatmaps)[1]){
      Heatmap_Title = sprintf("~ %s, all %d genes. CPM**0.5 ", GLM.model, x)
      png.mod(filename = sprintf('%s/%s.heatmap.all.%d.genes.png', out.dir, Analysis.name, x),
              units="in", width=hm.width, height=hm.height, pointsize=12, res=300)
    } else {
      Heatmap_Title = sprintf("~ %s, top %d DE genes. CPM**0.5 ", GLM.model, x)
      png.mod(filename = sprintf('%s/%s.heatmap.%d.top.genes.png', out.dir, Analysis.name, x),
              units="in", width=hm.width, height=hm.height, pointsize=12, res=300)
    }
    
    pheatmap(mat = data.matrix(pmin(sorted.norm.counts.for.heatmaps[1:x, ],limit))**0.5, scale = 'none', color = hmcols,
             cluster_cols = cluster.samples, border_color = if (x <= 70) "grey60" else NA,
             main = Heatmap_Title, fontsize = 7, fontsize_row = min(9, 9/x*50) * hm.height / 8,
             fontsize_col = min(22*min(0.5, 0.8/max.chars*10), 23*min(1, 30/dim(sorted.norm.counts.for.heatmaps)[2])),
             annotation_col = sample.annotations.clf, annotation_row = current.genes.annotation,
             annotation_colors = ann.colors)
    
    dev.off()
    # heatmap.2(main = Heatmap_Title, data.matrix(pmin(sorted.norm.counts.for.heatmaps[1:x,],limit))**0.5,Colv = FALSE,Rowv = T,scale = 'none',
    #           col = hmcols, key=TRUE, symkey=FALSE,
    #           density.info="none", trace="none",
    #           cexRow=0.55/x*50, cexCol=min(1, 0.8/max.chars*10),dendrogram='row')
    
    
    
    if (x == dim(sorted.norm.counts.for.heatmaps)[1]){
      Heatmap_Title = sprintf("~ %s, all %d genes. CPM z-score ", GLM.model, x)
      png.mod(filename = sprintf('%s/%s.heatmap.norm.all.%d.genes.png', out.dir, Analysis.name, x),
          units="in", width=hm.width, height = hm.height, pointsize=12, res=300)
    } else {
      Heatmap_Title = sprintf("~ %s, top %d DE genes. CPM z-score", GLM.model, x)
      png.mod(filename = sprintf('%s/%s.heatmap.norm.%d.top.genes.png', out.dir, Analysis.name, x),
          units="in", width=hm.width, height = hm.height, pointsize=12, res=300)
    }
    
    
    intra.norm = sorted.norm.counts.for.heatmaps[1:x,]
    intra.norm = intra.norm/rowSums(intra.norm)*ncol(intra.norm)
    if(ncol(intra.norm) > 70){
      intra.norm = t(apply(intra.norm,1,function(y){ 
        #print(quantile(y, 0.99))
        lim = quantile(y, 0.97)
        if(lim > 0){
          return(pmin(y, lim))
        } else {
          return(y)
        }
        }))
    }
    #dim(intra.norm2)
    
    par(cex.main=0.4)
    hmcols<-colorRampPalette(c("#7372b2","#78abf8",'white',"#ffb31f","#ff6521"))(512)
    pheatmap(mat = intra.norm, scale = 'row', color = hmcols, cluster_cols = cluster.samples,
             border_color = if (x <= 70) "grey60" else NA,
             main = Heatmap_Title, fontsize = 7, fontsize_row = min(9, 9/x*50) * hm.height / 8,
             fontsize_col = min(22*min(0.5, 0.8/max.chars*10), 23*min(1, 30/dim(sorted.norm.counts.for.heatmaps)[2])),
             annotation_col = sample.annotations.clf, annotation_row = current.genes.annotation,
             annotation_colors = ann.colors)
    
    # heatmap.2(main = Heatmap_Title, intra.norm[1:x,]**1,Colv = FALSE,Rowv = T,scale = 'row',
    #           col = hmcols, key=TRUE, symkey=FALSE,
    #           density.info="none", trace="none",
    #           cexRow=0.55/x*50,cexCol=min(1, 0.8/max.chars*10),dendrogram='row')
    dev.off()
  }
  Add.Completed.step(sprintf("%s.%s.%s.heatmaps",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
  Analysis.Data$Heatmaps.created = TRUE
  cat('\nCreating heatmaps completed.\n')
}

Get.LogFCs.for.Ensembl.IDs.vector <- function(Analysis.Data, genes.vector, PValue.thr = 1, logCPM.thr = 0, sort.logFCs = TRUE){
  return(
    sapply(genes.vector, function(gene.ID.array){
      gene.ID.array = gene.ID.array[!duplicated(gene.ID.array)]
      gene.ID.array = gene.ID.array[gene.ID.array %in% rownames(Analysis.Data$ResTable.gene.part)]
      sub.table = Analysis.Data$ResTable.gene.part[gene.ID.array,]
      sub.table = sub.table[sub.table$logCPM >= logCPM.thr,]
      sub.table$logFC[which(sub.table$PValue > PValue.thr)] = 0
      
      if(sort.logFCs)  return(sub.table$logFC[rev(order(sub.table$logFC))])
      return(sub.table$logFC)
    })
  )  
}


Get.complete.LogFCs.with.Entrez.IDs <- function(Startup.Data, Pars, Analysis.Data){
  #if(Pars$Species == 'dme')  gene_ID_type = 'flybasecgid_gene' else  gene_ID_type = 'entrezgene'
  gene_ID_type = 'entrezgene'
  EntrezGene.IDs = Startup.Data$General.maRt.table[!(Startup.Data$General.maRt.table[,gene_ID_type] %in% NA),]
  rownames(EntrezGene.IDs) = EntrezGene.IDs[,"ensembl_gene_id"]
  Entrez.Stats.table = Analysis.Data$ResTable.gene.part[,c('logFC','PValue','Score')]
  Entrez.Stats.table = Entrez.Stats.table[(rownames(Entrez.Stats.table) %in% rownames(EntrezGene.IDs)),]
  Entrez.Stats.table = cbind(Entrez.Stats.table,data.frame(EntrezGene.IDs[rownames(Entrez.Stats.table),c(gene_ID_type)]))
  colnames(Entrez.Stats.table)[4] = gene_ID_type
  
  my_geneList = Entrez.Stats.table$logFC
  names(my_geneList) = Entrez.Stats.table[,gene_ID_type]
  my_geneList = my_geneList[rev(order(my_geneList))]
  return(my_geneList)
}

Do.end.run = FALSE
Get.Entrez.IDs <- function(Startup.Data, Pars, Analysis.Data,
                           max.DE.genes = NULL, DE.type = 'both', min.abs.LogFC = 0, max.PValue = 1,
                           min.Score = 0, disable.filters = FALSE, use.Org.Db__instead.of__biomaRt.table = FALSE){

  if(use.Org.Db__instead.of__biomaRt.table && Pars$Species == 'dme'){
    cat('\nUsing org.Dm.eg.db is not compartible with enrichKEGG / enrichPathway / enrichGO functions. biomaRt tables will be used instead of org.Dm.eg.db\n')
    use.Org.Db__instead.of__biomaRt.table = FALSE
  }
  if(use.Org.Db__instead.of__biomaRt.table){
    if (Pars$Species %in% names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
      org.DB = Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes[Pars$Species]
    } else {
      msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s)\n',
                    Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                    paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                    paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
      cat(msg)
      warning(msg)
      return(NULL)
    }
    # if (Pars$Species == 'hsa')  org.DB = org.Hs.eg.db
    # if (Pars$Species == 'mmu')  org.DB = org.Mm.eg.db
    # if (Pars$Species == 'dme')  org.DB = org.Dm.eg.db
    
    #gn = rownames(Analysis.Data$ResTable.gene.part)[1:40]
    #sum.mod(Startup.Data$General.maRt.table[gn,gene_ID_type] %in% NA)
    #sum.mod(mapIds(org.DB, keys=gn, column="ENTREZID", keytype="ENSEMBL", multiVals="first") %in% NA)
    
    ## use the following to obtain all available keys in org.Hs.eg.db: keytypes(org.Hs.eg.db)
    ##
    
    #if(Pars$Species == 'dme')  gene_ID_type = 'flybasecgid_gene' else  gene_ID_type = 'entrezgene'
    
    Entrez.Stats.table = Analysis.Data$ResTable.gene.part[,c('logFC','PValue','Score')]
    Entrez.IDs = mapIds(org.DB, keys=rownames(Entrez.Stats.table), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
    Entrez.IDs = Entrez.IDs[!(Entrez.IDs %in% NA)]
    #cat(sprintf('\n%d of %d genes unmapped\n',sum.mod(Entrez.IDs %in% NA),length(Entrez.IDs)))
    Entrez.Stats.table = Entrez.Stats.table[(rownames(Entrez.Stats.table) %in% names(Entrez.IDs)),]
    Entrez.Stats.table = cbind(Entrez.Stats.table,data.frame(Entrez.IDs[rownames(Entrez.Stats.table)]))
    head(Entrez.Stats.table)

    gene_ID_type = 'entrezgene'
    colnames(Entrez.Stats.table)[4] = gene_ID_type
    
  } else {
    #if(Pars$Species == 'dme')  gene_ID_type = 'flybasecgid_gene' else  gene_ID_type = 'entrezgene'
    gene_ID_type = 'entrezgene'
    EntrezGene.IDs = Startup.Data$General.maRt.table[!(Startup.Data$General.maRt.table[,gene_ID_type] %in% NA),]
    rownames(EntrezGene.IDs) = EntrezGene.IDs[,"ensembl_gene_id"]
    Entrez.Stats.table = Analysis.Data$ResTable.gene.part[,c('logFC','PValue','Score')]
    Entrez.Stats.table = Entrez.Stats.table[(rownames(Entrez.Stats.table) %in% rownames(EntrezGene.IDs)),]
    Entrez.Stats.table = cbind(Entrez.Stats.table,data.frame(EntrezGene.IDs[rownames(Entrez.Stats.table),c(gene_ID_type)]))
    #head(data.frame(EntrezGene.IDs[rownames(Entrez.Stats.table),c("entrezgene")]))
    
    colnames(Entrez.Stats.table)[4] = gene_ID_type
  }
  

  
  if(!disable.filters){
    if (DE.type == 'de' | DE.type == 'up-down' | DE.type == 'upreg-downreg' | DE.type == 'both'){
      Passed = (abs(Entrez.Stats.table[,"logFC"]) > min.abs.LogFC) & (Entrez.Stats.table[,"PValue"] < max.PValue) & (Entrez.Stats.table[,"Score"] > min.Score)
    }
    if (DE.type == 'upreg' | DE.type == 'up'){
      Passed = (Entrez.Stats.table[,"logFC"] > 0) & (abs(Entrez.Stats.table[,"logFC"]) > min.abs.LogFC) & (Entrez.Stats.table[,"PValue"] < max.PValue) & (Entrez.Stats.table[,"Score"] > min.Score)
    }
    if (DE.type == 'downreg' | DE.type == 'down'){
      Passed = (Entrez.Stats.table[,"logFC"] < 0) & (abs(Entrez.Stats.table[,"logFC"]) > min.abs.LogFC) & (Entrez.Stats.table[,"PValue"] < max.PValue) & (Entrez.Stats.table[,"Score"] > min.Score)
    }
    if (max.DE.genes >= length(Passed[Passed])){
      Do.end.run <<- TRUE
    }
    
    max.DE.genes = min(length(Passed[Passed])-1,max.DE.genes)
    Entrez.Stats.table = Entrez.Stats.table[Passed,]
    min.Score.threshold.effective = max(min.Score, sort(Entrez.Stats.table[,"Score"],decreasing = T)[max.DE.genes+1])
    Entrez.Stats.table = Entrez.Stats.table[Entrez.Stats.table[,"Score"] > min.Score.threshold.effective,]
  }
  Entrez.IDs = levels(as.factor(Startup.Data$General.maRt.table[rownames(Entrez.Stats.table),gene_ID_type]))
  ## if(Pars$Species == 'dme')  Entrez.IDs = sprintf('Dmel_CG%s',Entrez.IDs)
  return(Entrez.IDs)
}


Get.Ensembl.IDs <- function(Startup.Data, Pars, Analysis.Data,
                           max.DE.genes = NULL, DE.type = 'both', min.abs.LogFC = 0, max.PValue = 1,
                           min.Score = 0, disable.filters = FALSE){
  
  Ensembl.Stats.table = Analysis.Data$ResTable.gene.part[,c('logFC','PValue','Score')]
  rownames(Ensembl.Stats.table) = rownames(Analysis.Data$ResTable.gene.part)
  if(!disable.filters){
    if (DE.type == 'de' | DE.type == 'up-down' | DE.type == 'upreg-downreg' | DE.type == 'both'){
      Passed = (abs(Ensembl.Stats.table[,"logFC"]) > min.abs.LogFC) & (Ensembl.Stats.table[,"PValue"] < max.PValue) & (Ensembl.Stats.table[,"Score"] > min.Score)
    }
    if (DE.type == 'upreg' | DE.type == 'up'){
      Passed = (Ensembl.Stats.table[,"logFC"] > 0) & (abs(Ensembl.Stats.table[,"logFC"]) > min.abs.LogFC) & (Ensembl.Stats.table[,"PValue"] < max.PValue) & (Ensembl.Stats.table[,"Score"] > min.Score)
    }
    if (DE.type == 'downreg' | DE.type == 'down'){
      Passed = (Ensembl.Stats.table[,"logFC"] < 0) & (abs(Ensembl.Stats.table[,"logFC"]) > min.abs.LogFC) & (Ensembl.Stats.table[,"PValue"] < max.PValue) & (Ensembl.Stats.table[,"Score"] > min.Score)
    }
    if (max.DE.genes >= length(Passed[Passed])){
      Do.end.run <<- TRUE
    }
    
    max.DE.genes = min(length(Passed[Passed])-1,max.DE.genes)
    Ensembl.Stats.table = Ensembl.Stats.table[Passed,]
    min.Score.threshold.effective = max(min.Score, sort(Ensembl.Stats.table[,"Score"],decreasing = T)[max.DE.genes+1])
    Ensembl.Stats.table = Ensembl.Stats.table[Ensembl.Stats.table[,"Score"] > min.Score.threshold.effective,]
  }
  
  return(rownames(Ensembl.Stats.table))
}


Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table = function(GSEA.summary.table, Startup.Data, Pars, gene.col, Ensembl.mode = FALSE){
  if(nrow(GSEA.summary.table)==0) return(GSEA.summary.table)
  if(!Ensembl.mode){
    #if(Pars$Species == 'dme')  gene_ID_type = 'flybasecgid_gene' else  gene_ID_type = 'entrezgene'
    gene_ID_type = 'entrezgene'
  } else {
    gene_ID_type = 'ensembl_gene_id'
  }
  
  tmp.maRt.table = Startup.Data$General.maRt.table[!(Startup.Data$General.maRt.table[,gene_ID_type] %in% NA),]
  tmp.maRt.table = tmp.maRt.table[!duplicated(tmp.maRt.table[,gene_ID_type]),]
  rownames(tmp.maRt.table) = sapply(tmp.maRt.table[,gene_ID_type],function(y) { as.character(y) })
  
  for(x in 1:dim(GSEA.summary.table)[1]){
    gene.IDs = strsplit(GSEA.summary.table[x,gene.col],'/',fixed = TRUE)[[1]]
    gene.IDs = gene.IDs[gene.IDs %in% rownames(tmp.maRt.table)]
    gene.Names = tmp.maRt.table[gene.IDs,'external_gene_name']
    GSEA.summary.table[x,gene.col] = paste(gene.Names, collapse = '/')
  }
  return(GSEA.summary.table)
}

Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table__ENSEMBL = function(GSEA.summary.table, Startup.Data, Pars, gene.col){
  if(nrow(GSEA.summary.table)==0) return(GSEA.summary.table)
  Ensembl.names__to__gene.symbols = Startup.Data$Ensembl.names__to__gene.symbols

  for(x in 1:dim(GSEA.summary.table)[1]){
    gene.IDs = strsplit(GSEA.summary.table[x,gene.col],'/',fixed = TRUE)[[1]]
    gene.IDs = gene.IDs[gene.IDs %in% names(Ensembl.names__to__gene.symbols)]
    gene.Names = Ensembl.names__to__gene.symbols[gene.IDs]
    GSEA.summary.table[x,gene.col] = paste(gene.Names, collapse = '/')
  }
  return(GSEA.summary.table)
}

Get.Gene.Ensembl.IDs__for__GO.term.IDs = function(term.list, Startup.Data, forced.parameters = NULL, remove.notfound.terms = TRUE){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if (Pars$Species %in% names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
    org.DB = Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes[Pars$Species]
  } else {
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s). Gene ID conversion will not be performed\n',
                  Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                  paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                  paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
    cat(msg)
    warning(msg)
    org.DB = NULL
    return(vector(mode = 'list'))
  }
  
  # term.list = c('GO:0042416','GO:0051971','GO:0022008')
  # use Startup.Data$GO.full.descriptions to get to know GO terms IDs

  eval(parse(text = sprintf("GO.ID__to__Entrez = as.list(%s::%sGO2ALLEGS)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  if(remove.notfound.terms)  term.list = term.list[term.list %in% names(GO.ID__to__Entrez)]
  eval(parse(text = sprintf("Entrez__to__Ensembl = as.list(%s::%sENSEMBL)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  #gene.list = sapply(Reactome.ID__to__Entrez[term.list],function(x) { unlist(mapIds(org.Mm.eg.db, keys = x, column="ENSEMBL", keytype="ENTREZID", multiVals="list")) })
  return(sapply(GO.ID__to__Entrez[term.list], function(x) { unlist(Entrez__to__Ensembl[x]) }))
}



Get.Gene.Ensembl.IDs__for__Reactome.pathway.IDs = function(pathway.list, Startup.Data, forced.parameters = NULL, remove.notfound.pathways = TRUE){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if (Pars$Species %in% names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
    org.DB = Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes[Pars$Species]
  } else {
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s). Gene ID conversion will not be performed\n',
                  Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                  paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                  paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
    cat(msg)
    warning(msg)
    org.DB = NULL
    #return(pathway.list)
  }
  
  # pathway.list = c('R-MMU-1236974',
  #             'R-MMU-70171',
  #             'R-MMU-70263',
  #             'R-MMU-1222556',
  #             'R-MMU-5663213',
  #             'R-MMU-2029482',
  #             'R-MMU-1236978',
  #             'R-MMU-69229',
  #             'left')
  
  Reactome.ID__to__Entrez = as.list(reactomePATHID2EXTID)
  # Reactome.ID__to__Reactome.Name = as.list(reactomePATHID2NAME)
  if(remove.notfound.pathways)  pathway.list = pathway.list[pathway.list %in% names(Reactome.ID__to__Entrez)]
  
  if(!is.null(org.DB)){
    eval(parse(text = sprintf("Entrez__to__Ensembl = as.list(%s::%sENSEMBL)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
    #gene.list = sapply(Reactome.ID__to__Entrez[pathway.list],function(x) { unlist(mapIds(org.Mm.eg.db, keys = x, column="ENSEMBL", keytype="ENTREZID", multiVals="list")) })
    return(sapply(Reactome.ID__to__Entrez[pathway.list], function(x) { unlist(Entrez__to__Ensembl[x]) }))
  } else {
    return(Reactome.ID__to__Entrez[pathway.list])
  }
}

Get.Gene.Ensembl.IDs__for__KEGG.pathway.IDs___not.using.clusterProfiler = function(pathway.list, Startup.Data, forced.parameters = NULL, remove.notfound.pathways = TRUE){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if (Pars$Species %in% names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
    org.DB = Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes[Pars$Species]
  } else {
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s). KEGG pathway ID conversion to gene names will not be performed\n',
                  Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                  paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                  paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
    cat(msg)
    warning(msg)
    org.DB = NULL
    return(pathway.list)
  }
  
  # pathway.list = c('mmu04142',
  #   'mmu00511',
  #   'mmu04915',
  #   'mmu04960',
  #   'mmu04150',
  #   'mmu04213',
  #   'mmu05221',
  #   'mmu04211')

  #names(KEGG.ID__to__Entrez)
  
  pathway.list.num = gsub("[[:alpha:] ]", "", pathway.list)
  if(remove.notfound.pathways)  pathway.list.num = pathway.list.num[pathway.list.num %in% names(KEGG.ID__to__Entrez)]
  eval(parse(text = sprintf("KEGG.ID__to__Entrez = as.list(%s::%sPATH2EG)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  eval(parse(text = sprintf("Entrez__to__Ensembl = as.list(%s::%sENSEMBL)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  # KEGG.ID__to__KEGG.Name = as.list(....)
  
  #gene.list = sapply(Reactome.ID__to__Entrez[pathway.list],function(x) { unlist(mapIds(org.Mm.eg.db, keys = x, column="ENSEMBL", keytype="ENTREZID", multiVals="list")) })
  
  res = sapply(KEGG.ID__to__Entrez[pathway.list.num], function(x) { unlist(Entrez__to__Ensembl[x]) })
  names(res) = sprintf("%s%s",gsub("[^[:alpha:] ]", "", pathway.list[1]),names(res))
  return(res)
}


Get.Gene.Ensembl.IDs__for__KEGG.pathway.IDs = function(pathway.list, Startup.Data, forced.parameters = NULL, remove.notfound.pathways = TRUE){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters

  if (!(Pars$Species %in% rownames(Startup.Data$DB.data$KEGG.code.table))){
    msg = sprintf('\nSpecies "%s" is not found among available KEGG species. Ensure that you provided adequate KEGG code (e.g. hsa, mmu, dme, ...)\n', Pars$Species)
    cat(msg)
    warning(msg)
  }

  if (Pars$Species %in% names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
    org.DB = Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes[Pars$Species]
  } else {
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s). Gene ID conversion will not be performed\n',
                  Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                  paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                  paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
    cat(msg)
    warning(msg)
    org.DB = NULL
    #return(pathway.list)
  }
  
  # Startup.Data$KEGG_DATA = KEGG_DATA
  # pathway.list = c('mmu04142',
  #   'mmu00511',
  #   'mmu04915',
  #   'mmu04960',
  #   'mmu04150',
  #   'mmu04213',
  #   'mmu05221',
  #   'mmu04211',
  #   'left')
  

  # KEGG.ID__to__KEGG.Name = KEGG_DATA$PATHID2NAME
  if(remove.notfound.pathways)  pathway.list = pathway.list[pathway.list %in% names(Startup.Data$KEGG_DATA$PATHID2EXTID)]
  
  if(!is.null(org.DB)){
    eval(parse(text = sprintf("Entrez__to__Ensembl = as.list(%s::%sENSEMBL)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
    #gene.list = sapply(Reactome.ID__to__Entrez[pathway.list],function(x) { unlist(mapIds(org.Mm.eg.db, keys = x, column="ENSEMBL", keytype="ENTREZID", multiVals="list")) })
    return(sapply(Startup.Data$KEGG_DATA$PATHID2EXTID[pathway.list], function(x) { unlist(Entrez__to__Ensembl[x]) %>% .[!is.na(.)] }))
  } else {
    return(Startup.Data$KEGG_DATA$PATHID2EXTID[pathway.list])
  }
}



# Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table.Ensembl = function(GSEA.summary.table, Startup.Data, Pars, gene.col){
#   if(Pars$Species == 'dme')  gene_ID_type = 'flybasecgid_gene' else  gene_ID_type = 'entrezgene'
#   tmp.maRt.table = Startup.Data$General.maRt.table[!(Startup.Data$General.maRt.table[,gene_ID_type] %in% NA),]
#   tmp.maRt.table = tmp.maRt.table[!duplicated(tmp.maRt.table[,gene_ID_type]),]
#   rownames(tmp.maRt.table) = sapply(tmp.maRt.table[,gene_ID_type],function(y) { as.character(y) })
#   
#   for(x in 1:dim(GSEA.summary.table)[1]){
#     gene.IDs = strsplit(GSEA.summary.table[x,gene.col],'/',fixed = TRUE)[[1]]
#     gene.IDs = gene.IDs[gene.IDs %in% rownames(tmp.maRt.table)]
#     gene.Names = tmp.maRt.table[gene.IDs,'external_gene_name']
#     GSEA.summary.table[x,gene.col] = paste(gene.Names, collapse = '/')
#   }
#   return(GSEA.summary.table)
# }
# 
# KEGG.classic.Enrichment = function(Startup.Data, Analysis.Data, bypass.if.completed = FALSE,
#                                    forced.parameters = NULL, additional.plots = TRUE){
#   if (!Analysis.Data$GLM.analysis.performed) stop('KEGG.classic.Enrichment can be processed only after GLM/DE analysis. It is not performed.')
#   if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
#   } else Pars = forced.parameters
#   
#   if (bypass.if.completed & 
#       Read.Completed.steps.status(sprintf("%s.%s.%s.KEGG.classic.Enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file) &
#       Read.Completed.steps.status(sprintf("%s.%s.%s.pathways.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
#     message(sprintf('\n KEGG.classic.Enrichment pathway enrichment and visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
#     return(Analysis.Data)
#   }
#   
#   ######################################################
#   ### KEGG classic enrichment with clusterProfiler
#   
#   Analysis.Data$KEGG.classic.Enrich.working.dir = sprintf("%s/KEGG Pathways - Classic Enrichment Analysis",Analysis.Data$current.GLM.results.dir)
#   dir.create(Analysis.Data$KEGG.classic.Enrich.working.dir,showWarnings = F)
#   
#   KEGG.classic.Enrichment__results__by.DE.type__by.max.DE.genes = vector(mode = 'list')
#   Gene.Lists__by.DE.type.and.max.DE.genes = vector(mode = 'list')
# 
#   Enriched.Pathways.Classic = c()
#   setwd(Analysis.Data$KEGG.classic.Enrich.working.dir)
#   
#   KEGG.classic.Enrichment.output.file.names = c()
#   Pars$XX.classic.Enrich__with.clusterProfiler___max.DE.genes.list = sort(Pars$XX.classic.Enrich__with.clusterProfiler___max.DE.genes.list,decreasing = F)
#   XX.classic.Enrich__with.clusterProfiler___mode = Pars$XX.classic.Enrich__with.clusterProfiler___mode.list[1]  #for debug
#   for (XX.classic.Enrich__with.clusterProfiler___mode in Pars$XX.classic.Enrich__with.clusterProfiler___mode.list) {
#     KEGG.classic.Enrichment__results__by.DE.type__by.max.DE.genes[[XX.classic.Enrich__with.clusterProfiler___mode]] = vector(mode = 'list')
#     Do.end.run <<- FALSE
#     XX.classic.Enrich__with.clusterProfiler___max.DE.genes = Pars$XX.classic.Enrich__with.clusterProfiler___max.DE.genes.list[1]  #for debug
#     for (XX.classic.Enrich__with.clusterProfiler___max.DE.genes in Pars$XX.classic.Enrich__with.clusterProfiler___max.DE.genes.list) {
#       if (Do.end.run){
#         KEGG.classic.Enrichment__results__by.DE.type__by.max.DE.genes[[XX.classic.Enrich__with.clusterProfiler___mode]][[XX.classic.Enrich__with.clusterProfiler___max.DE.genes]] = kk
#         #Gene.Lists__by.DE.type.and.max.DE.genes[[sprintf('top%d %s',XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode)]] = Entrez.IDs
#         next
#       }
# 
#       Entrez.IDs = Get.Entrez.IDs(Startup.Data, Pars, Analysis.Data,
#                                   max.DE.genes = XX.classic.Enrich__with.clusterProfiler___max.DE.genes,  DE.type = XX.classic.Enrich__with.clusterProfiler___mode,
#                                   min.abs.LogFC = Pars$XX.classic.Enrich__with.clusterProfiler___gene.min.abs.LogFC.threshold,
#                                   max.PValue = Pars$XX.classic.Enrich__with.clusterProfiler___gene.max.PValue.threshold,
#                                   min.Score = Pars$XX.classic.Enrich__with.clusterProfiler___gene.min.Score.threshold)
# 
#       cat(sprintf('\n~ %s:   Performing KEGG classic enrichment for %d %s genes',Analysis.Data$GLM.model,length(Entrez.IDs),XX.classic.Enrich__with.clusterProfiler___mode))
#       
#       kk <- enrichKEGG(Entrez.IDs, organism=Pars$Species, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2)
#       ##kk2 = enrichPathway(ids, organism=Pars$Species, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=1)
#       ## minGSSize = ???
#       #kk <- gseKEGG(geneList, organism=Pars$Species, pvalueCutoff=0.05) #, pAdjustMethod="BH",qvalueCutoff=1)
# 
#       # kk = kk2
#       # kk2 <- gseKEGG(geneList     = my_geneList,
#       #                organism     = 'mmu',
#       #                nPerm        = 1000,
#       #                minGSSize    = 120,
#       #                pvalueCutoff = 0.05,
#       #                verbose      = TRUE)      
#       # 
#       # ego3 <- gseGO(geneList     = my_geneList,
#       #               OrgDb        = org.Mm.eg.db,
#       #               ont          = "BP",
#       #               nPerm        = 1000,
#       #               minGSSize    = 100,
#       #               maxGSSize    = 500,
#       #               pvalueCutoff = 0.05,
#       #               verbose      = TRUE)      
#       # 
#       # kk=ego3
#       
#       if (typeof(kk)!='S4') {
#         cat(sprintf('\nNo enriched KEGG pathways have been found for %s, top %d %s genes',Analysis.Data$GLM.model,XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode))
#         next
#       }
#       
#       if (dim(as.data.frame(kk))[1] == 0){
#         cat(sprintf('\nNo enriched KEGG pathways have been found for %s, top %d %s genes',Analysis.Data$GLM.model,XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode))
#         next
#       }
#       
#       Gene.Lists__by.DE.type.and.max.DE.genes[[sprintf('top%d %s',XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode)]] = Entrez.IDs
#       KEGG.classic.Enrichment__results__by.DE.type__by.max.DE.genes[[XX.classic.Enrich__with.clusterProfiler___mode]][[XX.classic.Enrich__with.clusterProfiler___max.DE.genes]] = kk
#       
#       cat(sprintf('\n%d enriched KEGG pathways have been found for %s, top %d %s genes',dim(as.data.frame(kk))[1],Analysis.Data$GLM.model,XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode))
#       
#       if(additional.plots){
#         png.mod(filename = sprintf("[cneplot] %s, KEGG classic enrich. for top-%d %s genes.png",Analysis.Data$Analysis.name,XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode),
#             units="in", width=7, height=7, pointsize=3, res=300)
#         cnetplot(kk,showCategory=20, categorySize="geneNum")
#         dev.off()
#         
#         #par(mar=c(1,1,1,1))
#         png.mod(filename = sprintf("[dotplot] %s, KEGG classic enrich. for top-%d %s genes.png",Analysis.Data$Analysis.name,XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode),
#             units="in", width=12, height=7, pointsize=8, res=300)
#         print({
#           dotplot(kk,showCategory=30) + ggplot2::xlab(sprintf("Gene ratio; KEGG, top-%d %s genes",XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode))
#         })
#         dev.off()
#         
#         ### optionally, one can use:
#         ## p = dotplot(kk,showCategory=30)
#         ## ggplot2::ggsave(file="plot.png", plot=p)
#         
#         #par(mar=c(1,1,1,1))
#         png.mod(filename = sprintf("[barplot] %s, KEGG classic enrich. for top-%d %s genes.png",Analysis.Data$Analysis.name,XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode),
#             units="in", width=12, height=7, pointsize=8, res=300)
#         print({
#           barplot(kk,showCategory=30) + ggplot2::ylab(sprintf("Gene count (per KEGG pathway) among top-%d %s ones",XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode))
#         })
#         dev.off()
#       }
#       
#       KEGG.classic.Enrich.summary.table = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table(as.data.frame(kk), Startup.Data, Pars, gene.col = 'geneID')
# 
#       file.name = sprintf("%s, KEGG classic enrichment for top-%d %s genes.tsv",Analysis.Data$Analysis.name,XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode)
#       KEGG.classic.Enrichment.output.file.names = c(KEGG.classic.Enrichment.output.file.names,file.name)
#       write.table.mod(x = KEGG.classic.Enrich.summary.table, file = file.name,sep='\t')
#       
#       Enriched.Pathways.Classic = c(Enriched.Pathways.Classic, rownames(as.data.frame(kk))) # This is a cumulative list of all found enriched pathways
#       
#     }
#   }
#   ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes, fun = "enrichKEGG", qvalueCutoff=1, organism=Pars$Species, pAdjustMethod="BH")
#   png.mod(filename = sprintf("[MULTI-SET dotplot] %s, KEGG Enriched Pathways multi-set.png",Analysis.Data$Analysis.name),
#       units="in", width=14, height=8, pointsize=8, res=300)
#   print({
#     dotplot(ck, showCategory=30)
#   })
#   dev.off()
#   
#   save(KEGG.classic.Enrichment__results__by.DE.type__by.max.DE.genes, file='KEGG classic enrichment.db')
#   
#   if (Pars$Create.Excel.results){
#     CL = sprintf('%s "%s/KEGG.classic.enrichment.results.to.Excel.py"',Pars$python.bin,Pars$suppl.data.dir)
#     CL = gsub('/','\\',CL,fixed = T)
#     CL = paste(c(CL,sprintf('"%s, KEGG classic enrichment.xlsx"',Analysis.Data$Analysis.name),sprintf('"%s"',KEGG.classic.Enrichment.output.file.names)),collapse = ' ')
#     system(CL)
#     if (!file.exists(sprintf('%s, KEGG classic enrichment.xlsx',Analysis.Data$Analysis.name))){
#       #file.remove(output.file.names)
#       message(sprintf('Excel worksheet generation (KEGG classic enrichment results) for ~%s FAILED',Analysis.Data$Analysis.name)) }
#   }
#   
#   Add.Completed.step(sprintf("%s.%s.%s.KEGG.classic.Enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
#   Analysis.Data$Enriched.Pathways.Classic = levels(factor(Enriched.Pathways.Classic))
# 
#   Analysis.Data$KEGG.classic.Enrichment.performed = T
#   setwd(Analysis.Data$current.GLM.results.dir)
#   cat('\nKEGG Pathway enrichment (classic) completed.\n')
#   return(invisible(Analysis.Data))
# }
# 
# 


classic.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                              database = 'KEGG',
                              printed.name = 'KEGG Pathways',
                              GO.type = 'BP',
                              bypass.if.completed = FALSE,
                              forced.parameters = NULL, additional.plots = TRUE, make.simplify = TRUE, simplify.cutoff = 0.55, make.simplify.plots = TRUE){
  
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(make.simplify && !grepl('GO', database)){
    msg = '\nSimplifying result structure (e.g. reducing redundancy) works only for Gene Ontology\n'
    cat(msg)
    make.simplify = FALSE
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  if (!Analysis.Data$GLM.analysis.performed) stop(sprintf('%s.classic.Enrichment can be processed only after GLM/DE analysis. It is not performed.',database))
  # if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  # } else Pars = forced.parameters
  
  if(database == 'GO'){
    database.plus = sprintf('%s.%s', database, GO.type)
  } else {  database.plus = database  }
  
  if (bypass.if.completed & 
      Read.Completed.steps.status(sprintf("%s.%s.%s.%s.classic.Enrichment", database.plus, Analysis.Data$GLM.model, Startup.Data$Startup.hashmd5, Analysis.Data$step.hashmd5), Startup.Data$Completed.steps.file) &
      Read.Completed.steps.status(sprintf("%s.%s.%s.pathways.visualization", Analysis.Data$GLM.model, Startup.Data$Startup.hashmd5, Analysis.Data$step.hashmd5), Startup.Data$Completed.steps.file)){
    message(sprintf('\n %s.classic.Enrichment pathway enrichment and visualization for GLM model "%s" was previously completed. Bypassing...', database.plus, Analysis.Data$GLM.model))
    return(invisible(Analysis.Data))
  }
  
  ######################################################
  ### classic enrichment with clusterProfiler
  
  convert.ENSEMBL.ids.to.ENTREZID = FALSE
  
  wd = sprintf("%s/%s - classic Enrichment Analysis", Analysis.Data$current.GLM.results.dir, printed.name)
  eval(parse(text = sprintf('Analysis.Data$%s.classic.Enrich.working.dir = wd',database.plus)))
  dir.create(wd, showWarnings = F)
  setwd(wd)
  
  Enrichment__results__by.DE.type__by.max.DE.genes = vector(mode = 'list')
  Gene.Lists__by.DE.type.and.max.DE.genes = vector(mode = 'list')
  
  Enriched.Pathways.Classic = c()

  classic.Enrichment.output.file.names = c()
  simplified.classic.Enrichment.output.file.names = c()
  
  Create.DE.plots.from.clusterProfile.res = TRUE
  if(Create.DE.plots.from.clusterProfile.res){
    min.p.value_by.term.id.UP = vector(mode="list", length=0)
    min.FDR_by.term.id.UP = vector(mode="list", length=0)
    min.p.value_by.term.id.DOWN = vector(mode="list", length=0)
    min.FDR_by.term.id.DOWN = vector(mode="list", length=0)

    enrichment.scores_by.term.id.UP = vector(mode="list", length=0)
    enrichment.scores_by.term.id.DOWN = vector(mode="list", length=0)
    enrichment.score_by.term.id.final = vector(mode="list", length=0)
    enrichment.score_by.term.id.final.UP = vector(mode="list", length=0)
    enrichment.score_by.term.id.final.DOWN = vector(mode="list", length=0)
  }
  
  
  eval(parse(text = sprintf('max.DE.genes.list = Pars$%s.classic.Enrich__with.clusterProfiler___max.DE.genes.list', database)))
  eval(parse(text = sprintf('DE.mode.list = Pars$%s.classic.Enrich__with.clusterProfiler___mode.list', database)))
  eval(parse(text = sprintf('gene.min.abs.LogFC.threshold = Pars$%s.classic.Enrich__with.clusterProfiler___gene.min.abs.LogFC.threshold', database)))
  eval(parse(text = sprintf('gene.max.PValue.threshold = Pars$%s.classic.Enrich__with.clusterProfiler___gene.max.PValue.threshold', database)))
  eval(parse(text = sprintf('gene.min.Score.threshold = Pars$%s.classic.Enrich__with.clusterProfiler___gene.min.Score.threshold', database)))
  
  eval(parse(text = sprintf('generate.Summary.Plot = Pars$%s.classic.Enrich__with.clusterProfiler___generate.Summary.Plot', database)))

  eval(parse(text = sprintf('EEP___limit.Pvalues.to = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___limit.Pvalues.to', database)))
  eval(parse(text = sprintf('EEP___Max.Pvalue.threshold = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___Max.Pvalue.threshold__classic.test', database)))
  eval(parse(text = sprintf('EEP___Max.FDR.threshold = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___Max.FDR.threshold__classic.test', database)))
  eval(parse(text = sprintf('EEP___scoring.adj.list.range.optimal.size = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___scoring.adj.list.range.optimal.size', database)))
  eval(parse(text = sprintf('EEP___scoring.adj.list.range.start = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___scoring.adj.list.range.start', database)))
  eval(parse(text = sprintf('EEP___scoring.adj.list.range.end = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___scoring.adj.list.range.end', database)))

  max.DE.genes.list = sort(max.DE.genes.list, decreasing = FALSE)
  
  enrich.res = NULL
  DE.mode = DE.mode.list[1]  #for debug
  for (DE.mode in DE.mode.list) {
    Enrichment__results__by.DE.type__by.max.DE.genes[[DE.mode]] = vector(mode = 'list')
    Do.end.run <<- FALSE
    max.DE.genes = max.DE.genes.list[1]  #for debug
    for (max.DE.genes in max.DE.genes.list) {
      if (Do.end.run){
        Enrichment__results__by.DE.type__by.max.DE.genes[[DE.mode]][[max.DE.genes]] = enrich.res
        #Gene.Lists__by.DE.type.and.max.DE.genes[[sprintf('top%d %s',max.DE.genes,DE.mode)]] = Entrez.IDs
        next
      }

      if(database == 'Reactome'){
        ### organism can be either "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"
        if (Pars$Species %in% names(Startup.Data$DB.data$Reactome.names.by.KEGG.codes)){
          org = Startup.Data$DB.data$Reactome.names.by.KEGG.codes[Pars$Species]
        } else {
          msg = sprintf('\nSpecies "%s" (%s) is not found among available Reactome species. Allowed only: %s (%s)\n',
                        Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                        paste(Startup.Data$DB.data$Reactome.names.by.KEGG.codes, collapse = ', '),
                        paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Reactome.names.by.KEGG.codes)], collapse = ', '))
          cat(msg)
          warning(msg)
          return(invisible(Analysis.Data))
        }
        
        my_geneList = Get.Entrez.IDs(Startup.Data, Pars, Analysis.Data,
                                    max.DE.genes = max.DE.genes,  DE.type = DE.mode,
                                    min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                    max.PValue = gene.max.PValue.threshold,
                                    min.Score = gene.min.Score.threshold)
        
        if(is.null(my_geneList)) return(invisible(Analysis.Data))
        cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
        #### fun = "gse.res = gsePathway(my_geneList, organism = org, pAdjustMethod = 'BH', pvalueCutoff = FDR.Cutoff)"
        enrich.res = enrichPathway(my_geneList, organism=org, pvalueCutoff=1.00, pAdjustMethod="BH", qvalueCutoff=1)
      } else if (database == 'KEGG'){
        my_geneList = Get.Entrez.IDs(Startup.Data, Pars, Analysis.Data,
                                     max.DE.genes = max.DE.genes,  DE.type = DE.mode,
                                     min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                     max.PValue = gene.max.PValue.threshold,
                                     min.Score = gene.min.Score.threshold)
        if(is.null(my_geneList)) return(invisible(Analysis.Data))
        cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
        #### fun = "gse.res = gseKEGG(my_geneList, organism = Pars$Species, pAdjustMethod = 'BH', pvalueCutoff = FDR.Cutoff)"
        if(Pars$Species == 'dme'){
          enrich.res = enrichKEGG(my_geneList, organism=Pars$Species, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0, keyType = 'ncbi-geneid')
        } else {
          enrich.res = enrichKEGG(my_geneList, organism=Pars$Species, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)
        }
        # print(my_geneList)
        # stop()
      } else if (grepl('GO', database)){
        if (Pars$Species %in% names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
          org.DB = Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes[Pars$Species]
        } else {
          msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s)\n',
                        Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                        paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                        paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
          cat(msg)
          warning(msg)
          return(invisible(Analysis.Data))
        }
        
        if(convert.ENSEMBL.ids.to.ENTREZID){
          my_geneList = Get.Entrez.IDs(Startup.Data, Pars, Analysis.Data,
                                       max.DE.genes = max.DE.genes,  DE.type = DE.mode,
                                       min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                       max.PValue = gene.max.PValue.threshold,
                                       min.Score = gene.min.Score.threshold)
          cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
          #### fun = "gse.res = gseGO(geneList = my_geneList, keyType = 'ENTREZID', OrgDb = org.DB, ont = GO.type, pAdjustMethod='BH', pvalueCutoff = FDR.Cutoff)"
          enrich.res = enrichGO(geneList = my_geneList, OrgDb = org.DB, keyType = 'ENTREZID', ont = GO.type, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)
        } else {
          my_geneList = Get.Ensembl.IDs(Startup.Data, Pars, Analysis.Data,
                                      max.DE.genes = max.DE.genes, DE.type = DE.mode,
                                      min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                      max.PValue = gene.max.PValue.threshold,
                                      min.Score = gene.min.Score.threshold)
          cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
          #### fun = "gse.res = gseGO(geneList = my_geneList, keyType = 'ENSEMBL', OrgDb = org.DB, ont = GO.type, pAdjustMethod='BH', pvalueCutoff = FDR.Cutoff)"
          enrich.res = enrichGO(gene = my_geneList, OrgDb = org.DB, keyType = 'ENSEMBL', ont = GO.type, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)
        }
      } else {
        stop('Unknown database')
      }

      if (typeof(enrich.res)!='S4') {
        cat(sprintf('\nNo enriched %s have been found for "~ %s", top %d %s genes', printed.name, Analysis.Data$GLM.model, length(my_geneList), DE.mode))
        next
      }
      
      if (dim(as.data.frame(enrich.res))[1] == 0){
        cat(sprintf('\nNo enriched %s have been found for "~ %s", top %d %s genes', printed.name, Analysis.Data$GLM.model, length(my_geneList), DE.mode))
        next
      }
      
      
      if(Create.DE.plots.from.clusterProfile.res){
        if (grepl('up', DE.mode, ignore.case = TRUE)){
          n = 1
          for(n in 1:dim(enrich.res@result)[1]){
            term = enrich.res@result$ID[n]
            
            # looking fo minimal P-values
            current.p = enrich.res@result$pvalue[n]
            if(is.na(current.p) || is.null(current.p)) current.p = 1
            val = min.p.value_by.term.id.UP[[term]]
            if (is.null(val)) { min.p.value_by.term.id.UP[[term]] = current.p
            } else if (val > current.p) {
              min.p.value_by.term.id.UP[[term]] = current.p
            }
            
            # looking for minimal FDR
            current.FDR = enrich.res@result$p.adjust[n]
            if(is.na(current.FDR) || is.null(current.FDR)) current.FDR = 1
            val = min.FDR_by.term.id.UP[[term]]
            if (is.null(val)) { min.FDR_by.term.id.UP[[term]] = current.FDR
            } else if (val > current.FDR) {
              min.FDR_by.term.id.UP[[term]] = current.FDR
            }
            
            # Enrichment scores [stage I]; first, preparing individual score lists
            current.p = max(current.p, EEP___limit.Pvalues.to)
            current.score = 0
            if (current.p > EEP___Max.Pvalue.threshold){
              current.score = -1
            } else if(current.FDR > EEP___Max.FDR.threshold) {
              current.score = -1
            }
            
            if(current.score >= 0){
              current.score = max(0,log10(0.07)-log10(current.p))
              adjusting.multiplier = EEP___scoring.adj.list.range.optimal.size/(
                min(max(EEP___scoring.adj.list.range.start, max.DE.genes), EEP___scoring.adj.list.range.end))
              if (current.p > EEP___scoring.adj.list.range.optimal.size){
                adjusting.multiplier = adjusting.multiplier^0.4
              }
              current.score = current.score*adjusting.multiplier * 10
            }
            
            val = enrichment.scores_by.term.id.UP[[term]]
            if (is.null(val)) { enrichment.scores_by.term.id.UP[[term]] = c(current.score)
            } else { enrichment.scores_by.term.id.UP[[term]] = c(val, current.score) }
          }
        }
        if (grepl('down', DE.mode, ignore.case = TRUE)){
          n = 1
          for(n in 1:dim(enrich.res@result)[1]){
            term = enrich.res@result$ID[n]
            
            # looking fo minimal P-values
            current.p = enrich.res@result$pvalue[n]
            if(is.na(current.p) || is.null(current.p)) current.p = 1
            val = min.p.value_by.term.id.DOWN[[term]]
            if (is.null(val)) { min.p.value_by.term.id.DOWN[[term]] = current.p
            } else if (val > current.p) {
              min.p.value_by.term.id.DOWN[[term]] = current.p
            }
            
            # looking for minimal FDR
            current.FDR = enrich.res@result$p.adjust[n]
            if(is.na(current.FDR) || is.null(current.FDR)) current.FDR = 1
            val = min.FDR_by.term.id.DOWN[[term]]
            if (is.null(val)) { min.FDR_by.term.id.DOWN[[term]] = current.FDR
            } else if (val > current.FDR) {
              min.FDR_by.term.id.DOWN[[term]] = current.FDR
            }
            
            # Enrichment scores [stage I]; first, preparing individual score lists
            current.p = max(current.p, EEP___limit.Pvalues.to)
            current.score = 0
            if (current.p > EEP___Max.Pvalue.threshold){
              current.score = -1
            } else if(current.FDR > EEP___Max.FDR.threshold) {
              current.score = -1
            }
            
            if(current.score >= 0){
              current.score = max(0,log10(0.07)-log10(current.p))
              adjusting.multiplier = EEP___scoring.adj.list.range.optimal.size/(
                min(max(EEP___scoring.adj.list.range.start, max.DE.genes), EEP___scoring.adj.list.range.end))
              if (current.p > EEP___scoring.adj.list.range.optimal.size){
                adjusting.multiplier = adjusting.multiplier^0.4
              }
              current.score = current.score*adjusting.multiplier * 10
            }
            
            val = enrichment.scores_by.term.id.DOWN[[term]]
            if (is.null(val)) { enrichment.scores_by.term.id.DOWN[[term]] = c(current.score)
            } else { enrichment.scores_by.term.id.DOWN[[term]] = c(val, current.score) }
          }
        }
      }
      # enrich.res.complete = enrich.res
      
      eval(parse(text = sprintf('pvalue.lim = Pars$%s.classic.Enrich__with.clusterProfiler___term.max.Pvalue.threshold', database)))
      eval(parse(text = sprintf('FDR.lim = Pars$%s.classic.Enrich__with.clusterProfiler___term.max.FDR.threshold', database)))
      eval(parse(text = sprintf('qvalue.lim = Pars$%s.classic.Enrich__with.clusterProfiler___term.max.Qvalue.threshold', database)))
      keep = (enrich.res@result$pvalue <= pvalue.lim) & (enrich.res@result$p.adjust <= FDR.lim) & (enrich.res@result$qvalue <= qvalue.lim)
      enrich.res@result = enrich.res@result[keep, ]

      enrich.res@result = enrich.res@result[!duplicated(tolower(enrich.res@result$Description)), ]
      
      Gene.Lists__by.DE.type.and.max.DE.genes[[sprintf('top%d %s',max.DE.genes, DE.mode)]] = my_geneList
      Enrichment__results__by.DE.type__by.max.DE.genes[[DE.mode]][[max.DE.genes]] = enrich.res
      
      cat(sprintf('\n%d enriched %s have been found for "~ %s", top %d %s genes', dim(as.data.frame(enrich.res))[1], printed.name, Analysis.Data$GLM.model, length(my_geneList), DE.mode))

      enrich.res.converted = enrich.res

      if(grepl('GO',database) & !convert.ENSEMBL.ids.to.ENTREZID){
        enrich.res.converted@result = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table__ENSEMBL(enrich.res@result, Startup.Data, Pars, gene.col = 'geneID')
      } else {
        enrich.res.converted@result = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table(enrich.res@result, Startup.Data, Pars, gene.col = 'geneID')
      }
      

      simplified.enrich.res.converted = enrich.res.converted
      if(make.simplify & (nrow(enrich.res@result) > 0)){
        simplified.enrich.res.converted@result = simplified.enrich.res.converted@result[1:min(250, dim(enrich.res.converted@result)[1]),]
        simplified.enrich.res.converted = simplify(simplified.enrich.res.converted, cutoff = simplify.cutoff)
      }

      if(additional.plots & (nrow(enrich.res@result) > 0)){
        if(make.simplify)  dir.create('Reduced redund.',showWarnings = FALSE)
        res.list = vector(mode = 'list')
        res.list[[1]] = enrich.res.converted
        dirs.list = c('')
        if(make.simplify && make.simplify.plots){
          dirs.list = c('', 'Reduced redund./')
          res.list[[2]] = simplified.enrich.res.converted
        }
        
        for(t in 1:(1 + make.simplify)){
          if(nrow(as.data.frame(res.list[[t]])) < 2) next
          png.mod(filename = sprintf("%s[cneplot] %s classic enrich. for top-%d %s genes.png", dirs.list[t], printed.name, max.DE.genes, DE.mode),
              units="in", width=7, height=7, pointsize=3, res=300)
          cnetplot(res.list[[t]],showCategory=20, categorySize="geneNum")
          dev.off()
          
          png.mod(filename = sprintf("%s[dotplot] %s classic enrich. for top-%d %s genes.png", dirs.list[t], printed.name, max.DE.genes, DE.mode),
              units="in", width=12, height=7, pointsize=8, res=300)
          print({
            dotplot(res.list[[t]], showCategory=30) + ggplot2::xlab(sprintf("Gene ratio; %s, top-%d %s genes", printed.name, max.DE.genes, DE.mode)) #+ scale_colour_gradient(limits=c(0, 0.05), high = 'blue', low="red")
          })
          dev.off()
          
          ### optionally, one can use:
          ## p = dotplot(enrich.res.converted,showCategory=30)
          ## ggplot2::ggsave(file="plot.png", plot=p)
          
          #par(mar=c(1,1,1,1))
          png.mod(filename = sprintf("%s[barplot] %s classic enrich. for top-%d %s genes.png", dirs.list[t], printed.name, max.DE.genes, DE.mode),
              units="in", width=12, height=7, pointsize=8, res=300)
          print({
            barplot(res.list[[t]], showCategory=30) + ggplot2::ylab(sprintf("Gene count (per %s) among top-%d %s ones", substr(printed.name, 1, nchar(printed.name)-1), max.DE.genes, DE.mode))
          })
          dev.off()
        }
      }

      # head(as.data.frame(enrich.res.converted))
      # if(grepl('GO',database) & !convert.ENSEMBL.ids.to.ENTREZID){
      #   classic.Enrich.summary.table = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table__ENSEMBL(as.data.frame(enrich.res), Startup.Data, Pars, gene.col = 'geneID')
      # } else {
      #   classic.Enrich.summary.table = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table(as.data.frame(enrich.res), Startup.Data, Pars, gene.col = 'geneID')
      # }
      
      file.name = Verify.path(sprintf("%s, %s classic enrichment for top-%d %s genes.tsv", Analysis.Data$Analysis.name, printed.name, max.DE.genes, DE.mode))
      classic.Enrichment.output.file.names = c(classic.Enrichment.output.file.names,file.name)
      write.table.mod(x = as.data.frame(enrich.res.converted), file = file.name, sep='\t')
      
      if(make.simplify){
        file.name = Verify.path(sprintf("%s, %s classic enrichment for top-%d %s genes (RR).tsv", Analysis.Data$Analysis.name, printed.name, max.DE.genes, DE.mode))
        simplified.classic.Enrichment.output.file.names = c(simplified.classic.Enrichment.output.file.names, file.name)
        write.table.mod(x = as.data.frame(simplified.enrich.res.converted), file = file.name, sep='\t')
      }
      
      # save(enrich.res.complete,sprintf("%s, %s classic enrichment for top-%d %s genes.db", Analysis.Data$Analysis.name, printed.name, max.DE.genes, DE.mode))
      
      Enriched.Pathways.Classic = c(Enriched.Pathways.Classic, rownames(as.data.frame(enrich.res))) # This is a cumulative list of all found enriched pathways
      
    }
  }
  
  if(generate.Summary.Plot){
    tryCatch(expr = {
      if(database == 'Reactome'){
        org = Startup.Data$DB.data$Reactome.names.by.KEGG.codes[Pars$Species]
        ## enrichPathway(my_geneList, organism=org, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=1)
        ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
                             fun = "enrichPathway", organism=org,  pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1)
      } else if (database == 'KEGG') {
        ## enrichKEGG(my_geneList, organism=Pars$Species, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=1)
        if (Pars$Species == 'dme'){
          ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
                               fun = "enrichKEGG", organism=Pars$Species,  pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1,keyType = 'ncbi-geneid')
        } else {
          ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
                               fun = "enrichKEGG", organism=Pars$Species,  pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1)
        }
      } else if (grepl('GO', database)){
        if(convert.ENSEMBL.ids.to.ENTREZID){
          ## enrichGO(geneList = my_geneList, OrgDb = org.DB, keytype = 'ENTREZID', ont = GO.type, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=1)
          ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
                             fun = "enrichGO", OrgDb = org.DB, keyType = 'ENTREZID', ont = GO.type, pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1)
        } else {
          ## enrichGO(geneList = my_geneList, OrgDb = org.DB, keytype = 'ENSEMBL', ont = GO.type, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=1)
          ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
                               fun = "enrichGO", OrgDb = org.DB, keyType = 'ENSEMBL', ont = GO.type, pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1)
        }
      }
    }, error = function (err){
      warning('Creating summary plot has failed')
      ck <- NULL
    })
        
    if(!is.null(ck)){
      ck.mod = ck
      ck.mod@compareClusterResult = ck.mod@compareClusterResult[!duplicated(ck.mod@compareClusterResult$ID),]
      ids.of.duplicated.descriptions = ck.mod@compareClusterResult$ID[duplicated(tolower(ck.mod@compareClusterResult$Description))]
      ck@compareClusterResult = ck@compareClusterResult[!(ck@compareClusterResult$ID %in% ids.of.duplicated.descriptions),]
      
      max.term.name.length = 70
      ck@compareClusterResult$Description = sapply(ck@compareClusterResult$Description, function(x){ if (nchar(x) > max.term.name.length) sprintf("%s...", substr(x, 1, max.term.name.length)) else x })
    
      ### ck@compareClusterResult$Description
      
      font.size = min(14, max(4, 1500/dim(ck@compareClusterResult)[1]))
      font.size.for.min = min(14, max(4, 1000/dim(ck@compareClusterResult)[1]))
    
      png.mod(filename = sprintf("[MS dotplot] %s, Enriched %s multi.png",Analysis.Data$Analysis.name, printed.name),
          units="in", width=14, height=8, pointsize=8, res=300)
      print({
        suppressMessages(dotplot(ck, showCategory = 10) + theme_dose(font.size = font.size) + scale_colour_gradient(limits=c(0, 0.05), high = 'blue', low="red"))
      })
      dev.off()
      
      png.mod(filename = sprintf("[MULTI-SET dotplot] %s, Enriched %s multi-set, max. terms count.png",Analysis.Data$Analysis.name, printed.name),
          units="in", width=14, height=8, pointsize=8, res=300)
      print({
        suppressMessages(dotplot(ck, showCategory = 70) + theme_dose(font.size = font.size.for.min) + scale_colour_gradient(limits=c(0, 0.05), high = 'blue', low="red"))
      })
      dev.off()
    }
  }  
  eval(parse(text = sprintf('%s.classic.Enrichment__results__by.DE.type__by.max.DE.genes = Enrichment__results__by.DE.type__by.max.DE.genes', database.plus)))
  eval(parse(text = sprintf('saveRDS(%s.classic.Enrichment__results__by.DE.type__by.max.DE.genes, file="%s, %s classic enrichment.db")', database.plus, Analysis.Data$Analysis.name, database.plus)))
  
  
  
  if (Pars$Create.Excel.results){
    CL = sprintf('%s "%s/Classic.enrichment.results.to.Excel.py" %s', Pars$python.bin, Pars$suppl.data.dir, database)
    CL = gsub('/','\\',CL,fixed = TRUE)
    res.file.name = sprintf('%s, %s classic enrichment.xlsx', Analysis.Data$Analysis.name, printed.name)
    CL = paste(c(CL, sprintf('"%s"', Verify.path(res.file.name)), sprintf('"%s"', classic.Enrichment.output.file.names)),collapse = ' ')
    system(CL)
    if (!file.exists(Verify.path(res.file.name))){
      message(sprintf('Excel worksheet generation (%s classic enrichment results) for "~ %s" FAILED', printed.name, Analysis.Data$Analysis.name))
    } else {
      file.remove(classic.Enrichment.output.file.names)
      file.copy(Verify.path(res.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, res.file.name)), overwrite = TRUE)
    }
  }
  if (Pars$Create.Excel.results && make.simplify){
    CL = sprintf('%s "%s/Classic.enrichment.results.to.Excel.py" %s', Pars$python.bin, Pars$suppl.data.dir, database)
    CL = gsub('/','\\',CL,fixed = TRUE)
    res.file.name = sprintf('%s, %s classic enrichment (RR).xlsx', Analysis.Data$Analysis.name, printed.name)
    CL = paste(c(CL, sprintf('"%s"',Verify.path(res.file.name)), sprintf('"%s"', simplified.classic.Enrichment.output.file.names)),collapse = ' ')
    system(CL)
    if (!file.exists(Verify.path(res.file.name))){
      message(sprintf('Excel worksheet generation (%s classic enrichment results, RR) for "~ %s" FAILED', printed.name, Analysis.Data$Analysis.name))
    } else {
      file.remove(simplified.classic.Enrichment.output.file.names)
      file.copy(Verify.path(res.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, res.file.name)), overwrite = TRUE)
    }
  }
  
  if(Create.DE.plots.from.clusterProfile.res){
    all.DB.entries = c(names(enrichment.scores_by.term.id.UP),names(enrichment.scores_by.term.id.DOWN))
    all.DB.entries = all.DB.entries[!duplicated(all.DB.entries)]
    for (gt in all.DB.entries){
      score_list_DOWN = c()
      score_list_UP = c()
      if(!is.null(enrichment.scores_by.term.id.UP[[gt]])) score_list_UP = c(score_list_UP,enrichment.scores_by.term.id.UP[[gt]])
      if(!is.null(enrichment.scores_by.term.id.DOWN[[gt]])) score_list_DOWN = c(score_list_DOWN,enrichment.scores_by.term.id.DOWN[[gt]])
      
      score_list = c(score_list_UP, score_list_DOWN)
      score_list_UP = score_list_UP[score_list_UP > 0]
      score_list_DOWN = score_list_DOWN[score_list_DOWN > 0]
      score_list = c(score_list_UP, score_list_DOWN)
      
      eval(parse(text = sprintf('pow = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___sum.scores.with.power', database)))
      score.final = sum(score_list^pow)^(1/pow)
      score.final.UP = sum(score_list_UP^pow)^(1/pow)
      score.final.DOWN = sum(score_list_DOWN^pow)^(1/pow)
      
      enrichment.score_by.term.id.final[[gt]] = score.final
      enrichment.score_by.term.id.final.DOWN[[gt]] = score.final.DOWN
      enrichment.score_by.term.id.final.UP[[gt]] = score.final.UP
      
    }
    
    Analysis.Data$Enrich.Results.Summary.by.DB[[database]] = vector(mode = "list")
    Analysis.Data$Enrich.Results.Summary.by.DB[[database]]$min.p.value_by.term.id.UP = min.p.value_by.term.id.UP
    Analysis.Data$Enrich.Results.Summary.by.DB[[database]]$min.FDR_by.term.id.UP = min.FDR_by.term.id.UP
    Analysis.Data$Enrich.Results.Summary.by.DB[[database]]$min.p.value_by.term.id.DOWN = min.p.value_by.term.id.DOWN
    Analysis.Data$Enrich.Results.Summary.by.DB[[database]]$min.FDR_by.term.id.DOWN = min.FDR_by.term.id.DOWN
    
    Analysis.Data$Enrich.Results.Summary.by.DB[[database]]$enrichment.scores_by.term.id.UP = enrichment.scores_by.term.id.UP
    Analysis.Data$Enrich.Results.Summary.by.DB[[database]]$enrichment.scores_by.term.id.DOWN = enrichment.scores_by.term.id.DOWN
    Analysis.Data$Enrich.Results.Summary.by.DB[[database]]$enrichment.score_by.term.id.final = enrichment.score_by.term.id.final
    Analysis.Data$Enrich.Results.Summary.by.DB[[database]]$enrichment.score_by.term.id.final.UP = enrichment.score_by.term.id.final.UP
    Analysis.Data$Enrich.Results.Summary.by.DB[[database]]$enrichment.score_by.term.id.final.DOWN = enrichment.score_by.term.id.final.DOWN
    
    ## creating data.frame with Enrichment stats
    all.DB.entries = levels(factor(c(names(min.p.value_by.term.id.UP),names(min.p.value_by.term.id.DOWN))))
    all.DB.entries <- sort(all.DB.entries, decreasing = FALSE)
    
    Enrich.info.split = vector(mode="list")
    all.stats = c('min.p.value_by.term.id.UP','min.FDR_by.term.id.UP','min.p.value_by.term.id.DOWN','min.FDR_by.term.id.DOWN',
                  'enrichment.score_by.term.id.final',
                  'enrichment.score_by.term.id.final.UP','enrichment.score_by.term.id.final.DOWN')
    
    for (x in all.stats){
      eval(parse(text=sprintf("prep = %s[all.DB.entries]",x)))
      names(prep) = all.DB.entries
      eval(parse(text=sprintf('Enrich.info.split$%s = prep',x)))
      eval(parse(text=sprintf('Enrich.info.split$%s[sapply(Enrich.info.split$%s,is.null)] = 1',x,x)))
      #eval(parse(text=sprintf('Enrich.info.split$%s[which(is.null(Enrich.info.split$%s))] = 1',x,x)))
      eval(parse(text=sprintf('Enrich.info.split$%s = do.call(rbind,Enrich.info.split$%s)',x,x)))
    }
    
    Enrich.info.split = do.call(cbind, Enrich.info.split)
    colnames(Enrich.info.split) = all.stats
    
    write.table.mod(x=Enrich.info.split, file=sprintf('%s, %s classic Enrich. stats.tsv', Analysis.Data$Analysis.name, database.plus), sep='\t', quote = FALSE, na = "")
    saveRDS(Enrich.info.split, file = sprintf('%s, %s classic Enrich. stats.db', Analysis.Data$Analysis.name, database.plus))
  }
  saveRDS(Enriched.Pathways.Classic, file = sprintf('%s, %s enriched pathways (classic).db', Analysis.Data$Analysis.name, database.plus))
  
  Add.Completed.step(sprintf("%s.%s.%s.%s.classic.Enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5,database.plus),Startup.Data$Completed.steps.file)
  eval(parse(text = sprintf("Analysis.Data$%s.classic.Enriched.Pathways = levels(factor(Enriched.Pathways.Classic))", database.plus)))
  eval(parse(text = sprintf('Analysis.Data$%s.classic.Enrichment.performed = TRUE', database.plus)))
  
  
  setwd(Analysis.Data$current.GLM.results.dir)
  cat(sprintf('\n%s enrichment (classic) completed.\n', printed.name))
  return(invisible(Analysis.Data))
}



Reactome.classic.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                       bypass.if.completed = FALSE, forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = classic.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'Reactome',
                                    printed.name = 'Reactome Pathways',
                                    bypass.if.completed = bypass.if.completed,
                                    forced.parameters = forced.parameters, additional.plots = additional.plots, make.simplify = FALSE)
  
  return(invisible(Analysis.Data))
}

KEGG.classic.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                   bypass.if.completed = FALSE, forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = classic.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'KEGG',
                                    printed.name = 'KEGG Pathways',
                                    bypass.if.completed = bypass.if.completed,
                                    forced.parameters = forced.parameters, additional.plots = additional.plots, make.simplify = FALSE)
  
  return(invisible(Analysis.Data))
}

GO.classic.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                 GO.types = '{from.parameters}', bypass.if.completed = FALSE, forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  
  # if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  # } else Pars = forced.parameters
  
  if (GO.types[1] == '{from.parameters}'){
    GO.types = Pars$GO.Enrich__with.clusterProfiler___Ontology.list
  }
  
  for(GO.type in GO.types){
    Analysis.Data = classic.Enrichment(Startup.Data, Analysis.Data,
                                      database = sprintf('GO',GO.type),
                                      printed.name = sprintf('GO Terms (%s)', GO.type),
                                      GO.type = GO.type,
                                      bypass.if.completed = bypass.if.completed,
                                      forced.parameters = forced.parameters, additional.plots = additional.plots, make.simplify = TRUE)
  }
  return(invisible(Analysis.Data))
}




# 
# KEGG.trends.Enrichment = function(Startup.Data, Analysis.Data, bypass.if.completed = FALSE,
#                                   forced.parameters = NULL, additional.plots = TRUE){
#   if (!Analysis.Data$GLM.analysis.performed) stop('KEGG.trends.Enrichment can be processed only after GLM/DE analysis. It is not performed.')
#   if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
#   } else Pars = forced.parameters
#   
#   if (bypass.if.completed & 
#       Read.Completed.steps.status(sprintf("%s.%s.%s.KEGG.trends.Enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file) &
#       Read.Completed.steps.status(sprintf("%s.%s.%s.pathways.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
#     message(sprintf('\n KEGG.trends.Enrichment pathway enrichment and visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
#     return(Analysis.Data)
#   }
#   
#   ######################################################
#   ### KEGG trends enrichment with clusterProfiler
#   
#   Analysis.Data$KEGG.trends.Enrich.working.dir = sprintf("%s/KEGG Pathways - Trends Enrichment Analysis",Analysis.Data$current.GLM.results.dir)
#   dir.create(Analysis.Data$KEGG.trends.Enrich.working.dir,showWarnings = F)
#   
#   #Enriched.Pathways.Trended = c()
#   setwd(Analysis.Data$KEGG.trends.Enrich.working.dir)
#   
#   my_geneList = Get.complete.LogFCs.with.Entrez.IDs(Startup.Data, Pars, Analysis.Data)
#   
#   ### organism = "human", "mouse" and "yeast" 
#   gse.res = gseKEGG(my_geneList, organism=Pars$Species, pvalueCutoff=0.05)
#   
#   Enrichment.Successful = TRUE
#   if (typeof(gse.res)!='S4') {
#     cat(sprintf('\nNo enriched KEGG pathways have been found for %s (trends Enrich.)',Analysis.Data$GLM.model))
#     Enrichment.Successful = FALSE
#   }
#   
#   if(Enrichment.Successful){
#     if (dim(as.data.frame(gse.res))[1] == 0){
#       cat(sprintf('\nNo enriched KEGG pathways have been found for %s, top %d %s genes',Analysis.Data$GLM.model,XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode))
#       Enrichment.Successful = FALSE
#     }
#   }
#   save(gse.res, file='KEGG trends enrichment.db')
#   if(! Enrichment.Successful){
#     return()
#   }
#   cat(sprintf('\n%d enriched KEGG pathways have been found for %s (trends Enrich.)\n',dim(as.data.frame(gse.res))[1],Analysis.Data$GLM.model))
# 
#   KEGG.trends.Enrich.summary.table = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table(as.data.frame(gse.res), Startup.Data, Pars, gene.col = 'core_enrichment')
#   
#   png.mod(filename = sprintf("[enrichMap] %s, KEGG trends enrich.png",Analysis.Data$Analysis.name),
#       units="in", width=12, height=7, pointsize=8, res=300)
#   #print({
#     enrichMap(gse.res, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
#   #})
#   dev.off()
# 
#   png.mod(filename = sprintf("[cneplot] %s, KEGG trends enrich.png",Analysis.Data$Analysis.name),
#       units="in", width=7, height=7, pointsize=3, res=300)
#   #cnetplot(gse.res,showCategory=10, categorySize="geneNum", foldChange=my_geneList)
#   cnetplot(gse.res,showCategory=10, categorySize="geneNum")
#   dev.off()
#   
#   ## dotplot, splitted by UP- DOWN-regulated genes
#   #  Thanks crazyhottommy and GuangchuangYu.
#   #  taken from https://github.com/GuangchuangYu/DOSE/issues/20
#   gene_count<- gse.res@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
#   dot_df<- left_join(gse.res@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
#   png.mod(filename = sprintf("[dotplot] %s, KEGG trends enrich.png",Analysis.Data$Analysis.name),
#       units="in", width=12, height=7, pointsize=8, res=300)
#   print({
#     #dotplot(gse.res, showCategory=30) + ggplot2::xlab(sprintf("Gene ratio; KEGG trends enrichment"))
#     
#     ## single
#     # ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
#     #   geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     #   theme_bw(base_size = 14) +
#     #   scale_colour_gradient(limits=c(0, 0.05), low="red") +
#     #   ggtitle("KEGG pathway enrichment")
#     
#     # split by UP-DOWN.
#     dot_df$type = "UP-reg"
#     dot_df$type[dot_df$NES < 0] = "DOWN-reg"
#     ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
#       geom_point(aes(size = GeneRatio, color = pvalue)) +
#       theme_bw(base_size = 14) +
#       #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#       scale_colour_gradient(low="red") +
#       xlab("Gene ratio; KEGG pathway enrichment") + ylab(NULL) + facet_grid(.~type)
#     
#   })
#   dev.off()
# 
#   file.name = sprintf("%s, KEGG trends enrichment.tsv",Analysis.Data$Analysis.name)
#   write.table.mod(x = KEGG.trends.Enrich.summary.table, file = file.name,sep='\t')
#   
#   Enriched.Pathways.Trended = rownames(KEGG.trends.Enrich.summary.table)
#   
#   if(additional.plots){
#     gse.res.table = as.data.frame(gse.res)
#     dir.create('GSEA plots', showWarnings = FALSE)
#     n = 1
#     max_GSEA_plots_per_direction = 30
#     pathway.ID.list.UP = gse.res.table$ID[gse.res.table$NES > 0]
#     pathway.ID.list.UP = pathway.ID.list.UP[1:min(max_GSEA_plots_per_direction, length(pathway.ID.list.UP))]
#     pathway.ID.list.DOWN = gse.res.table$ID[gse.res.table$NES < 0]
#     pathway.ID.list.DOWN = pathway.ID.list.UP[1:min(max_GSEA_plots_per_direction, length(pathway.ID.list.DOWN))]
#     pathway.ID.list = c(pathway.ID.list.UP, pathway.ID.list.DOWN)
#     
#     for (pathway.ID in gse.res.table$ID[gse.res.table$NES > 0]){
#       cat(sprintf('\r Rendering GSEA plot %d of %d...',n,dim(gse.res.table)[1]))
#       pathway.name = gsub("[^[:alnum:] ]", "", gse.res.table[pathway.ID,'Description'])
#       pathway.pvalue = gse.res.table[pathway.ID,'pvalue']
#       
#       png.mod(filename = sprintf("GSEA plots/%s - %s, DE genes distrib.png", pathway.ID, pathway.name),
#           units="in", width=12, height=7, pointsize=8, res=300)
#       gseaplot(gse.res, geneSetID = pathway.ID)
#       dev.off()
#       n = n+1
#     }
#   }
#   
# 
#   if (Pars$Create.Excel.results){
#     CL = sprintf('%s "%s/KEGG.trends.enrichment.results.to.Excel.py"',Pars$python.bin,Pars$suppl.data.dir)
#     CL = gsub('/','\\',CL,fixed = T)
#     CL = paste(c(CL,sprintf('"%s, KEGG trends enrichment.xlsx"',Analysis.Data$Analysis.name),sprintf('"%s"',file.name)),collapse = ' ')
#     system(CL)
#     if (!file.exists(sprintf('%s, KEGG trends enrichment.xlsx',Analysis.Data$Analysis.name))){
#       #file.remove(output.file.names)
#       message(sprintf('Excel worksheet generation (KEGG trends enrichment results) for ~%s FAILED',Analysis.Data$Analysis.name)) }
#   }
#   
#   Add.Completed.step(sprintf("%s.%s.%s.KEGG.trends.Enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
#   Analysis.Data$Enriched.Pathways = levels(factor(Enriched.Pathways))
#   
#   Analysis.Data$KEGG.trends.Enrichment.performed = T
#   setwd(Analysis.Data$current.GLM.results.dir)
#   cat('\rKEGG Pathway enrichment (trends) completed.       \n')
#   return(invisible(Analysis.Data))
# }

# Reactome.trends.Enrichment = function(Startup.Data, Analysis.Data, bypass.if.completed = FALSE,
#                                   forced.parameters = NULL, additional.plots = TRUE){
#   if (!Analysis.Data$GLM.analysis.performed) stop('Reactome.trends.Enrichment can be processed only after GLM/DE analysis. It is not performed.')
#   if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
#   } else Pars = forced.parameters
#   
#   if (bypass.if.completed & 
#       Read.Completed.steps.status(sprintf("%s.%s.%s.Reactome.trends.Enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file) &
#       Read.Completed.steps.status(sprintf("%s.%s.%s.pathways.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
#     message(sprintf('\n Reactome.trends.Enrichment pathway enrichment and visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
#     return(Analysis.Data)
#   }
#   
#   ######################################################
#   ### Reactome trends enrichment with clusterProfiler
#   
#   Analysis.Data$Reactome.trends.Enrich.working.dir = sprintf("%s/Reactome Pathways - Trends Enrichment Analysis",Analysis.Data$current.GLM.results.dir)
#   dir.create(Analysis.Data$Reactome.trends.Enrich.working.dir,showWarnings = F)
#   
#   #Enriched.Pathways.Trended = c()
#   setwd(Analysis.Data$Reactome.trends.Enrich.working.dir)
#   
#   my_geneList = Get.complete.LogFCs.with.Entrez.IDs(Startup.Data, Pars, Analysis.Data)
#   
#   ## "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#   #if Pars$Species == 'hsa'
#   gse.res = gsePathway(my_geneList, organism='mouse', pAdjustMethod = 'BH', pvalueCutoff = 0.1)
#   # ego3 <- gseGO(geneList     = my_geneList,
#   #               OrgDb        = org.Mm.eg.db,
#   #               ont          = "BP",
#   #               nPerm        = 1000,
#   #               minGSSize    = 100,
#   #               maxGSSize    = 500,
#   #               pvalueCutoff = 0.05,
#   #               verbose      = TRUE)      
#   # gse.res = gseGO(my_geneList, organism='mouse', pAdjustMethod = 'BH', pvalueCutoff = 0.1)
# 
#   Enrichment.Successful = TRUE
#   if (typeof(gse.res)!='S4') {
#     cat(sprintf('\nNo enriched Reactome pathways have been found for %s (trends Enrich.)',Analysis.Data$GLM.model))
#     Enrichment.Successful = FALSE
#   }
#   
#   if(Enrichment.Successful){
#     if (dim(as.data.frame(gse.res))[1] == 0){
#       cat(sprintf('\nNo enriched Reactome pathways have been found for %s, top %d %s genes',Analysis.Data$GLM.model,XX.classic.Enrich__with.clusterProfiler___max.DE.genes,XX.classic.Enrich__with.clusterProfiler___mode))
#       Enrichment.Successful = FALSE
#     }
#   }
#   save(gse.res, file='Reactome trends enrichment.db')
#   if(! Enrichment.Successful){
#     return()
#   }
#   cat(sprintf('\n%d enriched Reactome pathways have been found for %s (trends Enrich.)\n',dim(as.data.frame(gse.res))[1],Analysis.Data$GLM.model))
#   
#   Reactome.trends.Enrich.summary.table = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table(as.data.frame(gse.res), Startup.Data, Pars, gene.col = 'core_enrichment')
#   
#   png.mod(filename = sprintf("[enrichMap] %s, Reactome trends enrich.png",Analysis.Data$Analysis.name),
#       units="in", width=12, height=7, pointsize=8, res=300)
#   #print({
#   enrichMap(gse.res, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
#   #})
#   dev.off()
#   
#   png.mod(filename = sprintf("[cneplot] %s, Reactome trends enrich.png",Analysis.Data$Analysis.name),
#       units="in", width=7, height=7, pointsize=3, res=300)
#   #cnetplot(gse.res,showCategory=10, categorySize="geneNum", foldChange=my_geneList)
#   cnetplot(gse.res,showCategory=10, categorySize="geneNum")
#   dev.off()
#   
#   ## dotplot, splitted by UP- DOWN-regulated genes
#   #  Thanks crazyhottommy and GuangchuangYu.
#   #  taken from https://github.com/GuangchuangYu/DOSE/issues/20
#   gene_count<- gse.res@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
#   dot_df<- left_join(gse.res@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
#   png.mod(filename = sprintf("[dotplot] %s, Reactome trends enrich.png",Analysis.Data$Analysis.name),
#       units="in", width=12, height=7, pointsize=8, res=300)
#   print({
#     #dotplot(gse.res, showCategory=30) + ggplot2::xlab(sprintf("Gene ratio; Reactome trends enrichment"))
#     
#     ## single
#     # ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
#     #   geom_point(aes(size = GeneRatio, color = p.adjust)) +
#     #   theme_bw(base_size = 14) +
#     #   scale_colour_gradient(limits=c(0, 0.05), low="red") +
#     #   ggtitle("Reactome pathway enrichment")
#     
#     # split by UP-DOWN.
#     dot_df$type = "UP-reg"
#     dot_df$type[dot_df$NES < 0] = "DOWN-reg"
#     ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
#       geom_point(aes(size = GeneRatio, color = pvalue)) +
#       theme_bw(base_size = 14) +
#       #scale_colour_gradient(limits=c(0, 0.10), low="red") +
#       scale_colour_gradient(low="red") +
#       xlab("Gene ratio; Reactome pathway enrichment") + ylab(NULL) + facet_grid(.~type)
#     
#   })
#   dev.off()
#   
#   file.name = sprintf("%s, Reactome trends enrichment.tsv",Analysis.Data$Analysis.name)
#   write.table.mod(x = Reactome.trends.Enrich.summary.table, file = file.name,sep='\t')
#   
#   Enriched.Pathways.Trended = rownames(Reactome.trends.Enrich.summary.table)
#   
#   if(additional.plots){
#     gse.res.table = as.data.frame(gse.res)
#     dir.create('GSEA plots', showWarnings = FALSE)
#     max_GSEA_plots_per_direction = 20
#     pathway.ID.list.UP = gse.res.table$ID[gse.res.table$NES > 0]
#     pathway.ID.list.UP = pathway.ID.list.UP[1:min(max_GSEA_plots_per_direction, length(pathway.ID.list.UP))]
#     pathway.ID.list.DOWN = gse.res.table$ID[gse.res.table$NES < 0]
#     pathway.ID.list.DOWN = pathway.ID.list.UP[1:min(max_GSEA_plots_per_direction, length(pathway.ID.list.DOWN))]
#     pathway.ID.list = c(pathway.ID.list.UP, pathway.ID.list.DOWN)
#     
#     n = 1
#     for (pathway.ID in pathway.ID.list){
#       cat(sprintf('\r Rendering GSEA plot %d of %d...',n,min(max_GSEA_plots_per_direction*2,dim(gse.res.table)[1])))
#       pathway.name = gsub("[^[:alnum:] ]", "", gse.res.table[pathway.ID,'Description'])
#       if(nchar(pathway.name) > 50) pathway.name = sprintf('%s...', substr(pathway.name,1,45))
#       pathway.pvalue = gse.res.table[pathway.ID,'pvalue']
#       
#       png.mod(filename = sprintf("GSEA plots/%s - %s, DE genes distrib.png", pathway.ID, pathway.name),
#           units="in", width=12, height=7, pointsize=8, res=150)
#       gseaplot(gse.res, geneSetID = pathway.ID)
#       dev.off()
#       n = n+1
#     }
#   }
#   
#   
#   if (Pars$Create.Excel.results){
#     CL = sprintf('%s "%s/KEGG.trends.enrichment.results.to.Excel.py"',Pars$python.bin,Pars$suppl.data.dir)
#     CL = gsub('/','\\',CL,fixed = T)
#     CL = paste(c(CL,sprintf('"%s, Reactome trends enrichment.xlsx"',Analysis.Data$Analysis.name),sprintf('"%s"',file.name)),collapse = ' ')
#     system(CL)
#     if (!file.exists(sprintf('%s, Reactome trends enrichment.xlsx',Analysis.Data$Analysis.name))){
#       #file.remove(output.file.names)
#       message(sprintf('Excel worksheet generation (Reactome trends enrichment results) for ~%s FAILED',Analysis.Data$Analysis.name)) }
#   }
#   
#   Add.Completed.step(sprintf("%s.%s.%s.Reactome.trends.Enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
#   Analysis.Data$Enriched.Pathways = levels(factor(Enriched.Pathways))
#   
#   Analysis.Data$Reactome.trends.Enrichment.performed = T
#   setwd(Analysis.Data$current.GLM.results.dir)
#   cat('\rReactome Pathway enrichment (trends) completed.       \n')
#   return(invisible(Analysis.Data))
# }


trends.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                             database = 'KEGG',
                             printed.name = 'KEGG Pathways',
                             GO.type = 'BP',
                             bypass.if.completed = FALSE,
                             forced.parameters = NULL, additional.plots = TRUE, make.simplify = FALSE, simplify.cutoff = 0.55, make.simplify.plots = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if(make.simplify && !grepl('GO', database)){
    msg = '\nSimplifying result structure (e.g. reducing redundancy) works only for Gene Ontology\n'
    cat(msg)
    make.simplify = FALSE
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  if (!Analysis.Data$GLM.analysis.performed) stop(sprintf('%s.trends.Enrichment can be processed only after GLM/DE analysis. It is not performed.',database))
  # if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  # } else Pars = forced.parameters
  
  if(database == 'GO'){
    database.plus = sprintf('%s.%s', database, GO.type)
  } else {  database.plus = database  }
  
  if (bypass.if.completed & 
      Read.Completed.steps.status(sprintf("%s.%s.%s.%s.trends.Enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5,database.plus),Startup.Data$Completed.steps.file) &
      Read.Completed.steps.status(sprintf("%s.%s.%s.pathways.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
    message(sprintf('\n %s.trends.Enrichment pathway enrichment and visualization for GLM model "~ %s" was previously completed. Bypassing...',database.plus,Analysis.Data$GLM.model))
    return(invisible(Analysis.Data))
  }
  
  ######################################################
  ### trends enrichment with clusterProfiler
  
  wd = sprintf("%s/%s - trends Enrichment Analysis", Analysis.Data$current.GLM.results.dir, printed.name)
  dir.create(wd,showWarnings = F)
  setwd(wd)
  eval(parse(text = sprintf('Analysis.Data$%s.trends.Enrich.working.dir = wd', database.plus)))

  convert.ENSEMBL.ids.to.ENTREZID = FALSE
  # FDR.Cutoff = 0.1
  # FDR.Cutoff.increment = 0.8
  # FDR.Cutoff.max = 0.8
  
  ## "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  if(database == 'Reactome'){
    ### organism can be either "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"
    if (Pars$Species %in% names(Startup.Data$DB.data$Reactome.names.by.KEGG.codes)){
      org = Startup.Data$DB.data$Reactome.names.by.KEGG.codes[Pars$Species]
    } else {
      msg = sprintf('\nSpecies "%s" (%s) is not found among available Reactome species. Allowed only: %s (%s)\n',
                    Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                    paste(Startup.Data$DB.data$Reactome.names.by.KEGG.codes, collapse = ', '),
                    paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Reactome.names.by.KEGG.codes)], collapse = ', '))
      cat(msg)
      warning(msg)
      return(invisible(Analysis.Data))
    }
    my_geneList = Get.complete.LogFCs.with.Entrez.IDs(Startup.Data, Pars, Analysis.Data)
    if(is.null(my_geneList)) return(invisible(Analysis.Data))
    fun = "gse.res = gsePathway(my_geneList, organism = org, pAdjustMethod = 'BH', pvalueCutoff = 1.0)"
  } else if (database == 'KEGG'){
    my_geneList = Get.complete.LogFCs.with.Entrez.IDs(Startup.Data, Pars, Analysis.Data)
    if(is.null(my_geneList)) return(invisible(Analysis.Data))
    if(Pars$Species == 'dme'){
      fun = "gse.res = gseKEGG(my_geneList, organism = Pars$Species, pAdjustMethod = 'BH', pvalueCutoff = 1.0, keyType = 'ncbi-geneid')"
    } else {
      fun = "gse.res = gseKEGG(my_geneList, organism = Pars$Species, pAdjustMethod = 'BH', pvalueCutoff = 1.0)"
    }

  } else if (grepl('GO', database)){
    if (Pars$Species %in% names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
      org.DB = Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes[Pars$Species]
    } else {
      msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s)\n',
                    Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                    paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                    paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
      cat(msg)
      warning(msg)
      return(invisible(Analysis.Data))
    }
    if(convert.ENSEMBL.ids.to.ENTREZID){
      my_geneList = Get.complete.LogFCs.with.Entrez.IDs(Startup.Data, Pars, Analysis.Data)
      fun = "gse.res = gseGO(geneList = my_geneList, keyType = 'ENTREZID', OrgDb = org.DB, ont = GO.type, pAdjustMethod='BH', pvalueCutoff = 1.0)"
      
    } else {
      my_geneList = Analysis.Data$ResTable.gene.part[,c('logFC')]
      names(my_geneList) = rownames(Analysis.Data$ResTable.gene.part)
      my_geneList = my_geneList[rev(order(my_geneList))]
      fun = "gse.res = gseGO(geneList = my_geneList, keyType = 'ENSEMBL', OrgDb = org.DB, ont = GO.type, pAdjustMethod='BH', pvalueCutoff = 1.0)"
    }
  } else {
    stop('Unknown database')
  }
  
  #gse.res = gseGO(geneList = my_geneList, keyType = 'ENSEMBL', OrgDb = org.DB, ont = GO.type, pAdjustMethod='BH', pvalueCutoff = FDR.Cutoff)
  eval(parse(text = fun))
  # repeat{
  #   eval(parse(text = fun))
  #   if (dim(as.data.frame(gse.res))[1] > 0 | FDR.Cutoff >= 1){
  #     break
  #   } else {
  #     if(FDR.Cutoff * (1 + FDR.Cutoff.increment) < FDR.Cutoff.max){
  #       new.FDR.Cutoff = min(1, FDR.Cutoff * (1 + FDR.Cutoff.increment))
  #       cat(sprintf('\nNo enriched %s have been found for "~ %s" (trends Enrich.) under specified FDR=%g cutoff. Increasing FDR cutoff to %g...',
  #                   printed.name, Analysis.Data$GLM.model, FDR.Cutoff, new.FDR.Cutoff))
  #       FDR.Cutoff = new.FDR.Cutoff
  #       
  #     } else {
  #       break
  #     }
  #   }
  # }
  
  Enrich.info = vector(mode="list")
  all.stats = c('min.p.value_by.term.id.UP','min.FDR_by.term.id.UP','min.p.value_by.term.id.DOWN','min.FDR_by.term.id.DOWN',
                'enrichment.score_by.term.id.final',
                'enrichment.score_by.term.id.final.UP','enrichment.score_by.term.id.final.DOWN')
  Enrich.info = data.frame(array(dim = c(dim(gse.res@result)[1],length(all.stats))))
  colnames(Enrich.info) = all.stats
  rownames(Enrich.info) = rownames(gse.res@result)
  
  Enrich.info$min.p.value_by.term.id.UP = rep(1, dim(Enrich.info)[1])
  Enrich.info$min.FDR_by.term.id.UP = rep(1, dim(Enrich.info)[1])
  Enrich.info$min.p.value_by.term.id.DOWN = rep(1, dim(Enrich.info)[1])
  Enrich.info$min.FDR_by.term.id.DOWN = rep(1, dim(Enrich.info)[1])
  Enrich.info$enrichment.score_by.term.id.final = rep(0, dim(Enrich.info)[1])
  Enrich.info$enrichment.score_by.term.id.final.UP = rep(0, dim(Enrich.info)[1])
  Enrich.info$enrichment.score_by.term.id.final.DOWN = rep(0, dim(Enrich.info)[1])
  
  
  eval(parse(text = sprintf('EEP___limit.Pvalues.to = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___limit.Pvalues.to', database)))
  eval(parse(text = sprintf('EEP___Max.Pvalue.threshold = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___Max.Pvalue.threshold__trends.test', database)))
  eval(parse(text = sprintf('EEP___Max.FDR.threshold = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___Max.FDR.threshold__trends.test', database)))

  x = 1
  for(x in 1:dim(gse.res@result)[1]){
    current.p = max(gse.res@result$pvalue[x], EEP___limit.Pvalues.to)
    current.FDR = gse.res@result$p.adjust[x]
    if (current.p > EEP___Max.Pvalue.threshold){
      current.score = -1
    } else if(current.FDR > EEP___Max.FDR.threshold) {
      current.score = -1
    } else {
      #current.score = 10 * abs(gse.res@result$NES[x])
      current.score = 50*max(0,log10(0.07)-log10(current.p))
    }
    Enrich.info[x, 'enrichment.score_by.term.id.final'] = current.score
    
    if(gse.res@result$enrichmentScore[x] > 0){
      Enrich.info[x, 'min.p.value_by.term.id.UP'] = gse.res@result$pvalue[x]
      Enrich.info[x, 'min.FDR_by.term.id.UP'] = gse.res@result$p.adjust[x]
      Enrich.info[x, 'enrichment.score_by.term.id.final.UP'] = current.score
    } else if (gse.res@result$enrichmentScore[x] < 0) {
      Enrich.info[x, 'min.p.value_by.term.id.DOWN'] = gse.res@result$pvalue[x]
      Enrich.info[x, 'min.FDR_by.term.id.DOWN'] = gse.res@result$p.adjust[x]
      Enrich.info[x, 'enrichment.score_by.term.id.final.DOWN'] = current.score
    }
  }
  
  Enrich.info.split = Enrich.info
  
  # for (x in all.stats){
  #   eval(parse(text=sprintf("prep = %s[all.DB.entries]",x)))
  #   names(prep) = all.DB.entries
  #   eval(parse(text=sprintf('Enrich.info.split$%s = prep',x)))
  #   eval(parse(text=sprintf('Enrich.info.split$%s[sapply(Enrich.info.split$%s,is.null)] = 1',x,x)))
  #   eval(parse(text=sprintf('Enrich.info.split$%s = do.call(rbind,Enrich.info.split$%s)',x,x)))
  # }
  # 
  # Enrich.info.split = do.call(cbind, Enrich.info.split)
  # colnames(Enrich.info.split) = all.stats
  
  write.table.mod(x=Enrich.info.split, file=sprintf('%s, %s trends Enrich. stats.tsv', Analysis.Data$Analysis.name, database.plus), sep='\t', quote = FALSE, na = "")
  saveRDS(Enrich.info.split, file = sprintf('%s, %s trends Enrich. stats.db', Analysis.Data$Analysis.name, database.plus))
  
  eval(parse(text = sprintf('pvalue.lim = Pars$%s.trends.Enrich__with.clusterProfiler___term.max.Pvalue.threshold', database)))
  eval(parse(text = sprintf('FDR.lim = Pars$%s.trends.Enrich__with.clusterProfiler___term.max.FDR.threshold', database)))
  eval(parse(text = sprintf('qvalue.lim = Pars$%s.trends.Enrich__with.clusterProfiler___term.max.Qvalue.threshold', database)))
  
  keep = (gse.res@result$pvalue <= pvalue.lim) & (gse.res@result$p.adjust <= FDR.lim) & (gse.res@result$qvalue <= qvalue.lim)
  gse.res@result = gse.res@result[keep, ]
  
  Enrichment.Successful = TRUE
  if (typeof(gse.res)!='S4') {
    cat(sprintf('\nNo enriched %s have been found for "~ %s" (trends Enrich.)',printed.name, Analysis.Data$GLM.model))
    Enrichment.Successful = FALSE
  }
  
  if(Enrichment.Successful){
    if (dim(as.data.frame(gse.res))[1] == 0){
      cat(sprintf('\nNo enriched %s have been found for "~ %s" (trends Enrich.)',printed.name, Analysis.Data$GLM.model))
      Enrichment.Successful = FALSE
    }
  }
  saveRDS(gse.res, file=sprintf('%s, %s trends enrichment.db', Analysis.Data$Analysis.name, database.plus))
  if(! Enrichment.Successful){
    return(invisible(Analysis.Data))
  }
  cat(sprintf('\n%d enriched %s have been found for "~ %s" (trends Enrich.)\n',dim(as.data.frame(gse.res))[1], printed.name, Analysis.Data$GLM.model))
  
  gse.res.converted = gse.res
  if(grepl('GO',database) & !convert.ENSEMBL.ids.to.ENTREZID){
    gse.res.converted@result = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table__ENSEMBL(gse.res@result, Startup.Data, Pars, gene.col = 'core_enrichment')
  } else {
    gse.res.converted@result = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table(gse.res@result, Startup.Data, Pars, gene.col = 'core_enrichment')
  }
  
  simplified.gse.res.converted = gse.res.converted
  if(make.simplify & (nrow(gse.res@result) > 0)){
    result.pos = gse.res@result[gse.res@result$NES > 0,]
    result.pos = result.pos[1:max(250, ncol(result.pos)),]
    temp.pos = new('enrichResult', ontology = gse.res@setType, organism = gse.res@organism, keytype = gse.res@keytype, result = result.pos)
    simplified.temp.pos = simplify(temp.pos, cutoff = simplify.cutoff)
    
    result.neg = gse.res@result[gse.res@result$NES < 0,]
    result.neg = result.neg[1:max(250, ncol(result.neg)),]
    temp.neg = new('enrichResult', ontology = gse.res@setType, organism = gse.res@organism, keytype = gse.res@keytype, result = result.neg)
    simplified.temp.neg = simplify(temp.neg, cutoff = simplify.cutoff)
    
    result.joint = rbind(simplified.temp.pos@result, simplified.temp.neg@result)
    result.joint = result.joint[order(result.joint$pvalue),]
    
    simplified.gse.res.converted = gse.res.converted
    simplified.gse.res.converted@result = result.joint
  }
  
  if(nrow(gse.res@result) > 0){
    if(make.simplify)  dir.create('Reduced redund.',showWarnings = FALSE)
    res.list = vector(mode = 'list')
    res.list[[1]] = gse.res.converted
    prefixes = c('')
    if(make.simplify && make.simplify.plots){
      prefixes = c('', 'Reduced redund. ')
      res.list[[2]] = simplified.gse.res.converted
    }
    
    for(t in 1:(1 + make.simplify)){
      png.mod(filename = sprintf("%s[enrichMap] %s, %s trends enrich.png", prefixes[[t]], Analysis.Data$Analysis.name, printed.name),
          units="in", width=12, height=7, pointsize=8, res=300)
      #print({
      enrichMap(res.list[[t]], layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
      #})
      dev.off()
      
      png.mod(filename = sprintf("%s[cneplot] %s, %s trends enrich.png", prefixes[[t]], Analysis.Data$Analysis.name, printed.name),
          units="in", width=7, height=7, pointsize=3, res=300)
      #cnetplot(gse.res,showCategory=10, categorySize="geneNum", foldChange=my_geneList)
      cnetplot(res.list[[t]],showCategory=10, categorySize="geneNum")
      dev.off()
      
      ## dotplot, splitted by UP- DOWN-regulated genes
      #  Thanks crazyhottommy and GuangchuangYu.
      #  taken from https://github.com/GuangchuangYu/DOSE/issues/20
      gene_count<- res.list[[t]]@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
      
      dot_df<- left_join(res.list[[t]]@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
      max.terms.count = 70
      max.term.name.length = 70
      if(dim(dot_df)[1] > max.terms.count) dot_df = dot_df[1:max.terms.count,]
      dot_df[,'Description'] = sapply(dot_df[,'Description'], function(x){ if (nchar(x) > max.term.name.length) sprintf("%s...", substr(x, 1, max.term.name.length)) else x })
      dot_df$type = "UP-reg"
      dot_df$type[dot_df$NES < 0] = "DOWN-reg"
      
      png.mod(filename = sprintf("%s[dotplot] %s, %s trends enrich.png", prefixes[[t]], Analysis.Data$Analysis.name, printed.name),
          units="in", width=12, height=7, pointsize=8, res=300)
      print({
        #dotplot(gse.res, showCategory=30) + ggplot2::xlab(sprintf("Gene ratio; Reactome trends enrichment"))
        
        ## single
        # ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
        #   geom_point(aes(size = GeneRatio, color = p.adjust)) +
        #   theme_bw(base_size = 14) +
        #   scale_colour_gradient(limits=c(0, 0.05), low="red") +
        #   ggtitle("Reactome pathway enrichment")
        
        # split by UP-DOWN.
        ggplot(dot_df[1:max.terms.count,], aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
          geom_point(aes(size = GeneRatio, color = pvalue)) +
          theme_bw(base_size = min(14, 14 * 40 / dim(dot_df)[1])) +
          #scale_colour_gradient(limits=c(0, 0.10), low="red") +
          #scale_colour_gradient(limits=c(0, 0.05), high = 'blue', low="red") +
          scale_colour_gradient(high = 'blue', low="red") +
          xlab(sprintf("Gene ratio; %s enrichment",printed.name)) + ylab(NULL) + facet_grid(.~type)
        
      })
      dev.off()
    }
  }
  
  enrich.file.name = Verify.path(sprintf("%s, %s trends enrichment.tsv",Analysis.Data$Analysis.name, printed.name))
  write.table.mod(x = as.data.frame(gse.res.converted), file = enrich.file.name,sep='\t')
  if(make.simplify){
    simplified.enrich.file.name = Verify.path(sprintf("%s, %s trends enrichment (RR).tsv",Analysis.Data$Analysis.name, printed.name))
    write.table.mod(x = as.data.frame(simplified.gse.res.converted), file = simplified.enrich.file.name, sep='\t')
  }
  
  Enriched.Pathways.Trended = rownames(as.data.frame(gse.res.converted))
  
  if(additional.plots & (nrow(gse.res@result) > 0)){
    gse.res.table = as.data.frame(gse.res)
    dir.create('GSEA plots', showWarnings = FALSE)
    max_GSEA_plots_per_direction = 10
    pathway.ID.list.UP = gse.res.table$ID[gse.res.table$NES > 0]
    pathway.ID.list.UP = pathway.ID.list.UP[1:min(max_GSEA_plots_per_direction, length(pathway.ID.list.UP))]
    pathway.ID.list.DOWN = gse.res.table$DOWN[gse.res.table$NES < 0]
    pathway.ID.list.DOWN = pathway.ID.list.DOWN[1:min(max_GSEA_plots_per_direction, length(pathway.ID.list.DOWN))]
    pathway.ID.list = c(pathway.ID.list.UP, pathway.ID.list.DOWN)
    
    n = 1
    for (pathway.ID in pathway.ID.list){
      cat(sprintf('\r Rendering GSEA plot %d of %d...',n,min(max_GSEA_plots_per_direction*2,dim(gse.res.table)[1])))
      pathway.name = gsub("[^[:alnum:] ]", "", gse.res.table[pathway.ID,'Description'])
      pathway.ID.str = gsub("[^[:alnum:] ]", "_", pathway.ID)
      if(nchar(pathway.name) > 50) pathway.name = sprintf('%s...', substr(pathway.name,1,45))
      pathway.pvalue = gse.res.table[pathway.ID,'pvalue']
      
      png.mod(filename = sprintf("GSEA plots/%s - %s, DE genes distrib.png", pathway.ID.str, pathway.name),
          units="in", width=12, height=7, pointsize=8, res=150)
      gseaplot(gse.res, geneSetID = pathway.ID)
      dev.off()
      n = n+1
    }
  }
  
  
  if (Pars$Create.Excel.results){
    CL = sprintf('%s "%s/Trends.enrichment.results.to.Excel.py" %s',Pars$python.bin, Pars$suppl.data.dir, database)
    CL = gsub('/','\\',CL,fixed = TRUE)
    res.file.name = sprintf('%s, %s trends enrichment.xlsx',Analysis.Data$Analysis.name, printed.name)
    CL = paste(c(CL, sprintf('"%s"',Verify.path(res.file.name)), sprintf('"%s"', enrich.file.name)), collapse = ' ')
    system(CL)
    if (!file.exists(Verify.path(res.file.name))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (%s trends enrichment results) for "~ %s" FAILED', printed.name, Analysis.Data$Analysis.name))
    } else {
      file.remove(enrich.file.name)
      file.copy(Verify.path(res.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, res.file.name)), overwrite = TRUE)
    }
  }
  
  if (Pars$Create.Excel.results && make.simplify){
    CL = sprintf('%s "%s/Trends.enrichment.results.to.Excel.py" %s',Pars$python.bin, Pars$suppl.data.dir, database)
    CL = gsub('/','\\',CL,fixed = T)
    res.file.name = sprintf('%s, %s trends enrichment (RR).xlsx',Analysis.Data$Analysis.name, printed.name)
    CL = paste(c(CL, sprintf('"%s"',Verify.path(res.file.name)), sprintf('"%s"', simplified.enrich.file.name)), collapse = ' ')
    system(CL)
    if (!file.exists(Verify.path(res.file.name))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (%s trends enrichment results) for "~ %s" FAILED', printed.name, Analysis.Data$Analysis.name))
    } else {
      file.remove(simplified.enrich.file.name)
      file.copy(Verify.path(res.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, res.file.name)), overwrite = TRUE)
    }
  }
  
  # sprintf('%s, %s trends Enrich. stats.db', Analysis.Data$Analysis.name, printed.name)
  saveRDS(Enriched.Pathways.Trended, file = sprintf('%s, %s enriched pathways (trends).db', Analysis.Data$Analysis.name, database.plus))
  
  Add.Completed.step(sprintf("%s.%s.%s.%s.trends.Enrichment", 
                             Analysis.Data$GLM.model, Startup.Data$Startup.hashmd5, Analysis.Data$step.hashmd5, database.plus), Startup.Data$Completed.steps.file)
  tmp = sprintf("Analysis.Data$%s.trends.Enriched.Pathways = levels(factor(Enriched.Pathways.Trended))", database.plus)
  eval(parse(text = tmp))
  
  tmp = sprintf("Analysis.Data$%s.trends.Enrichment.performed = TRUE", database.plus)
  eval(parse(text = tmp))
  
  setwd(Analysis.Data$current.GLM.results.dir)
  
  cat(sprintf('\r%s enrichment (trends) completed.       \n',printed.name))
  return(invisible(Analysis.Data))
}



Reactome.trends.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                      bypass.if.completed = FALSE, forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = trends.Enrichment(Startup.Data, Analysis.Data,
                             database = 'Reactome',
                             printed.name = 'Reactome Pathways',
                             bypass.if.completed = bypass.if.completed,
                             forced.parameters = forced.parameters, additional.plots = additional.plots,
                             make.simplify = FALSE)
  
  return(invisible(Analysis.Data))
}

KEGG.trends.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                  bypass.if.completed = FALSE, forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = trends.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'KEGG',
                                    printed.name = 'KEGG Pathways',
                                    bypass.if.completed = bypass.if.completed,
                                    forced.parameters = forced.parameters, additional.plots = additional.plots,
                                    make.simplify = FALSE)
  
  return(invisible(Analysis.Data))
}

GO.trends.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                GO.types = '{from.parameters}', bypass.if.completed = FALSE, forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  # if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  # } else Pars = forced.parameters
  
  if (GO.types[1] == '{from.parameters}'){
    GO.types = Pars$GO.Enrich__with.clusterProfiler___Ontology.list
  }
  
  for(GO.type in GO.types){
    Analysis.Data = trends.Enrichment(Startup.Data, Analysis.Data,
                                      database = sprintf('GO',GO.type),
                                      printed.name = sprintf('GO Terms (%s)', GO.type),
                                      GO.type = GO.type,
                                      bypass.if.completed = bypass.if.completed,
                                      forced.parameters = forced.parameters, additional.plots = additional.plots,
                                      make.simplify = FALSE)
  }
  return(invisible(Analysis.Data))
}




KEGG.Pathways.Visualization = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                       Add.enriched.pathways = TRUE,
                                       bypass.if.completed = FALSE, forced.parameters = NULL, forced.pathways = NULL,
                                       same.layer = 'adaptive', kegg.native = TRUE, write.de.info = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  if (!Analysis.Data$GLM.analysis.performed) stop('Pathway.Visualization can be processed only after GLM/DE analysis. It is not performed.')
  # if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  # } else Pars = forced.parameters
  
  if (bypass.if.completed & Read.Completed.steps.status(sprintf("%s.%s.%s.pathways.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
    message(sprintf('\nKEGG pathway visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
    return()
  }
  
  ######################################################
  ### Pathway visualization
  
  if(is.null(forced.pathways)){
    Pathways.to.Visualize = Startup.Data$Custom.Pathways.to.Visualize
    if(Add.enriched.pathways){
      if(!is.null(Analysis.Data$KEGG.trends.Enriched.Pathways))      Pathways.to.Visualize = c(Pathways.to.Visualize, Analysis.Data$KEGG.trends.Enriched.Pathways)
      if(!is.null(Analysis.Data$KEGG.classic.Enriched.Pathways))      Pathways.to.Visualize = c(Pathways.to.Visualize, Analysis.Data$KEGG.classic.Enriched.Pathways)
      
      file.name = sprintf("%s/KEGG Pathways - classic Enrichment Analysis/%s, KEGG enriched pathways (classic).db", Analysis.Data$current.GLM.results.dir, Analysis.Data$Analysis.name)
      if(file.exists(file.name)){
        ep = readRDS(file.name)
        Pathways.to.Visualize = c(Pathways.to.Visualize, ep)
      }
      
      file.name = sprintf("%s/KEGG Pathways - trends Enrichment Analysis/%s, KEGG enriched pathways (trends).db", Analysis.Data$current.GLM.results.dir, Analysis.Data$Analysis.name)
      if(file.exists(file.name)){
        ep = readRDS(file.name)
        Pathways.to.Visualize = c(Pathways.to.Visualize, ep)
      }
    }
    
    Pathways.to.Visualize = levels(as.factor(Pathways.to.Visualize))
    #Pathways.to.Visualize = c("hsa01100","mmu01100","dme01100","hsa00511","hsa00514","hsa00533","hsa04215")
    
    forbidden.pathways = c(Startup.Data$Omit.pathways,"01100","00511","00514","00533","04215")
    # fb = forbidden.pathways[1]
    for(fb in forbidden.pathways){
      Pathways.to.Visualize = Pathways.to.Visualize[!grepl(fb, Pathways.to.Visualize)]
    }
    
    #if(same.layer) 
    # Pathways.to.Visualize = Pathways.to.Visualize[!(Pathways.to.Visualize %in% Omit.pathways)]
  } else Pathways.to.Visualize = levels(as.factor(forced.pathways))
  
  #Pathways.to.Visualize = c('dme00190')
  if(same.layer == 'adaptive' & kegg.native){
    pathway.IDs.num = as.numeric(as.character(substr(Pathways.to.Visualize,start=nchar(Pars$Species) + 1,stop=100)))
    
    KEGG.Pathways.Visualization(Startup.Data, Analysis.Data = Analysis.Data, bypass.if.completed = bypass.if.completed,
                                     forced.parameters = forced.parameters, forced.pathways = Pathways.to.Visualize[pathway.IDs.num < 2000],
                                     same.layer = FALSE, kegg.native = TRUE, write.de.info = write.de.info)
    KEGG.Pathways.Visualization(Startup.Data, Analysis.Data = Analysis.Data, bypass.if.completed = bypass.if.completed,
                                     forced.parameters = forced.parameters, forced.pathways = Pathways.to.Visualize[pathway.IDs.num >= 2000],
                                     same.layer = TRUE, kegg.native = TRUE, write.de.info = write.de.info)
    Add.Completed.step(sprintf("%s.%s.%s.pathways.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
    Analysis.Data$Pathway.Visualization.performed = TRUE
    cat('\nPathway visualization completed.\n')
    
  } else {
    cat(sprintf("\n%s:   Performing pathway visualization...",Analysis.Data$GLM.model))
    Analysis.Data$PA.working.dir = sprintf("%s/KEGG Pathways - Visualization",Analysis.Data$current.GLM.results.dir)
    dir.create(Analysis.Data$PA.working.dir, showWarnings = FALSE)
    setwd(Analysis.Data$PA.working.dir)
    
    DE.info.file.names.by.pathway.id = vector(mode="list")
    
    if(Pars$Pathway.Visualization___CPM.aware){
      logFCs = Analysis.Data$ResTable.gene.part$logFC
      CPMs = 2^Analysis.Data$ResTable.gene.part$logCPM
      PValues = Analysis.Data$ResTable.gene.part$PValue
      logFCs [PValues > Pars$Pathway.Visualization___max.gene.PValue] = 0
      names(logFCs) = rownames(Analysis.Data$ResTable.gene.part)
      names(CPMs) = rownames(Analysis.Data$ResTable.gene.part)
      #logFCs = logFCs[tt$table$PValue < Pathway.Visualization___max.gene.PValue]
      #CPMs = CPMs[tt$table$PValue < Pathway.Visualization___max.gene.PValue]
      
      cat(sprintf("\nTotal %d genes has been selected for pathway analysis\n",sum.mod(abs(logFCs) > 0)))
      
      #Pathways.to.Visualize=c('hsa04110')
      pv.out <- suppressWarnings(pathviewmod(gene.data = logFCs,
                                             gene.idtype = "ENSEMBL",
                                             pathway.id = Pathways.to.Visualize, species = Pars$Species, same.layer = same.layer,
                                             out.suffix = Analysis.Data$Analysis.name, keys.align = "y", kegg.native = kegg.native,
                                             limit = list(gene = Pars$Pathway.Visualization___LogFC.limits), bins = 20,cpm.data = CPMs,
                                             write.de.info = write.de.info,cached.dir=sprintf('%s/KEGG.cache',Pars$suppl.data.dir)))
    
    } else {
      logFCs = Analysis.Data$ResTable.gene.part$logFC
      names(logFCs) = rownames(Analysis.Data$ResTable.gene.part)
      logFCs = logFCs[Analysis.Data$ResTable.gene.part$PValue < Pars$Pathway.Visualization___max.gene.PValue]
      
      cat(sprintf("\nTotal %d genes has been selected for the pathway analysis\n",sum.mod(abs(logFCs) > 0)))
      
      pv.out <- suppressWarnings(pathview(gene.data = logFCs,
                                          gene.idtype = "ENSEMBL",
                                          pathway.id = Pathways.to.Visualize, species = Pars$Species, same.layer = same.layer,
                                          out.suffix = Analysis.Data$Analysis.name, keys.align = "y", kegg.native = kegg.native,
                                          limit = list(gene = Pars$Pathway.Visualization___LogFC.limits), bins = list(gene = 12)))
    }
    
    temp = suppressWarnings(file.remove(sprintf("%s.png",Pathways.to.Visualize)))
    temp = suppressWarnings(file.remove(sprintf("%s.xml",Pathways.to.Visualize)))
    keep = Pathways.to.Visualize %in% rownames(Startup.Data$Pathway.names.table)
    ext = ""
    if (!same.layer) ext = " (nsl)"
    temp=suppressWarnings(file.remove(sprintf("%s, %s - %s%s.png",Pathways.to.Visualize[keep],Analysis.Data$Analysis.name,Startup.Data$Pathway.names.table[Pathways.to.Visualize[keep],"Pathway name"],ext)))
    temp=suppressWarnings(file.remove(sprintf("%s, %s - %s, DE info.txt",Pathways.to.Visualize[keep],Analysis.Data$Analysis.name,Startup.Data$Pathway.names.table[Pathways.to.Visualize[keep],"Pathway name"])))
    temp=suppressWarnings(file.remove(sprintf("%s, %s%s.png",Pathways.to.Visualize[!keep],Analysis.Data$Analysis.name,ext)))
    temp=suppressWarnings(file.remove(sprintf("%s, %s, DE info.txt",Pathways.to.Visualize[!keep],Analysis.Data$Analysis.name)))
    temp=file.rename(
      from = sprintf("%s.%s.png",Pathways.to.Visualize[keep],Analysis.Data$Analysis.name),
      to = sprintf("%s, %s - %s%s.png",Pathways.to.Visualize[keep],Analysis.Data$Analysis.name,Startup.Data$Pathway.names.table[Pathways.to.Visualize[keep],"Pathway name"],ext)
    )
    temp=file.rename(
      from = sprintf("%s.%s.png",Pathways.to.Visualize[!keep],Analysis.Data$Analysis.name),
      to = sprintf("%s, %s%s.png",Pathways.to.Visualize[!keep],Analysis.Data$Analysis.name,ext)
    )
    if (write.de.info){
      temp=file.rename(
        from = sprintf("%s DE info.txt",Pathways.to.Visualize[keep]),
        to = sprintf("%s, %s - %s, DE info.txt",Pathways.to.Visualize[keep],Analysis.Data$Analysis.name,Startup.Data$Pathway.names.table[Pathways.to.Visualize[keep],"Pathway name"])
      )
      temp=file.rename(
        from = sprintf("%s DE info.txt",Pathways.to.Visualize[!keep]),
        to = sprintf("%s, %s, DE info.txt",Pathways.to.Visualize[!keep],Analysis.Data$Analysis.name)
      )
    }
    for (pw in Pathways.to.Visualize[keep])  DE.info.file.names.by.pathway.id[[pw]] = sprintf("%s, %s - %s, DE info.txt",pw,Analysis.Data$Analysis.name,Startup.Data$Pathway.names.table[pw,"Pathway name"])
    for (pw in Pathways.to.Visualize[!keep])  DE.info.file.names.by.pathway.id[[pw]] = sprintf("%s, %s, DE info.txt",pw,Analysis.Data$Analysis.name)
    
    ### creating pathway-centric summary
    abs.logFC.ticks = c(0.1,0.25,0.5,1,1.5,2,3,4,5)
    if (write.de.info){
      DE.info.file.names.by.pathway.id = unlist(DE.info.file.names.by.pathway.id, recursive = F, use.names = TRUE)
      DE.info.file.names.by.pathway.id = DE.info.file.names.by.pathway.id[file.exists(DE.info.file.names.by.pathway.id)]
      pw.DE.info = data.frame(array(dim = c(length(DE.info.file.names.by.pathway.id),2+4*length(abs.logFC.ticks))))
      rownames(pw.DE.info) = names(DE.info.file.names.by.pathway.id)
      colnames(pw.DE.info) = c('name','nodes',sprintf('nodes with LogFC > %g',abs.logFC.ticks),sprintf('nodes with LogFC < %g',-abs.logFC.ticks),
                               sprintf('nodes with LogFC > %g, %%',abs.logFC.ticks),sprintf('nodes with LogFC < %g, %%',-abs.logFC.ticks))
      pw.DE.info[,'name'] = Startup.Data$Pathway.names.table[names(DE.info.file.names.by.pathway.id),"Pathway name"]
  
      pw = names(DE.info.file.names.by.pathway.id)[1]
      for (pw in names(DE.info.file.names.by.pathway.id)){
        tmp = read.table(DE.info.file.names.by.pathway.id[pw],sep='\t',header = T)
        pw.DE.info[pw,'nodes'] = dim(tmp)[1]
        tick = abs.logFC.ticks[1]
        for (tick in abs.logFC.ticks){
          pw.DE.info[pw,sprintf('nodes with LogFC > %g',tick)] = sum.mod(tmp[,'ge1'] > tick)
          pw.DE.info[pw,sprintf('nodes with LogFC < -%g',tick)] = sum.mod(tmp[,'ge1'] < -tick)
          pw.DE.info[pw,sprintf('nodes with LogFC > %g, %%',tick)] = sum.mod(tmp[,'ge1'] > tick)/dim(tmp)[1]*100
          pw.DE.info[pw,sprintf('nodes with LogFC < -%g, %%',tick)] = sum.mod(tmp[,'ge1'] < -tick)/dim(tmp)[1]*100
        }
      }
      write.table.mod(pw.DE.info,file=sprintf('%s, KEGG pathways DE summary.tsv',Analysis.Data$Analysis.name),sep='\t',quote = FALSE)
    }
    
    
    ## other options for gene.idtype: "SYMBOL"       "GENENAME"     "ENSEMBL"
    # "ENSEMBLPROT"  "UNIGENE"      "UNIPROT"      "ACCNUM"       "ENSEMBLTRANS" "REFSEQ"
    # "ENZYME"       "TAIR"  "PROSITE"      "ORF"
    
  }    
  setwd(Analysis.Data$current.GLM.results.dir)
  
}


topGO.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                            GO.types = c('BP','MF','CC'),
                            bypass.if.completed = FALSE, forced.parameters = NULL){
  for(GO.type in GO.types){
    topGO.Enrichment.single.ontology(Startup.Data, Analysis.Data = Analysis.Data, GLM.model = GLM.model, GLM.working.dir = GLM.working.dir,
                                     Analysis.Data.RDS.file = Analysis.Data.RDS.file, GO.type = GO.type,
                                     bypass.if.completed = bypass.if.completed, forced.parameters = forced.parameters)
  }
}

topGO.Enrichment.single.ontology = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                            GO.type = 'BP',
                            bypass.if.completed = FALSE, forced.parameters = NULL){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  # if (!Analysis.Data$GLM.analysis.performed) stop('topGO.Enrichment can be processed only after GLM/DE analysis. It is not performed.')
  # if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  # } else Pars = forced.parameters
  
  
  if (bypass.if.completed & 
      Read.Completed.steps.status(sprintf("%s.%s.%s.GO.expression.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file) &
      Read.Completed.steps.status(sprintf("%s.%s.%s.GO.enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
    message(sprintf('\nGO enrichment and GO-centric expression profiles visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
    return(Analysis.Data)
  }
  
  ######################################################
  ### GO Enrichment analysis
  
  Analysis.Data$GO.working.dir = sprintf("%s/GO Terms (%s) - topGO - Enrichment Analysis", Analysis.Data$current.GLM.results.dir, GO.type)
  dir.create(Analysis.Data$GO.working.dir, showWarnings = FALSE)
  setwd(Analysis.Data$GO.working.dir)
  
  min.p.value_by.term.id.UP = vector(mode="list", length=0)
  min.FDR_by.term.id.UP = vector(mode="list", length=0)
  min.p.value_by.term.id.DOWN = vector(mode="list", length=0)
  min.FDR_by.term.id.DOWN = vector(mode="list", length=0)
  
  min.p.value.elim_by.term.id.UP = vector(mode="list", length=0)
  min.FDR.elim_by.term.id.UP = vector(mode="list", length=0)
  min.p.value.elim_by.term.id.DOWN = vector(mode="list", length=0)
  min.FDR.elim_by.term.id.DOWN = vector(mode="list", length=0)
  
  ### additional features for the selecting top enriched GO terms across tests (up/down, 40,1000,200,... difexpr genes)
  enrichment.scores_by.term.id.UP = vector(mode="list", length=0)
  enrichment.scores_by.term.id.DOWN = vector(mode="list", length=0)
  enrichment.score_by.term.id.final = vector(mode="list", length=0)
  enrichment.score_by.term.id.final.UP = vector(mode="list", length=0)
  enrichment.score_by.term.id.final.DOWN = vector(mode="list", length=0)
  
  Pars$GO.Enrich__with.topGO___max.DE.genes.list = sort(Pars$GO.Enrich__with.topGO___max.DE.genes.list,decreasing = F)
  
  #all.GOdata <- new("topGOdata", ontology = "BP", allGenes = Universe.Genes, geneSel =  function(Score) TRUE,
  #              description = "Test", annot = annFUN.org, mapping="org.Mm.eg.db", ID="Ensembl")
  
  
  output.file.names = c()
  
  GO.Enrich__with.topGO___mode = Pars$GO.Enrich__with.topGO___mode.list[1]  # for debug
  for (GO.Enrich__with.topGO___mode in Pars$GO.Enrich__with.topGO___mode.list){
    Do.end.run = F
    GO.Enrich__with.topGO___max.DE.genes = Pars$GO.Enrich__with.topGO___max.DE.genes.list[1] # for debug, too
    for (GO.Enrich__with.topGO___max.DE.genes in Pars$GO.Enrich__with.topGO___max.DE.genes.list){
      GO.Enrich__with.topGO___max.DE.genes = min(GO.Enrich__with.topGO___max.DE.genes,dim(Analysis.Data$ResTable.gene.part)[1]-1)
      # if (Read.Completed.steps.status(sprintf("%s.%s.GO.top.%d.%s",GLM.model,step.hashmd5,GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode))) next
      if (Do.end.run){
        # Add.Completed.step(sprintf("%s.%s.GO.top.%d.%s",GLM.model,step.hashmd5,GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode))
        next
      }
      
      
      if (GO.Enrich__with.topGO___mode == 'de'){
        Passed = (abs(Analysis.Data$ResTable.gene.part[,"logFC"]) > Pars$GO.Enrich__with.topGO___gene.min.abs.LogFC.threshold) & (Analysis.Data$ResTable.gene.part[,"PValue"] < Pars$GO.Enrich__with.topGO___gene.max.PValue.threshold) & (Analysis.Data$ResTable.gene.part[,"Score"] > Pars$GO.Enrich__with.topGO___gene.min.Score.threshold)
      }
      if (GO.Enrich__with.topGO___mode == 'upreg'){
        Passed = (Analysis.Data$ResTable.gene.part[,"logFC"] > 0) & (abs(Analysis.Data$ResTable.gene.part[,"logFC"]) > Pars$GO.Enrich__with.topGO___gene.min.abs.LogFC.threshold) & (Analysis.Data$ResTable.gene.part[,"PValue"] < Pars$GO.Enrich__with.topGO___gene.max.PValue.threshold) & (Analysis.Data$ResTable.gene.part[,"Score"] > Pars$GO.Enrich__with.topGO___gene.min.Score.threshold)
      }
      if (GO.Enrich__with.topGO___mode == 'downreg'){
        Passed = (Analysis.Data$ResTable.gene.part[,"logFC"] < 0) & (abs(Analysis.Data$ResTable.gene.part[,"logFC"]) > Pars$GO.Enrich__with.topGO___gene.min.abs.LogFC.threshold) & (Analysis.Data$ResTable.gene.part[,"PValue"] < Pars$GO.Enrich__with.topGO___gene.max.PValue.threshold) & (Analysis.Data$ResTable.gene.part[,"Score"] > Pars$GO.Enrich__with.topGO___gene.min.Score.threshold)
      }
      if (sum.mod(Passed) <= GO.Enrich__with.topGO___max.DE.genes)  Do.end.run = T
      
      Universe.Genes = as.vector(Analysis.Data$ResTable.gene.part[,"Score"])
      Universe.Genes[!Passed] = 0
      names(Universe.Genes) = rownames(Analysis.Data$ResTable.gene.part)
      
      GO.Enrich__with.topGO___gene.min.Score.threshold.effective = sort(Universe.Genes,decreasing = T)[GO.Enrich__with.topGO___max.DE.genes+1]
      selected.genes.count = sum.mod(Universe.Genes > GO.Enrich__with.topGO___gene.min.Score.threshold.effective)

      
      cat(sprintf("\n~ %s:   Performing Gene Ontology GSEA for top %d %s genes...",Analysis.Data$GLM.model,sum.mod(Universe.Genes > GO.Enrich__with.topGO___gene.min.Score.threshold.effective),GO.Enrich__with.topGO___mode))
      
      topDiffGenes <- function(Score) Score > GO.Enrich__with.topGO___gene.min.Score.threshold.effective
      
      if (Pars$Species %in% names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
        org.DB = Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes[Pars$Species]
      } else {
        msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s)\n',
                      Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                      paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                      paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
        cat(msg)
        warning(msg)
        return(Analysis.Data)
      }
      
      # if (Pars$Species == 'hsa'){
      GOdata <- new("topGOdata", ontology = GO.type, allGenes = Universe.Genes, geneSel = topDiffGenes,
                    description = "Test", annot = annFUN.org, mapping=org.DB, ID="Ensembl")
      # } else if (Pars$Species == 'mmu'){
      #   GOdata <- new("topGOdata", ontology = "BP", allGenes = Universe.Genes, geneSel = topDiffGenes,
      #                 description = "Test", annot = annFUN.org, mapping="org.Mm.eg.db", ID="Ensembl")
      # } else if (Pars$Species == 'dme'){
      #   GOdata <- new("topGOdata", ontology = "BP", allGenes = Universe.Genes, geneSel = topDiffGenes,
      #                 description = "Test", annot = annFUN.org, mapping="org.Dm.eg.db", ID="Ensembl")
      # } else {
      #   stop(sprintf('Incorrect organism: %s', Pars$Species))
      # }
      Total.GO.terms = dim(termStat(GOdata))[1]
      
      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
      resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
      
      allRes <- GenTable(GOdata, P.Fisher = resultFisher,P.Fisher.elim = resultFisher.elim,orderBy = "P.Fisher", topNodes = min(Total.GO.terms,10000))
      #allRes <- GenTable(GOdata, classicFisher = resultFisher,orderBy = "classicFisher", topNodes = 10000)
      
      ## correcting non-numerical "< 1e-30"
      for (n in grep(pattern = '<',x = allRes[,"P.Fisher"])){
        allRes[n,"P.Fisher"] = "1e-30"
      }
      for (n in grep(pattern = '<',x = allRes[,"P.Fisher.elim"])){
        allRes[n,"P.Fisher.elim"] = "1e-30"
      }

      allRes[,"P.Fisher"] = as.numeric(allRes[,"P.Fisher"])
      allRes[,"P.Fisher.elim"] = as.numeric(allRes[,"P.Fisher.elim"])

      #retList = array(dim=c(dim(allRes)[1],1))
      FDR.Fisher = p.adjust(allRes$P.Fisher,method = 'BH')
      FDR.Fisher.elim = p.adjust(allRes$P.Fisher.elim,method = 'BH')
      
      allRes = cbind(allRes,FDR.Fisher,FDR.Fisher.elim)
      
      ## preparing some data for custom GO terms gene expression visualization
      if(GO.Enrich__with.topGO___mode == 'upreg'){
  		  i = 1   # for debug only
        for (i in 1:dim(allRes)[1]){
    			gt = allRes[i,'GO.ID']
    			# P-VALUES
    			val = min.p.value_by.term.id.UP[[gt]]
    			if (is.null(val)) { min.p.value_by.term.id.UP[[gt]] = allRes[i,'P.Fisher']
    			} else if (val > allRes[i,'P.Fisher']) {
    			  min.p.value_by.term.id.UP[[gt]] = allRes[i,'P.Fisher']
    			}
    			# FDRs
    			val = min.FDR_by.term.id.UP[[gt]]
    			if (is.null(val)) { min.FDR_by.term.id.UP[[gt]] = allRes[i,'FDR.Fisher']
    			} else if (val > allRes[i,'FDR.Fisher']) {
    			  min.FDR_by.term.id.UP[[gt]] = allRes[i,'FDR.Fisher']
    			}
    			
    			# P-VALUES, eliminative test
    			val = min.p.value.elim_by.term.id.UP[[gt]]
    			if (is.null(val)) { min.p.value.elim_by.term.id.UP[[gt]] = allRes[i,'P.Fisher.elim']
    			} else if (val > allRes[i,'P.Fisher.elim']) {
    			  min.p.value.elim_by.term.id.UP[[gt]] = allRes[i,'P.Fisher.elim']
    			}
    			# FDRs, eliminative test
    			val = min.FDR.elim_by.term.id.UP[[gt]]
    			if (is.null(val)) { min.FDR.elim_by.term.id.UP[[gt]] = allRes[i,'FDR.Fisher.elim']
    			} else if (val > allRes[i,'FDR.Fisher.elim']) {
    			  min.FDR.elim_by.term.id.UP[[gt]] = allRes[i,'FDR.Fisher.elim']
    			}
    			
    			# Enrichment scores [stage I]; first, preparing individual score lists
    			if (is.na(allRes[i,'P.Fisher'])) { current.p.value = 1
    			} else {  current.p.value = max(allRes[i,'P.Fisher'],Pars$topGO.Enriched.Expression.Profiles___limit.Pvalues.to)  }
    			
    			current.score = 0
    			if (current.p.value > Pars$topGO.Enriched.Expression.Profiles___Max.Pvalue.threshold){
    			  current.score = -1
    			} else if (is.na(allRes[i,'FDR.Fisher'])){ current.score = -1
    			} else if(allRes[i,'FDR.Fisher'] > Pars$topGO.Enriched.Expression.Profiles___Max.FDR.threshold) {
    			  current.score = -1  }
    			
    			if(current.score >= 0){
    			  current.score = max(0,log10(0.07)-log10(current.p.value))
    			  adjusting.multiplier = Pars$topGO.Enriched.Expression.Profiles___scoring.adj.list.range.optimal.size/(
    			    min(max(Pars$topGO.Enriched.Expression.Profiles___scoring.adj.list.range.start,
    			            GO.Enrich__with.topGO___max.DE.genes),Pars$topGO.Enriched.Expression.Profiles___scoring.adj.list.range.end))
    			  if (GO.Enrich__with.topGO___max.DE.genes > Pars$topGO.Enriched.Expression.Profiles___scoring.adj.list.range.optimal.size){
    			    adjusting.multiplier = adjusting.multiplier^0.4
    			  }
    			  current.score = current.score*adjusting.multiplier * 10
    			}
    			
    			val = enrichment.scores_by.term.id.UP[[gt]]
    			if (is.null(val)) { enrichment.scores_by.term.id.UP[[gt]] = c(current.score)
    			} else { enrichment.scores_by.term.id.UP[[gt]] = c(val,current.score) }
    			
  		  }
      } else if(GO.Enrich__with.topGO___mode == 'downreg'){
        for (i in 1:dim(allRes)[1]){
          gt = allRes[i,'GO.ID']
          # P-VALUES
          val = min.p.value_by.term.id.DOWN[[gt]]
          if (is.null(val)) { min.p.value_by.term.id.DOWN[[gt]] = allRes[i,'P.Fisher']
          } else if (val > allRes[i,'P.Fisher']) {
            min.p.value_by.term.id.DOWN[[gt]] = allRes[i,'P.Fisher']
          }
          # FDRs
          val = min.FDR_by.term.id.DOWN[[gt]]
          if (is.null(val)) { min.FDR_by.term.id.DOWN[[gt]] = allRes[i,'FDR.Fisher']
          } else if (val > allRes[i,'FDR.Fisher']) {
            min.FDR_by.term.id.DOWN[[gt]] = allRes[i,'FDR.Fisher']
          }

          # P-VALUES, eliminative test
          val = min.p.value.elim_by.term.id.DOWN[[gt]]
          if (is.null(val)) { min.p.value.elim_by.term.id.DOWN[[gt]] = allRes[i,'P.Fisher.elim']
          } else if (val > allRes[i,'P.Fisher.elim']) {
            min.p.value.elim_by.term.id.DOWN[[gt]] = allRes[i,'P.Fisher.elim']
          }
          # FDRs, eliminative test
          val = min.FDR.elim_by.term.id.DOWN[[gt]]
          if (is.null(val)) { min.FDR.elim_by.term.id.DOWN[[gt]] = allRes[i,'FDR.Fisher.elim']
          } else if (val > allRes[i,'FDR.Fisher.elim']) {
            min.FDR.elim_by.term.id.DOWN[[gt]] = allRes[i,'FDR.Fisher.elim']
          }

          # Enrichment scores [stage I]; first, preparing individual score lists
          if (is.na(allRes[i,'P.Fisher'])) { current.p.value = 1
          } else {
            current.p.value = max(allRes[i,'P.Fisher'],Pars$topGO.Enriched.Expression.Profiles___limit.Pvalues.to)  }
          
          current.score = 0
          if (current.p.value > Pars$topGO.Enriched.Expression.Profiles___Max.Pvalue.threshold){
            current.score = -1
          } else if (is.na(allRes[i,'FDR.Fisher'])){ current.score = -1
          } else if(allRes[i,'FDR.Fisher'] > Pars$topGO.Enriched.Expression.Profiles___Max.FDR.threshold) {
            current.score = -1  }
          
          if(current.score >= 0){
            current.score = max(0,log10(0.07)-log10(current.p.value))
            adjusting.multiplier = Pars$topGO.Enriched.Expression.Profiles___scoring.adj.list.range.optimal.size/(
              min(max(Pars$topGO.Enriched.Expression.Profiles___scoring.adj.list.range.start,
                      GO.Enrich__with.topGO___max.DE.genes),Pars$topGO.Enriched.Expression.Profiles___scoring.adj.list.range.end))
            if (GO.Enrich__with.topGO___max.DE.genes > Pars$topGO.Enriched.Expression.Profiles___scoring.adj.list.range.optimal.size){
              adjusting.multiplier = adjusting.multiplier^0.4
            }
            current.score = current.score*adjusting.multiplier*10
          }
          

          val = enrichment.scores_by.term.id.DOWN[[gt]]
          if (is.null(val)) { enrichment.scores_by.term.id.DOWN[[gt]] = c(current.score)
          } else { enrichment.scores_by.term.id.DOWN[[gt]] = c(val,current.score) }

        }
      } else stop('Unknown GSEA mode. upreg/downreg are only allowed')
      allRes = allRes[allRes$P.Fisher < 0.05,]
      
      #!!!!!!!
      
      #printGenes(object = GOdata, whichTerms = allRes$GO.ID, file = "test.txt",chip = "org.Mm.eg.")
      #GOdata, whichTerms = allRes$GO.ID ,chip = "org.Mm.eg.db", geneCutOff = 1000
      
      
      pValue.fisher <- score(resultFisher)
      for (N in c(5,10,20)) {
        png.mod(filename = sprintf('%s, GO GSEA for top-%s %s genes, top %d GO.terms.png',Analysis.Data$Analysis.name,GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode,N),
            units="in",width=8.23,height=6.3,pointsize=12,
            res=600)
        suppressWarnings(showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = N, useInfo = 'all'))
        dev.off()
      }
      
      
      pdf(file = sprintf('%s, GO GSEA for top-%d %s genes.pdf',Analysis.Data$Analysis.name,GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode),
          width=8.23, height=6.3,onefile = T)
      #graphAttrs$node$fontsize <- "14"
      for (N in c(5,10,20,50)) {
        suppressWarnings(showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = N, useInfo = 'all'))
      }
      dev.off()
      
      
      
      ############################################################
      #### Assigning genes names for each Significant Go Term
      
      whichTerms = allRes$GO.ID
      term.genes = genesInTerm(GOdata, allRes$GO.ID)
      
      retList = array(dim=c(dim(allRes)[1],1))
      colnames(retList) = c('genes')
      rownames(retList) = allRes$GO.ID
      
      ## adding to the results the comlete list of gene names ang LogFCs for enriched GO terms
      all.term.genes = genesInTerm(GOdata, allRes$GO.ID)
      
      all.GT.genes.list = array(dim=c(dim(allRes)[1],1))
      colnames(all.GT.genes.list) = c('All genes in term')
      rownames(all.GT.genes.list) = allRes$GO.ID
      
      all.GT.logFCs.list = array(dim=c(dim(allRes)[1],1))
      colnames(all.GT.logFCs.list) = c('All LogFCs')
      rownames(all.GT.logFCs.list) = allRes$GO.ID
      
      
      for (gt in whichTerms){
        affID <- term.genes[[gt]]
        gene.scores <- sort(geneScore(GOdata, affID, use.names = TRUE),decreasing = T)
        ## adding to the results the complete list of gene names ang LogFCs for enriched GO terms
        all.GT.logFCs = Analysis.Data$ResTable.gene.part[names(gene.scores),"logFC"]
        all.GT.genes = Startup.Data$General.maRt.table[names(gene.scores),"external_gene_name"]
        ## sorting LogFC
        if (Pars$topGO.Expression.Profiles___Sort.LogFCs){
          ord = rev(order(all.GT.logFCs))
          all.GT.genes = all.GT.genes[ord]
          all.GT.logFCs = all.GT.logFCs[ord]
        }
        all.GT.logFCs = paste(all.GT.logFCs, collapse = ', ')
        all.GT.genes = paste(all.GT.genes, collapse = ', ')
        
        all.GT.genes.list[gt,] = all.GT.genes
        all.GT.logFCs.list[gt,] = all.GT.logFCs
        
        ## restricting gene list
        gene.scores <- gene.scores[gene.scores > GO.Enrich__with.topGO___gene.min.Score.threshold.effective]
        GT.genes <- names(gene.scores)
        
        #### Converting EnsemblID to Gene names
        GT.genes <- Startup.Data$General.maRt.table[GT.genes,"external_gene_name"]
        
        
        #### we restrict the output to the number of genes
        if (length(GT.genes) > Pars$GO.Enrich__with.topGO___max.Genes.in.term.to.show){
          length(GT.genes) <- Pars$GO.Enrich__with.topGO___max.Genes.in.term.to.show
          GT.genes = c(GT.genes,'...')
        }
        
        GT.genes = paste(GT.genes, collapse = ', ')
        retList[gt,] = GT.genes
      }
      
      
      allRes = cbind(allRes,retList,Startup.Data$GO.full.descriptions[allRes[,'GO.ID'],c('full description')],all.GT.genes.list,all.GT.logFCs.list)
      
      
      colnames(allRes)[match("Annotated",colnames(allRes))] = 'Genes in genome (background)'
      colnames(allRes)[match("Significant",colnames(allRes))] = sprintf('Genes in top-%d %s',GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode)
      colnames(allRes)[match("Expected",colnames(allRes))] = sprintf('Expected genes in top-%d %s',GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode)
      
      filename = Verify.path(sprintf("%s, GO GSEA for top-%s %s genes.txt",Analysis.Data$Analysis.name,GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode))
      write.table.mod(allRes, file = filename,sep='\t',na = '')
      output.file.names = c(output.file.names,filename)
      
      #Add.Completed.step(sprintf("%s.%s.GO.top.%d.%s",Analysis.Data$GLM.model,step.hashmd5,GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode))
    }
  }

  
  all.GO.term.names = c(names(enrichment.scores_by.term.id.UP),names(enrichment.scores_by.term.id.DOWN))
  all.GO.term.names = all.GO.term.names[!duplicated(all.GO.term.names)]
  for (gt in all.GO.term.names){
    score_list_DOWN = c()
    score_list_UP = c()
    if(!is.null(enrichment.scores_by.term.id.UP[[gt]])) score_list_UP = c(score_list_UP,enrichment.scores_by.term.id.UP[[gt]])
    if(!is.null(enrichment.scores_by.term.id.DOWN[[gt]])) score_list_DOWN = c(score_list_DOWN,enrichment.scores_by.term.id.DOWN[[gt]])
    
    score_list = c(score_list_UP, score_list_DOWN)
    score_list_UP = score_list_UP[score_list_UP > 0]
    score_list_DOWN = score_list_DOWN[score_list_DOWN > 0]
    score_list = c(score_list_UP, score_list_DOWN)
    
    pow = Pars$topGO.Enriched.Expression.Profiles___sum.scores.with.power
    score.final = sum(score_list^pow)^(1/pow)
    score.final.UP = sum(score_list_UP^pow)^(1/pow)
    score.final.DOWN = sum(score_list_DOWN^pow)^(1/pow)
    
    enrichment.score_by.term.id.final[[gt]] = score.final
    enrichment.score_by.term.id.final.DOWN[[gt]] = score.final.DOWN
    enrichment.score_by.term.id.final.UP[[gt]] = score.final.UP
    
  }
  
    
  #Analysis.Data$output.file.names = output.file.names
  if (Pars$Create.Excel.results){
    CL = sprintf('%s "%s/GO.GSEA.results.to.Excel.py"',Pars$python.bin,Pars$suppl.data.dir)
    CL = gsub('/','\\',CL,fixed = T)
    file.base.name = sprintf('%s, GO Terms (%s) topGO enrichment.xlsx', Analysis.Data$Analysis.name, GO.type)
    CL = paste(c(CL, Verify.path(sprintf('"%s"', file.base.name)),sprintf('"%s"',output.file.names)),collapse = ' ')
    system(CL)
    if (!file.exists(Verify.path(file.base.name))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (GO GSEA results) for ~%s FAILED', Analysis.Data$Analysis.name)) }
    
    file.copy(Verify.path(file.base.name),
              Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, file.base.name)),
              overwrite = TRUE)
  }
  
  
  Analysis.Data$min.p.value_by.term.id.UP = min.p.value_by.term.id.UP
  Analysis.Data$min.FDR_by.term.id.UP = min.FDR_by.term.id.UP
  Analysis.Data$min.p.value_by.term.id.DOWN = min.p.value_by.term.id.DOWN
  Analysis.Data$min.FDR_by.term.id.DOWN = min.FDR_by.term.id.DOWN
  
  Analysis.Data$min.p.value.elim_by.term.id.UP = min.p.value.elim_by.term.id.UP
  Analysis.Data$min.FDR.elim_by.term.id.UP = min.FDR.elim_by.term.id.UP
  Analysis.Data$min.p.value.elim_by.term.id.DOWN = min.p.value.elim_by.term.id.DOWN
  Analysis.Data$min.FDR.elim_by.term.id.DOWN = min.FDR.elim_by.term.id.DOWN

  Analysis.Data$enrichment.score_by.term.id.final = enrichment.score_by.term.id.final
  Analysis.Data$enrichment.score_by.term.id.final.UP = enrichment.score_by.term.id.final.UP
  Analysis.Data$enrichment.score_by.term.id.final.DOWN = enrichment.score_by.term.id.final.DOWN
  
    
  Analysis.Data$GO.Enrichment.performed = T
  
  ## creating data.frame with GO stats
  all.GO.terms = levels(factor(c(names(min.p.value_by.term.id.UP),names(min.p.value_by.term.id.DOWN))))
  all.GO.terms <- sort(all.GO.terms,decreasing = F)

  GO.terms.GSEA.info.split = vector(mode="list")
  all.stats = c('min.p.value_by.term.id.UP','min.FDR_by.term.id.UP','min.p.value_by.term.id.DOWN','min.FDR_by.term.id.DOWN',
                'min.p.value.elim_by.term.id.UP','min.FDR.elim_by.term.id.UP','min.p.value.elim_by.term.id.DOWN',
                'min.FDR.elim_by.term.id.DOWN','enrichment.score_by.term.id.final',
                'enrichment.score_by.term.id.final.UP','enrichment.score_by.term.id.final.DOWN')
  for (x in all.stats){
    eval(parse(text=sprintf("prep = %s[all.GO.terms]",x)))
    names(prep) = all.GO.terms
    eval(parse(text=sprintf('GO.terms.GSEA.info.split$%s = prep',x)))
    eval(parse(text=sprintf('GO.terms.GSEA.info.split$%s[sapply(GO.terms.GSEA.info.split$%s,is.null)] = 1',x,x)))
    #eval(parse(text=sprintf('GO.terms.GSEA.info.split$%s[which(is.null(GO.terms.GSEA.info.split$%s))] = 1',x,x)))
    eval(parse(text=sprintf('GO.terms.GSEA.info.split$%s = do.call(rbind,GO.terms.GSEA.info.split$%s)',x,x)))
  }
  
  GO.terms.GSEA.info = do.call(cbind,GO.terms.GSEA.info.split)
  colnames(GO.terms.GSEA.info) = all.stats

  write.table.mod(x=GO.terms.GSEA.info,file=sprintf('%s, GO terms GSEA stats.tsv',Analysis.Data$Analysis.name),sep='\t',quote = F,na = "")

  setwd(Analysis.Data$current.GLM.results.dir)
  Add.Completed.step(sprintf("%s.%s.%s.GO.enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
  
  cat('\nGene Ontology enrichment completed.\n')
  return(invisible(Analysis.Data))
}

topGO.Expression.Profiles.Custom = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                            GO.type = 'BP', profile.type = 'custom',
                                            bypass.if.completed = FALSE, forced.parameters = NULL, forced.terms = NULL,
                                            out.dir = NULL, forced.Analysis.name = NULL,GO.stats.file.name = NULL, sparklines.axis.limits = 2, min.term.size = NULL){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  if(!Analysis.Data$GLM.analysis.performed) stop('GO.Expression.Profiles(...) must be run after GLM/DE analysis. Run Analyze.GLM(...) first.')
  
  # if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  # } else Pars = forced.parameters
  
  if (is.null(forced.Analysis.name)) { Analysis.name = Analysis.Data$Analysis.name
  } else Analysis.name = forced.Analysis.name
  
  if(is.null(out.dir) != is.null(forced.terms))  stop('GO.Expression.Profiles: both forced.terms and out.dir parameters should be both NULL or both have a value')
  
  if (bypass.if.completed & 
      Read.Completed.steps.status(sprintf("%s.%s.%s.GO.expression.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
    message(sprintf('\nGO-centric expression profiles visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
    return()
  }
  
  if (is.null(forced.terms)){
    topGO.Expression.Profiles___Custom.GO.terms = Startup.Data$topGO.Expression.Profiles___Custom.GO.terms
    GO.custom.terms.working.dir = sprintf("%s/GO Terms (%s) - topGO - custom terms DE", Analysis.Data$current.GLM.results.dir, GO.type)
  } else {
    topGO.Expression.Profiles___Custom.GO.terms = forced.terms
    GO.custom.terms.working.dir = out.dir
  }
  
  if (length(topGO.Expression.Profiles___Custom.GO.terms) == 0){
    message('No GO terms to create expression profiles')
    return (invisible())
  }

  dir.create(GO.custom.terms.working.dir, showWarnings = FALSE)
  setwd(GO.custom.terms.working.dir)
  
  
  ######################################################
  ### Creating expression profiles for custom go GO terms
  
  ## looking up if there is GO.terms.GSEA.info table with p, FDR for GO terms. It is produced with GO.Enrichment() method
  if (is.null(GO.stats.file.name)){
    tmp = sprintf('%s/GO Terms (%s) - topGO - Enrichment Analysis/%s, GO terms GSEA stats.tsv', Analysis.Data$current.GLM.results.dir, GO.type, Analysis.name)
    if(file.exists(tmp)) GO.stats.file.name = tmp
  }
    
  if (!is.null(GO.stats.file.name)){
    GO.terms.GSEA.info = read.table(file = GO.stats.file.name,header = T,sep = '\t')
    
    min.p.value_by.term.id.sel.UP = GO.terms.GSEA.info[topGO.Expression.Profiles___Custom.GO.terms,'min.p.value_by.term.id.UP']
    min.p.value_by.term.id.sel.UP[is.na(min.p.value_by.term.id.sel.UP)] = 1
    
    min.FDR_by.term.id.sel.UP = GO.terms.GSEA.info[topGO.Expression.Profiles___Custom.GO.terms,'min.FDR_by.term.id.UP']
    min.FDR_by.term.id.sel.UP[is.na(min.FDR_by.term.id.sel.UP)] = 1
    
    min.p.value_by.term.id.sel.DOWN = GO.terms.GSEA.info[topGO.Expression.Profiles___Custom.GO.terms,'min.p.value_by.term.id.DOWN']
    min.p.value_by.term.id.sel.DOWN[is.na(min.p.value_by.term.id.sel.DOWN)] = 1
    
    min.FDR_by.term.id.sel.DOWN = GO.terms.GSEA.info[topGO.Expression.Profiles___Custom.GO.terms,'min.FDR_by.term.id.DOWN']
    min.FDR_by.term.id.sel.DOWN[is.na(min.FDR_by.term.id.sel.DOWN)] = 1
  } else if (Analysis.Data$GO.Enrichment.performed) {
    min.p.value_by.term.id.sel.UP = data.matrix(Analysis.Data$min.p.value_by.term.id.UP[topGO.Expression.Profiles___Custom.GO.terms])
    rownames(min.p.value_by.term.id.sel.UP) = topGO.Expression.Profiles___Custom.GO.terms
    colnames(min.p.value_by.term.id.sel.UP) = c('min. p-value, upreg')
    
    min.FDR_by.term.id.sel.UP = data.matrix(Analysis.Data$min.FDR_by.term.id.UP[topGO.Expression.Profiles___Custom.GO.terms])
    rownames(min.FDR_by.term.id.sel.UP) = topGO.Expression.Profiles___Custom.GO.terms
    colnames(min.FDR_by.term.id.sel.UP) = c('min. FDR, upreg')
    
    min.p.value_by.term.id.sel.DOWN = data.matrix(Analysis.Data$min.p.value_by.term.id.DOWN[topGO.Expression.Profiles___Custom.GO.terms])
    rownames(min.p.value_by.term.id.sel.DOWN) = topGO.Expression.Profiles___Custom.GO.terms
    colnames(min.p.value_by.term.id.sel.DOWN) = c('min. p-value, downreg')
    
    min.FDR_by.term.id.sel.DOWN = data.matrix(Analysis.Data$min.FDR_by.term.id.DOWN[topGO.Expression.Profiles___Custom.GO.terms])
    rownames(min.FDR_by.term.id.sel.DOWN) = topGO.Expression.Profiles___Custom.GO.terms
    colnames(min.FDR_by.term.id.sel.DOWN) = c('min. FDR, downreg')
  } else {
    warning('GO enrichment analysis has not been performed yet. Run GO.Enrichment(...) or specify GO.stats.file.name argument to include GO stats in the results.')

    min.p.value_by.term.id.sel.UP = data.matrix(rep(1,length(topGO.Expression.Profiles___Custom.GO.terms)))
    rownames(min.p.value_by.term.id.sel.UP) = topGO.Expression.Profiles___Custom.GO.terms
    colnames(min.p.value_by.term.id.sel.UP) = c('min. p-value, upreg')
    
    min.FDR_by.term.id.sel.UP = data.matrix(rep(1,length(topGO.Expression.Profiles___Custom.GO.terms)))
    rownames(min.FDR_by.term.id.sel.UP) = topGO.Expression.Profiles___Custom.GO.terms
    colnames(min.FDR_by.term.id.sel.UP) = c('min. FDR, upreg')
    
    min.p.value_by.term.id.sel.DOWN = data.matrix(rep(1,length(topGO.Expression.Profiles___Custom.GO.terms)))
    rownames(min.p.value_by.term.id.sel.DOWN) = topGO.Expression.Profiles___Custom.GO.terms
    colnames(min.p.value_by.term.id.sel.DOWN) = c('min. p-value, downreg')
    
    min.FDR_by.term.id.sel.DOWN = data.matrix(rep(1,length(topGO.Expression.Profiles___Custom.GO.terms)))
    rownames(min.FDR_by.term.id.sel.DOWN) = topGO.Expression.Profiles___Custom.GO.terms
    colnames(min.FDR_by.term.id.sel.DOWN) = c('min. FDR, downreg')
  }
  
  
  Universe.Genes = as.vector(Analysis.Data$ResTable.gene.part[,"Score"])   ### this is only for topGOdata to be created properly
  names(Universe.Genes) = rownames(Analysis.Data$ResTable.gene.part)
  topDiffGenes <- function(Score) Score > 0           ### this is only for topGOdata to be created properly

  if (Pars$Species %in% names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
    org.DB = Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes[Pars$Species]
  } else {
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s)\n',
                  Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                  paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                  paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
    cat(msg)
    warning(msg)
    return(NULL)
  }
  
  
  # if (Pars$Species == 'hsa'){
  GOdata <- new("topGOdata", ontology = GO.type, allGenes = Universe.Genes, geneSel = topDiffGenes,
                description = "Test", annot = annFUN.org, mapping=org.DB, ID="Ensembl")
  # }
  # if (Pars$Species == 'mmu'){
  #   GOdata <- new("topGOdata", ontology = "BP", allGenes = Universe.Genes, geneSel = topDiffGenes,
  #                 description = "Test", annot = annFUN.org, mapping="org.Mm.eg.db", ID="Ensembl")
  # }
  # if (Pars$Species == 'dme'){
  #   GOdata <- new("topGOdata", ontology = "BP", allGenes = Universe.Genes, geneSel = topDiffGenes,
  #                 description = "Test", annot = annFUN.org, mapping="org.Dm.eg.db", ID="Ensembl")
  # }
  
  countCharOccurrences <- function(char, s) {
    s2 <- gsub(char,"",s)
    return (nchar(s) - nchar(s2))
  }
  
  
  output.file.names = c()
  first.run = TRUE
  topGO.Expression.Profiles___Max.P = Pars$topGO.Expression.Profiles___Max.gene.PValue.list[1] ## for debug
  for (topGO.Expression.Profiles___Max.P in Pars$topGO.Expression.Profiles___Max.gene.PValue.list){
    topGO.Expression.Profiles___Min.gene.logCPM = Pars$topGO.Expression.Profiles___Min.gene.logCPM.list[1] ## for debug
    for (topGO.Expression.Profiles___Min.gene.logCPM in Pars$topGO.Expression.Profiles___Min.gene.logCPM.list){
      L = length(topGO.Expression.Profiles___Custom.GO.terms)
      all.GT.genes.list = array(dim=c(L,1))
      colnames(all.GT.genes.list) = c('All genes in term')
      rownames(all.GT.genes.list) = topGO.Expression.Profiles___Custom.GO.terms
      
      all.GT.logFCs.list = array(dim=c(L,1))
      colnames(all.GT.logFCs.list) = c('All LogFCs')
      rownames(all.GT.logFCs.list) = topGO.Expression.Profiles___Custom.GO.terms
      
      all.GT.genes.counts.list = array(dim=c(L,1))
      colnames(all.GT.genes.counts.list) = c('Gene counts')
      rownames(all.GT.genes.counts.list) = topGO.Expression.Profiles___Custom.GO.terms
      
      term.genes = suppressWarnings(genesInTerm(GOdata, topGO.Expression.Profiles___Custom.GO.terms))
      #gt = 'GO:0016055'
      unavailable.terms = c()
      too.small.terms = c()
      for (gt in topGO.Expression.Profiles___Custom.GO.terms){
        if (is.null(term.genes[[gt]])){
    		  if (!(gt %in% unavailable.terms)){
    		    if (first.run){
      		    if (length(unavailable.terms) < 5){  cat(sprintf('\n Term %s is not found in DB',gt))
      		    } else if (length(unavailable.terms) == 5) cat('\n .... ')
    		    }
            unavailable.terms = c(unavailable.terms,gt)
    		  }
          next
        }
        
        if(!is.null(min.term.size)){
          if(length(term.genes[[gt]]) < min.term.size){
            if (!(gt %in% too.small.terms)){
              if (first.run){
                if (length(too.small.terms) < 5){  cat(sprintf('\n Term %s is too small (%d genes)', gt, length(term.genes[[gt]])))
                } else if (length(too.small.terms) == 5) cat('\n .... ')
              }
              too.small.terms = c(too.small.terms, gt)
            }
            next
          }
        }
        gene.names <- names(sort(geneScore(GOdata, term.genes[[gt]], use.names = TRUE),decreasing = T))
        all.GT.logFCs = Analysis.Data$ResTable.gene.part[gene.names,"logFC"]
        names(all.GT.logFCs) = gene.names
        all.GT.logFCs = all.GT.logFCs[Analysis.Data$ResTable.gene.part[gene.names,"logCPM"] > topGO.Expression.Profiles___Min.gene.logCPM]
        
        gene.names = names(all.GT.logFCs)
        all.GT.logFCs = all.GT.logFCs*(Analysis.Data$ResTable.gene.part[gene.names,"PValue"] < topGO.Expression.Profiles___Max.P)
        all.GT.genes = Startup.Data$General.maRt.table[gene.names,"external_gene_name"]
        ## sorting LogFC
        if (Pars$topGO.Expression.Profiles___Sort.LogFCs){
          ord = rev(order(all.GT.logFCs))
          all.GT.genes = all.GT.genes[ord]
          all.GT.logFCs = all.GT.logFCs[ord]
        }
        
        all.GT.genes.counts.list[gt,] = length(gene.names)
        
        all.GT.logFCs = paste(sprintf("%.2f",all.GT.logFCs), collapse = '\t')
        all.GT.genes = paste(all.GT.genes, collapse = '\t')
        
        all.GT.genes.list[gt,] = all.GT.genes
        all.GT.logFCs.list[gt,] = all.GT.logFCs
        
      }
	    if (first.run && length(unavailable.terms) > 0) cat(sprintf('\n\nTotal %d terms are not found\n',length(unavailable.terms)))
      if (first.run && length(too.small.terms) > 0) cat(sprintf('\n\nTotal %d terms are too small\n', length(too.small.terms)))
      
      first.run = FALSE
      #write.table.mod(cbind(data.matrix(topGO.Expression.Profiles___Custom.GO.terms),GO.full.descriptions[topGO.Expression.Profiles___Custom.GO.terms,c('name','full description')],all.GT.genes.list,all.GT.logFCs.list),file = sprintf('%s, custom GO terms gene-centric info.tsv',Analysis.name),sep='\t')
      file.name = Verify.path(sprintf('%s, GO-DE info, logCPM.gt.%g, P.lt.%g.tsv',Analysis.name,topGO.Expression.Profiles___Min.gene.logCPM,topGO.Expression.Profiles___Max.P))
      
      all.GT.genes.counts.list[is.na(all.GT.genes.counts.list),] = 0
      Custom.GO.table = cbind(data.matrix(topGO.Expression.Profiles___Custom.GO.terms),
                              Startup.Data$GO.full.descriptions[topGO.Expression.Profiles___Custom.GO.terms,c('name')],
                              min.p.value_by.term.id.sel.UP,
                              min.FDR_by.term.id.sel.UP,
                              min.p.value_by.term.id.sel.DOWN,
                              min.FDR_by.term.id.sel.DOWN,
                              all.GT.genes.counts.list,
                              all.GT.logFCs.list)
      
      Custom.GO.table = Custom.GO.table[!(rownames(Custom.GO.table) %in% unavailable.terms),]
      write.table.mod(Custom.GO.table, file = file.name ,sep='\t',quote = F)
      output.file.names = c(output.file.names,file.name)
    }
  }
  
  if (Pars$Create.Excel.results){
    CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces no ',Pars$python.bin,Pars$suppl.data.dir)
    CL = gsub('/','\\',CL,fixed = T)
    res.file.name.base = sprintf('%s, GO Terms (%s) - topGO - %s DE profiles.xlsx', Analysis.name, GO.type, profile.type)
    CL = paste(c(CL,sprintf('--out-excel "%s"', Verify.path(res.file.name.base)),'--in',sprintf('"%s"',output.file.names),
                 '--maximal-genes-count',as.character(Pars$topGO.Expression.Profiles___Maximal.genes.in.terms.to.visualize),
                 sprintf('--sparklines-axis-limit %g',sparklines.axis.limits)),
               collapse = ' ')
    system(CL)
    if (!file.exists(Verify.path(res.file.name.base))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (GO-centric DE profiles) for ~%s FAILED',Analysis.name))
    } else {
      file.copy(Verify.path(res.file.name.base), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, res.file.name.base)), overwrite = TRUE)
    }
    
    CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces yes',Pars$python.bin,Pars$suppl.data.dir)
    CL = gsub('/','\\',CL,fixed = T)
    res.file.name.base = sprintf('%s, GO Terms (%s) - topGO - %s DE profiles (with spaces).xlsx', Analysis.name, GO.type, profile.type)
    CL = paste(c(CL,sprintf('--out-excel "%s"', Verify.path(res.file.name.base)),'--in',sprintf('"%s"',output.file.names),
                 '--maximal-genes-count',as.character(Pars$topGO.Expression.Profiles___Maximal.genes.in.terms.to.visualize),
                 sprintf('--sparklines-axis-limit %g',sparklines.axis.limits)),collapse = ' ')
    system(CL)
    if (!file.exists(Verify.path(res.file.name.base))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (GO-centric DE profiles) for ~%s FAILED',Analysis.name))
    }
  }
  else {
    file.copy(Verify.path(res.file.name.base), sprintf('%s/%s', Verify.path(Analysis.Data$current.GLM.results.dir, res.file.name.base)), overwrite = TRUE)
  }
  
  Add.Completed.step(sprintf("%s.%s.%s.GO.expression.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
  
  
  if (is.null(out.dir)) setwd(Analysis.Data$current.GLM.results.dir)
  Analysis.Data$GO.Expression.Profiles.Created = T
  cat('\nCreating GO-centric gene expression profiles completed.\n')
}

topGO.Expression.Profiles.Enriched = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                              GO.type = 'BP',
                                              bypass.if.completed = FALSE, forced.parameters = NULL, out.dir = NULL, forced.Analysis.name = NULL,
                                              GO.stats.file.name = NULL,sparklines.axis.limits = 2){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  if(!Analysis.Data$GLM.analysis.performed) stop('GO.Expression.Profiles(...) must be run after GLM/DE analysis. Run Analyze.GLM(...) first.')
  
  # if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  # } else Pars = forced.parameters
  
  if (is.null(forced.Analysis.name)) { Analysis.name = Analysis.Data$Analysis.name
  } else Analysis.name = forced.Analysis.name
  
  if (bypass.if.completed & 
      Read.Completed.steps.status(sprintf("%s.%s.%s.GO.expression.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
    message(sprintf('\nGO-centric expression profiles visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
    return()
  }
  
  
  if(is.null(out.dir)){
    GO.enriched.terms.working.dir = sprintf("%s/GO Terms (%s) - topGO - Enriched terms DE", Analysis.Data$current.GLM.results.dir, GO.type)
  } else { GO.enriched.terms.working.dir = out.dir }
  
  # if (length(top.Scored.GO.terms) == 0){
  #   message('No GO terms to create expression profiles')
  #   return ()
  # }
  
  dir.create(GO.enriched.terms.working.dir, showWarnings = FALSE)
  # setwd(GO.custom.terms.working.dir)
  
  
  
  ## looking up if there is GO.terms.GSEA.info table with p, FDR for GO terms. It is produced with GO.Enrichment() method
  if (is.null(GO.stats.file.name)){
    tmp = sprintf('%s/GO Terms (%s) - topGO - Enrichment Analysis/%s, GO terms GSEA stats.tsv', Analysis.Data$current.GLM.results.dir, GO.type, Analysis.name)
    if(file.exists(tmp)) GO.stats.file.name = tmp
  }
  
  
  ### Loading scoring data for various GO terms
  if (!is.null(GO.stats.file.name)){
    GO.terms.GSEA.info = read.table(file = GO.stats.file.name,header = T,sep = '\t')
    enrichment.score_by.term.id.final = GO.terms.GSEA.info[,'enrichment.score_by.term.id.final']
    names(enrichment.score_by.term.id.final) = rownames(GO.terms.GSEA.info)
  } else if (Analysis.Data$GO.Enrichment.performed) {
    enrichment.score_by.term.id.final = data.matrix(Analysis.Data$enrichment.score_by.term.id.final)
  } else {
    stop('GO enrichment analysis has not been performed yet. Run GO.Enrichment(...) or specify GO.stats.file.name argument to include GO stats in the results.')
  }
  
  #Pars$topGO.Enriched.Expression.Profiles___desired.GO.terms.count__to.discard.score = 40
  if(sum.mod(enrichment.score_by.term.id.final > Pars$topGO.Enriched.Expression.Profiles___minimal.score) < Pars$topGO.Enriched.Expression.Profiles___desired.GO.terms.count__to.discard.score){
    adj.score.threshold = sort(enrichment.score_by.term.id.final,decreasing = TRUE)[min(length(enrichment.score_by.term.id.final),Pars$topGO.Enriched.Expression.Profiles___desired.GO.terms.count__to.discard.score)]
    print(sprintf('GO.Expression.Profiles.Enriched(...): score threshol was adjusted from %g to %g',Pars$topGO.Enriched.Expression.Profiles___minimal.score,adj.score.threshold))
  } else { adj.score.threshold = Pars$topGO.Enriched.Expression.Profiles___minimal.score }

  pre.passed.top.Scores.of.GO.terms = enrichment.score_by.term.id.final[enrichment.score_by.term.id.final > adj.score.threshold]
  top.Scores.of.GO.terms = sort(pre.passed.top.Scores.of.GO.terms,decreasing = TRUE)[1:min(length(pre.passed.top.Scores.of.GO.terms),Pars$topGO.Enriched.Expression.Profiles___Max.terms.to.visualize)]
  top.Scores.of.GO.terms = enrichment.score_by.term.id.final[names(enrichment.score_by.term.id.final) %in% names(top.Scores.of.GO.terms)]
  
  ### sorting GO terms with one of three criteria
  
  if (Pars$topGO.Enriched.Expression.Profiles___sort.terms.by == 'score'){
    top.Scored.GO.terms = names(sort(top.Scores.of.GO.terms,decreasing = TRUE))
    
  } else if (Pars$topGO.Enriched.Expression.Profiles___sort.terms.by == 'id'){
    top.Scored.GO.terms = sort(names(top.Scores.of.GO.terms))
    
  } else if (Pars$topGO.Enriched.Expression.Profiles___sort.terms.by == 'name'){
    top.Scored.GO.terms.named_ids = names(top.Scores.of.GO.terms)
    top.Scored.GO.terms.full.names = Startup.Data$GO.full.descriptions[top.Scored.GO.terms.named_ids,c('name')]
    
    names(top.Scored.GO.terms.named_ids) = top.Scored.GO.terms.full.names
    NA.coords = names(top.Scored.GO.terms.named_ids) %in% NA
    NA.count = 1
    for (x in 1:length(NA.coords)){
      if (NA.coords[x]){
        names(top.Scored.GO.terms.named_ids)[x] = sprintf('unknown %d',NA.count)
        NA.count = NA.count+1
      }
    }
    top.Scored.GO.terms.named_ids = top.Scored.GO.terms.named_ids[sort(names(top.Scored.GO.terms.named_ids))]
    top.Scored.GO.terms = top.Scored.GO.terms.named_ids
    names(top.Scored.GO.terms) = NULL
  } else if (Pars$topGO.Enriched.Expression.Profiles___sort.terms.by == 'none'){
    top.Scored.GO.terms = names(top.Scores.of.GO.terms)
  } else {
    stop('Unknown sort mode')
  }
  
  topGO.Expression.Profiles.Custom(Startup.Data, Analysis.Data = Analysis.Data, GO.type = GO.type,
                                   profile.type = 'enriched', forced.parameters = forced.parameters,
                                   forced.terms = top.Scored.GO.terms, out.dir = GO.enriched.terms.working.dir,
                                   forced.Analysis.name = forced.Analysis.name, GO.stats.file.name = GO.stats.file.name,
                                   sparklines.axis.limits = sparklines.axis.limits, min.term.size = Pars$topGO.Enriched.Expression.Profiles___minimal.term.size)
}


Summarize.GLM.results = function(Startup.Data, GLM.models = NULL, bypass.if.completed = FALSE, forced.parameters = NULL, out.file.name = NULL){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if (!Pars$Create.Excel.results){
    warning('Create.Excel.results parameter is turned OFF. "Summarize.GLM.results" will not run')
    return()
  }
  
  if(is.null(GLM.models)) GLM.models = Pars$Models.to.Test
  file.names = c()
  for (GLM.model in GLM.models){
    GLM.model = Delete.prefix.in.model(GLM.model)
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    f = sprintf("%s/~ %s, results/%s_%s_combined.txt",Pars$results.dir,Analysis.name,Analysis.name,Pars$DE.package)
    f = gsub("/","\\",f,fixed = T)
    f = Verify.path(f)
    file.names = append(file.names,f)
  }
  if (sum.mod(!file.exists(file.names))>0){
    warning(sprintf('The following files with DE/GLM info do not exist: %s',toString(file.names[!file.exists(file.names)])))
  }

  
  
  
  ### running non-merging txt > Excel conversion
  CL = sprintf('%s "%s/DE.results.to.Excel.py"',Pars$python.bin, Pars$suppl.data.dir)
  CL = gsub('/','\\',CL,fixed = T)
  #if(is.null(out.dir)) out.dir = Pars$results.dir
  if(is.null(out.file.name)){
    out.file.name = sprintf('%s/Summary expression info.xlsx', Pars$results.dir)
    out.file.name.ff = sprintf('%s/Summary expression info, joint.xlsx', Pars$results.dir)
    out.file.name.lite = sprintf('%s/Summary expression info, joint, lite.xlsx', Pars$results.dir)
    out.file.name.lite.logfconly = sprintf('%s/Summary expression info, joint, lite, LogFC table.xlsx', Pars$results.dir)
  } else {
    out.file.name.ff = sprintf("%s, joint.xlsx", gsub('.{5}$', '', out.file.name))
    out.file.name.lite = sprintf("%s, joint lite.xlsx", gsub('.{5}$', '', out.file.name))
    out.file.name.lite.logfconly = sprintf("%s, joint lite, LogFC table.xlsx", gsub('.{5}$', '', out.file.name))
  }
  CL = paste(c(CL,sprintf('"%s"',Verify.path(out.file.name)),sprintf('"%s"',file.names)),collapse = ' ')
  system(CL)
  if (!file.exists(Verify.path(out.file.name))){
    #file.remove(output.file.names)
    stop(sprintf('Excel worksheet generation FAILED [non-joint], %s',out.file.name))
  }
  
  
  ### running merging txt > Excel conversion
  CL_base = sprintf('%s "%s/DE.results.to.Excel.joint.py"', Pars$python.bin, Pars$suppl.data.dir)
  CL_base = gsub('/','\\', CL_base, fixed = TRUE)
  
  CL.ff = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs_%s.tsv" --lr no --score no --logcpm yes --logcpm-common yes --spearmanp no --pearsonr no --pearsonp no',
                  CL_base, paste(sprintf('"%s"',file.names),collapse=' '), Verify.path(out.file.name.ff), Pars$results.dir, Pars$DE.package)
  CL.lite = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs_%s.tsv" --lr no --score no --logcpm yes --logcpm-common yes --spearmanp no --pearsonr no --pearsonp no --simple-logfc-format-mode yes',
                  CL_base, paste(sprintf('"%s"',file.names),collapse=' '), Verify.path(out.file.name.lite), Pars$results.dir, Pars$DE.package)
  CL.lite.logfconly = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs_%s.tsv" --p no --fdr no --lr no --score no --logcpm no --logcpm-common yes --spearmanr no --spearmanp no --pearsonr no --pearsonp no --simple-logfc-format-mode yes',
                    CL_base, paste(sprintf('"%s"',file.names),collapse=' '), Verify.path(out.file.name.lite.logfconly), Pars$results.dir, Pars$DE.package)
  
  system(CL.ff)
  if (!file.exists(Verify.path(out.file.name.ff))){
    stop(sprintf('Excel worksheet generation FAILED [joint, full format], %s',out.file.name.ff)) }
  system(CL.lite)
  if (!file.exists(Verify.path(out.file.name.lite))){
    stop(sprintf('Excel worksheet generation FAILED [joint, lite format], %s',out.file.name.lite)) }
  system(CL.lite.logfconly)
  if (!file.exists(Verify.path(out.file.name.lite.logfconly))){
    stop(sprintf('Excel worksheet generation FAILED [joint, lite format, LogFC table], %s',out.file.name.lite.logfconly)) }
}

Summarize.inDetails.GLM.results = function(Startup.Data, GLM.models = NULL, bypass.if.completed = FALSE, forced.parameters = NULL,
                                           out.dir = NULL,
                                           GO.terms.inDetails = NULL, add.biotypes = c('lincRNA + antisense_RNA + macro_lincRNA')){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(GLM.models)) GLM.models = Pars$Models.to.Test
  if(is.null(GO.terms.inDetails))  GO.terms.inDetails = Startup.Data$GO.terms.inDetails
  
  if(length(GO.terms.inDetails) == 0){
    cat('\nNo inDetails groups are declared - see "GO terms (inDetails)" sheet\n')
    return(invisible())
  }
  
  all.group.names = c(names(GO.terms.inDetails), add.biotypes)
  all.group.Name_s = gsub("[^[:alnum:] ]", "", all.group.names)
  if(is.null(out.dir))  out.dir = sprintf('%s/inDetails DE summary', Pars$results.dir)
  dir.create(out.dir, showWarnings = FALSE)
  
  for(group.Name_s in all.group.Name_s){
    all.inDetails.DE.file.names = c()
    for(GLM.model in GLM.models){
      Analysis.name  = Cleanup.Model.Name(GLM.model)
      GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
      inDetails.main.wd = sprintf("%s/inDetails", GLM.working.dir)
      
      inDetails.DE.file.name = sprintf("%s/%s/%s - %s_%s_combined.txt",
                                      inDetails.main.wd, group.Name_s, group.Name_s, Analysis.name, Pars$DE.package)
      inDetails.DE.file.name = Verify.path(inDetails.DE.file.name)
      all.inDetails.DE.file.names = c(all.inDetails.DE.file.names, inDetails.DE.file.name)
    }
    keep = file.exists(all.inDetails.DE.file.names)
    if(sum(keep) != length(keep)){
      cat(sprintf('\n[%s]: %d of %d DE analysis results (*.tsv) are not found\n', group.Name_s, length(keep) - sum(keep), length(keep)))
      next
    } else {
      cat(sprintf('\n[%s]: all %d DE analysis results are present', group.Name_s, length(keep)))
    }
    
    all.inDetails.DE.file.names = all.inDetails.DE.file.names[keep]
    all.inDetails.DE.file.names = gsub("/", "\\", all.inDetails.DE.file.names, fixed = TRUE)
    
    out.file.name = Verify.path(sprintf("%s/%s - Summary expression info.xlsx", out.dir, group.Name_s))
    
    ### running non-merging txt > Excel conversion
    CL = sprintf('%s "%s/DE.results.to.Excel.py"',Pars$python.bin,Pars$suppl.data.dir)
    CL = gsub('/','\\',CL,fixed = T)
    CL = paste(c(CL,sprintf('"%s"',out.file.name),sprintf('"%s"',all.inDetails.DE.file.names)),collapse = ' ')
    system(CL)
    if (!file.exists(out.file.name)){
      stop(sprintf('Excel worksheet generation FAILED [non-joint], %s',out.file.name)) }

        
    ### running merging txt > Excel conversion
    out.file.name.ff = Verify.path(sprintf("%s, joint.xlsx",gsub('.{5}$', '', out.file.name)))
    out.file.name.lite = Verify.path(sprintf("%s, joint lite.xlsx",gsub('.{5}$', '', out.file.name)))
    out.file.name.lite.logfconly = Verify.path(sprintf("%s, joint lite, LogFC table.xlsx",gsub('.{5}$', '', out.file.name)))
    CL_base = sprintf('%s "%s/DE.results.to.Excel.joint.py"',Pars$python.bin,Pars$suppl.data.dir)
    CL_base = gsub('/','\\',CL_base,fixed = T)
    CL.ff = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs_%s.tsv" --lr no --score no --logcpm no --logcpm-common yes --spearmanp no --pearsonr no --pearsonp no --exclude-genes-without-de-info yes',
                    CL_base, paste(sprintf('"%s"',all.inDetails.DE.file.names),collapse=' '), out.file.name.ff, Pars$results.dir, Pars$DE.package)
    CL.lite = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs_%s.tsv" --lr no --score no --logcpm no --logcpm-common yes --spearmanp no --pearsonr no --pearsonp no --simple-logfc-format-mode yes --exclude-genes-without-de-info yes',
                      CL_base, paste(sprintf('"%s"',all.inDetails.DE.file.names),collapse=' '), out.file.name.lite, Pars$results.dir, Pars$DE.package)
    CL.lite.logfconly = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs_%s.tsv" --p no --fdr no --lr no --score no --logcpm no --logcpm-common yes --spearmanr no --spearmanp no --pearsonr no --pearsonp no --simple-logfc-format-mode yes --exclude-genes-without-de-info yes',
                                CL_base, paste(sprintf('"%s"',all.inDetails.DE.file.names),collapse=' '), out.file.name.lite.logfconly, Pars$results.dir, Pars$DE.package)
    
    system(CL.ff)
    if (!file.exists(out.file.name.ff)){
      stop(sprintf('Excel worksheet generation FAILED [joint, full format], %s',out.file.name.ff)) }
    system(CL.lite)
    if (!file.exists(out.file.name.lite)){
      stop(sprintf('Excel worksheet generation FAILED [joint, lite format], %s',out.file.name.lite)) }
    system(CL.lite.logfconly)
    if (!file.exists(out.file.name.lite.logfconly)){
      stop(sprintf('Excel worksheet generation FAILED [joint, lite format, LogFC table], %s',out.file.name.lite.logfconly)) }
  }
}


Summarize.topGO.Expression.Profiles.Custom = function(Startup.Data, GLM.models = NULL, GO.type = 'BP',
                                                      bypass.if.completed = FALSE, forced.parameters = NULL,
                                                      out.dir = NULL, out.file.prefix = ''){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  if (length(Startup.Data$topGO.Expression.Profiles___Custom.GO.terms) == 0) return()

  if (!Pars$Create.Excel.results){
    warning('Create.Excel.results parameter is turned OFF. "Summarize.topGO.Expression.Profiles.Custom" will not run')
    return()
  }
  
  if(is.null(GLM.models)){
    GLM.models = Pars$Models.to.Test
  }
  
  if(is.null(out.dir)) out.dir = Pars$results.dir
  
  file.names = c()
  for (GLM.model in GLM.models){
    GLM.model = Delete.prefix.in.model(GLM.model)
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    f = sprintf("%s/~ %s, results/GO Terms (%s) - topGO - custom terms DE/%s, GO-DE info, logCPM.gt.%g, P.lt.%g.tsv",
                Pars$results.dir, Analysis.name, GO.type, Analysis.name,
                Pars$topGO.Expression.Profiles___Min.gene.logCPM.for.summary.across.models, Pars$topGO.Expression.Profiles___Max.gene.PValue.for.summary.across.models)
    f = gsub("/","\\", f, fixed = TRUE)
    f = Verify.path(f)
    file.names = append(file.names, f)
  }
  if (sum.mod(!file.exists(file.names))>0){
    cat(sprintf('[Summarize.topGO.Expression.Profiles.Custom] Total %d files with custom GO-terms expression profiles do not exist: %s\n',sum.mod(!file.exists(file.names)),toString(file.names[!file.exists(file.names)])))
    file.names = file.names[file.exists(file.names)]
    if(length(file.names) == 0){
      cat('No files found\n')
      return(invisible())
    }
  } else {
    cat(sprintf('[Summarize.topGO.Expression.Profiles.Custom] All %d files with custom GO-terms expression profiles exist\n', length(file.names)))
  }
  #setwd(Pars$results.dir)
  
  
  CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces no ',Pars$python.bin,Pars$suppl.data.dir)
  out.file.name = sprintf('%s/%sSummary topGO (%s) DE profiles.xlsx', out.dir, out.file.prefix, GO.type)
  CL = paste(c(CL,sprintf('--out-excel "%s"',Verify.path(out.file.name)),'--in',sprintf('"%s"',file.names),
               '--maximal-genes-count',as.character(Pars$topGO.Expression.Profiles___Maximal.genes.in.terms.to.visualize)),collapse = ' ')
  CL = gsub('/','\\',CL,fixed = TRUE)
  system(CL)
  if (!file.exists(Verify.path(out.file.name)))    stop(sprintf('[Summarize.topGO.Expression.Profiles.Custom] Excel worksheet generation FAILED'))
  
  CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces yes ',Pars$python.bin,Pars$suppl.data.dir)
  out.file.name = sprintf('%s/%sSummary topGO (%s) DE profiles (with spaces).xlsx', out.dir, out.file.prefix, GO.type)
  CL = paste(c(CL,sprintf('--out-excel "%s"',Verify.path(out.file.name)),'--in',sprintf('"%s"',file.names),
               '--maximal-genes-count',as.character(Pars$topGO.Expression.Profiles___Maximal.genes.in.terms.to.visualize)),collapse = ' ')
  CL = gsub('/','\\',CL,fixed = TRUE)
  system(CL)
  if (!file.exists(Verify.path(out.file.name)))    stop(sprintf('[Summarize.topGO.Expression.Profiles.Custom] Excel worksheet generation FAILED'))
}


Summarize.clusterProfiler.Expression.Profiles.Custom = function(Startup.Data, GLM.models = NULL, database = 'GO', printed.name = 'GO Terms',
                                                      bypass.if.completed = FALSE, forced.parameters = NULL, out.dir = NULL, out.file.prefix = ""){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  # if (length(Startup.Data$topGO.Expression.Profiles___Custom.GO.terms) == 0) return()
  
  if (!Pars$Create.Excel.results){
    warning('Create.Excel.results parameter is turned OFF. "Summarize.clusterProfiler.Expression.Profiles.Custom" will not run')
    return()
  }
  
  if(is.null(GLM.models)){
    GLM.models = Pars$Models.to.Test
  }
  
  eval(parse(text = sprintf('PValue.thr = Pars$%s.clusterProfiler.Expression.Profiles___Max.gene.PValue.for.summary.across.models', database)))
  eval(parse(text = sprintf('logCPM.thr = Pars$%s.clusterProfiler.Expression.Profiles___Min.gene.logCPM.for.summary.across.models', database)))
  eval(parse(text = sprintf('max.Genes = Pars$%s.clusterProfiler.Expression.Profiles___Maximal.genes.in.terms.to.visualize', database)))
  
  file.names = c()
  for (GLM.model in GLM.models){
    GLM.model = Delete.prefix.in.model(GLM.model)
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    f = sprintf("%s/~ %s, results/%s - custom terms DE/%s, %s-DE info, logCPM.gt.%g, P.lt.%g.tsv",
                Pars$results.dir, Analysis.name, printed.name,
                Analysis.name, database, logCPM.thr, PValue.thr)
    f = gsub("/","\\", f, fixed = TRUE)
    f = Verify.path(f)
    file.names = append(file.names, f)
  }
  if (sum.mod(!file.exists(file.names))>0){
    warning(sprintf('Total %d files with custom GO-terms expression profiles do not exist: %s',sum.mod(!file.exists(file.names)),toString(file.names[!file.exists(file.names)])))
  }
  
  if(is.null(out.dir)) out.dir = Pars$results.dir
  
  file.names = file.names[file.exists(file.names)]
  if(length(file.names) == 0){
    cat('No files found\n')
    return(invisible())
  }
  
  CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces no --database %s ',Pars$python.bin,Pars$suppl.data.dir,database)
  CL = gsub('/','\\',CL,fixed = TRUE)
  out.file.name = sprintf('%s/%sSummary %s DE profiles.xlsx', out.dir, out.file.prefix, database)
  CL = paste(c(CL,sprintf('--out-excel "%s"',Verify.path(out.file.name)),'--in',sprintf('"%s"',file.names),
               '--maximal-genes-count', as.character(max.Genes)),collapse = ' ')
  system(CL)
  if (!file.exists(Verify.path(out.file.name)))    stop(sprintf('Excel worksheet generation FAILED'))
  
  CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces yes --database %s ',Pars$python.bin,Pars$suppl.data.dir,database)
  CL = gsub('/','\\',CL,fixed = TRUE)
  out.file.name = sprintf('%s/%sSummary %s DE profiles (with spaces).xlsx', out.dir, out.file.prefix, database)
  CL = paste(c(CL,sprintf('--out-excel "%s"',Verify.path(out.file.name)),'--in',sprintf('"%s"',file.names),
               '--maximal-genes-count', as.character(max.Genes)),collapse = ' ')
  system(CL)
  if (!file.exists(Verify.path(out.file.name)))    stop(sprintf('Excel worksheet generation FAILED'))
}


Summarize.GO.Expression.Profiles.Custom = function(Startup.Data, GLM.models = NULL, bypass.if.completed = FALSE, forced.parameters = NULL,
                                                   out.dir = NULL, out.file.prefix = ''){
  Summarize.clusterProfiler.Expression.Profiles.Custom(Startup.Data = Startup.Data, GLM.models = GLM.models,
                                                       database = 'GO', printed.name = 'GO Terms', out.dir = out.dir, out.file.prefix = out.file.prefix,
                                                       bypass.if.completed = bypass.if.completed, forced.parameters = forced.parameters)
}

Summarize.KEGG.Expression.Profiles.Custom = function(Startup.Data, GLM.models = NULL, bypass.if.completed = FALSE, forced.parameters = NULL,
                                                     out.dir = NULL, out.file.prefix = ""){
  Summarize.clusterProfiler.Expression.Profiles.Custom(Startup.Data = Startup.Data, GLM.models = GLM.models,
                                                       database = 'KEGG', printed.name = 'KEGG Pathways', out.dir = out.dir, out.file.prefix = out.file.prefix,
                                                       bypass.if.completed = bypass.if.completed, forced.parameters = forced.parameters)
}

Summarize.Reactome.Expression.Profiles.Custom = function(Startup.Data, GLM.models = NULL, bypass.if.completed = FALSE, forced.parameters = NULL,
                                                         out.dir = NULL, out.file.prefix = ''){
  Summarize.clusterProfiler.Expression.Profiles.Custom(Startup.Data = Startup.Data, GLM.models = GLM.models,
                                                       database = 'Reactome', printed.name = 'Reactome Pathways', out.dir = out.dir, out.file.prefix = out.file.prefix,
                                                       bypass.if.completed = bypass.if.completed, forced.parameters = forced.parameters)
}



Summarize.KEGG.DE.info = function(Startup.Data, GLM.models = NULL, bypass.if.completed = FALSE, forced.parameters = NULL,
                                  out.dir = NULL, out.file.prefix = ''){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(GLM.models)) GLM.models = Startup.Data$Pars$Models.to.Test
  
  if(is.null(out.dir)) out.dir = Pars$results.dir
  
  ## Integrating
  #setwd(Pars$results.dir)
  all.Analysis.names = sapply(GLM.models, Cleanup.Model.Name)
  # all.Analysis.names  = gsub("\\*","(m)",gsub("\\:","(d)", GLM.models))
  all.GLM.results.dir = sprintf("%s/~ %s, results",Pars$results.dir,all.Analysis.names)
  all.PA.working.dirs = sprintf("%s/KEGG Pathways - Visualization",all.GLM.results.dir)
  all.src.files = sprintf('%s/%s, KEGG pathways DE summary.tsv',all.PA.working.dirs, all.Analysis.names)
  
  #KEGG.de.summary.dir = sprintf("%s/%sKEGG DE summary", out.dir, out.file.prefix)
  #dir.create(KEGG.de.summary.dir, showWarnings = FALSE)
  
  existing.all.src.files = all.src.files[file.exists(all.src.files)]
  if(!all(file.exists(all.src.files))){
    cat(sprintf('[Summarize.KEGG.DE.info] %d out of %d KEGG DE summary files do not exist\n', sum.mod(!file.exists(all.src.files)), length(all.src.files)))
  } else {
    cat(sprintf('[Summarize.KEGG.DE.info] All %d KEGG DE summary files are present\n', length(all.src.files)))
  }
  
  
  #setwd(KEGG.de.summary.dir)
  CL=sprintf('%s "%s\\Integrate.KEGG.pathway.summaries.py" --input-file-names "%s" --out-tsv-prefix KEGG.de.summary_ --out-excel "%s\\%sKEGG DE summary.xlsx" --bar-max-value-override 100',
             Pars$python.bin, Pars$suppl.data.dir, paste(existing.all.src.files, collapse = ';'), out.dir, out.file.prefix)
  system(CL)
  
  #file.copy(from = 'KEGG.DE.summary.xlsx', to = sprintf('%s/KEGG.DE.summary.xlsx',Pars$results.dir),overwrite = TRUE)
  CL=sprintf('%s "%s\\Integrate.KEGG.pathway.summaries.py" --input-file-names "%s" --out-tsv-prefix KEGG.de.summary_ --out-excel "%s\\%sKEGG DE summary combined.xlsx" --bar-max-value-override 100 --combine-node-count-and-perc yes',
             Pars$python.bin, Pars$suppl.data.dir, paste(existing.all.src.files, collapse = ';'), out.dir, out.file.prefix)
  system(CL)
  #file.copy(from = 'KEGG.DE.summary.combined.xlsx', to = sprintf('%s/KEGG.DE.summary.combined.xlsx',Pars$results.dir),overwrite = TRUE)
  #setwd(Pars$results.dir)
  
  #"KEGG.DE.summary.to.Excel.py"
}



Process.ext.DE.data = function(Startup.Data, ext.DE.data.file, custom.analysis.name = NULL, bypass.if.completed = F, forced.parameters = NULL){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  setwd(Pars$results.dir)
  Analysis.Data = new("RTransAnalysisData")
  if(!is.null(custom.analysis.name)){
    Analysis.Data$GLM.model = custom.analysis.name
    Analysis.name = custom.analysis.name
  } else {
    Analysis.Data$GLM.model = ext.DE.data.file
    Analysis.name = ext.DE.data.file
  }
  Analysis.Data$Analysis.name = Analysis.name
  
  current.GLM.results.dir = sprintf("%s/~ %s, results",Pars$results.dir,Analysis.name)
  Analysis.Data$current.GLM.results.dir = current.GLM.results.dir
  dir.create(current.GLM.results.dir, showWarnings = FALSE)
  setwd(current.GLM.results.dir)
  
  Analysis.Data$step.hashmd5 = digest(ext.DE.data.file)
  
  ResTable.gene.part = Startup.Data$ext.DE.data[[ext.DE.data.file]][,c('gene','score','logfc','pvalue','fdr','logcpm')]
  colnames(ResTable.gene.part) = c('Gene name','Score','logFC','PValue','FDR','logCPM')
  rownames(ResTable.gene.part) = Startup.Data$ext.DE.data[[ext.DE.data.file]][,'ensembl_gene_id']

  Mart.table = Startup.Data$General.maRt.table[rownames(ResTable.gene.part),c("external_gene_name","gene_biotype","description")]
  mart.absent = !(rownames(ResTable.gene.part) %in% rownames(Startup.Data$General.maRt.table))
  for (n in 1:length(mart.absent)){
    if (mart.absent[n]){
      rownames(Mart.table)[n] = rownames(ResTable.gene.part)[n]
    }
  }
  
  #sum(grepl("^NA",rownames(Mart.table)))
  
  if (Pars$Use.Info.Table){
    Mart.table = cbind(Mart.table,data.frame(Startup.Data$Info.table[rownames(ResTable.gene.part),"RefSeq_Summary"]))
    colnames(Mart.table) = c("Gene name","Biotype","description","RefSeq_Summary")
  }  else{
    colnames(Mart.table) = c("Gene name","Biotype","description")
  }

  ResTable.gene.part = cbind(Mart.table,ResTable.gene.part[,-1])

  output.file.name = sprintf("%s_external_data.txt",Analysis.Data$Analysis.name)
  write.table.mod(ResTable.gene.part, file = output.file.name,sep='\t',na = '',col.names = colnames(ResTable.gene.part))
  
  if (Pars$Create.Excel.results){
    CL = sprintf('%s "%s/DE.results.to.Excel.py"',Pars$python.bin,Pars$suppl.data.dir)
    CL = gsub('/','\\',CL,fixed = T)
    CL = paste(c(CL,sprintf('"%s, DE results.xlsx"',Analysis.Data$Analysis.name),sprintf('"%s"',c(output.file.name))),collapse = ' ')
    system(CL)
    if (!file.exists(sprintf('%s, DE results.xlsx',Analysis.Data$Analysis.name))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (DE results) for ~%s FAILED',Analysis.Data$Analysis.name)) }
    
  }
  
  Analysis.Data$ResTable.gene.part = ResTable.gene.part
  saveRDS(ResTable.gene.part,file='DE.results.db')

  rm(ResTable.gene.part)

  Analysis.Data$GLM.analysis.performed = T
  cat('\nProcessing completed.\n')
  return(Analysis.Data)
}

#OrgDb = org.Mm.eg.db
#GO_DATA <- clusterProfiler::get_GO_data(OrgDb, ont = 'BP', keyType = 'ENSEMBL')

# KEGG_DATA <- prepare_KEGG(species, "KEGG", keyType)
# Reactome_DATA <- get_Reactome_DATA(organism)

clusterProfiler.Expression.Profiles = function(Startup.Data, DB.entries, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                              database = 'KEGG',
                              printed.name = 'KEGG Pathways',
                              # GO.type = 'BP',
                              profile.type = 'custom',
                              bypass.if.completed = FALSE,
                              forced.parameters = NULL, out.dir = NULL,
                              forced.Analysis.name = NULL, stats.file.name = NULL, sparklines.axis.limits = 2){
  
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  

  # if(database == 'GO'){
  #   database.plus = sprintf('%s.%s', database, GO.type)
  # } else {  database.plus = database  }

  # if(is.null(out.dir) != is.null(forced.DB.entries))  stop('clusterProfiler.Expression.Profiles.Custom: both forced.terms and out.dir parameters should be both NULL or both have a value')

  if(is.null(out.dir)){
    wd = sprintf("%s/%s - custom terms DE", Analysis.Data$current.GLM.results.dir, printed.name)
  } else {
    wd = out.dir
  }
  dir.create(wd, showWarnings = FALSE)
  setwd(wd)
  
  DB.entries = DB.entries[!duplicated(DB.entries)]

  if (length(DB.entries) == 0){
    message(sprintf('No %s to create expression profiles', printed.name))
    return(invisible())
  }
  
  if(database == 'GO'){
    Genes.per.DB.entry = Get.Gene.Ensembl.IDs__for__GO.term.IDs(DB.entries, Startup.Data, forced.parameters = Pars, remove.notfound.terms = TRUE)
  } else if(database == 'KEGG'){
    # Startup.Data$KEGG_DATA <- prepare_KEGG(species = Pars$Species, keyType = 'kegg')
    Genes.per.DB.entry = Get.Gene.Ensembl.IDs__for__KEGG.pathway.IDs(DB.entries, Startup.Data, forced.parameters = Pars, remove.notfound.pathways = TRUE)
  } else if(database == 'Reactome'){
    Genes.per.DB.entry = Get.Gene.Ensembl.IDs__for__Reactome.pathway.IDs(DB.entries, Startup.Data, forced.parameters = Pars, remove.notfound.pathways = TRUE)
  } else {
    stop('Unknown database')
  }
  
  not.found.DB.entries = DB.entries[!(DB.entries %in% names(Genes.per.DB.entry))]
  DB.entries = DB.entries[DB.entries %in% names(Genes.per.DB.entry)]
  
  if(length(not.found.DB.entries) > 10){
    cat(sprintf('\nThe following %s were not found among "%s -> genes" conversion database: %s.... (total %d entries). These %s will be excluded\n', printed.name,
                printed.name, paste(not.found.DB.entries[1:10], collapse = ', '), length(not.found.DB.entries), printed.name))
  } else if(length(not.found.DB.entries) > 0){
    cat(sprintf('\nThe following %s were not found among "%s -> genes" conversion database: %s (total %d entries). These %s will be excluded\n', printed.name,
                printed.name, paste(not.found.DB.entries, collapse = ', '), length(not.found.DB.entries), printed.name))
  }
  
  if (length(DB.entries) == 0){
    message(sprintf('No %s to create expression profiles', printed.name))
    return(invisible())
  }
  
  min.p.value_by.DB.entry.id.sel.UP = rep(1,length(DB.entries))
  names(min.p.value_by.DB.entry.id.sel.UP) = DB.entries

  min.FDR_by.DB.entry.id.sel.UP = rep(1,length(DB.entries))
  names(min.FDR_by.DB.entry.id.sel.UP) = DB.entries
  
  min.p.value_by.DB.entry.id.sel.DOWN = rep(1,length(DB.entries))
  names(min.p.value_by.DB.entry.id.sel.DOWN) = DB.entries
  
  min.FDR_by.DB.entry.id.sel.DOWN = data.matrix(rep(1,length(DB.entries)))
  names(min.FDR_by.DB.entry.id.sel.DOWN) = DB.entries
  
  ## looking up if there is GO.terms.GSEA.info table with p, FDR for GO terms. It is produced with GO.Enrichment() method
  if (is.null(stats.file.name)){
    # tmp = sprintf('%s/%s - classic Enrichment Analysis/%s, GO terms GSEA stats.tsv', Analysis.Data$current.GLM.results.dir, printed.name, Analysis.Data$Analysis.name)
    # if(file.exists(tmp)) stats.file.name = tmp
  # }
    Enrich.info.split.list = vector(mode = 'list')
    enrich.test.type = 'classic'
    for(enrich.test.type in c('classic', 'trends')){
      Enrich.info.split__is.found = FALSE
      if(database == 'GO'){
        all.stats = c('min.p.value_by.term.id.UP','min.FDR_by.term.id.UP','min.p.value_by.term.id.DOWN','min.FDR_by.term.id.DOWN',
                      'enrichment.score_by.term.id.final',
                      'enrichment.score_by.term.id.final.UP','enrichment.score_by.term.id.final.DOWN')
        Enrich.info.split =  data.frame(array(dim=c(0,length(all.stats))))
        colnames(Enrich.info.split) = all.stats
        for (GO.type in c('BP', 'MF', 'CC')){
          file.name = sprintf("%s/GO Terms (%s) - %s Enrichment Analysis/%s, GO.%s %s Enrich. stats.db",
                              Analysis.Data$current.GLM.results.dir, GO.type, enrich.test.type, Analysis.name, GO.type, enrich.test.type)
          if(file.exists(file.name)){
            Enrich.info.split = rbind(Enrich.info.split, readRDS(file = file.name))
          }
        }
        if(dim(Enrich.info.split)[1] > 0)  Enrich.info.split__is.found = TRUE
      } else {
        if (database == 'KEGG'){
          file.name = sprintf('%s/KEGG Pathways - %s Enrichment Analysis/%s, KEGG %s Enrich. stats.db',
                              Analysis.Data$current.GLM.results.dir, enrich.test.type, Analysis.name, enrich.test.type)
        } else if (database == 'Reactome'){
          file.name = sprintf('%s/Reactome Pathways - %s Enrichment Analysis/%s, Reactome %s Enrich. stats.db',
                              Analysis.Data$current.GLM.results.dir, enrich.test.type, Analysis.name, enrich.test.type)
        }
        if(file.exists(file.name)){
          Enrich.info.split = readRDS(file = file.name)
          Enrich.info.split__is.found = TRUE
        }
      }
      if(!Enrich.info.split__is.found){
        Enrich.info.split =  data.frame(array(dim=c(0,length(all.stats))))
        msg = sprintf('\n Results of %s %s enrichemnt are not found for analysis "~ %s" in the directory "%s"\n', enrich.test.type, printed.name, Analysis.name, dirname(file.name))
        cat(msg)
        # warning(msg)
      } else {
        Enrich.info.split.list[[enrich.test.type]] = Enrich.info.split
      }
    }

    not.found.DB.entries = c()
    if (('classic' %in% names(Enrich.info.split.list)) && ('trends' %in% names(Enrich.info.split.list))){
      e = DB.entries[1]
      for(e in DB.entries){
        if( (e %in% rownames(Enrich.info.split.list[['classic']])) & (e %in% rownames(Enrich.info.split.list[['trends']])) ){
          min.p.value_by.DB.entry.id.sel.UP[e] = min(Enrich.info.split.list[['classic']][e,'min.p.value_by.term.id.UP'],
                                                     Enrich.info.split.list[['trends']][e,'min.p.value_by.term.id.UP'])
          min.FDR_by.DB.entry.id.sel.UP[e] = min(Enrich.info.split.list[['classic']][e,'min.FDR_by.term.id.UP'],
                                                     Enrich.info.split.list[['trends']][e,'min.FDR_by.term.id.UP'])
          min.p.value_by.DB.entry.id.sel.DOWN[e] = min(Enrich.info.split.list[['classic']][e,'min.p.value_by.term.id.DOWN'],
                                                     Enrich.info.split.list[['trends']][e,'min.p.value_by.term.id.DOWN'])
          min.FDR_by.DB.entry.id.sel.DOWN[e] = min(Enrich.info.split.list[['classic']][e,'min.FDR_by.term.id.DOWN'],
                                                 Enrich.info.split.list[['trends']][e,'min.FDR_by.term.id.DOWN'])
          
        } else if( (e %in% rownames(Enrich.info.split.list[['classic']])) & !(e %in% rownames(Enrich.info.split.list[['trends']])) ){
          min.p.value_by.DB.entry.id.sel.UP[e] = Enrich.info.split.list[['classic']][e,'min.p.value_by.term.id.UP']
          min.FDR_by.DB.entry.id.sel.UP[e] = Enrich.info.split.list[['classic']][e,'min.FDR_by.term.id.UP']
          min.p.value_by.DB.entry.id.sel.DOWN[e] = Enrich.info.split.list[['classic']][e,'min.p.value_by.term.id.DOWN']
          min.FDR_by.DB.entry.id.sel.DOWN[e] = Enrich.info.split.list[['classic']][e,'min.FDR_by.term.id.DOWN']
          
        } else if( !(e %in% rownames(Enrich.info.split.list[['classic']])) & (e %in% rownames(Enrich.info.split.list[['trends']])) ){
          min.p.value_by.DB.entry.id.sel.UP[e] = Enrich.info.split.list[['trends']][e,'min.p.value_by.term.id.UP']
          min.FDR_by.DB.entry.id.sel.UP[e] = Enrich.info.split.list[['trends']][e,'min.FDR_by.term.id.UP']
          min.p.value_by.DB.entry.id.sel.DOWN[e] = Enrich.info.split.list[['trends']][e,'min.p.value_by.term.id.DOWN']
          min.FDR_by.DB.entry.id.sel.DOWN[e] = Enrich.info.split.list[['trends']][e,'min.FDR_by.term.id.DOWN']
          
        } else {
          not.found.DB.entries = c(not.found.DB.entries, e)
        }
      }
    } else if(('classic' %in% names(Enrich.info.split.list)) && !('trends' %in% names(Enrich.info.split.list))){
      e = DB.entries[1]
      for(e in DB.entries){
        if( e %in% rownames(Enrich.info.split.list[['classic']])){
          min.p.value_by.DB.entry.id.sel.UP[e] = Enrich.info.split.list[['classic']][e,'min.p.value_by.term.id.UP']
          min.FDR_by.DB.entry.id.sel.UP[e] = Enrich.info.split.list[['classic']][e,'min.FDR_by.term.id.UP']
          min.p.value_by.DB.entry.id.sel.DOWN[e] = Enrich.info.split.list[['classic']][e,'min.p.value_by.term.id.DOWN']
          min.FDR_by.DB.entry.id.sel.DOWN[e] = Enrich.info.split.list[['classic']][e,'min.FDR_by.term.id.DOWN']
        } else {
          not.found.DB.entries = c(not.found.DB.entries, e)
        }
      }
    } else if(!('classic' %in% names(Enrich.info.split.list)) && ('trends' %in% names(Enrich.info.split.list))){
      e = DB.entries[1]
      for(e in DB.entries){
        if( e %in% rownames(Enrich.info.split.list[['trends']])){
          min.p.value_by.DB.entry.id.sel.UP[e] = Enrich.info.split.list[['trends']][e,'min.p.value_by.term.id.UP']
          min.FDR_by.DB.entry.id.sel.UP[e] = Enrich.info.split.list[['trends']][e,'min.FDR_by.term.id.UP']
          min.p.value_by.DB.entry.id.sel.DOWN[e] = Enrich.info.split.list[['trends']][e,'min.p.value_by.term.id.DOWN']
          min.FDR_by.DB.entry.id.sel.DOWN[e] = Enrich.info.split.list[['trends']][e,'min.FDR_by.term.id.DOWN']
        } else {
          not.found.DB.entries = c(not.found.DB.entries, e)
        }
      }
    } else {
      msg = sprintf("\n%s Enrichment analysis seems to be not performed. P-values, FDRs will be set as 1\n", printed.name)
      cat(msg)
      warning(msg)
      # Enrich.info.split = NULL
    }
    
    if(length(not.found.DB.entries) > 10){
      cat(sprintf('\nThe following %s were not found among enrichment tests results: %s.... (total %d entries out of %d). p-values and FDR will be set as 1\n', printed.name,
                  paste(not.found.DB.entries[1:10], collapse = ', '), length(not.found.DB.entries), length(DB.entries)))
    } else if(length(not.found.DB.entries) > 0){
      cat(sprintf('\nThe following %s were not found among enrichment tests results: %s (total %d entries out of %d). p-values and FDR will be set as 1\n', printed.name,
                  paste(not.found.DB.entries, collapse = ', '), length(not.found.DB.entries), length(DB.entries)))
    }
  } else {
    Enrich.info = read.table(file = stats.file.name,header = T,sep = '\t')
    
    min.p.value_by.DB.entry.id.sel.UP = Enrich.info[DB.entries,'min.p.value_by.term.id.UP']
    min.p.value_by.DB.entry.id.sel.UP[is.na(min.p.value_by.term.id.sel.UP)] = 1
    
    min.FDR_by.DB.entry.id.sel.UP = Enrich.info[DB.entries,'min.FDR_by.term.id.UP']
    min.FDR_by.DB.entry.id.sel.UP[is.na(min.FDR_by.term.id.sel.UP)] = 1
    
    min.p.value_by.DB.entry.id.sel.DOWN = Enrich.info[DB.entries,'min.p.value_by.term.id.DOWN']
    min.p.value_by.DB.entry.id.sel.DOWN[is.na(min.p.value_by.term.id.sel.DOWN)] = 1
    
    min.FDR_by.DB.entry.id.sel.DOWN = Enrich.info[DB.entries,'min.FDR_by.term.id.DOWN']
    min.FDR_by.DB.entry.id.sel.DOWN[is.na(min.FDR_by.term.id.sel.DOWN)] = 1
    
  }
  
  if(database == 'GO'){
    DB.entries.names = Startup.Data$GO.full.descriptions[DB.entries, c('name')]
    DB.entries.descriptions = Startup.Data$GO.full.descriptions[DB.entries, c('full description')]
  } else if (database == 'KEGG'){
    DB.entries.names = Startup.Data$KEGG_DATA$PATHID2NAME[DB.entries]
    DB.entries.descriptions = NULL
  } else if (database == 'Reactome'){
    DB.entries.names = unlist(as.list(reactome.db::reactomePATHID2NAME[DB.entries]))
    DB.entries.descriptions = NULL
  } else {
    stop('Unknown database')
  }
  
  eval(parse(text = sprintf('PValue.thr.list = Pars$%s.clusterProfiler.Expression.Profiles___Max.gene.PValue.list', database)))
  eval(parse(text = sprintf('logCPM.thr.list = Pars$%s.clusterProfiler.Expression.Profiles___Min.gene.logCPM.list', database)))
  
  output.file.names = c()

  PValue.thr = PValue.thr.list[1]
  logCPM.thr = logCPM.thr.list[1]
  for(PValue.thr in PValue.thr.list){
    for(logCPM.thr in logCPM.thr.list){
      all.LogFCs.per.DB.entry = Get.LogFCs.for.Ensembl.IDs.vector(Analysis.Data, Genes.per.DB.entry, PValue.thr = PValue.thr, logCPM.thr = logCPM.thr, sort.logFCs = TRUE)
      all.GT.logFCs.list = unlist(lapply(all.LogFCs.per.DB.entry, function(x){ paste(sprintf("%.2f", x), collapse = '\t')}))
      #all.GT.genes.counts.list = sapply(Genes.per.DB.entry[DB.entries], length)
      all.GT.genes.counts.list = sapply(all.LogFCs.per.DB.entry, length)
      
      DE.profile.table = as.data.frame(cbind(DB.entries,
                              DB.entries.names,
                              min.p.value_by.DB.entry.id.sel.UP,
                              min.FDR_by.DB.entry.id.sel.UP,
                              min.p.value_by.DB.entry.id.sel.DOWN,
                              min.FDR_by.DB.entry.id.sel.DOWN,
                              all.GT.genes.counts.list,
                              all.GT.logFCs.list))
      
      DE.profile.table = DE.profile.table[which(sapply(all.LogFCs.per.DB.entry, length) > 0),]
      file.name = Verify.path(sprintf('%s, %s-DE info, logCPM.gt.%g, P.lt.%g.tsv', Analysis.name, database, logCPM.thr, PValue.thr))
      write.table.mod(DE.profile.table,file = file.name ,sep='\t',quote = FALSE)
      output.file.names = c(output.file.names, file.name)
    }
  }
  
  if (Pars$Create.Excel.results){
    CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces no --database %s',Pars$python.bin,Pars$suppl.data.dir,database)
    CL = gsub('/','\\',CL,fixed = T)
    base.res.file.name = sprintf('%s, %s %s DE profiles.xlsx', Analysis.name, printed.name, profile.type)
    CL = paste(c(CL,sprintf('--out-excel "%s"', Verify.path(base.res.file.name)),'--in',sprintf('"%s"',output.file.names),
                 '--maximal-genes-count',as.character(Pars$topGO.Expression.Profiles___Maximal.genes.in.terms.to.visualize),
                 sprintf('--sparklines-axis-limit %g',sparklines.axis.limits)),
               collapse = ' ')
    system(CL)
    if (!file.exists(Verify.path(base.res.file.name))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (%s-centric DE profiles) for ~%s FAILED', database, Analysis.name))
    } else {
      file.copy(Verify.path(base.res.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, base.res.file.name)))
    }
    
    CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces yes --database %s',Pars$python.bin,Pars$suppl.data.dir,database)
    CL = gsub('/','\\',CL,fixed = T)
    base.res.file.name = sprintf('%s, %s %s DE profiles (with spaces).xlsx', Analysis.name, printed.name, profile.type)
    CL = paste(c(CL,sprintf('--out-excel "%s"', Verify.path(base.res.file.name)),'--in',sprintf('"%s"',output.file.names),
                 '--maximal-genes-count',as.character(Pars$topGO.Expression.Profiles___Maximal.genes.in.terms.to.visualize),
                 sprintf('--sparklines-axis-limit %g',sparklines.axis.limits)),collapse = ' ')
    system(CL)
    if (!file.exists(Verify.path(base.res.file.name))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (%s-centric DE profiles) for ~%s FAILED', database, Analysis.name))
    } else {
      file.copy(Verify.path(base.res.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, base.res.file.name)))
    }
  }
  #saveRDS()
  if (is.null(out.dir)) setwd(Analysis.Data$current.GLM.results.dir)
  # load
  return(invisible(Analysis.Data))
}

GO.Expression.Profiles.Custom =  function(
  Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
  bypass.if.completed = FALSE, forced.parameters = NULL, out.dir = NULL,
  forced.Analysis.name = NULL, stats.file.name = NULL, sparklines.axis.limits = 2){
    
    DB.entries = Startup.Data$GO.clusterProfiler.Expression.Profiles___Custom.DB.entries
    clusterProfiler.Expression.Profiles(Startup.Data, DB.entries, Analysis.Data = Analysis.Data, GLM.model = GLM.model,
                                        GLM.working.dir = GLM.working.dir, Analysis.Data.RDS.file = Analysis.Data.RDS.file,
                                        database = 'GO', printed.name = 'GO Terms', profile.type = 'custom',
                                        bypass.if.completed = bypass.if.completed,
                                        forced.parameters = forced.parameters, out.dir = out.dir,
                                        forced.Analysis.name = forced.Analysis.name, stats.file.name = stats.file.name,
                                        sparklines.axis.limits = sparklines.axis.limits)
}


KEGG.Expression.Profiles.Custom =  function(
  Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
  bypass.if.completed = FALSE, forced.parameters = NULL, out.dir = NULL,
  forced.Analysis.name = NULL, stats.file.name = NULL, sparklines.axis.limits = 2){
  
  DB.entries = Startup.Data$KEGG.clusterProfiler.Expression.Profiles___Custom.DB.entries
  clusterProfiler.Expression.Profiles(Startup.Data, DB.entries, Analysis.Data = Analysis.Data, GLM.model = GLM.model,
                                      GLM.working.dir = GLM.working.dir, Analysis.Data.RDS.file = Analysis.Data.RDS.file,
                                      database = 'KEGG', printed.name = 'KEGG Pathways', profile.type = 'custom',
                                      bypass.if.completed = bypass.if.completed,
                                      forced.parameters = forced.parameters, out.dir = out.dir,
                                      forced.Analysis.name = forced.Analysis.name, stats.file.name = stats.file.name,
                                      sparklines.axis.limits = sparklines.axis.limits)
}


Reactome.Expression.Profiles.Custom =  function(
  Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
  bypass.if.completed = FALSE, forced.parameters = NULL, out.dir = NULL,
  forced.Analysis.name = NULL, stats.file.name = NULL, sparklines.axis.limits = 2){
  
  DB.entries = Startup.Data$Reactome.clusterProfiler.Expression.Profiles___Custom.DB.entries
  clusterProfiler.Expression.Profiles(Startup.Data, DB.entries, Analysis.Data = Analysis.Data, GLM.model = GLM.model,
                                      GLM.working.dir = GLM.working.dir, Analysis.Data.RDS.file = Analysis.Data.RDS.file,
                                      database = 'Reactome', printed.name = 'Reactome Pathways', profile.type = 'custom',
                                      bypass.if.completed = bypass.if.completed,
                                      forced.parameters = forced.parameters, out.dir = out.dir,
                                      forced.Analysis.name = forced.Analysis.name, stats.file.name = stats.file.name,
                                      sparklines.axis.limits = sparklines.axis.limits)
}


inDetails = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                     forced.parameters = NULL, GO.terms.inDetails = NULL, out.dir = NULL, cluster.samples = FALSE,
                     add.biotypes = c('lincRNA + antisense_RNA + macro_lincRNA'), sort.by.predictor = NULL){
  
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.db', GLM.working.dir, Analysis.name)
    if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  if(is.null(GO.terms.inDetails))  GO.terms.inDetails = Startup.Data$GO.terms.inDetails
  
  if(length(GO.terms.inDetails) == 0){
    cat('\nNo inDetails groups are declared - see "GO terms (inDetails)" sheet\n')
    return(invisible())
  }
  
  if(is.null(out.dir)){
    main.wd = sprintf("%s/inDetails", Analysis.Data$current.GLM.results.dir)
  } else {
    main.wd = out.dir
  }
  dir.create(main.wd, showWarnings = FALSE)
  setwd(main.wd)
  
  all.GO.terms = unlist(GO.terms.inDetails)
  all.GO.terms = all.GO.terms[!duplicated(all.GO.terms)]
  pre__genes.per.GO.term = Get.Gene.Ensembl.IDs__for__GO.term.IDs(term.list = all.GO.terms, Startup.Data = Startup.Data, forced.parameters = forced.parameters)
  
  all__group.Types = vector(mode = 'list')
  all__genes.of.interest = vector(mode = 'list')
  for(group.Name in names(GO.terms.inDetails)){
    genes.per.GO.term = pre__genes.per.GO.term[GO.terms.inDetails[[group.Name]]]
    genes.of.interest = unlist(genes.per.GO.term)
    if(length(genes.of.interest) == 0){
      cat(sprintf('No genes were found for group "%s"\n',group.Name))
      next
    } else if(length(genes.of.interest) < 2){
      cat(sprintf('Group "%s" is bypassed since only 1 gene is found within\n',group.Name))
      next
    }
    genes.of.interest = genes.of.interest[!duplicated(genes.of.interest)]
    all__genes.of.interest[[group.Name]] = genes.of.interest
    all__group.Types[[group.Name]] = 'GO'
  }
  
  biotypes = add.biotypes[1]
  for(biotypes in add.biotypes){
    group.Name = biotypes
    #biotypes= 'lincRNA + antisense_RNA + macro_lincRNA'
    biotypes %<>%  gsub(" ", "", .) %>% strsplit(., split = '+', fixed = TRUE) %>% unlist
    keep = Analysis.Data$ResTable.gene.part$Biotype %in% biotypes
    if(sum(keep) == 0){
      cat(sprintf('No genes were found for additional biotype-group "%s"\n',group.Name))
      next
    } else if(sum(keep) < 2){
      cat(sprintf('Group "%s" is bypassed since only 1 gene is found within\n',group.Name))
      next
    }
    all__genes.of.interest[[group.Name]] = rownames(Analysis.Data$ResTable.gene.part)[keep]
    all__group.Types[[group.Name]] = 'by_biotype'
  }

  
  all.output.file.names = c()
  
  group.Name = names(all__group.Types)[1]
  for(group.Name in names(all__group.Types)){
    genes.of.interest = all__genes.of.interest[[group.Name]]
    if(all__group.Types[[group.Name]] == 'GO'){
      cat(sprintf('inDetails analysis of group "%s" (%d terms, %d genes)...\n', group.Name, length(GO.terms.inDetails[[group.Name]]), length(genes.of.interest)))
    } else {
      cat(sprintf('inDetails analysis of group "%s" (%d genes)...\n', group.Name, length(GO.terms.inDetails[[group.Name]]), length(genes.of.interest)))
    }
    if(length(genes.of.interest) < 2){
      cat('The analysis is bypassed since the genes count is less than 2\n')
      next
    }
    
    group.Name_s = gsub("[^[:alnum:] ]", "", group.Name)
    dir.create(group.Name_s, showWarnings = FALSE)
    setwd(group.Name_s)
    
    if('edgeR_d' %in% names(Analysis.Data)){
      schema = Analysis.Data$schema
      Current.Predictors.list = Analysis.Data$Current.Predictors.list
      d = Analysis.Data$edgeR_d
      if(sum(rownames(d$counts) %in% genes.of.interest) == 0){
        cat(sprintf('There are no genes that passed DE analysis in "%s" group\n', group.Name))
        next
      } else if(sum(rownames(d$counts) %in% genes.of.interest) == 1){
        temp <- d$counts[rownames(d$counts) %in% genes.of.interest,] %>% as.data.frame %>% t
        rownames(temp) = rownames(d$counts)[rownames(d$counts) %in% genes.of.interest]
        d$counts = temp
      } else {
        d$counts = d$counts[rownames(d$counts) %in% genes.of.interest,]
      }

      Main.Component = Analysis.Data$Main.Component
      Main.Predictor.name = Analysis.Data$Main.Predictor.name
      
      plot.subtitle = ''
      if (length(levels(factor(Main.Component))) == 2) plot.subtitle = sprintf('color indicates %s (quantity; blue-orange)',Main.Predictor.name)
      if (length(levels(factor(Main.Component))) == 3) plot.subtitle = sprintf('color indicates %s (quantity; blue-green-orange)',Main.Predictor.name)
      if (length(levels(factor(Main.Component))) >= 4) plot.subtitle = sprintf('color indicates %s (quantity; blue-green-orange gradient)',Main.Predictor.name)
      
      col.variety <- colorRampPalette(c("#0b00c4","#02babc","#25ad00","#cfcd00","#ff8400"))(512)
      Main.Component.mod = Main.Component - min(Main.Component)
      cl = col.variety[1 + as.integer(Main.Component.mod/max(Main.Component.mod)*511)]
      
      cpms = cpm(d)
      if(dim(d$counts)[1] > 2){
        png.mod(filename = sprintf('CPM density, normalized (%s).png',Pars$RNA.Seq.norm.method),
            units="in",width=8.23,height=6.3,pointsize=12,
            res=600)
        plot(density(log(cpms[,1])),col=cl[1],xlim=c(-4,10),main = sprintf('Log2(CPM) density, normalized (%s)',Pars$RNA.Seq.norm.method), sub = plot.subtitle)
        for (x in 2:dim(cpms)[2]){
          lines(density(log(cpms[,x])),col=cl[x])
        }
        dev.off()
      }
      
      
      col.variety<-colorRampPalette(c("#0b00c4","#02babc","#25ad00","#cfcd00","#ff8400"))(length(levels(factor(Main.Component))))
      
      plot.title = sprintf('MDS plot, %s',Main.Predictor.name)
      plot.subtitle = ''
      if (length(levels(factor(Main.Component))) == 2) plot.subtitle = sprintf('color indicates %s (blue-orange)',Main.Predictor.name)
      if (length(levels(factor(Main.Component))) == 3) plot.subtitle = sprintf('color indicates %s (blue-green-orange)',Main.Predictor.name)
      if (length(levels(factor(Main.Component))) >= 4) plot.subtitle = sprintf('color indicates %s (blue-green-orange gradient)',Main.Predictor.name)
      
      
      if(dim(d$counts)[2] > 2 && dim(d$counts)[1] > 1){
        tryCatch(expr = {
          plotMDS(d,dim.plot = c(1,2), col = col.variety[factor(Main.Component)])
          title(main = plot.title, sub = plot.subtitle)
          Create.MDS <<- TRUE
        }, error = function (err){
          warning('MDS plots creating failed')
          Create.MDS <<- FALSE
        })
      } else {
        Create.MDS <<- FALSE
        msg = sprintf('MDS plots will not be created bacause only %d samples and %d genes are present', dim(d$counts)[2], dim(d$counts)[1])
        cat(sprintf('\n%s\n', msg))
        #warning(msg)
      }
      
      
      if(Create.MDS){
        pdf(file = 'MDS.dim.1-2.pdf',
            width=18.23,height=16.3,pointsize=12)
        MDS.info = plotMDS(d,dim.plot = c(1,2), col = col.variety[factor(Main.Component)])
        two.dim.MDS.info = cbind(MDS.info$x,MDS.info$y)
        title(main = plot.title, sub = plot.subtitle)
        dev.off()
        
        if(dim(d)[2] > 3){
          pdf(file = 'MDS.dim.2-3.pdf',
              width=18.23,height=16.3,pointsize=12)
          MDS.info = plotMDS(d,dim.plot = c(2,3), col = col.variety[factor(Main.Component)])
          title(main = plot.title, sub = plot.subtitle)
          dev.off()
          three.dim.MDS.info = cbind(two.dim.MDS.info,MDS.info$y)
          colnames(three.dim.MDS.info) = c('MDS dim1','MDS dim2','MDS dim3')
          write.table.mod(x = three.dim.MDS.info,file = 'MDS info.tsv',sep = '\t')
        }
        
        png.mod(filename = 'MDS.dim.1-2.png',units="in",res=600,
            width=8.23,height=6.3,pointsize=12)
        plotMDS(d,dim.plot = c(1,2), col = col.variety[factor(Main.Component)])
        title(main = plot.title, sub = plot.subtitle)
        dev.off()
        
        if(dim(d)[2] > 3){
          png.mod(filename = 'MDS.dim.2-3.png',units="in",res=600,
              width=8.23,height=6.3,pointsize=12)
          a = plotMDS(d,dim.plot = c(2,3), col = col.variety[factor(Main.Component)])
          title(main = plot.title, sub = plot.subtitle)
          dev.off()
        }
      }
      
      ## Calculating distance matrices
      C.Create.Distance.matrices = Pars$Create.Distance.matrices
      if (Pars$Create.Distance.matrices == '{auto}'){
        C.Create.Distance.matrices = (dim(cpms)[2] <= 50)
      }
      
      if (C.Create.Distance.matrices){
        if (Pars$Distance.method %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")){
          dist.matrix = data.frame(dist(cpms,Pars$Distance.method))
        } else {
          N = dim(cpms)[2]
          dist.matrix = data.frame(array(dim = c(N,N)))
          colnames(dist.matrix) = colnames(cpms)
          rownames(dist.matrix) = colnames(cpms)
          
          #F = colnames(cpms)[1] # this is for debug only
          ready = 0
          for (F in colnames(cpms)){
            #S = colnames(cpms)[2]
            for (S in colnames(cpms)){
              if (Pars$Distance.method == 'canberra.weighted') dist.matrix[F,S] = Weighted.Canberra.dist(cpms[,F],cpms[,S],Pars$Min.CPM.to.include.in.dist,Pars$Canberra.weight.power)
              if (Pars$Distance.method == '1-cor') dist.matrix[F,S] = 1 - suppressWarnings(cor.test(cpms[,F],cpms[,S],method = 'spearman')$estimate)
              ready = ready+1
              cat(sprintf('\rCalculating dist matrices. Completed %.2f%%',ready/N/N*100))
              
              #Weighted.Canberra.dist(c(1,2,3,4,5,6,7,8,9,10),c(1,20,30,4,5,6,7,8,9,10),min.value.to.calc.dist = 0.99)
            }
          }
        }
        
        tmp=as.dist(dist.matrix)
        dend=hclust(tmp)
        dend=as.dendrogram(dend)
        
        dend.colors = col.variety[factor(Main.Component)]
        names(dend.colors) = schema$'Sample names'
        labels_colors(dend) = dend.colors[labels(dend)]
        png.mod(filename = sprintf('Clustering, dist as %s.png',Pars$Distance.method),units="in",res=600,
            width=8.23,height=6.3,pointsize=12)
        plot(dend,main = sprintf('%s, clustering dendrogram [%s]',GLM.model,Pars$Distance.method),
             sub = sprintf('distances are calculated with "%s" method',Pars$Distance.method))
        dev.off()
        ##hang=-1
        
        
        #Current.Predictors.list = c('Age','Genotype','Gender')
        info.block = data.frame(as.data.frame(schema)[colnames(cpms),Current.Predictors.list])
        colnames(info.block) = Current.Predictors.list
        info.block.mod = cbind(array(dim = c(dim(info.block)[2],dim(info.block)[2])),t(info.block))
        colnames(info.block.mod) = colnames(cbind(info.block,dist.matrix))
        class(info.block.mod) <- "numeric"
        write.table.mod(x=rbind(info.block.mod,
                            cbind(info.block,dist.matrix)),
                    file='distance matrix.tsv',sep="\t",na = "")
      }
      
      
    }
    
    ResTable.gene.part = Analysis.Data$ResTable.gene.part[rownames(Analysis.Data$ResTable.gene.part) %in% genes.of.interest,]
    ResTable.sorted.by.Score = rbind(Analysis.Data$ResTable.pred.part, ResTable.gene.part)
    output.file.name = sprintf("%s - %s_%s_combined.txt", group.Name_s, Analysis.name, Pars$DE.package)
    all.output.file.names = c(all.output.file.names, output.file.name)
    
    write.table.mod(ResTable.sorted.by.Score, file = output.file.name,sep='\t',na = '',col.names = colnames(ResTable.sorted.by.Score))
    if (Pars$Create.Excel.results){
      CL = sprintf('%s "%s/DE.results.to.Excel.py"', Pars$python.bin, Pars$suppl.data.dir)
      CL = gsub('/','\\',CL,fixed = T)
      CL = paste(c(CL, Verify.path(sprintf('"%s [%s], DE results.xlsx"', Analysis.name, group.Name_s)), sprintf('"%s"',Verify.path(output.file.name))),collapse = ' ')
      system(CL)
      if (!file.exists(Verify.path(sprintf('%s [%s], DE results.xlsx', Analysis.name, group.Name_s)))){
        #file.remove(output.file.names)
        message(sprintf('Excel worksheet generation (DE results) for ~%s FAILED', Analysis.name))
      }
    }
    
    Create.Heatmaps(Startup.Data, Analysis.Data = Analysis.Data, bypass.if.completed = FALSE,
                    limit.to.genes = genes.of.interest, out.dir = '.', max.Genes.counts.list = c(50, 100, 500, 2000, 5000),
                    cluster.samples = cluster.samples, sort.by.predictor = sort.by.predictor)
    setwd(main.wd)
  }
  
  # all.output.file.names = c()
  # for(name in names(all__group.Types)){
  #   if(length(all__genes.of.interest[[group.Name]]) < 2) next
  #   
  #   all.output.file.names = c(all.output.file.names,
  #        sprintf("%s/%s - %s_%s_combined.txt", gsub("[^[:alnum:] ]", "", name),
  #                gsub("[^[:alnum:] ]", "", name), Analysis.name, Pars$DE.package))
  # }
  
  
  all.group.Name_s = gsub("[^[:alnum:] ]", "", names(all__group.Types))
  all.output.file.names = sprintf("%s/%s - %s_%s_combined.txt", all.group.Name_s, all.group.Name_s, Analysis.name, Pars$DE.package)
  
  all.output.file.names = sapply(all.output.file.names,Verify.path)
  
  CL = sprintf('%s "%s/DE.results.to.Excel.py"',Pars$python.bin,Pars$suppl.data.dir)
  CL = gsub('/','\\',CL,fixed = T)
  out.file.name = sprintf('inDetails - %s, DE results.xlsx',Analysis.name)
  CL = paste(c(CL,sprintf('"%s"', Verify.path(out.file.name)),sprintf('"%s"', all.output.file.names)),collapse = ' ')
  system(CL)
  if (!file.exists(Verify.path(out.file.name))){
    #file.remove(output.file.names)
    stop(sprintf('Excel worksheet generation FAILED [non-joint], %s',out.file.name))
  }
  file.copy(Verify.path(out.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, out.file.name)), overwrite = TRUE)
  
  if (is.null(out.dir)) setwd(Analysis.Data$current.GLM.results.dir)
}


Collect.topGO.enrichment.result.files <- function(Startup.Data, GLM.models = NULL, bypass.if.completed = FALSE, forced.parameters = NULL,
                                                  out.dir = NULL, GO.types = c('BP','MF','CC')){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if (!Pars$Create.Excel.results){
    warning('Create.Excel.results parameter is turned OFF. "Collect.topGO.enrichment.result.files" will not run')
    return()
  }
  
  if(is.null(out.dir)){
    out.dir = sprintf('%s/topGO Enrichment test - collected results', Pars$results.dir)
  }
  dir.create(out.dir, showWarnings = FALSE)
  
  GO.type = GO.types[1]
  for(GO.type in GO.types){
    res.files = c()
    not.found.res.files = c()
    dest.files = c()
    GLM.model = GLM.models[1]
    for(GLM.model in GLM.models){
      Analysis.name  = Cleanup.Model.Name(GLM.model)
      current.GLM.results.dir = sprintf("%s/~ %s, results",Pars$results.dir, Analysis.name)
      file.base.name = sprintf('%s, GO Terms (%s) topGO enrichment.xlsx', Analysis.name, GO.type)
      res.file = sprintf('%s/%s', current.GLM.results.dir, file.base.name)
      
      if(file.exists(res.file)){
        res.files = c(res.files, res.file)
      } else if (file.exists(Verify.path(res.file))){
        res.files = c(res.files, Verify.path(res.file))
      } else {
        not.found.res.files = c(not.found.res.files, res.file)
        next
      }
      dest.file = Verify.path(sprintf('%s/%s', out.dir, file.base.name))
      dest.files = c(dest.files, dest.file)
      # out.dir
    }
    if(length(res.files) == 0){
      cat(sprintf('No topGO results are found for "%s" ontology type\n', GO.type))
      next
    }
    if(length(not.found.res.files) > 0){
      cat(sprintf('The following topGO (%s) results are not found: %s\n', GO.type, toString(not.found.res.files)))
    } else {
      cat(sprintf('All topGO (%s) results are present\n', GO.type))
    }
    file.copy(res.files, dest.files, overwrite = TRUE)
    #dest.files
  }
}



Collect.Enrich.results <- function(Startup.Data, GLM.models = NULL, forced.parameters = NULL, out.dir = NULL, database = 'GO', printed.name = 'GO terms classic enrichment',
                                  GO.types = c('BP','MF','CC'), file.name.skeleton = 'GO Terms (%s) - classic Enrichment Analysis/%s, GO Terms (%s) classic enrichment.xlsx',
                                  bypass.if.completed = FALSE){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if (!Pars$Create.Excel.results){
    warning('Create.Excel.results parameter is turned OFF. "Collect.Enrich.results" will not run')
    return()
  }
  
  if(is.null(out.dir)){
    out.dir = Pars$results.dir
  }
  dir.create(out.dir, showWarnings = FALSE)
  
  if(length(strsplit(file.name.skeleton, split = '/', fixed = TRUE)[[1]]) != 2){
    stop(sprintf('Collect.Enrich.results: incorrect skeleton %s', file.name.skeleton))
  }
  
  src.dir.name.skeleton = strsplit(file.name.skeleton, split = '/', fixed = TRUE)[[1]][1]
  src.file.name.skeleton = strsplit(file.name.skeleton, split = '/', fixed = TRUE)[[1]][2]
  
  if(database == 'GO' | database == 'topGO'){
    GO.type = GO.types[1]
    for(GO.type in GO.types){
      src.files = c()
      not.found.src.files = c()
      dest.files = c()
      GLM.model = GLM.models[1]
      for(GLM.model in GLM.models){
        Analysis.name  = Cleanup.Model.Name(GLM.model)
        current.GLM.results.dir = sprintf("%s/~ %s, results",Pars$results.dir, Analysis.name)
        #file.base.name = sprintf('%s, GO Terms (%s) topGO enrichment.xlsx', Analysis.name, GO.type)
        
        src.dir.name = sprintf(src.dir.name.skeleton, GO.type)
        src.file.name = sprintf(src.file.name.skeleton, Analysis.name, GO.type)
        
        src.file.full.name = sprintf('%s/%s/%s', current.GLM.results.dir, src.dir.name, src.file.name)

        if(file.exists(src.file.full.name)){
          src.files = c(src.files, src.file.full.name)
        } else if (file.exists(Verify.path(src.file.full.name))){
          src.files = c(src.files, Verify.path(src.file.full.name))
        } else {
          not.found.src.files = c(not.found.src.files, src.file.full.name)
          next
        }
        
        dest.file = Verify.path(sprintf('%s/~ %s', out.dir, src.file.name))
        dest.files = c(dest.files, dest.file)
        # out.dir
      }
      if(length(src.files) == 0){
        cat(sprintf('No %s results are found for "%s" ontology type\n', printed.name, GO.type))
        next
      }
      if(length(not.found.src.files) > 0){
        cat(sprintf('The following %s (%s) results are not found: %s\n', printed.name, GO.type, toString(not.found.src.files)))
      } else {
        cat(sprintf('All %s (%s) results are present\n', printed.name, GO.type))
      }
      file.copy(src.files, dest.files, overwrite = TRUE)
      #dest.files
    }
  } else {
    src.files = c()
    not.found.src.files = c()
    dest.files = c()
    GLM.model = GLM.models[1]
    for(GLM.model in GLM.models){
      Analysis.name  = Cleanup.Model.Name(GLM.model)
      current.GLM.results.dir = sprintf("%s/~ %s, results",Pars$results.dir, Analysis.name)
      #file.base.name = sprintf('%s, GO Terms (%s) topGO enrichment.xlsx', Analysis.name, GO.type)
      
      src.dir.name = src.dir.name.skeleton
      src.file.name = sprintf(src.file.name.skeleton, Analysis.name)
      
      src.file.full.name = sprintf('%s/%s/%s', current.GLM.results.dir, src.dir.name, src.file.name)
      
      if(file.exists(src.file.full.name)){
        src.files = c(src.files, src.file.full.name)
      } else if (file.exists(Verify.path(src.file.full.name))){
        src.files = c(src.files, Verify.path(src.file.full.name))
      } else {
        not.found.src.files = c(not.found.src.files, src.file.full.name)
        next
      }
      
      dest.file = Verify.path(sprintf('%s/~ %s', out.dir, src.file.name))
      dest.files = c(dest.files, dest.file)
      # out.dir
    }
    if(length(src.files) == 0){
      cat(sprintf('No %s results are found\n', printed.name))
    } else {
      if(length(not.found.src.files) > 0){
        cat(sprintf('The following %s results are not found: %s\n', printed.name, toString(not.found.src.files)))
      } else {
        cat(sprintf('All %s results are present\n', printed.name))
      }
      file.copy(src.files, dest.files, overwrite = TRUE)
    }
    #dest.files
  }
  return(invisible(NULL))
}

Summarize.results.by.Analyses.groups = function(Startup.Data, bypass.if.completed = FALSE, forced.parameters = NULL, out.dir = NULL){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if (!Pars$Create.Excel.results){
    warning('Create.Excel.results parameter is turned OFF. "Summarize.results.by.Analyses.groups" will not run')
    return()
  }
  
  if(is.null(out.dir)){
    out.dir = Pars$results.dir
  }
  dir.create(out.dir, showWarnings = FALSE)
  
  group.name = names(Pars$Analyses.groups)[1]
  for(group.name in names(Pars$Analyses.groups)){
    GLM.models = Pars$Analyses.groups[[group.name]]
    current.group.out.dir = sprintf('%s/Essential results - %s', out.dir, group.name)
    dir.create(current.group.out.dir, showWarnings = FALSE)
    
    Collect.Enrich.results(Startup.Data, GLM.models, forced.parameters = Pars, out.dir = current.group.out.dir,
                           database = 'topGO', printed.name = 'GO terms - topGO - classic enrichment', GO.types = c('BP','MF','CC'),
                           file.name.skeleton = 'GO Terms (%s) - topGO - Enrichment Analysis/%s, GO Terms (%s) topGO enrichment.xlsx')
    Collect.Enrich.results(Startup.Data, GLM.models, forced.parameters = Pars, out.dir = current.group.out.dir,
                           database = 'GO', printed.name = 'GO terms classic enrichment', GO.types = c('BP','MF','CC'),
                           file.name.skeleton = 'GO Terms (%s) - classic Enrichment Analysis/%s, GO Terms (%s) classic enrichment.xlsx')
    Collect.Enrich.results(Startup.Data, GLM.models, forced.parameters = Pars, out.dir = current.group.out.dir,
                           database = 'KEGG', printed.name = 'KEGG Pathways classic enrichment',
                           file.name.skeleton = 'KEGG Pathways - classic Enrichment Analysis/%s, KEGG Pathways classic enrichment.xlsx')
    Collect.Enrich.results(Startup.Data, GLM.models, forced.parameters = Pars, out.dir = current.group.out.dir,
                           database = 'Reactome', printed.name = 'Reactome Pathways classic enrichment',
                           file.name.skeleton = 'Reactome Pathways - classic Enrichment Analysis/%s, Reactome Pathways classic enrichment.xlsx')
    
    Summarize.GLM.results(Startup.Data, GLM.models = GLM.models, forced.parameters = Pars,
                          out.file.name = sprintf('%s/%s - Summary expression info.xlsx', current.group.out.dir, group.name))
    inDetails.dir = sprintf('%s/%s - inDetails DE summary', current.group.out.dir, group.name)
    dir.create(inDetails.dir, showWarnings = FALSE)
    Summarize.inDetails.GLM.results(Startup.Data, GLM.models = GLM.models, out.dir = inDetails.dir)
    
    Summarize.KEGG.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.models, out.dir = current.group.out.dir,
                                              out.file.prefix = sprintf('%s - ', group.name))
    Summarize.KEGG.DE.info(Startup.Data, GLM.models = GLM.models, out.dir = current.group.out.dir,
                           out.file.prefix = sprintf('%s - ', group.name))
    Summarize.Reactome.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.models, out.dir = current.group.out.dir,
                                                  out.file.prefix = sprintf('%s - ', group.name))
    #Summarize.GO.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.models)
    Summarize.topGO.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.models, out.dir = current.group.out.dir,
                                               out.file.prefix = sprintf('%s - ', group.name))
    
  }
}
