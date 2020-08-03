# source('/mnt/raid/illumina/geo/progs/RTrans/RTrans_base/RTrans.methods 1.51.R')
setClass("RTransParameters",
         representation("list")
)

setClass("RTransStartupData",
         representation("list")
)

RTransAnalysisData.prototype = vector(mode="list")
RTransAnalysisData.prototype$GLM.analysis.performed = FALSE
RTransAnalysisData.prototype$KEGG.classic.Enrichment.performed = FALSE
RTransAnalysisData.prototype$GO.Enrichment.performed = FALSE
RTransAnalysisData.prototype$Heatmaps.created = FALSE
RTransAnalysisData.prototype$Pathway.Visualization.performed = FALSE
RTransAnalysisData.prototype$GO.Expression.Profiles.Created = FALSE
RTransAnalysisData.prototype$Enriched.Pathways.Classic = c()
RTransAnalysisData.prototype$Enriched.Pathways.Trended = c()
RTransAnalysisData.prototype$Enrich.Results.Summary.by.DB = vector(mode="list")

setClass("RTransAnalysisData",
         representation("list"),
         prototype = prototype(RTransAnalysisData.prototype)
)

RTrans.WGCNA.AnalysisData.prototype = vector(mode="list")
RTrans.WGCNA.AnalysisData.prototype$WGCNA.performed = FALSE
setClass("RTrans.WGCNA.AnalysisData",
         representation("list"),
         prototype = prototype(RTrans.WGCNA.AnalysisData.prototype)
)

TimeStamp = function(){
  gsub(':', '-', format(Sys.time(), "%a %b %d %X %Y"), fixed = TRUE)
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

sum.mod = function(x){
  x=x[!is.na(x)]; return(length(x[x]))
}

rel_sd_trim = function(x) sd_trim(x,trim = trimming,const = F)/mean(x)+0.001

# src_array = c(10,0,0,0,0,0,0,0,0,0,0,0,0,0,50,0,0,0,10)
# desired_length = 10

DeDup = function(x) x[!duplicated(x)]

DeDup.na.rm = function(x){
  x = x[!is.na(x)]
  x[!duplicated(x)]
}

Get.LogFC.col.name = function(ResTable.gene.part__col.names, preferred_columns = c('TBA LogFC', 'delta LogFC')){
  if(!is.null(preferred_columns)){
    if(any(preferred_columns %in% ResTable.gene.part__col.names)){
      LogFC.col.name = preferred_columns[which(preferred_columns %in% ResTable.gene.part__col.names)[1]]
    } else {
      LogFC.col.name = Get.LogFC.col.name__internal(ResTable.gene.part__col.names)
      # cat(sprintf('[Get.LogFC.col.name]   Column %s is not found in ResTable.gene.part__col.names. Using LogFC variant "%s"...\n', preferred_columns, LogFC.col.name))
    }
  } else {
    LogFC.col.name = et.LogFC.col.name__internal(ResTable.gene.part__col.names)
  }
  return(LogFC.col.name)
}

Get.LogFC.col.name__internal = function(ResTable.gene.part__col.names){
  if('TBA LogFC' %in% ResTable.gene.part__col.names){
    LogFC.col.name = 'TBA LogFC'
  } else if('trimmed LogFC' %in% ResTable.gene.part__col.names){
    LogFC.col.name = 'trimmed LogFC'
  } else if('trimmed logFC' %in% ResTable.gene.part__col.names){
    LogFC.col.name = 'trimmed logFC'
  } else if('LogFC' %in% ResTable.gene.part__col.names){
    LogFC.col.name = 'LogFC'
  } else if('logFC' %in% ResTable.gene.part__col.names){
    LogFC.col.name = 'logFC'
  } else {
    LogFC.col.name = 'not_found_LogFC_colname'
  }
  return(LogFC.col.name)
}


Get.PValue.col.name = function(Main.test, ResTable.gene.part__col.names = NULL){
  if(!is.null(ResTable.gene.part__col.names)){
    if('TBA p-value' %in% ResTable.gene.part__col.names) return('TBA p-value')
  }

  if(Main.test == 'ET'){
    PValue.col = 'p (ex. test)'
  } else if(Main.test == 'LR'){
    PValue.col = 'p (LR test)'
  } else if(Main.test == 'QL'){
    PValue.col = 'p (QLF test)'
  } else if(Main.test == 'LR paired'){
    PValue.col = 'p (LR test paired)'
  } else if(Main.test == 'QL paired'){
    PValue.col = 'p (QLF test paired)'
  } else if(Main.test == 't-test (FC)'){
    PValue.col = 'p (t-test, FC)'
  } else if(Main.test == 'Mann-Wh. (FC)'){
    PValue.col = 'p (Mann-Wh., FC)'
  } else if(Main.test == 'Lin.Mod. (FC)'){
    PValue.col = 'p (Lin.Mod., FC)'
  } else {
    PValue.col = 'not_found_p_colname'
  }

  if(!is.null(ResTable.gene.part__col.names)){
    if(PValue.col %in% ResTable.gene.part__col.names){
      PValue.col = PValue.col
    } else if ('PValue' %in% ResTable.gene.part__col.names){
      PValue.col = 'PValue'
    } else if ('pvalue' %in% ResTable.gene.part__col.names){
      PValue.col = 'pvalue'
    } else if ('p-value' %in% ResTable.gene.part__col.names){
      PValue.col = 'p-value'
    } else if ('P' %in% ResTable.gene.part__col.names){
      PValue.col = 'P'
    } else if ('p' %in% ResTable.gene.part__col.names){
      PValue.col = 'p'
    } else {
      stop('Cannot find P-value column name')
    }
  }
  return(PValue.col)
}

Get.NP.PValue.col.name = function(ResTable.gene.part__col.names){
  if ('p (Mann-Wh., FC)' %in% ResTable.gene.part__col.names){
    PValue.col = 'p (Mann-Wh., FC)'
  } else if ('p (Mann-Wh.)' %in% ResTable.gene.part__col.names){
    PValue.col = 'p (Mann-Wh.)'
  } else if ('p (Wilcoxon paired)' %in% ResTable.gene.part__col.names){
    PValue.col = 'p (Wilcoxon paired)'
  } else if ('p (Mann-Wh.)' %in% ResTable.gene.part__col.names){
    PValue.col = 'p (Mann-Wh.)'
  } else if ('p (Krusk-W.)' %in% ResTable.gene.part__col.names){
    PValue.col = 'p (Krusk-W.)'
  } else {
    cat('Cannot find NP P-value column name\n')
    PValue.col = NA
  }
  return(PValue.col)
}


Get.FDR.col.name = function(Main.test, ResTable.gene.part__col.names = NULL, look.for.paired.test.res = FALSE){
  if(!is.null(ResTable.gene.part__col.names)){
    if('TBA p-value' %in% ResTable.gene.part__col.names) return('TBA p-value')
  }

  if(Main.test == 'ET'){
    FDR.col = 'FDR (ex. test)'
  } else if(Main.test == 'LR paired' | (look.for.paired.test.res & Main.test == 'LR')){
    FDR.col = 'FDR (LR test paired)'
  } else if(Main.test == 'QL paired' | (look.for.paired.test.res & Main.test == 'QL')){
    FDR.col = 'FDR (QLF test paired)'
  } else if(Main.test == 'LR'){
    FDR.col = 'FDR (LR test)'
  } else if(Main.test == 'QL'){
    FDR.col = 'FDR (QLF test)'
  } else if(Main.test == 't-test (FC)'){
    FDR.col = 'FDR (t-test, FC)'
  } else if(Main.test == 'Mann-Wh. (FC)'){
    FDR.col = 'FDR (Mann-Wh., FC)'
  } else if(Main.test == 'Lin.Mod. (FC)'){
    FDR.col = 'FDR (Lin.Mod., FC)'
  } else {
    FDR.col = 'not_found_FDR_colname'
  }

  if(!is.null(ResTable.gene.part__col.names)){
    if(FDR.col %in% ResTable.gene.part__col.names){
      FDR.col = FDR.col
    } else if ('FDR' %in% ResTable.gene.part__col.names){
      FDR.col = 'FDR'
    } else if ('fdr' %in% ResTable.gene.part__col.names){
      FDR.col = 'fdr'
    } else {
      stop('Cannot find FDR column name')
    }
  }
  return(FDR.col)
}

Get.LogCPM.col.name = function(ResTable.gene.part__col.names = NULL){
  if ('LogCPM' %in% ResTable.gene.part__col.names){
    LogCPM.col = 'LogCPM'
  } else if ('logCPM' %in% ResTable.gene.part__col.names){
    LogCPM.col = 'logCPM'
  } else {
    stop('Cannot find LogCPM column name')
  }
  return(LogCPM.col)
}

Convert..Sample.Cols..to.Character = function(schema.like){
  if('File names' %in% colnames(schema.like))
    schema.like[,'File names'] = as.character(schema.like[,'File names'])
  if('Sample names' %in% colnames(schema.like))
    schema.like[,'Sample names'] = as.character(schema.like[,'Sample names'])
  return(schema.like)
}

Read.phys.chem.data <- function(Excel.file, sheet.name = NULL, wd = NULL, convert.to.numeric..if.possible = TRUE,
  convert..to.numeric..if.levels.count..is.greater.than = 4){

  if(!is.null(wd)) Excel.file = sprintf('%s/%s', wd, Excel.file)
  if(!file.exists(Excel.file))  stop(sprintf('File "%s" is not found', Excel.file))
  if(is.null(sheet.name))   sheet.name = excel_sheets(Excel.file)[1]
  if(!(sheet.name %in% excel_sheets(Excel.file))) stop(sprintf('sheet "%s" is not found in Excel workbook "%s"', sheet.name, Excel.file))

  phys.chem.data = suppressMessages(readxl::read_excel(Excel.file, sheet=sheet.name, col_names = TRUE))
  phys.chem.data = data.frame(phys.chem.data, check.names = FALSE)
  phys.chem.data = Convert..Sample.Cols..to.Character(phys.chem.data)

  if('File names' %in% colnames(phys.chem.data)){
    phys.chem.data = phys.chem.data[!startsWith(phys.chem.data[, 'File names'], '#'), ]
    if(any(dim(phys.chem.data) == 0))  return(NULL)
  }

  if('Sample names' %in% colnames(phys.chem.data)){
    phys.chem.data = phys.chem.data[!startsWith(phys.chem.data[, 'Sample names'], '#'), ]
    if(any(dim(phys.chem.data) == 0))  return(NULL)
  }

  for(col in colnames(phys.chem.data)){
    x = phys.chem.data[,col]
    x[x %in% c('NA', 'N/A', 'n/a')] = NA
    #MR.StartupData$phys.chem.data[,col] = x
    if(convert.to.numeric..if.possible){
      #which.is..non.NA = which(!is.na(x))
      rest.x.is.numeric = suppressWarnings(all(!is.na(as.numeric(x[!is.na(x)]))))
      if(rest.x.is.numeric){
        x.levels.count = suppressWarnings(length(DeDup.na.rm(x)))
        if(x.levels.count > convert..to.numeric..if.levels.count..is.greater.than){
          x = suppressWarnings(as.numeric(x))
        }
      }
    }
    phys.chem.data[,col] = x
  }


  commented.lines__v1 = sapply(phys.chem.data[,1], FUN = function (x) substring(x, first=1, last=1)) %in% '#'
  commented.lines__v2 = gsub(' ', '', phys.chem.data[,1], fixed = TRUE) %in% ''
  phys.chem.data = subset(phys.chem.data, subset = !(commented.lines__v1 | commented.lines__v2 | is.na(phys.chem.data[,1])))

  if('Sample names' %in% colnames(phys.chem.data)){
    rownames(phys.chem.data) = phys.chem.data[ , which(colnames(phys.chem.data) %in% 'Sample names') ]
    phys.chem.data = subset(phys.chem.data, select = !(colnames(phys.chem.data) %in% 'Sample names'))
  } else {
    rownames(phys.chem.data) = phys.chem.data[,1]
    phys.chem.data = subset(phys.chem.data, select = c(TRUE, rep(FALSE, dim(phys.chem.data)[2] - 1)))
  }
  
  if(any(dim(phys.chem.data) == 0))  return(NULL)
  
  phys.chem.data = subset(phys.chem.data, select = !(colnames(phys.chem.data) %in% 'File names'))
  if(any(dim(phys.chem.data) == 0))  return(NULL)

  phys.chem.data = subset(phys.chem.data, select = !startsWith(colnames(phys.chem.data), '#'))
  if(any(dim(phys.chem.data) == 0))  return(NULL)

  rn = rownames(phys.chem.data)
  if(any(startsWith(rn, '#')))
    phys.chem.data = subset(phys.chem.data, subset = rn[!startsWith(rn, '#')])
  if(any(dim(phys.chem.data) == 0))  return(NULL)

  return(phys.chem.data)
}


lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


ResizeArray <- function (src_array, desired_length){
  initial_length = length(src_array)
  if(initial_length == desired_length){
    return(src_array)
  }
  
  res_array = rep(0.0, desired_length)
  part_size =  initial_length / desired_length
  part_number = 0
  for(part_number in 0:(desired_length-1)){
    start_coord = part_number*part_size
    end_coord = start_coord + part_size
    values = c()
    weights = c()
    
    start_fragment_weight = floor(min(start_coord + 1, end_coord)) - start_coord
    if(start_fragment_weight > 0){
      values = c(values, src_array[floor(start_coord) + 1])
      weights = c(weights, start_fragment_weight)
    }
    
    end_fragment_weight = end_coord - floor(max(end_coord, start_coord+1))
    if(end_fragment_weight > 0){
      values = c(values, src_array[min(initial_length-1, floor(end_coord)) + 1])
      weights = c(weights, end_fragment_weight)
    }
    
    
    if(start_fragment_weight <= 0 && end_fragment_weight <= 0){
      if(floor(start_coord) != floor(end_coord)){
        print('smth strange...')
      }
      values = c(values, src_array[floor(start_coord) + 1])
      weights = c(weights, 1)
    }
    
    for(x in (floor(start_coord)+1):(floor(end_coord) - 1)){
      values = c(values, src_array[x + 1])
      weights = c(weights, 1)
    }
    
    
    # acc = 0
    # for(x in 1:(length(values))){
    #   acc = acc + values[x] * weights[x]
    # }
    # 
    # res_array[part_number + 1] = acc/sum(weights)
    res_array[part_number + 1] = weighted.mean(values, weights)
  }
  return(res_array)
}


Smoothed_Trimmed_Mean <- function (array, low.trim = 0.1, high.trim = 0.1){
  if(low.trim + high.trim > 1){
    stop(sprintf('[Smoothed_Trimmed_Mean]  Very high trim, %g and %g', low.trim, high.trim))
  }
  
  from_percentile = low.trim*100
  to_percentile = 100 - high.trim*100
  
  weights = ResizeArray(c(
    rep(0, from_percentile),
    rep(1, to_percentile - from_percentile),
    rep(0 , 100 - to_percentile)),
    length(array))
  
  #sorted_array = array[order(array)]
  # return sum([weights[n]*sorted_array[n] for n in range(len(sorted_array))])/sum(weights)
  return(weighted.mean(array[order(array)], weights))
  #return numpy.average(sorted_array, weights=weights)
}

# Smoothed_Trimmed_StDev <- function (array, from_percentile, to_percentile,ForcedAverage = None){
#   if(low.trim + high.trim > 1){
#     stop(sprintf('[Smoothed_Trimmed_StDev]  Very high trim, %g and %g', low.trim, high.trim))
#   }
#   
#   from_percentile = low.trim*100
#   to_percentile = 100 - high.trim*100
#   
#   
#   weights = ResizeArray(c(
#     rep(0, from_percentile),
#     rep(1, to_percentile - from_percentile),
#     rep(0 , 100 - to_percentile)),
#     length(array))
#   
#   sorted_array = sorted(array)
#     # return sum([weights[n]*sorted_array[n] for n in range(len(sorted_array))])/sum(weights)
#     # print(sorted_array)
#     # print(weights)
#   return weighted_std(sorted_array, weights=weights,ForcedAverage = ForcedAverage)
# }

Smoothed_Trimmed_LogFC <- function(vector1, vector2, trim.low = 0.1, trim.high = 0.1, const.add = 0){
  return(log2((Smoothed_Trimmed_Mean(vector2, trim.low, trim.high) + const.add) / (Smoothed_Trimmed_Mean(vector1, trim.low, trim.high) + const.add)))
}

Smoothed_Trimmed_LogFC <- function(vector1, vector2, trim.low = 0.1, trim.high = 0.1, const.add = 0){
  #return(log2((Smoothed_Trimmed_Mean(vector2, trim.low, trim.high) + const.add) / (Smoothed_Trimmed_Mean(vector1, trim.low, trim.high) + const.add)))
  return(log2((mean(vector2, trim = max(trim.low, trim.high)) + const.add) / (mean(vector1, trim = max(trim.low, trim.high) + const.add))))
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

png.mod <- function(filename, width = 480, height = 480, pointsize = 12, res = NA, units="in", ...){
  png(filename = Verify.path(filename), width = width, height = height, pointsize = pointsize, res = res, units = units, ...)
}

write.table.mod <- function(x, file = "", append = FALSE, quote = TRUE, sep = " ", ...){
  write.table(x = x, file = Verify.path(file), append = append, quote = quote, sep = sep, ...)
}

Choose.Separator = function(csv.file.name, verbose = TRUE){
  con = file(csv.file.name, "r")
  first_line = readLines(con, n=1)
  close(con)
  possible.separators = c('\t', ',', ';')
  possible.separators..occurence = c(0, 0, 0)
  for(n_sep in 1:length(possible.separators)){
    possible.separators..occurence[n_sep] = str_count(first_line, possible.separators[n_sep])
  }
  sep = possible.separators[which.max(possible.separators..occurence)]
  if(verbose)
    cat(sprintf('Separator "%s" will be used for file "%s"\n', toString(sep), csv.file.name))
}


Check.availability.of.paired.test = function(sample.names, verbose = TRUE, verbose.OK = TRUE){
  #sample.names = c('17N', '17T', '19N', '19T', '15N', '21T')
  
  if(any(nchar(sample.names) <= 1)){
    if(verbose)  cat(sprintf('paired test is not available because not all sample names have > 1 character length: %s\n',
                             toString(all.sample.names[nchar(all.sample.names) <= 1])))
    return(FALSE)
  }
  
  suffixes = substr(sample.names, nchar(sample.names), nchar(sample.names))
  suffixes..dedup = sort(DeDup(suffixes))
  
  if(length(suffixes..dedup) != 2){
    if(verbose)  cat(sprintf('paired test is not available because there %d suffix variants: %s\nOnly 2 suffix variants (e.g. N, T) are permitted\n',
                             length(suffixes..dedup), toString(suffixes..dedup)))
    return(FALSE)
  }
  
  if(!all(grepl('[[:alpha:]]', suffixes..dedup))){
    if(verbose)  cat(sprintf('paired test is not available. Suffixes must be alphabetical, e.g. N, T. There %d non-alphabetical suffix variant(s): %s\n',
                             sum.mod(!grepl('[[:alpha:]]', suffixes..dedup)), toString(suffixes..dedup[!grepl('[[:alpha:]]', suffixes..dedup)])))
    return(FALSE)
  }
  
  samples.N.base.names = sample.names %>% .[endsWith(., suffixes..dedup[1])] %>% substr(., 1, nchar(.) - 1) %>% sort
  samples.T.base.names = sample.names %>% .[endsWith(., suffixes..dedup[2])] %>% substr(., 1, nchar(.) - 1) %>% sort
  
  if(length(samples.N.base.names) != length(samples.T.base.names)){
    if(verbose)  cat(sprintf('paired test is not available. The number of %s- and %s- samples differs: %d and %d\n',
                             suffixes..dedup[1], suffixes..dedup[2],
                             length(samples.N.base.names), length(samples.T.base.names)))
    return(FALSE)
  }
  
  diff.samples = c(setdiff(samples.T.base.names, samples.N.base.names), setdiff(samples.N.base.names, samples.T.base.names))
  
  if(length(diff.samples) > 0){
    if(verbose)  cat(sprintf('paired test is not available. Some samples are present in only one subset (%s or %s): %s\n',
                             suffixes..dedup[1], suffixes..dedup[2],
                             toString(diff.samples)))
    return(FALSE)
  }

  if(length(samples.N.base.names) <= 1){
    if(verbose)  cat(sprintf('paired test is not available. Too few samples (%d samples, %d pairs)\n',
                             2*length(samples.N.base.names), length(samples.N.base.names)))
    return(FALSE)
  }
  
  if(verbose.OK)  cat(sprintf('paired test is available. Sample names are self-consistent\n'))
  
  return(TRUE)
}

Reorder.sample.names..for.paired.test = function(sample.names, return.order = TRUE){
  if(!Check.availability.of.paired.test(sample.names, verbose.OK = FALSE))
    stop('Reorder.sample.names..for.paired.test cannot run because Check.availability.of.paired.test has failed')
  
  suffixes = substr(sample.names, nchar(sample.names), nchar(sample.names))
  suffixes..dedup = sort(DeDup(suffixes))
  
  samples..base.names = sample.names %>% .[endsWith(., suffixes..dedup[1])] %>% substr(., 1, nchar(.) - 1) %>% sort
  reordered..samples = c(sprintf("%s%s", samples..base.names, suffixes..dedup[1]), sprintf("%s%s", samples..base.names, suffixes..dedup[2]))
  
  if(return.order){
    return(match(reordered..samples, sample.names))
  } else {
    return(reordered..samples)
  }
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

install.packages.std <- function(pkg.list, force.re.install=FALSE){
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
  #for (pkg in pkg.list) eval(parse(text = sprintf("suppressPackageStartupMessages(library(\"%s\"))",pkg)))
}

load.packages.std.bio <- function(pkg.list){
  for (pkg in pkg.list) eval(parse(text = sprintf("suppressPackageStartupMessages(library(\"%s\"))",pkg)))
}

Check..R.version..gt.3.5 = function(){
  return(as.numeric(as.character(version$major)) > 3 | as.numeric(as.character(strsplit(version$minor, split = '.', fixed = TRUE)[[1]][1])) > 4)
}

install.packages.bio <- function(pkg.list, force.re.install=FALSE){
  if (any(!(pkg.list %in% rownames(installed.packages()))) | force.re.install){
    # if(Check..R.version..gt.3.5()){
      if(!('BiocManager' %in% rownames(installed.packages())))
        install.packages('BiocManager')
      if(force.re.install){
        pkg.list..to.install = pkg.list
      } else {
        pkg.list..to.install = pkg.list[!(pkg.list %in% rownames(installed.packages()))]
      }
      BiocManager::install(pkg.list..to.install)
    # } else {      
    #   source("https://bioconductor.org/biocLite.R")
    #   for (pkg in pkg.list){
    #     if(force.re.install){
    #       eval(parse(text = sprintf("biocLite(\"%s\", dep = TRUE)", pkg)))
    #     } else {
    #       if (!(pkg %in% rownames(installed.packages()))) eval(parse(text = sprintf("biocLite(\"%s\")", pkg)))
    #     }
    #   }
    # }
  }
  #for (pkg in pkg.list) eval(parse(text = sprintf("suppressPackageStartupMessages(library(\"%s\"))",pkg)))
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
          final.pred = 
          final.values = as.data.frame(predictor.values[,colnames(res.array)[keep]])
          colnames(final.values) = colnames(res.array)[keep]
          rownames(final.values) = rownames(predictor.values)
          return(final.values)
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


Split.Predictor.Name = function(predictor.name){
  splitted = strsplit(predictor.name ,' vs ')[[1]]
  if(length(splitted) != 2){
    splitted = strsplit(predictor.name ,' vs\\. ')[[1]]
    if(length(splitted) != 2){
      return(NULL)
    }
  }
  incount = str_count(splitted[2], '\\[')
  outcount = str_count(splitted[2], '\\]')
  if(incount == 1 & outcount == 1){
    # if()
    splitted2 = c(splitted[1], strsplit(splitted[2], '\\[')[[1]])
    if(str_count(splitted2[3], '\\]') == 1){
      splitted2[3] = gsub(']', '', splitted2[3], fixed = TRUE)
    } else  return(c(trimws(splitted), NA))
    return(trimws(splitted2))
  } else return(c(trimws(splitted), NA))
}

Get..Predictor.values..to..group.names = function(Predictor.values, GLM.model){
  Predictor.values.reordered = Predictor.values[order(Predictor.values)]
  Predictor.levels = Predictor.values.reordered[!duplicated(Predictor.values.reordered)]
  
  Predictor.values..to..group.names = vector(mode = 'list')
  for(ln in 1:length(Predictor.levels)){
    Predictor.values..to..group.names[[toString(Predictor.levels[ln])]] = Predictor.levels[ln]
  }

  if(!is.null(Split.Predictor.Name(GLM.model))){    
    splitted.GLM.model = Split.Predictor.Name(GLM.model)
    if(length(Predictor.levels) == length(splitted.GLM.model) - 1){
      Predictor.values..to..group.names = vector(mode = 'list')
      for(ln in 1:length(Predictor.levels)){
        if(!is.na(tail(splitted.GLM.model, n = 1))){
          Predictor.values..to..group.names[[toString(Predictor.levels[ln])]] = sprintf('%s [%s]', splitted.GLM.model[length(splitted.GLM.model) - ln], splitted.GLM.model[length(splitted.GLM.model)])
        } else {
          Predictor.values..to..group.names[[toString(Predictor.levels[ln])]] = splitted.GLM.model[length(splitted.GLM.model) - ln]
        }
      }
    }
  }
  return(Predictor.values..to..group.names)
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
  Default.Parameters[['Perform.DE.analysis']] = TRUE
  Default.Parameters[['Create.summary.reports']] = TRUE
  Default.Parameters[['Create.Venn_Euler.diagrams']] = TRUE
  Default.Parameters[['Create.heatmaps']] = TRUE
  Default.Parameters[['Perform.inDetails.analysis']] = TRUE
  Default.Parameters[['Create.heatmaps.for.inDetails.analyses']] = TRUE
  Default.Parameters[['Perform.KEGG.pathway.visualization']] = TRUE

  Default.Parameters[['Perform.classic.GO.enrichment']] = TRUE
  Default.Parameters[['Perform.trends.GO.enrichment']] = TRUE
  Default.Parameters[['GO.enrichment.namespaces']] = c('BP', 'CC', 'MF')
  Default.Parameters[['Perform.topGO.enrichment']] = TRUE
  Default.Parameters[['topGO.enrichment.namespaces']] = 'BP'
  Default.Parameters[['Perform.KEGG.classic.enrichment']] = TRUE
  Default.Parameters[['Perform.KEGG.trends.enrichment']] = TRUE
  Default.Parameters[['Perform.Reactome.classic.enrichment']] = TRUE
  Default.Parameters[['Perform.Reactome.trends.enrichment']] = TRUE
  Default.Parameters[['Perform.Disease.Ontology.classic.enrichment']] = TRUE
  Default.Parameters[['Perform.Disease.Ontology.trends.enrichment']] = TRUE
  Default.Parameters[['Perform.WikiPathways.classic.enrichment']] = TRUE
  # Default.Parameters[['Perform.WikiPathways.trends.enrichment']] = TRUE
  Default.Parameters[['Perform.DisGeNET.classic.enrichment']] = TRUE
  Default.Parameters[['Perform.DisGeNET.trends.enrichment']] = TRUE
  Default.Parameters[['Perform.Network.of.Cancer.Genes.classic.enrichment']] = TRUE
  Default.Parameters[['Perform.Network.of.Cancer.Genes.trends.enrichment']] = TRUE

  Default.Parameters[['Create.GO_centric.expression.profiles']] = TRUE
  Default.Parameters[['Create.GO_centric.expression.profiles.for.enriched.terms']] = FALSE
  Default.Parameters[['Create.topGO_centric.expression.profiles']] = TRUE
  Default.Parameters[['topGO_centric.expression.profiles__namespace']] = 'BP'
  Default.Parameters[['Create.topGO_centric.expression.profiles.for.enriched.terms']] = FALSE
  Default.Parameters[['topGO_centric.expression.profiles.for.enriched.terms__namespace']] = 'BP'
  Default.Parameters[['Create.KEGG_centric.expression.profiles']] = TRUE
  Default.Parameters[['Create.Reactome_centric.expression.profiles']] = TRUE

  Default.Parameters[['Stop.if.an.error.occurs']] = FALSE

  ###########################################################################

  Default.Parameters[['DE.package']] = 'edgeR'
  Default.Parameters[['Models.to.Test']] = '{auto}'
  Default.Parameters[['Species']] = 'hsa'
  Default.Parameters[['add.MDS.dims.as.predictors']] = FALSE
  Default.Parameters[['Normalize.predictors']] = TRUE
  Default.Parameters[['results.dir']] = '{script_dir}' 
  Default.Parameters[['counts.dir']] = ''
  Default.Parameters[['counts.suffix']] = ''
  Default.Parameters[['read.counts.table']] = ''
  Default.Parameters[['suppl.data.dir']] = '{script_dir}/../RTrans_main_scripts'
  Default.Parameters[['Create.Excel.results']] = TRUE
  Default.Parameters[['Max.genes.in.Excel.results']] = '{auto}'
  Default.Parameters[['Cluster.samples.in.Excel.results']] = TRUE
  Default.Parameters[['re.Cluster.samples.in.Details.Excel.results']] = TRUE
  Default.Parameters[['Clustering.method']] = 'complete'
  Default.Parameters[['Min.group.size.to..cluster']] = 4
  Default.Parameters[['Keep.original..in.Details..gene.order']] = FALSE
  Default.Parameters[['classic.paired.test']] = '{auto}'
  Default.Parameters[['FC.associations.paired.test']] = '{auto}'
  Default.Parameters[['paired.LogFCs.calculating..CPM.const.add']] = 0.25
  
  Default.Parameters[['use.N.top.Genes.to.Cluster.samples']] = 500
  Default.Parameters[['Create.CPM.profiles.in.Excel.results']] = FALSE
  Default.Parameters[['Create.rel.LogCPM.profiles.in.Excel.results']] = TRUE
  Default.Parameters[['Include.Kruskal.Wallis.FDR']] = TRUE
  Default.Parameters[['Include.Mann.Whitney.FDR']] = TRUE
  Default.Parameters[['Include.ttest.FDR.for.FC.association.paired.test']] = TRUE
  Default.Parameters[['Include.Wilcoxon.FDR.for.FC.association.paired.test']] = TRUE
  
  Default.Parameters[['Include.Spearman.corr']] = TRUE
  Default.Parameters[['Include.Pearson.corr']] = TRUE
  Default.Parameters[['Include.Corr.PValues']] = TRUE
  Default.Parameters[['Mann.Wh..Kr.Wall..zero.level.CPM']] = 0.5
  Default.Parameters[['Include.DeltaFreq.data']] = '{auto}'
  Default.Parameters[['DeltaFreq..delta.log.CPM.level']] = 0.5
  Default.Parameters[['DeltaFreq..CPM.const.add']] = 1
  
  Default.Parameters[['Include.GO.annotation.in.Excel.results']] = TRUE
  Default.Parameters[['Include.KEGG.annotation.in.Excel.results']] = TRUE
  Default.Parameters[['Include.Reactome.annotation.in.Excel.results']] = TRUE
  Default.Parameters[['Include.Anno.entry.names.in.Excel.results']] = TRUE
  Default.Parameters[['Include.annotation..top.genes.limit']] = 0

  Default.Parameters[['Bypass.Completed.steps']] = TRUE
  Default.Parameters[['RNA.Seq.norm.method']] = 'TMM'
  Default.Parameters[['min.samples.with.sufficient.CPM']] = '{auto}'
  Default.Parameters[['sufficient.CPM.to.analyze']] = '{auto}'
  Default.Parameters[['sufficient.CPM..autoadjust..multiplier']] = 1.0
  Default.Parameters[['sufficient.raw.read.counts.to.analyze']] = 0
  Default.Parameters[['sufficient.norm.read.counts.to.analyze']] = 0
  
  Default.Parameters[['use.quasi.likelihood.test']] = TRUE
  Default.Parameters[['use.likelihood.ratio.test']] = TRUE
  Default.Parameters[['use.exact.test.in.binary.predictors']] = TRUE
  Default.Parameters[['use.quasi.likelihood.test..paired']] = TRUE
  Default.Parameters[['use.likelihood.ratio.test..paired']] = FALSE
  Default.Parameters[['Main.test']] = 'QL'
  Default.Parameters[['Genes.sorting.criterium']] = 'default'   # also may be QL, LR, ET, Mann-Wh., Mann, Wilcoxon, Wilcox, np, Kruskal, Spearman, Pearson, default
  Default.Parameters[['Score.calculating.method']] = '{auto}'   # also may be standard, non-parametric, complex
  Default.Parameters[['PValue.score.weigth']] = 1
  Default.Parameters[['NP.PValue.score.weigth']] = 1
  Default.Parameters[['DeltaFreq.score.weigth']] = 1.6
  Default.Parameters[['LogFC.score.weigth']] = 0.75
  Default.Parameters[['LogCPM.score.weigth']] = 0.85
  Default.Parameters[['use.Mann.Whitney.criterium.for.FC.association.paired.test']] = TRUE
  Default.Parameters[['use.ttest.for.FC.association.paired.test']] = TRUE
  Default.Parameters[['use.Linear.Models.for.FC.association.paired.test']] = FALSE
  Default.Parameters[['Preferred.FC.association.test']] = 't-test (FC)'   # also may be 't-test (FC)'  'Mann-Wh. (FC)'  'Lin.Mod. (FC)'
  
  Default.Parameters[['add.Trimmed.LogFC.if.Samples.N.is.greater.than']] = 8
  Default.Parameters[['highest.expression.trimming.percents']]= 10
  Default.Parameters[['lowest.expression.trimming.percents']] = 10
  Default.Parameters[['trimmed.LogFC.const.add']] = 0.5
  Default.Parameters[['Max.Levels.to.split.Sparkline.groups']] = 4

  Default.Parameters[['Merge.samples.in.sample.setup']] = TRUE
  Default.Parameters[['Merge.samples..separator']] = ','

  Default.Parameters[['report.norm.read.counts.instead.of.CPM']] = FALSE
  Default.Parameters[['white.background.in.Excel.reports']] = TRUE
  Default.Parameters[['colored.LogFC.cells.in.Excel.reports']] = FALSE
  Default.Parameters[['bidirectional.bars.in.LogFC.cells.in.Excel.reports']] = TRUE
  Default.Parameters[['Include.small.heatmaps.in.Excel.reports']] = TRUE
  Default.Parameters[['freeze.header.in.Excel.reports']] = TRUE
  Default.Parameters[['freeze.gene.names.in.Excel.reports']] = FALSE
  
  Default.Parameters[['inDetails.additional.biotypes']] = list(lncRNA = c('lncRNA', 'lincRNA', 'antisense_RNA', 'macro_lincRNA'))
  Default.Parameters[['create.inDetails.Euler.Venn.diagrams']] = TRUE

  Default.Parameters[['disable.biomaRt']] = '{auto}'


  ### Pipeline for assembled transcriptomes

  # Default.Parameters[['Assembled.Transcripomes__transcript.lengths.by.gene.filename']] = ''
  # Default.Parameters[['Assembled.Transcripomes__GO.annotations.filename']] = ''
  # Default.Parameters[['Assembled.Transcripomes__KEGG.annotations.filename']] = ''
  # Default.Parameters[['Assembled.Transcripomes__BLAST.annotations.filename']] = ''

  # Default.Parameters[['Verbose']] = TRUE
  
  
  ### bias adjustment

  Default.Parameters[['Evaluate.bias.factors.Associations']] = TRUE
  Default.Parameters[['Bias.factors.to.analyze']] = c('Avg. Transcript Length', 'Expression level')
  Default.Parameters[['gene.bias.factors.file']] = ''
  Default.Parameters[['Draw.bias.associations.plots']] = TRUE
  Default.Parameters[['Draw.CPM..to..bias.factor.density.plots']] = TRUE

  Default.Parameters[['Adjust.Transcript.length.bias__in.Read.counts']] = FALSE
  Default.Parameters[['Adjust.Transcript.length.bias__in.Read.counts__exclude.bins']] = TRUE
  Default.Parameters[['Adjust.Transcript.length.bias__in.Read.counts__StDev.ratio__in.bin.to.exclude']] = 1.6
  Default.Parameters[['Adjust.Transcript.length.bias__in.Read.counts__StDev.abs.value__in.bin.to.exclude']] = 0.007
  Default.Parameters[['Adjust.Transcript.length.bias__in.Read.counts__max.bins.to.exclude__count']] = 1
  Default.Parameters[['Adjust.Transcript.length.bias__in.Read.counts__max.bins.to.exclude__percentage']] = 25

  Default.Parameters[['Adjust.Transcript.length.bias__in.LogFC']] = FALSE
  Default.Parameters[['Adjust.Transcript.length.bias__in.LogFC__mode']] = 'before GLM'

  Default.Parameters[['Adjust.Expression.level.bias__in.Read.counts']] = FALSE
  Default.Parameters[['Adjust.Expression.level.bias__in.Read.counts__exclude.bins']] = FALSE
  Default.Parameters[['Adjust.Expression.level.bias__in.Read.counts__StDev.ratio__in.bin.to.exclude']] = 1.6
  Default.Parameters[['Adjust.Expression.level.bias__in.Read.counts__StDev.abs.value__in.bin.to.exclude']] = 0.007
  Default.Parameters[['Adjust.Expression.level.bias__in.Read.counts__max.bins.to.exclude__count']] = 1
  Default.Parameters[['Adjust.Expression.level.bias__in.Read.counts__max.bins.to.exclude__percentage']] = 25

  Default.Parameters[['minimal.Transcript.length']] = 0
  Default.Parameters[['include.Transcript.length__quantile.statistics']] = TRUE
  Default.Parameters[['include.Expression.level__quantile.statistics']] = TRUE

  # Randomization test
  Default.Parameters[['Perform.Randomization.test']] = FALSE
  Default.Parameters[['Randomization.test__permutations.count']] = 100
  Default.Parameters[['Randomization.test__Score.threshold']] = 0
  Default.Parameters[['Randomization.test__PValue.threshold']] = 1
  Default.Parameters[['Randomization.test__FDR.threshold']] = 0.05
  Default.Parameters[['Randomization.test__abs.LogFC.threshold']] = 0
  Default.Parameters[['Randomization.test__LogCPM.threshold']] = 2
  Default.Parameters[['Randomization.test__Mann.Whitney.PValue.threshold']] = 1
  Default.Parameters[['Randomization.test__Mann.Whitney.FDR.threshold']] = 1
  Default.Parameters[['Randomization.test__deltaFreq.positive.threshold']] = 0
  Default.Parameters[['Randomization.test__deltaFreq.negative.threshold']] = 0
  Default.Parameters[['Randomization.test__max.deltaFreq.threshold']] = 0
  # Default.Parameters[['Randomization.test__delta.log.CPM.level']] = 0.5
  


  
  # Calculating distance matrices
  Default.Parameters[['Create.Distance.matrices']] = '{auto}'
  Default.Parameters[['Distance.method']] = 'canberra.weighted'
  Default.Parameters[['Minkowski.power']] = 0.5
  Default.Parameters[['Canberra.weight.power']] = 0.4
  Default.Parameters[['Min.CPM.to.include.in.dist']] = 25
  
  #Heatmaps
  Default.Parameters[['Create.Heatmaps']] =  TRUE
  Default.Parameters[['Top.genes.to.include.in.heatmaps.list']] = c(50,500,5000,25000)
  Default.Parameters[['Top.genes.to.include.inDetails.heatmaps.list']] = c(50, 100, 500, 2000, 5000)
  Default.Parameters[['Heatmaps.PValue.threshold']] = 1
  Default.Parameters[['Heatmaps.NP.PValue.threshold']] = 1
  Default.Parameters[['Heatmaps.FDR.threshold']] = 1
  Default.Parameters[['Heatmaps.abs.LogFC.threshold']] = 0
  Default.Parameters[['Heatmaps.LogCPM.threshold']] = -100
  Default.Parameters[['Create.heatmaps..with.Src.sample.order']] = TRUE
  Default.Parameters[['Create.heatmaps..with.Sorted.samples']] = TRUE
  Default.Parameters[['Create.heatmaps..with.Clustered.samples']] = TRUE
  Default.Parameters[['Create.heatmaps..with.re.Clustered.samples.inDetails']] = TRUE
  Default.Parameters[['Create.heatmaps..with.All.samples']] = TRUE
  Default.Parameters[['Create.heatmaps..CPM']] = TRUE
  Default.Parameters[['Create.heatmaps..ZScore']] = TRUE
  Default.Parameters[['Create.heatmaps..Log.rel.CPM']] = TRUE
  Default.Parameters[['Create.heatmaps..Log.rel.CPM..using.ext.samples']] = TRUE
  Default.Parameters[['Create.heatmaps..Log.rel.CPM..log.range']] = 2
  Default.Parameters[['Create.heatmaps..Log.rel.CPM..const.CPM.add']] = 1
  Default.Parameters[['Create.heatmaps..Log.rel.CPM..include.ref.samples.with.suffix']] = ""
  Default.Parameters[['Discretize.sample.annotations']] = TRUE
  Default.Parameters[['Max.Levels.to.discretize']] = 4
  
  # Gene ontology gene set enrichments analysis (GSEA)
  Default.Parameters[['Perform.GO.Enrich__with.topGO']] = TRUE
  Default.Parameters[['GO.Enrich__with.topGO___mode.list']] = c('upreg', 'downreg')
  Default.Parameters[['GO.Enrich__with.topGO___max.DE.genes.list']] = c(40, 80, 250, 500, 1000, 2000)
  Default.Parameters[['GO.Enrich__with.topGO___gene.min.Score.threshold']] = 0
  Default.Parameters[['GO.Enrich__with.topGO___gene.max.PValue.threshold']] = 0.05
  Default.Parameters[['GO.Enrich__with.topGO___gene.max.NP.PValue.threshold']] = 1.00
  Default.Parameters[['GO.Enrich__with.topGO___gene.min.LogCPM.threshold']] = -100
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
  Default.Parameters[['topGO.Enriched.Expression.Profiles___Max.PValue.threshold']] = 0.01
  Default.Parameters[['topGO.Enriched.Expression.Profiles___Max.FDR.threshold']] = 1
  Default.Parameters[['topGO.Enriched.Expression.Profiles___scoring.adjustment.to.gene.list.size']] = T
  Default.Parameters[['topGO.Enriched.Expression.Profiles___scoring.adj.list.range.optimal.size']] = 100
  Default.Parameters[['topGO.Enriched.Expression.Profiles___scoring.adj.list.range.start']] = 40
  Default.Parameters[['topGO.Enriched.Expression.Profiles___scoring.adj.list.range.end']] = 350
  Default.Parameters[['topGO.Enriched.Expression.Profiles___limit.PValues.to']] = 1e-30
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
  Default.Parameters[['Perform.Disease_Ont.classic.Enrich__with.clusterProfiler']] = TRUE
  Default.Parameters[['Perform.NCG.classic.Enrich__with.clusterProfiler']] = TRUE
  Default.Parameters[['Perform.DisGeNET.classic.Enrich__with.clusterProfiler']] = TRUE
  Default.Parameters[['Perform.WikiPathways.classic.Enrich__with.clusterProfiler']] = TRUE
  
  Default.Parameters[['classic.Enrich__with.clusterProfiler___mode.list']] = c('upreg', 'downreg')
  Default.Parameters[['classic.Enrich__with.clusterProfiler___max.DE.genes.list']] = c(40, 80, 250, 500)
  Default.Parameters[['classic.Enrich__with.clusterProfiler___gene.min.Score.threshold']] = 0
  Default.Parameters[['classic.Enrich__with.clusterProfiler___gene.max.PValue.threshold']] = 0.05
  Default.Parameters[['classic.Enrich__with.clusterProfiler___gene.max.NP.PValue.threshold']] = 1.00
  Default.Parameters[['classic.Enrich__with.clusterProfiler___gene.min.LogCPM.threshold']] = -100
  Default.Parameters[['classic.Enrich__with.clusterProfiler___term.max.PValue.threshold']] = 0.05
  Default.Parameters[['classic.Enrich__with.clusterProfiler___term.max.FDR.threshold']] = 1.0
  Default.Parameters[['classic.Enrich__with.clusterProfiler___gene.min.abs.LogFC.threshold']] = 0.3
  Default.Parameters[['classic.Enrich__with.clusterProfiler___term.max.Qvalue.threshold']] = 1.0
  Default.Parameters[['classic.Enrich__with.clusterProfiler___Render.Summary.Plot']] = TRUE
  Default.Parameters[['classic.Enrich__with.clusterProfiler___Render.Bar.plots']] = TRUE
  Default.Parameters[['classic.Enrich__with.clusterProfiler___Render.Dot.plots']] = TRUE
  Default.Parameters[['classic.Enrich__with.clusterProfiler___Render.CNE.plots']] = TRUE
  Default.Parameters[['classic.Enrich__with.clusterProfiler___Render.Heat.plots']] = TRUE
  
  Default.Parameters[['trends.Enrich__with.clusterProfiler___term.max.PValue.threshold']] = 0.05
  Default.Parameters[['trends.Enrich__with.clusterProfiler___term.max.FDR.threshold']] = 1.0
  Default.Parameters[['trends.Enrich__with.clusterProfiler___term.max.Qvalue.threshold']] = 1.0
  Default.Parameters[['trends.Enrich__with.clusterProfiler___Render.EnrichMap']] = TRUE
  Default.Parameters[['trends.Enrich__with.clusterProfiler___Render.Ridge.plot']] = TRUE
  
  # Parameters of creating gene expression profiles for Enriched GO terms + Custom GO terms (see 'GO terms (DE profiles)' sheet)
  Default.Parameters[['clusterProfiler.Expression.Profiles___Max.gene.PValue.list']] = c(0.01, 0.05, 1.00)
  Default.Parameters[['clusterProfiler.Expression.Profiles___Min.gene.logCPM.list']] = c(0, 3, 5, 7)
  Default.Parameters[['clusterProfiler.Expression.Profiles___Max.gene.PValue.for.summary.across.models']] = 1
  Default.Parameters[['clusterProfiler.Expression.Profiles___Min.gene.logCPM.for.summary.across.models']] = 3 
  Default.Parameters[['clusterProfiler.Expression.Profiles___Sort.LogFCs']] = TRUE
  Default.Parameters[['clusterProfiler.Expression.Profiles___Maximal.genes.in.terms.to.visualize']] = 120
  Default.Parameters[['clusterProfiler.Expression.Profiles___remove.redundant.entries']] = TRUE
  
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.terms.to.visualize']] = 300
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___balance.between.Classic.and.Trends.Enrich']] = 0.5
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.PValue.threshold__classic.test']] = 0.01
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.PValue.threshold__trends.test']] = 0.05
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.NP.PValue.threshold__classic.test']] = 1.00
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.NP.PValue.threshold__trends.test']] = 1.00
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.FDR.threshold__classic.test']] = 1
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___Max.FDR.threshold__trends.test']] = 1
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___scoring.adjustment.to.gene.list.size']] = TRUE
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___scoring.adj.list.range.optimal.size']] = 100
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___scoring.adj.list.range.start']] = 40
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___scoring.adj.list.range.end']] = 350
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___limit.PValues.to']] = 1e-30
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___sum.scores.with.power']] = 2.7
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___minimal.score']] = 40
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___desired.DB.entries.count__to.discard.score']] = 40
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___minimal.DB.entry.size']] = 5
  Default.Parameters[['clusterProfiler.Enriched.Expression.Profiles___sort.terms.by']] = 'name'
  Default.Parameters[['WikiPathways.GMT.file']] = '{auto}'
  
  
  for (db in c('GO', 'KEGG', 'Reactome', 'Disease_Ont', 'NCG', 'DisGeNET', 'WikiPathways')){
    for (item in names(Default.Parameters)[startsWith(names(Default.Parameters), 'classic.Enrich__with.clusterProfiler___')]){
      Default.Parameters[[sprintf('%s.%s',db,item)]] = Default.Parameters[[item]]
    }
  }
  
  for (db in c('GO','KEGG','Reactome', 'Disease_Ont', 'NCG', 'DisGeNET')){
    for (item in names(Default.Parameters)[startsWith(names(Default.Parameters), 'trends.Enrich__with.clusterProfiler___')]){
      Default.Parameters[[sprintf('%s.%s',db,item)]] = Default.Parameters[[item]]
    }
  }
  
  for (db in c('GO','KEGG','Reactome', 'Disease_Ont', 'NCG', 'DisGeNET', 'WikiPathways')){
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
  Default.Parameters[['Pathway.Visualization___CPM.aware']] = TRUE
  Default.Parameters[['Pathway.Visualization___max.gene.PValue']] = 0.05
  Default.Parameters[['Pathway.Visualization___LogFC.limits']] = 3
  
  #misc
  Default.Parameters[['Use.Info.Table']] = TRUE

  # WGCNA
  Default.Parameters[['sufficient.CPM.to.analyze.for.WGCNA.samples']] = '{auto}'
  Default.Parameters[['min.samples.with.sufficient.CPM.for.WGCNA.samples']] = '{auto}'
  Default.Parameters[['sufficient.CPM..autoadjust..multiplier.for.WGCNA.samples']] = 1.5

  Default.Parameters[['Create.CPM.density.plots.for.WGCNA.samples']] = TRUE
  Default.Parameters[['Create.MDS.plots.for.WGCNA.samples']] = TRUE
  Default.Parameters[['Create.Distance.matrices.for.WGCNA.samples']] = TRUE
  Default.Parameters[['Draw.bias.associations.plots..for.WGCNA.samples']] = FALSE
  Default.Parameters[['Draw.CPM..to..bias.factor.density.plots..for.WGCNA.samples']] = FALSE
  Default.Parameters[['WGCNA.fit.desired.slope']] = 0.05
  Default.Parameters[['WGCNA.threads']] = '{auto}'
  Default.Parameters[['WGCNA.power']] = '{auto}'

  # Euler-Venn diagrams
  Default.Parameters[['Venn.diagrams___types.list']] = c('Euler', 'Venn', 'Euler_NP', 'Venn_NP')
  Default.Parameters[['Venn.diagrams___threshold.names.list']] = c('light filt.', 'med. filt.', 'hard filt.')
  Default.Parameters[['Venn.diagrams___PValue.higher.bounds.list']] = c(0.05, 0.05, 0.05)
  Default.Parameters[['Venn.diagrams___PValue.lower.bounds.list']] = c(0.01, 0.01, 0.01)
  Default.Parameters[['Venn.diagrams___FDR.higher.bounds.list']] = c(1, 0.10, 0.01)
  Default.Parameters[['Venn.diagrams___FDR.lower.bounds.list']] = c(1, 0.05, 0.005)
  Default.Parameters[['Venn.diagrams___NP.PValue.higher.bounds.list']] = c(0.07, 0.03, 0.005)
  Default.Parameters[['Venn.diagrams___NP.PValue.lower.bounds.list']] = c(0.05, 0.01, 0.001)
  Default.Parameters[['Venn.diagrams___abs.LogFC.higher.bounds.list']] = c(0.5, 0.6, 1.0)
  Default.Parameters[['Venn.diagrams___abs.LogFC.lower.bounds.list']] = c(0.3, 0.4, 0.5)
  Default.Parameters[['Venn.diagrams___LogCPM.higher.bounds.list']] = c(1.5, 1.5, 1.5)
  Default.Parameters[['Venn.diagrams___LogCPM.lower.bounds.list']] = c(1, 1, 1)
  Default.Parameters[['Venn.diagrams___Include.NP.diagrams.if.avg.samples.N.per.group.exceeds']] = 8
  Default.Parameters[['Venn.diagrams___image.resolution']] = 300

  
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


Text.to.str.vector = function(x, sep = c('\\','/',' ',';')){
  for (symbol in sep) x = gsub(symbol, ',', x, fixed = TRUE)
  x = Delete.duplicated.symbols.in.text(x, ',')
  return(trimws(strsplit(x, ',')[[1]]))
}

Text.to.int.vector = function(x){
  for (symbol in c('\\','/',' ',';')) x = gsub(symbol,',',x,fixed=T)
  x = Delete.duplicated.symbols.in.text(x,',')
  return(as.numeric(as.character(strsplit(x,',')[[1]])))
}

Process.Parameters.table = function(Pars, par.table, show.warnings = TRUE, base.dir = NULL, script.dir = NULL, current.dir = NULL, suppress.base.dir.errors = FALSE){
  if(!is.null(base.dir))     if(endsWith(base.dir, '/'))     base.dir = substr(base.dir, 1, nchar(base.dir) - 1)
  if(!is.null(script.dir))   if(endsWith(script.dir, '/'))   script.dir = substr(script.dir, 1, nchar(script.dir) - 1)
  if(!is.null(current.dir))  if(endsWith(current.dir, '/'))  current.dir = substr(current.dir, 1, nchar(current.dir) - 1)

  if(is.null(current.dir)) current.dir = getwd()
  par.table = Eliminate.non.informative.lines(par.table)
  par.table[,1] = sapply(par.table[,1],function(x) gsub('"',"",gsub('"','',x)))
  par.table[,2] = sapply(par.table[,2],function(x) gsub('"',"",gsub('"','',x)))
  par.table[,1] = sapply(par.table[,1],function(x) gsub("'","",gsub('"','',x)))
  par.table[,2] = sapply(par.table[,2],function(x) gsub("'","",gsub('"','',x)))
  par.table[,1] = sapply(par.table[,1],function(x) gsub(' ', ".", x, fixed = TRUE))
  par.table[,1] = sapply(par.table[,1],function(x) gsub('-', "_", x, fixed = TRUE))

  i = 1
  for (i in 1:dim(par.table)[1]){
    var = par.table[i,1]
    val = par.table[i,2]
    if (!(var %in% names(Pars))){
      message(sprintf('Warning. Unknown parameter "%s". Omitting',var))
      next
    }
    if (val == '{from basic parameters}'){
      next
    } else if (var == 'DE.package' | var == 'Method'){
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

    } else if (var == 'Main.test'){
      val = casefold(val)
      if(val %in% c('lr', 'lr-test', 'lr test', 'likelihood ratio', 'lrt', 'lrtfit', 'glmfit')){
        Pars[[var]] = 'LR'
      } else if (val %in% c('f', 'f-test', 'f test', 'quasi likelihood ratio', 'quasi', 'ql', 'qlf', 'qlf test', 'qlfit', 'qlglmfit', 'ql test', 'ql-test')){
        Pars[[var]] = 'QL'
      } else if(val %in% c('exact', 'exact test', 'et')){
        Pars[[var]] = 'ET'
      } else if(val %in% c('lr paired', 'lr-test paired', 'lr test paired', 'likelihood ratio paired', 'lrt paired', 'lrtfit paired', 'glmfit paired')){
        Pars[[var]] = 'LR paired'
      } else if (val %in% c('f paired', 'f-test paired', 'f test paired', 'quasi likelihood ratio paired', 'quasi paired', 'ql paired', 'qlf paired', 'qlf test paired', 'qlfit paired', 'qlglmfit paired', 'ql test paired', 'ql-test paired')){
        Pars[[var]] = 'QL paired'
      } else {
        stop('Incorrect parameter "Main.test". Should be either "LR", "QL", "Exact test". Please check "Parameters" worksheet')
      }

    } else if (var %in% c('GO.enrichment.namespaces', 'GO.enrichment.namespaces')){
      val = Text.to.str.vector(val, sep = c('; ', ' ;', ', ', ' ,', ';', ','))
      val = toupper(val)
      incorrect.vals = val[!(val %in% c('BP', 'CC', 'MF'))]
      if(length(incorrect.vals) > 0)
        stop(sprintf('Incorrect parameter %s ("%s"). Only BP, MF and CC namespaces are allowed', var, toString(val)))
      Pars[[var]] = val

    } else if (var == 'Genes.sorting.criterium'){
      val = casefold(val)
      if(val %in% c('default', 'auto', '{auto}', 'score')){
        Pars[[var]] = 'default'
      } else if(val %in% c('lr', 'lr-test', 'lr test', 'likelihood ratio', 'lrt', 'lrtfit', 'glmfit')){
        Pars[[var]] = 'LR test'
      } else if (val %in% c('f', 'f-test', 'f test', 'quasi likelihood ratio', 'quasi', 'ql', 'qlf', 'qlf test', 'qlfit', 'qlglmfit', 'ql test', 'ql-test')){
        Pars[[var]] = 'QLF test'
      } else if(val %in% c('exact', 'exact test', 'et')){
        Pars[[var]] = 'ET test'
      } else if(val %in% c('lr paired', 'lr-test paired', 'lr test paired', 'likelihood ratio paired', 'lrt paired', 'lrtfit paired', 'glmfit paired')){
        Pars[[var]] = 'LR test paired'
      } else if (val %in% c('f paired', 'f-test paired', 'f test paired', 'quasi likelihood ratio paired', 'quasi paired', 'ql paired', 'qlf paired', 'qlf test paired', 'qlfit paired', 'qlglmfit paired', 'ql test paired', 'ql-test paired')){
        Pars[[var]] = 'QLF test paired'
      } else if(val %in% c('mann-wh', 'mann-whitney', 'mann', 'mann-wh.')){
        Pars[[var]] = 'Mann-Wh.'
      } else if(val %in% c('mann-wh paired', 'mann-whitney paired', 'mann paired', 'mann-whitney','mann-wh. paired')){
        Pars[[var]] = 'Mann-Wh., FC'
      } else if(val %in% c('wilcox', 'wilcoxon', 'wilcox paired', 'wilcoxon paired')){
        Pars[[var]] = 'Wilcoxon paired'
      } else if(val %in% c('kruskal', 'kruskal wallis', 'kruskal-wallis')){
        Pars[[var]] = 'Krusk-W.'
      } else if(val %in% c('np', 'non parametric', 'non-parametric')){
        Pars[[var]] = 'NP'
      } else if(val %in% c('student', 't', 'ttest', 't-test', 't test', 'student\'s t-test', 'student t-test', 'student\'s t test', 'student t test')){
        Pars[[var]] = 't-test'
      } else if(val %in% c('student paired', 't paired', 'ttest paired', 't-test paired', 't test paired', 'student\'s t-test paired', 'student t-test paired', 'student\'s t test paired', 'student t test paired')){
        Pars[[var]] = 't-test, FC'
      } else if(val %in% c('spearman')){
        Pars[[var]] = 'Spearman'
      } else if(val %in% c('spearman paired', 'spearman fc')){
        Pars[[var]] = 'Spearman, FC'
      } else if(val %in% c('pearson')){
        Pars[[var]] = 'Pearson'
      } else if(val %in% c('pearson paired', 'pearson fc')){
        Pars[[var]] = 'Pearson, FC'
      } else if(val %in% c('delta', 'delta freq', 'delta-freq', 'deltafreq')){
        Pars[[var]] = 'DeltaFreq'
      } else {
        stop(sprintf('Incorrect parameter "Genes.sorting.criterium". It is set as "%s". Should be either "LR", "QL", "Exact test", "LR paired", "QL paired", "Wilcoxon", "Mann-Whitney", "Spearman", "Pearson", "t-test", "non-parametric", or "default", "score". Please check "Parameters" worksheet', val))
      }

    } else if (var == 'Preferred.FC.association.test'){
      val = casefold(val)
      if(val %in% c('t', 't-test', 't-test (fc)', 't test','t test (fc)', 'student', 'student t-test')){
        Pars[[var]] = 't-test (FC)'   # also may be 't-test (FC)'  'Mann-Wh. (FC)'  'Lin.Mod. (FC)'
      } else if (val %in% c('mann-wh', 'mann-wh.', 'mann-whitney', 'mann-wh (fc)', 'mann-wh. (fc)', 'mann-whitney (fc)')){
        Pars[[var]] = 'Mann-Wh. (FC)'
      } else if(val %in% c('lin.mod. (FC)', 'lin.mod.', 'lm (FC)', 'lm')){
        Pars[[var]] = 'Lin.Mod. (FC)'
      }

    } else if (var == 'Clustering.method'){
      allowed.clustering.methods = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
      # val = casefold(val)
      if(!(val %in% allowed.clustering.methods))
        stop(sprintf('Incorrect Clustering.method %s. Allowed: %s', val, toString(allowed.clustering.methods)))
      Pars[[var]] = val
    } else if (var %in% c('counts.dir','results.dir','suppl.data.dir','read.counts.table','counts.suffix')){
      Pars[[var]] = gsub('\\','/',val,fixed=T)
    } else if('.bias__in.LogFC__mode' %in% var){
      if(val %in% c('before GLM' , 'after GLM')){
        Pars[[var]] = val
      } else  stop(sprintf('Incorrect parameter %s = %s', var, val))
    } else if(var == 'Bias.factors.to.analyze'){
      val = Text.to.str.vector(val, sep = c('; ', ' ;', ', ', ' ,', ';', ','))
      Pars[[var]] = val
    } else if(var == 'inDetails.additional.biotypes'){
      split1 = Text.to.str.vector(val, sep = c(';', ','))
      split2 = list()
      for(s in split1){
        split_s = Text.to.str.vector(s, sep = c(':'))
        if(length(split_s) != 2)
          stop(sprintf('Incorrect format of "inDetails.additional.biotypes" parameter: "%s". Should be like this:\nlncRNA: lncRNA + lincRNA + antisense_RNA + macro_lincRNA, pseudo: processed_pseudogene + polymorphic_pseudogene', val))
        split2[[split_s[1]]] = Text.to.str.vector(split_s[2], sep = c('+'))
      }
      Pars[[var]] = split2
    } else if (val %in% c('{auto}', 'auto')){
      next
    } else if (var == 'Distance.method'){
      val = casefold(val)
      if (!(val %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski",'1-cor','canberra.weighted'))){
        stop(sprintf('Incorrect parameter "Distance.method = %s". These values are allowed: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski","1-cor", "canberra.weighted" (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)', val))
      }
      Pars[[var]] = val
    } else if (var == 'RNA.Seq.norm.method'){
      if (!(val %in% c("TMM", "upperquartile", "none", "RLE", 'Reference genes', 'RG'))){
        stop(sprintf('Incorrect parameter "RNA.Seq.norm.method = %s". These values are allowed: "TMM","upperquartile","none","RLE" (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)', val))
      }
      if(val == 'Reference genes')  val = 'RG'
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
        stop(sprintf('Incorrect parameter "topGO.Enriched.Expression.Profiles___sort.terms.by = %s". These values are allowed: "name","id","score","none" (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)', val))
      }
      # if (val == 'p') val = 'pvalue'
      Pars[[var]] = val
    } else if (var %in% c('min.samples.with.sufficient.CPM', 'sufficient.CPM.to.analyze', 'Max.genes.in.Excel.results')){
      if (val %in% c('{auto}' , 'auto', 'unlimited')) {
        Pars[[var]] = val
      } else {
        val = as.numeric(as.character(val))
        if (is.na(var)){
          stop(sprintf('Incorrect parameter "%s". Only numeric values or "{auto}" are allowed (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)',var))
        }
        Pars[[var]] = val
      }
    } else if (var %in% c('Venn.diagrams___threshold.names.list', 'Venn.diagrams___types.list')) {
      val = Text.to.str.vector(val, sep = c(';', ','))
      Pars[[var]] = val
    } else if (endswith(var,'ntology.list')) {
      val = toupper(val)
      val = Text.to.str.vector(val)
      if (length(val) == 0 || any(!(val %in% c('BP', 'MF', 'CC')))){
        stop(sprintf('Incorrect parameter "%s". Only "BP", "MF" or "CC" ontologies are allowed (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)',var))
      }  
      Pars[[var]] = val
    } else if (endswith(var,'_mode.list')) {
      val = casefold(val)
      val = Text.to.str.vector(val)
      if (length(val) == 0 || any(!(val %in% c('upreg', 'downreg', 'up', 'down')))){
        stop(sprintf('Incorrect parameter "%s". Only "upreg" or "downreg" values are allowed (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)',var))
      }  
      Pars[[var]] = val
    } else if (endswith(var,'.list')) {
      val = Text.to.int.vector(val)
      if(any(is.na(val))){
        stop(sprintf('Incorrect parameter "%s". Only numeric values are allowed, e.g. "20,40,80,500" (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)',var))
      }
      Pars[[var]] = val
    } else if (typeof(Pars[[var]]) == "logical"){
      val = casefold(val)
      if (!(val %in% c('true','t','f','false','yes','y','n','no','on','off'))){
        stop(sprintf('Incorrect parameter "%s". Only logical values are allowed, e.g. True/False (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)',var))
      }
      val = (substring(val,1,1) == 't' || substring(val,1,1) == 'y' || val == 'on')
      Pars[[var]] = val
    } else if (typeof(Pars[[var]]) %in% c("numeric","double","float")){
      val = as.numeric(as.character(val))
      if (is.na(var)){
        stop(sprintf('Incorrect parameter "%s". Only numeric values are allowed (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)',var))
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

  for(var in names(Pars)){
    if(casefold(toString(Pars[[var]])) %in% c('true', 'yes', 'on')){
      Pars[[var]] = TRUE
    } else if(casefold(toString(Pars[[var]])) %in% c('false', 'no', 'off')){
      Pars[[var]] = FALSE
    }
  }

  Venn.thr.list.length = length(Pars$Venn.diagrams___threshold.names.list)
  for(att in c('Venn.diagrams___PValue.lower.bounds.list',
    'Venn.diagrams___FDR.higher.bounds.list',
    'Venn.diagrams___FDR.lower.bounds.list',
    'Venn.diagrams___NP.PValue.higher.bounds.list',
    'Venn.diagrams___NP.PValue.lower.bounds.list',
    'Venn.diagrams___abs.LogFC.higher.bounds.list',
    'Venn.diagrams___abs.LogFC.lower.bounds.list',
    'Venn.diagrams___LogCPM.higher.bounds.list',
    'Venn.diagrams___LogCPM.lower.bounds.list')){
    if(length(Pars[[att]]) != Venn.thr.list.length){
      stop(sprintf('The length of %s (%d) should be equal to the length of Venn.diagrams___threshold.names.list (%d)', att, length(Pars[[att]]), Venn.thr.list.length))
    }
  }

  if(!Pars$Evaluate.bias.factors.Associations)   Pars$Bias.factors.to.analyze = c()
  
  if(Pars$Adjust.Transcript.length.bias__in.Read.counts)      Pars$Evaluate.bias.factors.Associations = TRUE
  if(Pars$Adjust.Transcript.length.bias__in.LogFC)            Pars$Evaluate.bias.factors.Associations = TRUE
  if(Pars$Adjust.Expression.level.bias__in.Read.counts)       Pars$Evaluate.bias.factors.Associations = TRUE

  if (Pars$Adjust.Transcript.length.bias__in.Read.counts | Pars$Adjust.Transcript.length.bias__in.LogFC){
    if(!('Avg. Transcript Length' %in% Pars$Bias.factors.to.analyze))  Pars$Bias.factors.to.analyze = c(Pars$Bias.factors.to.analyze, 'Avg. Transcript Length')
  }

  if(Pars$Adjust.Expression.level.bias__in.Read.counts){
   if(!('Expression level' %in% Pars$Bias.factors.to.analyze))  Pars$Bias.factors.to.analyze = c(Pars$Bias.factors.to.analyze, 'Expression level') 
  }


  if(Pars$suppl.data.dir %in% c('{auto}', 'auto', '')){
    if(is.null(base.dir) && !suppress.base.dir.errors)  stop('Please specify base.dir argument because Pars$suppl.data.dir is set as "auto"')
    Pars$suppl.data.dir = base.dir
    if(!file.exists(base.dir) && !suppress.base.dir.errors)
      stop(sprintf('Directory "%s" (base.dir) does not exist', base.dir))
  }

  if(endsWith(Pars$suppl.data.dir, '/'))  Pars$suppl.data.dir = substr(Pars$suppl.data.dir, 1, nchar(Pars$suppl.data.dir) - 1)
  if(endsWith(Pars$results.dir, '/'))  Pars$results.dir = substr(Pars$results.dir, 1, nchar(Pars$results.dir) - 1)
  if(!is.null(Pars$counts.dir))  if(endsWith(Pars$counts.dir, '/'))  Pars$counts.dir = substr(Pars$counts.dir, 1, nchar(Pars$counts.dir) - 1)
  
  if (Pars$Main.test == 'ET'){
    Pars$use.exact.test.in.binary.predictors = TRUE
  } else if (Pars$Main.test == 'QL'){
    Pars$use.quasi.likelihood.test = TRUE
  } else if (Pars$Main.test == 'LR'){
    Pars$use.likelihood.ratio.test = TRUE
  } else if(Pars$Main.test == 'QL paired'){
    Pars$use.quasi.likelihood.test..paired = TRUE
  } else if(Pars$Main.test == 'LR paired'){
    Pars$use.likelihood.ratio.test..paired = TRUE
  }
  
  # stop(Pars$counts.dir)
  
  if(!is.null(script.dir)){
    Pars$suppl.data.dir = gsub('\\{script_dir\\}', script.dir, Pars$suppl.data.dir)
    Pars$results.dir = gsub('\\{script_dir\\}', script.dir, Pars$results.dir)
    Pars$counts.dir = gsub('\\{script_dir\\}', script.dir, Pars$counts.dir)
    Pars$read.counts.table = gsub('\\{script_dir\\}', script.dir, Pars$read.counts.table)
    Pars$WikiPathways.GMT.file = gsub('\\{script_dir\\}', script.dir, Pars$WikiPathways.GMT.file)
  }

  if(!is.null(base.dir)){
    Pars$suppl.data.dir = gsub('\\{base_dir\\}', base.dir, Pars$suppl.data.dir)
    Pars$results.dir = gsub('\\{base_dir\\}', base.dir, Pars$results.dir)
    Pars$counts.dir = gsub('\\{base_dir\\}', base.dir, Pars$counts.dir)
    Pars$read.counts.table = gsub('\\{base_dir\\}', base.dir, Pars$read.counts.table)
    Pars$WikiPathways.GMT.file = gsub('\\{base_dir\\}', base.dir, Pars$WikiPathways.GMT.file)
  }

  for(var in c('counts.dir', 'results.dir', 'read.counts.table')){
    if(Pars[[var]] != '' & Pars[[var]] != 'auto' & Pars[[var]] != '{auto}' & !(substr(Pars[[var]], 1, 1) %in% c('/', '\\')) & !(substr(Pars[[var]], 2, 2) %in% c(':'))){
      Pars[[var]] = file.path(current.dir, Pars[[var]])
    }
  }

  if(Pars$WikiPathways.GMT.file %in% c('auto', '{auto}')){
    DB.data = Read.Organisms.Info(Organisms.info.file.name = sprintf('%s/Organisms.info.xlsx', Pars$suppl.data.dir), script.dir = script.dir)
    if(!(Pars$Species %in% names(DB.data$Taxons.by.KEGG.codes))){
      warning(sprintf('Unknown species "%s". There is no such KEGG code in "Organisms.info.xlsx". Hence, cannot locate corresponding WikiPathways GMT-file for this species. Please specify it manually (WikiPathways.GMT.file parameter - see Advanced parameters sheet).  You can download GMT from http://data.wikipathways.org/current/gmt/', Pars$Species))
      Pars$WikiPathways.GMT.file = NULL
    } else {
      possible.files = list.files(path = Pars$suppl.data.dir, pattern = sprintf('wikipathways-[0123456789]*-gmt-%s',
        gsub(' ', '_', DB.data$Taxons.by.KEGG.codes[Pars$Species], fixed = TRUE)), full.names=TRUE)
      if(length(possible.files) == 0){
        warning(sprintf('Cannot locate corresponding WikiPathways GMT-file for this species "%s" (%s). Please specify it manually (WikiPathways.GMT.file parameter - see Advanced parameters sheet). You can download GMT from http://data.wikipathways.org/current/gmt/',
          DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species))
        Pars$WikiPathways.GMT.file = NULL
      } else {
        Pars$WikiPathways.GMT.file = sort(possible.files, decreasing = TRUE)[1]
      }
    }

  } else if(!file.exists(Pars$WikiPathways.GMT.file)){
    msg = sprintf('WikiPathways GMT file "%s" does not exist. Please download it from http://data.wikipathways.org/current/gmt/\n', Pars$WikiPathways.GMT.file)
    cat(msg)
    warning(msg)
    Pars$WikiPathways.GMT.file = NULL
  }

  if(Pars$counts.dir == '' & Pars$read.counts.table == '')    stop('Please specify either counts.dir+counts.suffix or read.counts.table parameter')
  if(Pars$counts.dir != '' & Pars$read.counts.table != '')    stop('Please specify either counts.dir+counts.suffix OR read.counts.table parameter, but not both of them')
  if(Pars$counts.dir != '' & Pars$counts.suffix == '')    stop('Please specify counts.suffix parameter correctly')
  #Pars$Use.Info.Table = T
  
  # looking for unfilled parameters
  if(show.warnings){
    for (n in names(Pars)[!(names(Pars) %in% par.table[,1])]){
      if(!(n %in% c('counts.dir', 'counts.suffix', 'read.counts.table', 'gene.bias.factors.file')))
        message(sprintf('Parameter %s is missing. Setting to defaults (%s)', n, toString(Pars[[n]])))
    }
  }
  return(Pars)
}

Get.script.path = function(emit.warning = TRUE){
  tryCatch(expr = {
    funr::get_script_path()          
  }, error = function (err){
    tryCatch(expr = {
      dirname(rstudioapi::getActiveDocumentContext()$path)
    }, error = function (err){
      if(emit.warning) warning('Cannot locate script directory')
      getwd()
    })
  })
}

Read.Parameters <- function(Parameters.xlsx.file.name = NULL, prepare.for.ext.DE.data = FALSE, script.dir = NULL, show.warnings..on.parameters = TRUE, base.dir = NULL, current.dir = NULL){
  #!!!! as.numeric(as.character(testdata$x)), not as.numeric(testdata$x) 
  
  if(!is.null(base.dir))     if(endsWith(base.dir, '/'))     base.dir = substr(base.dir, 1, nchar(base.dir) - 1)
  if(!is.null(script.dir))   if(endsWith(script.dir, '/'))   script.dir = substr(script.dir, 1, nchar(script.dir) - 1)
  if(!is.null(current.dir))  if(endsWith(current.dir, '/'))  current.dir = substr(current.dir, 1, nchar(current.dir) - 1)

  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")
  if (!('rstudioapi' %in% rownames(installed.packages()))) install.packages("rstudioapi")
  
  if(is.null(script.dir))  script.dir = Get.script.path()
  if(is.null(current.dir))  current.dir = getwd()

  if (is.null(Parameters.xlsx.file.name)){
    Parameters.xlsx.file.name = sprintf("%s/RTrans.parameters.xlsx", script.dir)
    if (file.exists(Parameters.xlsx.file.name)){
      cat(sprintf('\nFile with parameters found: %s\n',Parameters.xlsx.file.name))
    } else stop(sprintf('Parameters.xlsx.file.name argument is not specified. No parameters *.xlsx file has been found in the default location '))
  }
  
  if(!is.null(script.dir))  Parameters.xlsx.file.name = gsub('\\{script_dir\\}', script.dir, Parameters.xlsx.file.name)
  
  
  if (!file.exists(Parameters.xlsx.file.name)) stop(sprintf('file %s does not exist',Parameters.xlsx.file.name))
  if (substring(Parameters.xlsx.file.name,nchar(Parameters.xlsx.file.name)-4,nchar(Parameters.xlsx.file.name)) != '.xlsx') stop('Please fill up input parameters in Excel workbook (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)')
  
  
  # reading main worksheet 'Parameters'
  
  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")

  available.sheets = readxl::excel_sheets(Parameters.xlsx.file.name)
  if('Basic parameters' %in% available.sheets && 'Advanced parameters' %in% available.sheets){
    par.table..basic = data.frame(suppressMessages(suppressWarnings(readxl::read_excel(Parameters.xlsx.file.name, sheet='Basic parameters', col_names = FALSE))))[,1:2]
    par.table..adv = data.frame(suppressMessages(suppressWarnings(readxl::read_excel(Parameters.xlsx.file.name, sheet='Advanced parameters', col_names = FALSE))))[,1:2]
    par.table = rbind(par.table..basic, par.table..adv)

  } else if('Parameters' %in% available.sheets){
    par.table = data.frame(suppressMessages(suppressWarnings(readxl::read_excel(Parameters.xlsx.file.name, sheet='Parameters', col_names = FALSE))))[,1:2]

  } else {
    stop(sprintf('Sheets "Basic parameters", "Advanced parameters" are found in workbook %s', Parameters.xlsx.file.name))
  }
  colnames(par.table) = c('variable','value')
  
  Pars = Get.Default.Parameters()
  
  schema = tryCatch(expr = {suppressMessages(suppressWarnings(readxl::read_excel(Parameters.xlsx.file.name,sheet='Sample setup',col_names = TRUE)))
  }, error = function (err){
    stop(sprintf('Sheet "Sample setup" is not found in workbook %s',Parameters.xlsx.file.name))
  })
  
  available.preds = colnames(schema)[-1]
  available.preds_to_abbreviations = sprintf('pred_%.4d',1:length(available.preds))
  names(available.preds_to_abbreviations) = available.preds
  abbreviations_to_available.preds = available.preds
  names(abbreviations_to_available.preds) = sprintf('pred_%.4d',1:length(available.preds))
  
  
  Pars = Process.Parameters.table(Pars, par.table, show.warnings = show.warnings..on.parameters, base.dir = base.dir, script.dir = script.dir, current.dir = current.dir)


  # testing "Models.to.test" parameter
  if(!prepare.for.ext.DE.data){
    if (!('Models.to.Test' %in% names(Pars))){
      stop('Mandatory parameter "Models.to.Test" is missing (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)')
    }
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
  cat(sprintf('\ngene list hashmd5 for biomaRt: %s (contains %d genes)\n', hashmd5, length(complete.gene.list)))
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

    if(!Pars$disable.biomaRt){
      cat('\nQuerying biomaRt for gene info and saving this data to disk... Sometimes this operation may freeze or interrupt, depending on the BioMart web-service load. \n')
      mart = tryCatch(expr = {
        useMart("ensembl", dataset=dataset.name) #, host="www.ensembl.org"
      }, error = function (err){
        xmart = useMart("ensembl") #, host="www.ensembl.org"
        stop(sprintf('\nThe dataset %s is not available from biomaRt. Available datasets: %s', dataset.name, toString(listDatasets(xmart))))
      })
    }
    # if (Pars$Species == 'mmu') mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl") #, host="www.ensembl.org"
    # if (Pars$Species == 'dme') mart <- useMart("ensembl", dataset="dmelanogaster_gene_ensembl") #, host="www.ensembl.org"
    #curl.handle = getCurlHandle()
    
    needed.attributes = c("ensembl_gene_id","external_gene_name", "description","gene_biotype","entrezgene_id")
    if (Pars$Species == 'dme')   needed.attributes = c(needed.attributes,'flybasename_gene','flybase_gene_id')
    
    if(Pars$disable.biomaRt){
      General.maRt.table = data.frame(array(dim = c(0, length(needed.attributes))))
      colnames(General.maRt.table) = needed.attributes
    } else if(use.official.gene.symbol){
      General.maRt.table = getBM(attributes=needed.attributes, filters="external_gene_name", values=complete.gene.list, mart=mart)
      General.maRt.table = General.maRt.table[!(duplicated(General.maRt.table[,"external_gene_name"])),]
      rownames(General.maRt.table) = General.maRt.table[,"external_gene_name"]
    } else {
      General.maRt.table = getBM(attributes=needed.attributes,filters="ensembl_gene_id",values=complete.gene.list, mart=mart)
      ##write.table(listAttributes(mart = mart),"ttt.tsv",sep = '\t')  ## listing all the attributes
      General.maRt.table = General.maRt.table[!(duplicated(General.maRt.table[,"ensembl_gene_id"])),]
      rownames(General.maRt.table) = General.maRt.table[,"ensembl_gene_id"]
    }
    General.maRt.table[,"description"] = gsub(pattern = "\\[Source:.*", replacement = "", x = General.maRt.table[,"description"], ignore.case = TRUE, perl = FALSE)
    write.table.mod(x = General.maRt.table, file = General.maRt.table_file.name, sep = '\t')
  }
  # org.Dm.eg.db::
  return(General.maRt.table)
}


# Organisms.info.file.name = 'RTrans.base/Organisms.info.xlsx'
Read.Organisms.Info <- function(Organisms.info.file.name = NULL, script.dir = NULL){

  DB.data = vector(mode = "list")
  
  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")
  if (!('rstudioapi' %in% rownames(installed.packages()))) install.packages("rstudioapi")
  
  if(is.null(script.dir))  script.dir = Get.script.path()

  if (is.null(Organisms.info.file.name)){
    Organisms.info.file.name = sprintf("%s/../RTrans_main_scripts/Organisms.info.xlsx", script.dir)
    if (file.exists(Organisms.info.file.name)){
      cat(sprintf('\nFile with organisms info found: %s\n',Organisms.info.file.name))
    } else stop(sprintf('Organisms.info.file.name argument is not specified. No organism info *.xlsx file has been found in the default location '))
  }
  
  if(!is.null(script.dir))  Organisms.info.file.name = gsub('\\{script_dir\\}',script.dir,Organisms.info.file.name)
  
  
  if (!file.exists(Organisms.info.file.name)) stop(sprintf('file %s does not exist',Organisms.info.file.name))
  if (substring(Organisms.info.file.name,nchar(Organisms.info.file.name)-4,nchar(Organisms.info.file.name)) != '.xlsx') stop('Please fill up organisms info in Excel workbook')
  
  
  # reading main worksheet 'KEGG organism codes'
  
  tryCatch(expr = { KEGG.code.table = suppressMessages(suppressWarnings(readxl::read_excel(Organisms.info.file.name,sheet='KEGG organism codes',col_names = TRUE)))
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
  
  tryCatch(expr = { Bioconductor.code.table = suppressMessages(suppressWarnings(readxl::read_excel(Organisms.info.file.name,sheet='Bioconductor DBs',col_names = TRUE)))
  }, error = function (err){
    stop(sprintf('Sheet "Bioconductor DBs" is not found in workbook %s',Organisms.info.file.name))
  })
  
  Bioconductor.code.table = as.data.frame(Bioconductor.code.table)
  Bioconductor.code.table = Bioconductor.code.table[!(duplicated(Bioconductor.code.table[,'KEGG organism code'])),]
  Bioconductor.Org.DBs.by.KEGG.codes = Bioconductor.code.table[,'Bioconductor DB']
  names(Bioconductor.Org.DBs.by.KEGG.codes) = Bioconductor.code.table[,'KEGG organism code']
  
  DB.data[['Bioconductor.Org.DBs.by.KEGG.codes']] = Bioconductor.Org.DBs.by.KEGG.codes
  
  tryCatch(expr = { Reactome.code.table = suppressMessages(suppressWarnings(readxl::read_excel(Organisms.info.file.name,sheet='Reactome names',col_names = TRUE)))
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

Read.DB.entries.from.Excel = function(file.name, sheet.name, stop..if.sheet.is.not.found = TRUE){
  #file.name = 'C:/Users/gskra/OneDrive/Moskalev flies-2016/RTrans.parameters.xlsx'
  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")
  if (substring(file.name,nchar(file.name)-4,nchar(file.name)) != '.xlsx') stop('Please fill up input parameters in Excel workbook (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)')
  if (!file.exists(file.name)) stop(sprintf('file %s does not exist',file.name))
  
  # reading worksheet
  par.table = tryCatch(expr = { suppressMessages(suppressWarnings(readxl::read_excel(file.name, sheet=sheet.name, col_names = FALSE)))
  }, error = function (err){
    if(stop..if.sheet.is.not.found){
      stop(sprintf('Sheet "%s" is not found in workbook %s', sheet.name, file.name))
    } else {
      return(NULL)
    }
  })
  
  if(is.null(par.table))  return(c())

  par.table = data.frame(par.table)
  par.table = par.table[,1]
  par.table = par.table[!is.na(par.table)]
  commented.lines = sapply(X = par.table,FUN = function (x) substring(x,first=1,last=1)) %in% '#'
  par.table = par.table[!commented.lines]
  return(par.table)
}

Read.Reference.Genes = function(file.name){
  Read.DB.entries.from.Excel(file.name, 'Reference genes', stop..if.sheet.is.not.found = FALSE)
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
  if (substring(file.name,nchar(file.name)-4,nchar(file.name)) != '.xlsx') stop('Please fill up input parameters in Excel workbook (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)')
  if (!file.exists(file.name)) stop(sprintf('file %s does not exist',file.name))
  
  # reading worksheet
  tryCatch(expr = { par.table = suppressMessages(suppressWarnings(readxl::read_excel(file.name,sheet=sheet.name,col_names = FALSE)))
  }, error = function (err){
    stop(sprintf('Sheet "%s" is not found in workbook %s', sheet.name, file.name))
  })
  
  par.table = data.frame(par.table)
  par.table = par.table[,1]
  par.table = par.table[!is.na(par.table)]
  commented.lines = sapply(X = par.table,FUN = function (x) substring(x,first=1,last=1)) %in% '#'
  par.table = par.table[!commented.lines]
  par.table = gsub("[^0-9]", "", par.table)
  final.list = as.numeric(as.character(par.table))
  if (any(is.na(final.list)))  message(sprintf('Please provide only numerical KEGG id (e.g. 5133 instead of hsa05133) in "%s" sheet. Such values:  %s', sheet.name, toString(par.table[is.na(final.list)])))
  final.list = final.list[!is.na(final.list)]
  
  return(sprintf("%s%.5d",Species,final.list))
}

Read.KEGG.pathways.to.omit <- function(file.name,Species){
  #file.name = 'C:/Users/gskra/OneDrive/Moskalev flies-2016/RTrans.parameters.xlsx'
  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")
  if (substring(file.name,nchar(file.name)-4,nchar(file.name)) != '.xlsx') stop('Please fill up input parameters in Excel workbook (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)')
  if (!file.exists(file.name)) stop(sprintf('file %s does not exist',file.name))
  
  # reading worksheet 'Omit KEGG pathways'
  tryCatch(expr = { par.table = suppressMessages(suppressWarnings(readxl::read_excel(file.name,sheet='Omit KEGG pathways',col_names = FALSE)))
  }, error = function (err){
    stop(sprintf('Sheet "Omit KEGG pathways" is not found in workbook %s',file.name))
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

Convert.GO.term.names_to_GO.term.IDs = function(terms, GO.term.ID.by.name, remove.redundant = TRUE, prefix = '', messages = TRUE){
  not_found_terms_count = 0
  terms__converted = c()
  for (x in terms){
    x.up = toupper(x)
    if (!(substring(x.up,1,3) %in% 'GO:')){
      if (x.up %in% names(GO.term.ID.by.name)){
        x = GO.term.ID.by.name[x.up]
      } else {
        if(messages)
          message(sprintf('\n%sWarning. Custom GO term "%s" is not found', prefix, x))
        not_found_terms_count = not_found_terms_count + 1
      }
    }
    if (!(remove.redundant) | (!(x %in% terms__converted))){
      terms__converted = c(terms__converted, x)
    }
  }
  if (not_found_terms_count > 0 & messages)   message(sprintf('\n%sTotal custom GO %d terms were not found', prefix, not_found_terms_count))
  return(terms__converted)
}

# Parameters.file.name = 'C:/Users/gskra/Documents/NA_5FU-2/RTrans.parameters for NA_5FU-3.xlsx'
# sheet.name = 'GO terms (inDetails)'

Read.DB.entries.from.Excel__MultiSet = function(file.name, sheet.name){
  #file.name = 'C:/Users/gskra/OneDrive/Moskalev flies-2016/RTrans.parameters.xlsx'
  if (!('readxl' %in% rownames(installed.packages()))) install.packages("readxl")
  if (substring(file.name,nchar(file.name)-4,nchar(file.name)) != '.xlsx') stop('Please fill up input parameters in Excel workbook (see an example in ../RTrans_main_scripts/Parameters.example.xlsx)')
  if (!file.exists(file.name)) stop(sprintf('file %s does not exist',file.name))
  
  # reading worksheet
  tryCatch(expr = {   DB.entry.table = suppressMessages(suppressWarnings(readxl::read_excel(file.name,sheet=sheet.name, col_names = TRUE)))
  }, error = function (err){
    stop(sprintf('Sheet "%s" is not found in workbook %s', sheet.name, file.name))
  })
  
  DB.entry.table = as.data.frame(DB.entry.table)
  DB.entry.table = DB.entry.table[, apply(DB.entry.table, 2, function(col) (!all(is.na(col)))), drop = FALSE]
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

Write.Density.table = function(cpms, tsv.file.name = 'CPM density, non-normalized.tsv', n_bins = 512, log = TRUE, from = -4, to = +10){
  density.table = Get.Density.table(cpms, n_bins, log, from, to)
  write.table(density.table, tsv.file.name, sep= '\t', quote = FALSE)
}

Get.Density.table = function(cpms, n_bins = 512, log = TRUE, from = -4, to = +10){
  if(log) cpms = log(cpms)
  density.table = data.frame(array(dim = c(ncol(cpms), n_bins)))
  first.sample.density = density(cpms[,1], n = n_bins, from = from, to = to)
  colnames(density.table) = first.sample.density$x
  rownames(density.table) = colnames(cpms)
  for(s in colnames(cpms)){
    density.table[s,] = density(cpms[,s], n = n_bins, from = from, to = to)$y
  }
  return(density.table)
}

Get.Density.table.avg.StDev = function(cpms, n_bins = 512, log = TRUE, from = -4, to = +10){
  density.table = Get.Density.table(cpms, n_bins, log, from, to)
  return(mean(apply(density.table, 2, sd)))
}

Render.Density.Plots__color.by.condition = function(cpms, cl, out.dir = NULL, plot.title = 'Log2(CPM) density, non-normalized', plot.subtitle = '',
  png.file.name = 'CPM density, non-normalized, color by condition.png'){
  if(is.null(out.dir)) out.dir = getwd()
  
  png.mod(filename = sprintf('%s/%s', out.dir, png.file.name),
      units="in",width=8.23,height=6.3,pointsize=12,
      res=600)
  plot(density(log(cpms[,1]), n = 512, from = -4, to = +10),col=cl[1],xlim=c(-4,10), main = plot.title, sub = plot.subtitle)
  for (x in 2:dim(cpms)[2]){
    lines(density(log(cpms[,x]), n = 512, from = -4, to = +10),col=cl[x])
  }
  dev.off()
}

Render.Density.Plots__color.by.sample = function(cpms, out.dir = NULL, plot.title = 'Log2(CPM) density, non-normalized', plot.subtitle = '',
  png.file.name = 'CPM density, non-normalized, color by sample'){
  if(is.null(out.dir)) out.dir = getwd()

  pallete = c("#F46D43", "#66C2A5", "#cd8845", "#3288BD", "#a8bf32", "#5E4FA2", "#D53E4F", "#d6d639", "#8ed384", "#9E0142", "#ebba2f")
  density.cols = colorRampPalette(pallete)(dim(cpms)[2])
  png.mod(filename = sprintf('%s/%s', out.dir, png.file.name),
          units="in",width=8.23,height=6.3,pointsize=12,
          res=600)
  lty.n = 1
  lty.count = 4
  plot(density(log(cpms[,1])),col=density.cols[1], xlim=c(-4,10), main = plot.title, lty = lty.n, sub = plot.subtitle)
  for (x in 2:dim(cpms)[2]){
    lty.n = (lty.n +1) %% lty.count 
    lines(density(log(cpms[,x])),col=density.cols[x], lty = lty.n)
  }
  
  legend.x.pos = (par("usr")[2] - par("usr")[1])*0.73 + par("usr")[1]
  legend.y.pos = (par("usr")[4] - par("usr")[3])*0.93 + par("usr")[3]
  legend(legend.x.pos, legend.y.pos, legend = colnames(cpms), col = density.cols, cex = 0.8, lty = 1:4)
  dev.off()
}

Normalize..Using.Reference.Genes = function(d, Reference.genes){
  raw.counts.table = d$counts
  absent.Reference.genes = Reference.genes[!(Reference.genes %in% rownames(raw.counts.table))]
  if(length(absent.Reference.genes) > 0)
    stop(sprintf('%d of %d are not found within expression dataset %d genes', length(absent.Reference.genes), length(Reference.genes), nrow(raw.counts.table)))
  
  RG.norm.factors = colSums(raw.counts.table[Reference.genes, , drop = FALSE])
  RG.norm.factors = RG.norm.factors / mean(RG.norm.factors)
  cat('\n Reference genes-based normalization factors:\n')
  print(RG.norm.factors)
  
  norm.CPMs.table = t(t(raw.counts.table) / RG.norm.factors)
  norm.CPMs.table = round(norm.CPMs.table)

  cs = colSums(norm.CPMs.table)
  dummy.row = t(as.data.frame(max(cs) - cs))
  rownames(dummy.row) = 'dummy'
  norm.CPMs.table = rbind(norm.CPMs.table, dummy.row)
  d = DGEList(counts = norm.CPMs.table)
  d = calcNormFactors(d, method = 'none') #,"TMM","upperquartile","none"
  d$counts = d$counts[!(rownames(d$counts) %in% 'dummy'),,drop = FALSE]
  return(d)
}

RTrans..pipeline = function(GLM.model.IDs = 'auto', Startup.Data = NULL, Startup.Data..RDS.file = NULL, Parameters.xlsx.file.name = NULL,
  script.dir = NULL, base.dir = NULL, current.dir = NULL, write.CPM.tables = TRUE,
  successful.steps.file = 'Successful.steps.txt', failed.steps.file = 'Failed.steps.txt', do_Summarize = TRUE,
  force.re.install.packages = FALSE, emit.package.message = TRUE){

  if(!is.null(base.dir))     if(endsWith(base.dir, '/'))     base.dir = substr(base.dir, 1, nchar(base.dir) - 1)
  if(!is.null(script.dir))   if(endsWith(script.dir, '/'))   script.dir = substr(script.dir, 1, nchar(script.dir) - 1)
  if(!is.null(current.dir))  if(endsWith(current.dir, '/'))  current.dir = substr(current.dir, 1, nchar(current.dir) - 1)

  if(is.null(current.dir)) current.dir = getwd()

  if(!is.null(Startup.Data..RDS.file) & !is.null(Startup.Data))
    stop('Both arguments Startup.Data..RDS.file and Startup.Data should not be specified. Only one is allowed')

  if(!is.null(Startup.Data..RDS.file)){
    tryCatch(expr = {
      Pre.load.RTrans.packages(Parameters.xlsx.file.name, script.dir = script.dir, base.dir = base.dir, current.dir = current.dir, force.re.install = force.re.install.packages, emit.message = emit.package.message)
      Startup.Data <<- readRDS(Startup.Data..RDS.file)
      if(length(GLM.model.IDs) == 0 & do_Summarize){
        msg = sprintf('[summarizing] Pre.load.RTrans.packages(...) successful')
      } else
        msg = sprintf('[models %s] Pre.load.RTrans.packages(...) successful', toString(GLM.model.IDs))
      cat(msg); cat('\n')
      write(msg, successful.steps.file, append = TRUE)
    }, error = function (err){
      if(length(GLM.model.IDs) == 0 & do_Summarize){
        msg = sprintf('[summarizing] Pre.load.RTrans.packages(...) FAILED with error: %s', err)
      } else
        msg = sprintf('[models %s] Pre.load.RTrans.packages(...) FAILED with error: %s', toString(GLM.model.IDs), err)
      write(msg, failed.steps.file, append = TRUE)
      cat(msg); cat('\n')
      stop(err)
    })
  } else if(is.null(Startup.Data)){
    tryCatch(expr = {
      Startup.Data <<- Prepare(Parameters.xlsx.file.name = Parameters.xlsx.file.name, write.CPM.tables = write.CPM.tables, script.dir = script.dir, base.dir = base.dir, current.dir = current.dir, emit.package.message = emit.package.message)
      if(length(GLM.model.IDs) == 0 & do_Summarize){
        msg = sprintf('[summarizing] Prepare(...) successful')
      } else
        msg = sprintf('[models %s] Prepare(...) successful', toString(GLM.model.IDs))
      cat(msg); cat('\n')
      write(msg, successful.steps.file, append = TRUE)
    }, error = function (err){
      if(length(GLM.model.IDs) == 0 & do_Summarize){
        msg = sprintf('[summarizing] Prepare(...) FAILED with error: %s', err)
      } else
        msg = sprintf('[models %s] Prepare(...) FAILED with error: %s', toString(GLM.model.IDs), err)
      write(msg, failed.steps.file, append = TRUE)
      cat(msg); cat('\n')
      stop(err)
    })
  }

  all.GLM.models = Startup.Data$Pars$Models.to.Test
  if(toString(GLM.model.IDs) == 'auto')
    GLM.model.IDs = 1:length(all.GLM.models)

  # GLM.model.IDs..numeric = c()
  # for(ID in GLM.model.IDs){
  #   if(as.character(ID) %in% all.GLM.models){
  #     GLM.model.IDs..numeric = c(GLM.model.IDs..numeric, which(all.GLM.models %in% as.character(ID)))
  #   }
  # }

  # GLM.model.IDs..letter = as.character(GLM.model.IDs)
  # GLM.model.IDs..letter %in%

  # GLM.model.IDs..letter = GLM.model.IDs[is.na(as.numeric(as.character(GLM.model.IDs)))]
  # GLM.model.IDs..numeric = as.numeric(as.character(GLM.model.IDs))[is.na(as.numeric(as.character(GLM.model.IDs)))]

  # do.Stop = FALSE
  # msg1 = ''
  # msg2 = ''

  if(any(GLM.model.IDs > length(all.GLM.models))){
    stop(sprintf('Incorrect GLM.model.IDs %s. The total number of models available is %d',
      toString(GLM.model.IDs..numeric[GLM.model.IDs > length(all.GLM.models)]), length(all.GLM.models)))
    # do.Stop = TRUE
  }

  # if(any(!(GLM.model.IDs..letter %in% all.GLM.models))){
  #   msg2 = sprintf('Unknown GLM.model.IDs %s. The available models: %s',
  #     toString(GLM.model.IDs..letter[!(GLM.model.IDs..letter %in% all.GLM.models)]), toString(all.GLM.models))
  #   do.Stop = TRUE
  # }

  # if(do.Stop){
  #   cat(sprintf('%s\n%s\n', msg1, msg2))
  #   stop(sprintf('%s\n%s\n', msg1, msg2))
  # }


  for(GLM.model.n in GLM.model.IDs){
    GLM.model = all.GLM.models[GLM.model.n]

    Individual.Pars = NULL
    if(GLM.model %in% colnames(Startup.Data$Ind.Par.table)){
      tmp = subset(Startup.Data$Ind.Par.table, select = GLM.model, subset = !is.na(Startup.Data$Ind.Par.table[,GLM.model]))
      par.table = cbind(as.data.frame(rownames(tmp)), tmp)
      colnames(par.table) = c('variable','value')
      cat(sprintf('Found %d individual parameters:\n', nrow(par.table)))
      for(x in 1:nrow(par.table)){
        cat(sprintf('      %s:    %s\n', par.table[x,1], par.table[x,2]))
      }
      Pars = Process.Parameters.table(Pars, par.table, show.warnings = FALSE, base.dir = base.dir, script.dir = script.dir, current.dir = current.dir)
    }
    
    log.prefix = sprintf('[model %d - "%s"] ', GLM.model.n, GLM.model)

    if(Pars$Perform.DE.analysis){
      tryCatch(expr = {
        Analyze.GLM(Startup.Data, GLM.model, script.dir = script.dir, Create.Excel.results = TRUE, verbose = FALSE)
        msg = sprintf('%s Analyze.GLM(...) successful', log.prefix)
        cat(msg); cat('\n')
        write(msg, successful.steps.file, append = TRUE)
        print(successful.steps.file)
      }, error = function (err){
        msg = sprintf('%s Analyze.GLM(...) FAILED with error: %s', log.prefix, err)
        cat(msg); cat('\n')
        write(msg, failed.steps.file, append = TRUE)
        stop(err)
      })
    }

    if(Pars$Create.heatmaps){
      tryCatch(expr = {
        Create.Heatmaps(Startup.Data, GLM.model = GLM.model)
        msg = sprintf('%s Create.Heatmaps(...) successful', log.prefix)
        cat(msg); cat('\n')
        write(msg, successful.steps.file, append = TRUE)
      }, error = function (err){
        msg = sprintf('%s Create.Heatmaps(...) FAILED with error: %s', log.prefix, err)
        write(msg, failed.steps.file, append = TRUE)
        if(Pars$Stop.if.an.error.occurs) { stop(msg)
        } else warning(msg)
      })
    }

    if(Pars$Perform.inDetails.analysis){
      tryCatch(expr = {
        inDetails(Startup.Data, GLM.model = GLM.model, verbose = FALSE, do.create.Heatmaps = Pars$Create.heatmaps.for.inDetails.analyses,
        	add.biotypes = Pars$inDetails.additional.biotypes)
        msg = sprintf('%s inDetails(...) successful', log.prefix)
        cat(msg); cat('\n')
        write(msg, successful.steps.file, append = TRUE)
      }, error = function (err){
        msg = sprintf('%s inDetails(...) FAILED with error: %s', log.prefix, err)
        write(msg, failed.steps.file, append = TRUE)
        if(Pars$Stop.if.an.error.occurs) { stop(msg)
        } else warning(msg)
      })
    }

    if(!Startup.Data$miR.mode){
      if(Pars$Perform.Disease.Ontology.classic.enrichment){
        tryCatch(expr = {
          Disease_Ont.classic.Enrichment(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s Disease_Ont.classic.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s Disease_Ont.classic.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Perform.Disease.Ontology.trends.enrichment){
        tryCatch(expr = {
          Disease_Ont.trends.Enrichment(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s Disease_Ont.trends.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s Disease_Ont.trends.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Perform.Network.of.Cancer.Genes.classic.enrichment){
        tryCatch(expr = {
          NCG.classic.Enrichment(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s NCG.classic.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s NCG.classic.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Perform.Network.of.Cancer.Genes.trends.enrichment){
        tryCatch(expr = {
          NCG.trends.Enrichment(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s NCG.trends.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s NCG.trends.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }
      
      if(Pars$Perform.DisGeNET.classic.enrichment){
        tryCatch(expr = {
          DisGeNET.classic.Enrichment(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s DisGeNET.classic.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s DisGeNET.classic.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Perform.DisGeNET.trends.enrichment){
        tryCatch(expr = {
          DisGeNET.trends.Enrichment(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s DisGeNET.trends.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s DisGeNET.trends.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Perform.WikiPathways.classic.enrichment){
        tryCatch(expr = {
          WikiPathways.classic.Enrichment(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s WikiPathways.classic.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s WikiPathways.classic.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Perform.KEGG.classic.enrichment){
        tryCatch(expr = {
          KEGG.classic.Enrichment(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s KEGG.classic.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s KEGG.classic.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Perform.KEGG.trends.enrichment){
        tryCatch(expr = {
          KEGG.trends.Enrichment(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s KEGG.trends.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s KEGG.trends.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Create.KEGG_centric.expression.profiles){
        tryCatch(expr = {
          KEGG.Expression.Profiles.Custom(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s KEGG.Expression.Profiles.Custom(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s KEGG.Expression.Profiles.Custom(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }
      
      if(Pars$Perform.KEGG.pathway.visualization){
        tryCatch(expr = {
          KEGG.Pathways.Visualization(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s KEGG.Pathways.Visualization(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s KEGG.Pathways.Visualization(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Perform.Reactome.classic.enrichment){
        tryCatch(expr = {
          Reactome.classic.Enrichment(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s Reactome.classic.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s Reactome.classic.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Perform.Reactome.trends.enrichment){
        tryCatch(expr = {
          Reactome.trends.Enrichment(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s Reactome.trends.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s Reactome.trends.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Create.Reactome_centric.expression.profiles){
        tryCatch(expr = {
          Reactome.Expression.Profiles.Custom(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s Reactome.Expression.Profiles.Custom(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s Reactome.Expression.Profiles.Custom(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }
      

      if(Pars$Perform.topGO.enrichment){
        tryCatch(expr = {
          topGO.Enrichment(Startup.Data, GLM.model = GLM.model, GO.types = Pars$topGO.enrichment.namespaces)
          msg = sprintf('%s topGO.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s topGO.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Create.topGO_centric.expression.profiles){
        tryCatch(expr = {
          topGO.Expression.Profiles.Custom(Startup.Data, GLM.model = GLM.model, GO.type = Pars$topGO_centric.expression.profiles__namespace)
          msg = sprintf('%s topGO.Expression.Profiles.Custom(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s topGO.Expression.Profiles.Custom(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Create.topGO_centric.expression.profiles.for.enriched.terms){
        tryCatch(expr = {
          topGO.Expression.Profiles.Enriched(Startup.Data, GLM.model = GLM.model, GO.type = Pars$topGO_centric.expression.profiles.for.enriched.terms__namespace)
          msg = sprintf('%s topGO.Expression.Profiles.Enriched(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s topGO.Expression.Profiles.Enriched(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Perform.trends.GO.enrichment){
        tryCatch(expr = {
          GO.trends.Enrichment(Startup.Data, GLM.model = GLM.model, GO.types = c('BP','CC','MF'))
          msg = sprintf('%s GO.trends.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s GO.trends.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Perform.classic.GO.enrichment){
        tryCatch(expr = {
          GO.classic.Enrichment(Startup.Data, GLM.model = GLM.model, GO.types = c('BP','CC','MF'))
          msg = sprintf('%s GO.classic.Enrichment(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s GO.classic.Enrichment(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }

      if(Pars$Create.GO_centric.expression.profiles){
        tryCatch(expr = {
          GO.Expression.Profiles.Custom(Startup.Data, GLM.model = GLM.model)
          msg = sprintf('%s GO.Expression.Profiles.Custom(...) successful', log.prefix)
          cat(msg); cat('\n')
          write(msg, successful.steps.file, append = TRUE)
        }, error = function (err){
          msg = sprintf('%s GO.Expression.Profiles.Custom(...) FAILED with error: %s', log.prefix, err)
          write(msg, failed.steps.file, append = TRUE)
          if(Pars$Stop.if.an.error.occurs) { stop(msg)
          } else warning(msg)
        })
      }
    }
  }

  if(do_Summarize & Pars$Create.summary.reports){
    log.prefix = sprintf('[summarizing] ')
    tryCatch(expr = {
      Summarize.results.by.Analyses.groups(Startup.Data, joint.GLM.results = TRUE, include.with.spaces = TRUE)
      msg = sprintf('%s Summarize.results.by.Analyses.groups(...) successful', log.prefix)
      cat(msg); cat('\n')
      write(msg, successful.steps.file, append = TRUE)
    }, error = function (err){
      msg = sprintf('%s Summarize.results.by.Analyses.groups(...) FAILED with error: %s', log.prefix, err)
      write(msg, failed.steps.file, append = TRUE)
      if(Pars$Stop.if.an.error.occurs) { stop(msg)
      } else warning(msg)
    })
  }
}


Check.Python.availability = function(){
  if(suppressWarnings(system('python3 -c "import xlsxwriter, numpy"',ignore.stdout = T, ignore.stderr = T)) == 0){
    python.bin = "python3"
  } else if (suppressWarnings(system('python -c "import xlsxwriter, numpy"',ignore.stdout = T, ignore.stderr = T)) == 0){
    python.bin = "python"
  } else {
    stop('Python3 with xlsxwriter and numpy packages is required for creating Excel reports. If python is already installed, open command line and type"python -m pip install xlsxwriter numpy" (or "python3 -m pip install xlsxwriter numpy")')
    #Pars$Create.Excel.results = FALSE
  }
  python.bin
}    



Pre.load.RTrans.packages = function(Parameters.xlsx.file.name = NULL, script.dir = NULL, base.dir = NULL, current.dir = NULL, force.re.install = FALSE, emit.message = TRUE){
  
  if(!is.null(base.dir))     if(endsWith(base.dir, '/'))     base.dir = substr(base.dir, 1, nchar(base.dir) - 1)
  if(!is.null(script.dir))   if(endsWith(script.dir, '/'))   script.dir = substr(script.dir, 1, nchar(script.dir) - 1)
  if(!is.null(current.dir))  if(endsWith(current.dir, '/'))  current.dir = substr(current.dir, 1, nchar(current.dir) - 1)

  if(is.null(current.dir)) current.dir = getwd()
  if(is.null(script.dir)){
    tryCatch(expr = {
      script.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
    }, error = function (err){
      script.dir = getwd()
      warning(sprintf('rstudioapi failed to locate script directory. Assuming scrpit dir as "%s"', script.dir))
    })
  }

  #ext.DE.data.files=c('ext.DE.data.example.txt')
  if (is.null(Parameters.xlsx.file.name)){
    Parameters.xlsx.file.name = sprintf("%s/RTrans.parameters.xlsx", script.dir)
    if (file.exists(Parameters.xlsx.file.name)){
      cat(sprintf('\nFile with parameters found: %s\n',Parameters.xlsx.file.name))
    } else stop(sprintf('Parameters.xlsx.file.name argument is not specified. No parameters *.xlsx file has been found in the default location '))
  }
  if (!is.null(script.dir))  Parameters.xlsx.file.name = gsub('\\{script_dir\\}',script.dir,Parameters.xlsx.file.name)
  # if (!is.null(script.dir) & !is.null(forced.maRt.table)) forced.maRt.table = gsub('\\{script_dir\\}',script.dir,forced.maRt.table)
  
  Pars = Read.Parameters(Parameters.xlsx.file.name, prepare.for.ext.DE.data = FALSE, script.dir = script.dir, show.warnings..on.parameters = FALSE, base.dir = base.dir, current.dir = current.dir)
  DB.data = Read.Organisms.Info(Organisms.info.file.name = sprintf('%s/Organisms.info.xlsx', Pars$suppl.data.dir), script.dir = script.dir)
  
  Pars$Species = tolower(Pars$Species)
  if (!(Pars$Species %in% DB.data$KEGG.codes.by.alias)){
    if(Pars$Species %in% names(DB.data$KEGG.codes.by.alias)){
      Pars$Species = DB.data$KEGG.codes.by.alias[Pars$Species]
    } else {
      msg = sprintf('Organism "%s" is not found among available taxons and KEGG abbreviations. Standard enrichment analyses will be disabled.', Pars$Species)
      cat(msg)
      warning(msg)
      Pars$Species = 'other'
    }
  }
  
  if(Pars$Species %in% names(DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
    org.DB.name = DB.data$Bioconductor.Org.DBs.by.KEGG.codes[[Pars$Species]]
  }

  pathviewmod.Installed.successfully = Load.RTrans.packages__internal(Pars, org.DB.name, force.re.install, emit.message = emit.message)
}


Load.RTrans.packages__internal = function(Pars, org.DB.name, force.re.install, emit.message = TRUE){
  if(emit.message)  cat('\nLoading/installing packages...\n')

  std.pkg.list = c('gplots','chemometrics','ggdendro','digest',"lme4","readxl","dendextend",'ggplot2','dplyr','stringr','forcats', 'pheatmap',
    'magrittr', 'utils', 'RColorBrewer', 'foreach', 'doParallel', 'stringr', 'parallel', 'funr', 'pryr', 'fs', 'eulerr', 'venn')
  
  bio.pkg.list.fp = c()
  if(!is.null(org.DB.name))  bio.pkg.list.fp = org.DB.name
  bio.pkg.list.fp = c(bio.pkg.list.fp, 'GO.db', "ReactomePA")
  bio.pkg.list.sp = c("Rgraphviz","edgeR","DESeq2","FactoMineR","biomaRt","topGO","clusterProfiler", 'GSEABase', "DOSE", "reactome.db",
    'AnnotationDbi', 'GOSemSim', 'WGCNA', 'goseq')
  bio.pkg.list = c(bio.pkg.list.fp, bio.pkg.list.sp)

  
  install.packages.std(std.pkg.list)
  install.packages.bio(bio.pkg.list)

  # if(force.re.install){
  #   install.packages(std.pkg.list, dep = TRUE)

  #   if(Check..R.version..gt.3.5()){
  #     if(!('BiocManager' %in% rownames(installed.packages())))
  #       install.packages('BiocManager')
  #     BiocManager::install(bio.pkg.list, dep = TRUE)
  #   } else {
  #     source("https://bioconductor.org/biocLite.R")
  #     biocLite(bio.pkg.list, dep = TRUE)
  #   }
  #   update.packages(ask = FALSE)
  # } else {
  #   std.pkg.list.to.install = std.pkg.list[!(std.pkg.list %in% installed.packages())]
  #   if(length(std.pkg.list.to.install > 0)) install.packages(std.pkg.list.to.install)
    
  #   bio.pkg.list.to.install = bio.pkg.list[!(bio.pkg.list %in% installed.packages())]
  #   if(length(bio.pkg.list.to.install > 0)){
  #     if(Check..R.version..gt.3.5()){
  #       if(!('BiocManager' %in% rownames(installed.packages())))
  #         install.packages('BiocManager')
  #       BiocManager::install(bio.pkg.list)
  #     } else {
  #       source("https://bioconductor.org/biocLite.R")
  #       biocLite(bio.pkg.list.to.install)
  #     }
  #   }
  # }
  
  for(lib in bio.pkg.list.fp)  eval(parse(text = sprintf('suppressPackageStartupMessages(library(%s))',lib)))
  for(lib in c(std.pkg.list, bio.pkg.list.sp))  eval(parse(text = sprintf('suppressPackageStartupMessages(library(%s))',lib)))

  # return()
  
  if (Pars$Pathway.Visualization___CPM.aware){
    pathviewmod.Installed.successfully = TRUE
    if (!('pathviewmod' %in% rownames(installed.packages()))){
      pathviewmod.Installed.successfully = FALSE
      tryCatch(expr = {
        install.packages(sprintf('%s/pathview_mod',Pars$suppl.data.dir), repos = NULL, type="source")
        pathviewmod.Installed.successfully <<- TRUE
      }, error = function (err){
        pathviewmod.Installed.successfully <<- FALSE
      })
    }
    
    tryCatch(expr = { suppressPackageStartupMessages(library('pathviewmod'))
    }, error = function (err){
      pathviewmod.Installed.successfully <<- FALSE
    })
    
    if (!pathviewmod.Installed.successfully | !('pathviewmod' %in% rownames(installed.packages()))){
      message('Warning. RTrans is provided with a modified, CPM-aware version of "pathview". It gives much better results for the analysis of massive gene expression profiles changes, e.g. age-associated. RTrans was unable to install the modified package, "pathviewmod" to your system. CPM-aware pathway visualization will be disabled.')
      Pars$Pathway.Visualization___CPM.aware = FALSE
    }
  }
  
  if (!Pars$Pathway.Visualization___CPM.aware){
    install.packages.bio(c("pathview"))
    suppressPackageStartupMessages(library('pathview'))
  }
  return(pathviewmod.Installed.successfully)
}
# 
# Parameters.xlsx.file.name = '/mnt/raid/illumina/geo/RStudioServer/SRP155256/RTrans.parameters.xlsx'
# forced.maRt.table = NULL
# force.re.install = FALSE
# disable.excel = FALSE
# ext.DE.data.files = NULL
# 
# ext.DE.data.species = NULL
# ext.DE.data.PValue.if.absent = 0.01
# ext.DE.data.logCPM.if.absent = 5.0
# ext.DE.data.Score.if.absent = 'auto'
# 
# write.CPM.tables = TRUE
# script.dir = '/mnt/raid/illumina/geo/progs/RTrans/RTrans_base/'
# base.dir = '/mnt/raid/illumina/geo/progs/RTrans/RTrans_base/'
# current.dir = '/mnt/raid/illumina/geo/RStudioServer/SRP155256'
# emit.package.message = TRUE

Prepare = function(Parameters.xlsx.file.name = NULL, forced.maRt.table = NULL, force.re.install = FALSE, disable.excel = FALSE, ext.DE.data.files = NULL,
                   ext.DE.data.species = NULL, ext.DE.data.PValue.if.absent = 0.01, ext.DE.data.logCPM.if.absent = 5.0, ext.DE.data.Score.if.absent = 'auto',
                   write.CPM.tables = TRUE, script.dir = NULL, base.dir = NULL, current.dir = NULL, emit.package.message = TRUE){
  #Parameters.xlsx.file.name  ='C:/Users/gskra/OneDrive/Evgeniev/RTrans.parameters.xlsx'
  if(!is.null(base.dir))     if(endsWith(base.dir, '/'))     base.dir = substr(base.dir, 1, nchar(base.dir) - 1)
  if(!is.null(script.dir))   if(endsWith(script.dir, '/'))   script.dir = substr(script.dir, 1, nchar(script.dir) - 1)
  if(!is.null(current.dir))  if(endsWith(current.dir, '/'))  current.dir = substr(current.dir, 1, nchar(current.dir) - 1)

  if (!('rstudioapi' %in% rownames(installed.packages()))) install.packages("rstudioapi")
  if (!('magrittr' %in% rownames(installed.packages()))) install.packages("magrittr")
  library(magrittr)

  if(is.null(current.dir)) current.dir = getwd()

  if(is.null(script.dir)){
    tryCatch(expr = {
      script.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
    }, error = function (err){
      script.dir = getwd()
      warning(sprintf('rstudioapi failed to locate script directory. Assuming scrpit dir as "%s"', script.dir))
    })
  }

  #ext.DE.data.files=c('ext.DE.data.example.txt')
  if (is.null(Parameters.xlsx.file.name)){
    Parameters.xlsx.file.name = sprintf("%s/RTrans.parameters.xlsx", script.dir)
    if (file.exists(Parameters.xlsx.file.name)){
      cat(sprintf('\nFile with parameters found: %s\n',Parameters.xlsx.file.name))
    } else stop(sprintf('Parameters.xlsx.file.name argument is not specified. No parameters *.xlsx file has been found in the default location '))
  }
  if (!is.null(script.dir))  Parameters.xlsx.file.name = gsub('\\{script_dir\\}',script.dir,Parameters.xlsx.file.name)
  if (!is.null(script.dir) & !is.null(forced.maRt.table)) forced.maRt.table = gsub('\\{script_dir\\}',script.dir,forced.maRt.table)
  
  #Startup.Data = new('RTransData')
  Pars = Read.Parameters(Parameters.xlsx.file.name, !is.null(ext.DE.data.files), script.dir = script.dir, base.dir = base.dir, current.dir =current.dir)
  if(!is.null(ext.DE.data.species)){
    Pars$Species = ext.DE.data.species
  }
  
  if(Pars$Cluster.samples.in.Excel.results < Pars$Create.heatmaps..with.Clustered.samples){
    msg = sprintf('You turned %s Cluster.samples.in.Excel.results but turned %s Create.heatmaps..with.Clustered.samples. Clustered heatmaps cannot be created if the sample clustering is turned FALSE. It will be turned TRUE\n',
                  toString(Pars$Cluster.samples.in.Excel.results), toString(Pars$Create.heatmaps..with.Clustered.samples))
    cat(msg)
    Pars$Cluster.samples.in.Excel.results = TRUE
  }

  DB.data = Read.Organisms.Info(Organisms.info.file.name = sprintf('%s/Organisms.info.xlsx', Pars$suppl.data.dir), script.dir = script.dir)
  
  Startup.Data = new('RTransStartupData')
  Startup.Data$DB.data = DB.data
  
  Pars$Species = tolower(Pars$Species)
  if (!(Pars$Species %in% DB.data$KEGG.codes.by.alias)){
    if(Pars$Species %in% names(DB.data$KEGG.codes.by.alias)){
      Pars$Species = DB.data$KEGG.codes.by.alias[Pars$Species]
    } else {
      msg = sprintf('Organism "%s" is not found among available taxons and KEGG abbreviations. Standard enrichment analyses will be disabled.', Pars$Species)
      cat(msg)
      warning(msg)
      Pars$Species = 'other'
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

  pathviewmod.Installed.successfully = Load.RTrans.packages__internal(Pars, Startup.Data$org.DB.name, force.re.install, emit.message = emit.package.message)
  
  cat('\nChecking pre-processing input expression data...\n')
  
  Pars$RefSeq.Info.File = sprintf('%s/%s.Ensembl.Genes.Info.txt',Pars$suppl.data.dir,Pars$Species)
  Pars$Pathway.names.table.File = sprintf('%s/KEGG.pathway.names.tsv',Pars$suppl.data.dir)
  Pars$GO.full.descriptions.File = sprintf('%s/Gene Ontology terms names and descriptions.tsv',Pars$suppl.data.dir)
  
  
  if (Pars$Use.Info.Table & file.exists(Pars$RefSeq.Info.File)) {
    Startup.Data$Info.table = read.table(Pars$RefSeq.Info.File,stringsAsFactors = FALSE, blank.lines.skip = FALSE,sep = '\t')
  } else {
    Pars$Use.Info.Table = FALSE
    Startup.Data$Info.table = NA
  }
  
  ## Loading individual parameters (if a sheet is present)
  Ind.Par.table = NULL
  all.sheets = readxl::excel_sheets(Parameters.xlsx.file.name)
  all.Ind.Par.sheets = c('Individual parameters', 'Ind. parameters', 'Ind. par.', 'Overriding parameters', 'overriding parameters', 'Forced parameters', 'forced parameters', 'OVP')
  if(sum.mod(all.sheets %in% all.Ind.Par.sheets) > 1){
    stop('Ambigous forced parameters sheet name')
  } else if(sum.mod(all.sheets %in% all.Ind.Par.sheets) == 1){
    Ind.Par.sheet = all.sheets[all.sheets %in% all.Ind.Par.sheets][1]
    Ind.Par.table = suppressMessages(suppressWarnings(readxl::read_excel(Parameters.xlsx.file.name, sheet = Ind.Par.sheet, col_names = TRUE)))
    Ind.Par.table %<>% as.data.frame
    Ind.Par.table = subset(Ind.Par.table, subset = !is.na(Ind.Par.table[,1]))
    if(any(duplicated(Ind.Par.table[,1]))){
      stop(sprintf('Duplicated entries are found in %s sheet: %s', Ind.Par.sheet, toString(Ind.Par.table[,1][duplicated(Ind.Par.table[,1])])))
    }
    Ind.Par.table = Ind.Par.table[!startsWith(Ind.Par.table[,1], '#'),, drop = FALSE]

    Ind.Par.table[,1] = sapply(Ind.Par.table[,1],function(x) gsub('"',"",gsub('"','',x)))
    Ind.Par.table[,1] = sapply(Ind.Par.table[,1],function(x) gsub("'","",gsub('"','',x)))
    Ind.Par.table[,1] = sapply(Ind.Par.table[,1],function(x) gsub(' ', ".", x, fixed = TRUE))
    Ind.Par.table[,1] = sapply(Ind.Par.table[,1],function(x) gsub('-', "_", x, fixed = TRUE))
    rownames(Ind.Par.table) = Ind.Par.table[,1]

    if(dim(Ind.Par.table)[2] == 2){    Ind.Par.table = subset(Ind.Par.table, select = c(FALSE, TRUE))
    } else  Ind.Par.table = Ind.Par.table[,-1]
    
    unknown.OVPs = rownames(Ind.Par.table)[!(rownames(Ind.Par.table) %in% names(Pars))]
    if(length(unknown.OVPs) > 0)   stop(sprintf('Unknown individual parameters are found: %s', toString(unknown.OVPs)))
    
    commented.cols = sapply(colnames(Ind.Par.table), FUN = function (x) substring(x,first=1,last=1)) %in% '#'
    Ind.Par.table = subset(Ind.Par.table, select = !commented.cols)
    commented.lines = sapply(Ind.Par.table[,1], FUN = function (x) substring(x,first=1,last=1)) %in% '#'
    Ind.Par.table = subset(Ind.Par.table, subset = !commented.lines)
    Ind.Par.table = subset(Ind.Par.table, select = !apply(Ind.Par.table, 2, function(x) all(is.na(x))))
    cat(sprintf('There are %d individual parameters for %d models: %s\n', nrow(Ind.Par.table), ncol(Ind.Par.table), toString(colnames(Ind.Par.table))))
  }


  WGCNA.schema = tryCatch(expr = { suppressMessages(suppressWarnings(readxl::read_excel(Parameters.xlsx.file.name, sheet='WGCNA groups',col_names = TRUE)))
  }, error = function (err){
    cat(sprintf('\nSheet "WGCNA groups" is not found in workbook %s. WGCNA analysis will not be performed\n\n',Parameters.xlsx.file.name))
    return(NULL)
  })
  
  if(!is.null(WGCNA.schema)){
    WGCNA.schema %<>% as.data.frame
    commented.cols = sapply(colnames(WGCNA.schema), FUN = function (x) substring(x,first=1,last=1)) %in% '#'
    WGCNA.schema = WGCNA.schema[,!commented.cols]
    commented.lines = sapply(WGCNA.schema[,1], FUN = function (x) substring(x,first=1,last=1)) %in% '#'
    WGCNA.schema = WGCNA.schema[!commented.lines,]
    WGCNA.schema = WGCNA.schema[,!apply(WGCNA.schema, 2, function(x) all(is.na(x)))]
    if (!("Sample names" %in% colnames(WGCNA.schema))) stop(sprintf('Column "Sample names" is not found in the sheet "WGCNA groups", workbook %s', Parameters.xlsx.file.name))
    rownames(WGCNA.schema) = WGCNA.schema$'Sample names'
  }

  Heatmaps.annotation = NULL
  if('Heatmaps annotation' %in% readxl::excel_sheets(Parameters.xlsx.file.name))
    Heatmaps.annotation = Read.phys.chem.data(Parameters.xlsx.file.name, sheet.name = 'Heatmaps annotation')


  Reference.genes = Read.Reference.Genes(Parameters.xlsx.file.name)
  if(length(Reference.genes) > 0)
    cat(sprintf('%d reference gene(s) are specified\n', length(Reference.genes)))
  
  if(length(Reference.genes) == 0 & Pars$RNA.Seq.norm.method == 'RG')
    stop('Reference genes must be specified in the "Reference genes" worksheet in order to perform RG-based normalization')


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
    cat(sprintf('%d Reactome pathways with unexpected species ID ("%s") have been detected. Species ID will be replaced with "%s"\n', R.Sp.incorrect, paste(levels(as.factor(R.Sp[Pars$Species != R.Sp])), collapse = ', '), Pars$Species))
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
  
  
  inDetails.entries = Read.DB.entries.from.Excel__MultiSet(Parameters.xlsx.file.name, 'inDetails')
  if(length(inDetails.entries) > 0){
    ## tryiung to assigned genes by GO ID
    inDetails.entries = lapply(inDetails.entries, function(x){
      Convert.GO.term.names_to_GO.term.IDs(terms = x, GO.term.ID.by.name = GO.term.ID.by.name,
                                           remove.redundant = Pars$topGO.Expression.Profiles___remove.redundant.GO.terms,
                                           prefix = '[GO terms (inDetails)]: ', messages = FALSE)
    })
    rm(GO.term.ID.by.name)
    Startup.Data$inDetails.entries = inDetails.entries

    # print('--------------')
    # print(inDetails.entries)
    
    res = Get.Gene.Ensembl.IDs__for__GO.term.IDs.pack(inDetails.entries, Pars$Species, DB.data)
    Startup.Data$inDetails.GO.genes = res[['genes.vector']]

    # print('--------------')
    # print(res[['not.found.terms.vector']])

    # print('--------------')
    ## tryiung to assigned genes by KEGG ID
    res = Get.Gene.Ensembl.IDs__for__KEGG.pathway.IDs.pack(res[['not.found.terms.vector']], Pars$Species, DB.data, Pars$suppl.data.dir)
    Startup.Data$inDetails.KEGG.genes = res[['genes.vector']]

    ## tryiung to assigned genes by Reactome ID
    # res = Get.Gene.Ensembl.IDs__for__Reactome.pathway.IDs.pack(res[['not.found.pathway.list']], Pars$Species, DB.data)
    # Startup.Data$inDetails.Reactome.genes = res[['genes.vector']]

    Startup.Data$inDetails.Standalone.genes = res[['not.found.pathway.list']]
  } else {
    Startup.Data$inDetails.entries = vector(mode = 'list')
    Startup.Data$inDetails.GO.genes = vector(mode = 'list')
    Startup.Data$inDetails.KEGG.genes = vector(mode = 'list')
    Startup.Data$inDetails.Reactome.genes = vector(mode = 'list')
    Startup.Data$inDetails.Standalone.genes = vector(mode = 'list')
  }

  # Startup.Data$GO.terms.inDetails = GO.terms.inDetails
  

  Completed.steps.file = sprintf("%s/Completed.steps.list",Pars$results.dir)
  
  setwd(Pars$results.dir)
  
  ##### Pre-creating complete gene list and retrieving biomaRt infos for these genes
  Startup.Data$use.ext.DE.data = FALSE
  
  
  if(is.null(ext.DE.data.files)){
    schema = tryCatch(expr = { suppressMessages(suppressWarnings(readxl::read_excel(Parameters.xlsx.file.name,sheet='Sample setup',col_names = TRUE)))
    }, error = function (err){
      stop(sprintf('Sheet "Sample setup" is not found in workbook %s',Parameters.xlsx.file.name))
    })
    schema %<>% as.data.frame
    commented.cols = sapply(colnames(schema), FUN = function (x) substring(x,first=1,last=1)) %in% '#'
    schema = schema[,!commented.cols]
    commented.lines = sapply(schema[,1], FUN = function (x) substring(x,first=1,last=1)) %in% '#'
    schema = schema[!commented.lines,]
    schema = schema[,!apply(schema, 2, function(x) all(is.na(x)))]
    if (!("Sample names" %in% colnames(schema))) stop(sprintf('Column "Sample names" is not found in the workbook %s',Parameters.xlsx.file.name))
    rownames(schema) = schema$'Sample names'

    all.src..sample.names = schema$'Sample names'
    all.split..sample.names = all.src..sample.names
    if(sum(duplicated(all.src..sample.names)) > 0){
      stop(sprintf('Duplicated sample names are not allowed: %s', all.src..sample.names[duplicated(all.src..sample.names)]))
    }
    
    if(!is.null(Heatmaps.annotation)){
      mSamples = rownames(schema)
      wSamples = rownames(Heatmaps.annotation)

      present.only.in.main.Schema = mSamples[!(mSamples %in% wSamples)]
      present.only.in.Heatmaps.annotation = wSamples[!(wSamples %in% mSamples)]
      present.in.both.Schemas = intersect(mSamples, wSamples)
      present.in.any.Schema = union(mSamples, wSamples)
      if(length(present.in.both.Schemas) == 0){
        msg = sprintf('Warning! WGCNA Heatmaps annotation ("Heatmaps annotation" sheet) completely does not match the main schema ("Sample setup" sheet) in file %s. Custom heatmaps annotation will be bypassed', Parameters.xlsx.file.name)
        cat(msg)
        cat('\n')
        warning(msg)
        Heatmaps.annotation = NULL
      } else{
        if(length(present.only.in.main.Schema) > 0){
          msg = sprintf('%d of %d samples are absent in "Heatmaps annotation" sheet: %s\n',
            length(present.only.in.main.Schema), length(present.in.any.Schema), length(present.only.in.main.Schema))
          stop(msg)
        }
        if(length(present.only.in.Heatmaps.annotation) > 0)
          cat(sprintf('%d of %d samples present only in "Heatmaps annotation" sheet but are absent in "Sample setup" sheet. They will be discarded. Only %d common samples will be used for WGCNA\n',
            length(present.only.in.Heatmaps.annotation), length(present.in.any.Schema), length(present.in.both.Schemas)))
        if(length(present.only.in.Heatmaps.annotation) == 0 & length(present.only.in.main.Schema) == 0)
          cat(sprintf('"Heatmaps annotation" sheet and "Sample setup" sheet perfectly match to each other\n'))
      }
    }


    if(!is.null(WGCNA.schema)){
      mSamples = rownames(schema)
      wSamples = rownames(WGCNA.schema)

      present.only.in.main.Schema = mSamples[!(mSamples %in% wSamples)]
      present.only.in.WGCNA.Schema = wSamples[!(wSamples %in% mSamples)]
      present.in.both.Schemas = intersect(mSamples, wSamples)
      present.in.any.Schema = union(mSamples, wSamples)
      if(length(present.in.both.Schemas) == 0){
        msg = sprintf('Warning! WGCNA schema ("WGCNA groups" sheet) completely does not match the main schema ("Sample setup" sheet) in file %s. WGCNA analysis will be disabled', Parameters.xlsx.file.name)
        cat(msg)
        cat('\n')
        warning(msg)
        WGCNA.schema = NULL
      } else{
        if(length(present.only.in.main.Schema) > 0)
          cat(sprintf('%d of %d samples present only in main schema ("Sample setup" sheet) but are absent in "WGCNA groups" sheet. They will be discarded. Only %d common samples will be used for WGCNA\n',
            length(present.only.in.main.Schema), length(present.in.any.Schema), length(present.in.both.Schemas)))
        if(length(present.only.in.WGCNA.Schema) > 0)
          cat(sprintf('%d of %d samples present only in WGCNA schema ("WGCNA groups" sheet) but are absent in "Sample setup" sheet. They will be discarded. Only %d common samples will be used for WGCNA\n',
            length(present.only.in.WGCNA.Schema), length(present.in.any.Schema), length(present.in.both.Schemas)))
        if(length(present.only.in.WGCNA.Schema) == 0 & length(present.only.in.main.Schema) == 0)
          cat(sprintf('"WGCNA groups" sheet and "Sample setup" sheet perfectly match to each other\n'))
      }
    }
    
    #colnames(schema) = replace.spec.symbols(colnames(schema))
  
    #if (!("Sample names" %in% colnames(schema))) stop(sprintf('Column "Sample names" is not found in the workbook %s',Parameters.xlsx.file.name))
    
    if(is.null(WGCNA.schema))  all.src..sample.names..in.WGCNA = WGCNA.schema$'Sample names'

    some.samples..should.be.merged = FALSE
    some.samples..should.be.merged..in.WGCNA = FALSE
    if(Pars$Merge.samples.in.sample.setup){
      some.samples..should.be.merged = (length(grep(Pars$Merge.samples..separator, all.src..sample.names)) > 0)
      all.split..sample.names = all.src..sample.names
      if(some.samples..should.be.merged){
        for(x in 1:10)
          all.split..sample.names = gsub(sprintf('%s ', Pars$Merge.samples..separator), Pars$Merge.samples..separator, all.split..sample.names)
        all.split..sample.names = unlist(strsplit(all.split..sample.names, split = ',')) %>% DeDup.na.rm
      } 

      if(is.null(WGCNA.schema)){
        some.samples..should.be.merged..in.WGCNA = (length(grep(Pars$Merge.samples..separator, all.src..sample.names..in.WGCNA)) > 0)
        all.split..sample.names..in.WGCNA = all.src..sample.names..in.WGCNA
        if(some.samples..should.be.merged..in.WGCNA){
          for(x in 1:10)
            all.split..sample.names..in.WGCNA = gsub(sprintf('%s ', Pars$Merge.samples..separator), Pars$Merge.samples..separator, all.split..sample.names..in.WGCNA)
          all.split..sample.names..in.WGCNA = unlist(strsplit(all.split..sample.names..in.WGCNA, split = ',')) %>% DeDup.na.rm
        }
      }
    }
    
    if(Pars$counts.dir != ''){
      if(!dir.exists(Pars$counts.dir))
        stop(sprintf('Directory "%s" does not exist', Pars$counts.dir))
      

      file.names = sprintf("%s/%s%s", Pars$counts.dir, all.split..sample.names, Pars$counts.suffix)
      #file.names = sprintf("%s/%s%s", Pars$counts.dir, schema$'Sample names', Pars$counts.suffix)
      # print()
      if (any(!file.exists(file.names))){
        #warning(sprintf('Warning; The following files do not exist: %s',toString(file.names[!file.exists(file.names)])))
        msg = sprintf('\nTotal %d of %d *.counts files do not exist in the directory "%s": %s\n',
                    sum.mod(!file.exists(file.names)), length(file.names), Pars$counts.dir, toString(
                      sapply(strsplit(file.names[!file.exists(file.names)], split='/'), function(x) { tail(x,1)})
                      ))
        if(some.samples..should.be.merged)
          stop(sprintf('Sample merging is activated. All the samples must be present. %s', msg))

        cat(msg)
        if (sum.mod(file.exists(file.names)) == 0) stop(sprintf('No *%s files are found in the directory "%s". Please fix "Sample setup" sheet', Pars$counts.suffix, Pars$counts.dir))
        
        schema = schema[file.exists(file.names),]
        all.split..sample.names = all.split..sample.names[file.exists(file.names)]
        rownames(schema) = schema$'Sample names'
        file.names = file.names[file.exists(file.names)]
      }
      
      samples_counts = suppressMessages(readDGE(file.names))

      if(!is.null(WGCNA.schema)){
        WGCNA.file.names = sprintf("%s/%s%s", Pars$counts.dir, all.split..sample.names..in.WGCNA, Pars$counts.suffix)
        if (any(!file.exists(WGCNA.file.names))){
          msg = sprintf('\nWGCNA: Total %d of %d *.counts files declared in WGCNA sheet do not exist in the directory "%s": %s\n',
                      sum.mod(!file.exists(WGCNA.file.names)), length(WGCNA.file.names), Pars$counts.dir, toString(
                        sapply(strsplit(WGCNA.file.names[!file.exists(WGCNA.file.names)], split='/'), function(x) { tail(x,1)})
                        ))
          if(some.samples..should.be.merged..in.WGCNA)
            stop(msg)
          cat(msg)
          if (sum.mod(file.exists(WGCNA.file.names)) == 0){
            msg = sprintf('WGCNA: No *%s files are found in the directory "%s". Please fix "WGCNA groups" sheet. WGCNA will be disabled', Pars$counts.suffix, Pars$counts.dir)
            cat(msg)
            warning(msg)
            WGCNA.schema = NULL
          } else {
            WGCNA.schema = WGCNA.schema[file.exists(WGCNA.file.names),]
            all.split..sample.names..in.WGCNA = all.split..sample.names..in.WGCNA[file.exists(WGCNA.file.names)]
            rownames(WGCNA.schema) = WGCNA.schema$'Sample names'
            WGCNA.file.names = WGCNA.file.names[file.exists(WGCNA.file.names)]
          }
        }
      }
    
    } else if (Pars$read.counts.table != ''){
      sep = Choose.Separator(Pars$read.counts.table)
      rt = read.table(Pars$read.counts.table, stringsAsFactors = FALSE, sep = sep, header = TRUE, check.names = FALSE)
      src.available.samples = colnames(rt)
      needed.samples = colnames(rt)[colnames(rt) %in% all.split..sample.names]
      rt = rt[, needed.samples, drop = FALSE]
      absent.samples = all.split..sample.names[!(all.split..sample.names %in% colnames(rt))]

      if(length(absent.samples) > 0){
        msg = sprintf('The following samples are absent in file %s (total %d of %d): %s\n. Available samples %s\n', Pars$read.counts.table, length(absent.samples),
          length(all.split..sample.names), toString(absent.samples), toString(src.available.samples))
        if(some.samples..should.be.merged)
          stop(msg)
        cat(msg)

        all.split..sample.names = all.split..sample.names[all.split..sample.names %in% colnames(rt)]
        schema = schema[schema$'Sample names' %in% colnames(rt), ]
        rownames(schema) = schema$'Sample names'
        rt = rt[, schema$'Sample names']
      }
      if(length(absent.samples) == length(all.split..sample.names))  stop(sprintf('All samples are absent in the expression data file %s', Pars$read.counts.table))

      samples_counts = DGEList(rt)

      if(!is.null(WGCNA.schema)){
        # sep = Choose.Separator(Pars$read.counts.table)
        rt = read.table(Pars$read.counts.table, stringsAsFactors = FALSE, sep = sep, header = TRUE, check.names = FALSE)
        src.available.samples = colnames(rt)
        needed.samples = colnames(rt)[colnames(rt) %in% all.split..sample.names..in.WGCNA]
        rt = rt[, needed.samples, drop = FALSE]
        absent.samples = all.split..sample.names..in.WGCNA[!(all.split..sample.names..in.WGCNA %in% colnames(rt))]

        if(length(absent.samples) > 0){
          mas = sprintf('WGCNA groups: the following samples are absent in file %s (total %d of %d): %s\n. Available samples %s\n', Pars$read.counts.table, length(absent.samples),
            length(all.split..sample.names..in.WGCNA), toString(absent.samples), toString(src.available.samples))
          if(some.samples..should.be.merged..in.WGCNA)
            stop(msg)
          cat(msg)

          all.split..sample.names..in.WGCNA = all.split..sample.names..in.WGCNA[all.split..sample.names..in.WGCNA %in% colnames(rt)]
          WGCNA.schema = WGCNA.schema[WGCNA.schema$'Sample names' %in% colnames(rt), ]
          rownames(WGCNA.schema) = WGCNA.schema$'Sample names'
          # rt = rt[, schema$'Sample names']
        }
        if(length(absent.samples) == length(all.split..sample.names..in.WGCNA)){
          msg = sprintf('\n WGCNA groups:   ALL samples declared in "WGCNA groups" are ABSENT in file %s. WGCNA will be disabled.', Pars$read.counts.table)
          cat(msg)
          warning(msg)
          WGCNA.schema = NULL
        } else {
          WGCNA.schema = WGCNA.schema[WGCNA.schema$'Sample names' %in% colnames(rt), ]
          rownames(WGCNA.schema) = WGCNA.schema$'Sample names'
        }
        #rt = rt[, WGCNA.schema$'Sample names']
        #samples_counts = DGEList(rt)

      }

    } else stop('Both Pars$read.counts.table and Pars$counts.dir are not specified')
    
    miR.ratio = sum.mod(grepl('mir-', rownames(samples_counts), ignore.case = TRUE)) / nrow(samples_counts)
    if(miR.ratio > 0.7){
      cat('microRNA mode activated\n')
      miR.mode = TRUE
    } else miR.mode = FALSE

    if(Pars$disable.biomaRt %in% c('auto', '{auto}'))  Pars$disable.biomaRt = miR.mode

    if(length(Reference.genes) > 0){
      absent.Reference.genes = Reference.genes[!(Reference.genes %in% rownames(samples_counts))]
      if(length(absent.Reference.genes) == 0){
        cat(sprintf('All %d reference gene(s) are found in the expression dataset\n', length(Reference.genes)))
      } else {
        msg.text = sprintf('%d of %d reference genes are NOT FOUND in the expression datasets. They will be excluded.', length(absent.Reference.genes), length(Reference.genes))
        cat(msg.text)
        cat('\n')
        warning(msg.text)
        Reference.genes = Reference.genes[Reference.genes %in% rownames(samples_counts)]
        if(length(Reference.genes) == 0 & Pars$RNA.Seq.norm.method == 'RG')
          stop('All reference genes are NOT found in the expression datasets')
      }
      
    }

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
    
    Analyses.groups = Read.DB.entries.from.Excel__MultiSet(Parameters.xlsx.file.name, 'Analyses groups')
    Analyses.groups %>% unlist %>% .[!duplicated(.)] -> models.from.Analyses.groups

    all.preds.from.Analyses.groups = c()
    conv.models = sapply(models.from.Analyses.groups, function(x){ r = Replace.Preds.with.Abbrebiations(x, Pars) })
    for (tmp.Pred.set in conv.models)  all.preds.from.Analyses.groups = append(all.preds.from.Analyses.groups, Extract.predictors(tmp.Pred.set))
    all.preds.from.Analyses.groups = levels(factor(all.preds.from.Analyses.groups))
    all.preds.from.Analyses.groups = sapply(all.preds.from.Analyses.groups, function(x){ r = Revert.Preds.from.Abbrebiations(x,Pars) })
    not.found.preds = all.preds.from.Analyses.groups[!(all.preds.from.Analyses.groups %in% colnames(schema))]
    if(length(not.found.preds) > 0){
      cat(sprintf('The following predictors from "Analyses groups" sheet are not found in "Sample Setup" sheet:\n'))
      for(pred in not.found.preds) cat(sprintf('     %s\n', pred))
      cat(sprintf('Total %d predictor(s). Please check "Sample Setup" row names or GLm models listed in "Analyses groups" sheet. These predictors will be removed now.', length(not.found.preds)))
    }
    Analyses.groups = lapply(Analyses.groups, function(x) x[!(x %in% not.found.preds)])
    Analyses.groups.to.remove.names = names(Analyses.groups)[lapply(Analyses.groups, length) %in% 0]
    for(name in Analyses.groups.to.remove.names){
      Analyses.groups[[name]] <- NULL
      cat(sprintf('Analysis group "%s" was removed since it is empty\n', name))
    }

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

    noint = rownames(samples_counts$counts) %in% c("__no_feature","__ambiguous",'__too_low_aQual','__not_aligned','__alignment_not_unique')
    complete.gene.list = sort(rownames(samples_counts$counts)[!noint],decreasing = FALSE)
    complete.Ensembl.gene.list = complete.gene.list
    General.maRt.table = Get.BiomaRt.table(complete.gene.list, Pars, DB.data, forced.maRt.table = NULL)

  } else {
    Analyses.groups = vector(mode = 'list')
    Startup.Data$use.ext.DE.data = TRUE
    Pars$Models.to.Test = ext.DE.data.files
    Startup.Data$ext.DE.data.files = ext.DE.data.files
    Startup.Data$ext.DE.data = vector(mode = "list")
    
    complete.Ensembl.gene.list = c()
    for(f in ext.DE.data.files){
      DE.data = read.table(f,stringsAsFactors = FALSE,sep = '\t',header = T,check.names = FALSE)
      colnames(DE.data) = tolower(colnames(DE.data))
      all.split..sample.names = colnames(DE.data)
      all.src..sample.names = colnames(DE.data)
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
 
  Gene.Symbols__to__Ensembl.IDs = tmp.maRt.table[,'ensembl_gene_id']
  names.with.dup = sapply(tmp.maRt.table[,'external_gene_name'],function(y) { as.character(y) })
  keep = !duplicated(names.with.dup)
  Gene.Symbols__to__Ensembl.IDs = Gene.Symbols__to__Ensembl.IDs[keep]  
  names(Gene.Symbols__to__Ensembl.IDs) = names.with.dup[keep]

  uppercase.names.with.dup = toupper(names(Gene.Symbols__to__Ensembl.IDs))
  keep = !duplicated(uppercase.names.with.dup)
  Uppercase.Gene.Symbols__to__Ensembl.IDs = Gene.Symbols__to__Ensembl.IDs[keep]
  names(Uppercase.Gene.Symbols__to__Ensembl.IDs) = uppercase.names.with.dup[keep]

  tmp.maRt.table = tmp.maRt.table[!duplicated(tmp.maRt.table[,'ensembl_gene_id']),]
  Ensembl.IDs__to__gene.symbols = tmp.maRt.table[,'external_gene_name']
  names(Ensembl.IDs__to__gene.symbols) = sapply(tmp.maRt.table[,'ensembl_gene_id'],function(y) { as.character(y) })

  cumulative__genes.found.by.ID = c()
  cumulative__genes.not.found.by.ID = c()
  cumulative__genes.found.by.name = c()
  cumulative__genes.not.found.by.name = c()
  for(group.name in names(Startup.Data$inDetails.Standalone.genes)){
    genes = Startup.Data$inDetails.Standalone.genes[[group.name]]
    if(is.null(genes))  next
    genes.found.by.ID = genes[genes %in% complete.Ensembl.gene.list]
    genes.not.found.by.ID = genes[!(genes %in% complete.Ensembl.gene.list)]
    genes.not.found.by.ID = toupper(genes.not.found.by.ID)
    genes.found.by.name = genes.not.found.by.ID[genes.not.found.by.ID %in% names(Uppercase.Gene.Symbols__to__Ensembl.IDs)]
    genes.not.found.by.name = genes.not.found.by.ID[!(genes.not.found.by.ID %in% names(Uppercase.Gene.Symbols__to__Ensembl.IDs))]
    genes = c(genes.found.by.ID, Uppercase.Gene.Symbols__to__Ensembl.IDs[genes.found.by.name])
    genes = genes[!duplicated(genes)]
    Startup.Data$inDetails.Standalone.genes[[group.name]] = genes

    cumulative__genes.found.by.ID = c(cumulative__genes.found.by.ID, genes.found.by.ID)
    cumulative__genes.not.found.by.ID = c(cumulative__genes.not.found.by.ID, genes.not.found.by.ID)
    cumulative__genes.found.by.name = c(cumulative__genes.found.by.name, genes.found.by.name)
    cumulative__genes.not.found.by.name = c(cumulative__genes.not.found.by.name)
  }
  cumulative__genes.found.by.ID = cumulative__genes.found.by.ID[!duplicated(cumulative__genes.found.by.ID)]
  cumulative__genes.not.found.by.ID = cumulative__genes.found.by.ID[!duplicated(cumulative__genes.not.found.by.ID)]
  cumulative__genes.found.by.name = cumulative__genes.found.by.ID[!duplicated(cumulative__genes.found.by.name)]
  cumulative__genes.not.found.by.name = cumulative__genes.found.by.ID[!duplicated(cumulative__genes.not.found.by.name)]
  cat(sprintf('\ninDetails:\n%d genes found from Gene Ontology terms\n%d genes found from KEGG pathway IDs\nAdditionally, %d genes were found by Ensembl ID, %d were found by Gene Name, %d genes were not found:\n%s\n',
    Startup.Data$inDetails.GO.genes %>% unlist %>% .[!duplicated(.)] %>% length,
    Startup.Data$inDetails.KEGG.genes %>% unlist %>% .[!duplicated(.)] %>% length,
    length(cumulative__genes.found.by.ID), length(cumulative__genes.found.by.name), length(cumulative__genes.not.found.by.name), toString(cumulative__genes.not.found.by.name)))
  
  
  
  names(Pars$Models.to.Test) = NULL
  
  Pars$Analyses.groups = Analyses.groups

  
  #####
  ##### Adding MDS coordinates of samples as predictors, if this option (add.MDS.dims.as.predictors) is et TRUE
  
  if ((Pars$add.MDS.dims.as.predictors | Pars$DE.package=='edger') & is.null(ext.DE.data.files)){
    counts = samples_counts$counts
    # excluding meta-tags and low expression genes
    meta.tags = rownames(counts)[which(startsWith(rownames(counts), '__'))]
    if(length(meta.tags) > 0){
      noint = rownames(counts) %in% meta.tags
      counts = counts[!noint,]
    }
    cpms = cpm(counts)
    if (Pars$min.samples.with.sufficient.CPM == '{auto}') {
      min.samples.with.sufficient.CPM = Calc.min.high.CPM.samples.from.total.samples.count(dim(schema)[1])
    } else {
      min.samples.with.sufficient.CPM = Pars$min.samples.with.sufficient.CPM
      if(min.samples.with.sufficient.CPM > ncol(counts))
        stop(sprintf('min.samples.with.sufficient.CPM argument (%d) is greater than the number of samples (%d)', min.samples.with.sufficient.CPM, ncol(counts)))
    }
    
    if(Pars$sufficient.CPM.to.analyze %in% c('{auto}', 'auto')){
      min.CPM.single.sample = max(1, 3e+7 / mean(colSums(counts)) / 1.25 ) * Pars$sufficient.CPM..autoadjust..multiplier
    } else {
      min.CPM.single.sample = Pars$sufficient.CPM.to.analyze
    }
    # cat(sprintf('\nSetting CPM limit %.1f for at least %d samples', min.CPM.single.sample, min.samples.with.sufficient.CPM))
    
    keep.by.CPM = rowSums(cpms >= min.CPM.single.sample) >= min.samples.with.sufficient.CPM
    keep.by.raw.read.count = rowSums(counts >= Pars$sufficient.raw.read.counts.to.analyze) >= min.samples.with.sufficient.CPM

    counts = counts[keep.by.CPM & keep.by.raw.read.count, ]
    colnames(counts) = all.split..sample.names
    
    #schema$'Sample names'
    d = DGEList(counts=counts)
    if(Pars$RNA.Seq.norm.method == 'RG'){
      d = Normalize..Using.Reference.Genes(d, Reference.genes)
    } else {
      d = calcNormFactors(d, method = Pars$RNA.Seq.norm.method) #,"TMM","upperquartile","none"
    }
    
    #colnames(d$counts)
    
    if(some.samples..should.be.merged){
      comma.separated..sample.names = schema$'Sample names'
      for(x in 1:10)
        comma.separated..sample.names = gsub(sprintf('%s ', Pars$Merge.samples..separator), Pars$Merge.samples..separator, comma.separated..sample.names)
      
      merged.counts..list = vector(mode = 'list')
      sample.names..list = strsplit(comma.separated..sample.names, split = ',')
      n = 1
      for(n in 1:length(comma.separated..sample.names)){
        if(length(sample.names..list[[n]]) == 1){
          merged.counts..list[[schema$'Sample names'[n]]] = d$counts[, sample.names..list[[n]], drop = TRUE]
            
        } else {
          cat(sprintf('merging samples %s into "%s" ....  \n', paste(sample.names..list[[n]], collapse = ' + '), schema$'Sample names'[n]))
          lib.sizes = d$samples[sample.names..list[[n]], 'lib.size', drop = TRUE]
          norm.factors = d$samples[sample.names..list[[n]], 'norm.factors', drop = TRUE]
          merged.counts..list[[schema$'Sample names'[n]]] = 
            apply(d$counts[, sample.names..list[[n]], drop = FALSE], 1, function(x) round(sum(x / norm.factors / lib.sizes * sum(lib.sizes)) / length(x)))
        }
      }
      
      d = DGEList(counts = as.data.frame(merged.counts..list, check.names = FALSE))
      if(Pars$RNA.Seq.norm.method == 'RG'){
        d = Normalize..Using.Reference.Genes(d, Reference.genes)
      } else {
        d = calcNormFactors(d, method = Pars$RNA.Seq.norm.method) #,"TMM","upperquartile","none"
      }
      cat('\nMerging completed.\n')
          
    }
    
    if(write.CPM.tables){
      write.table.mod(cpm(d, normalized.lib.sizes = TRUE), sprintf('CPMs.tsv'), sep='\t', quote = FALSE)
      saveRDS(cpm(d, normalized.lib.sizes = TRUE), file = 'CPMs.rds')
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

  absent.Ind.Par.Models = colnames(Ind.Par.table)[!(colnames(Ind.Par.table) %in% Pars$Models.to.Test)]
  if(length(absent.Ind.Par.Models) > 0)   cat(sprintf('%d of %d models with specified individual parameters are absent: %s', length(absent.Ind.Par.Models), ncol(Ind.Par.table), toString(absent.Ind.Par.Models)))


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
    norm.CPMs = cpm(counts(dds, normalized = TRUE))
    if(Pars$report.norm.read.counts.instead.of.CPM){
      lib.sizes = d$samples[,'lib.size']
      norm.factors = d$samples[,'norm.factors']
      norm.counts = t(t(d$counts) / lib.sizes / norm.factors * mean(lib.sizes))
    }

    write.table.mod(norm.CPMs,sprintf('CPMs.tsv'),sep='\t',quote = F)
    
  }
  
  # testing for the presence of Python
  if (disable.excel) Pars$Create.Excel.results = FALSE
  if (Pars$Create.Excel.results & !disable.excel){
    Pars$python.bin = Check.Python.availability()
  }


  all.preds = colnames(schema)[!(colnames(schema) %in% 'Sample names')]
  Startup.Data$predictor.values = suppressWarnings(data.matrix(schema)[, all.preds])
  if (length(all.preds) == 1){
    Startup.Data$predictor.values = as.data.frame(Startup.Data$predictor.values)
    colnames(Startup.Data$predictor.values) = c(all.preds)
  }
  Startup.Data$predictor.values.nr = Eliminate.Redundant.Predictors.Info(Startup.Data$predictor.values)

  # print(Startup.Data$predictor.values.nr)

  Startup.Data$Pars = Pars
  Startup.Data$schema = schema
  Startup.Data$WGCNA.schema = WGCNA.schema
  Startup.Data$WGCNA.groups = colnames(WGCNA.schema)[-1]
  Startup.Data$General.maRt.table = General.maRt.table
  Startup.Data$Ensembl.IDs__to__gene.symbols = Ensembl.IDs__to__gene.symbols
  Startup.Data$topGO.Expression.Profiles___Custom.GO.terms = topGO.Expression.Profiles___Custom.GO.terms
  Startup.Data$GO.clusterProfiler.Expression.Profiles___Custom.DB.entries = topGO.Expression.Profiles___Custom.GO.terms
  Startup.Data$KEGG.clusterProfiler.Expression.Profiles___Custom.DB.entries = KEGG.clusterProfiler.Expression.Profiles___Custom.DB.entries
  Startup.Data$Reactome.clusterProfiler.Expression.Profiles___Custom.DB.entries = Reactome.clusterProfiler.Expression.Profiles___Custom.DB.entries
  Startup.Data$Custom.Pathways.to.Visualize = Custom.Pathways.to.Visualize
  Startup.Data$Omit.pathways = Omit.pathways
  Startup.Data$Completed.steps.file = Completed.steps.file
  Startup.Data$Parameters.xlsx.file.name = Parameters.xlsx.file.name
  Startup.Data$Ind.Par.table = Ind.Par.table
  Startup.Data$Reference.genes = Reference.genes
  Startup.Data$Heatmaps.annotation = Heatmaps.annotation
  Startup.Data$miR.mode = miR.mode
  
  Startup.Data$Startup.hashmd5 = digest(toString(c(Startup.Data$Pars,Startup.Data$General.maRt.table,Startup.Data$topGO.Expression.Profiles___Custom.GO.terms,Startup.Data$Custom.Pathways.to.Visualize,Startup.Data$Omit.pathways)))
  saveRDS(Startup.Data, file = 'Startup.Data.rds')
  
  cat('\nPreparing completed.\n')
  return(Startup.Data)
  
  #####
}


Create.MDS.plots = function(d, Main.Component, Main.Predictor.name, out.dir = NULL, prefix = ''){
  if(is.null(out.dir))  out.dir = getwd()
  #col.variety = c("darkgreen","blue","red","darkgoldenrod2","chartreuse2","darkorchid2","black","gray59","lightsteelblue1","orange4","sienna3","bisque3")
  ### gradient: blue-teal-green-yellow-orange
  if(is.null(Main.Predictor.name)){
    plot.title = sprintf('MDS plot')
  } else plot.title = sprintf('MDS plot, %s', Main.Predictor.name)
  
  plot.subtitle = ''
  
  if(!is.null(Main.Component)){
    if (length(levels(factor(Main.Component))) == 2) plot.subtitle = sprintf('color indicates %s (blue-orange)',Main.Predictor.name)
    if (length(levels(factor(Main.Component))) == 3) plot.subtitle = sprintf('color indicates %s (blue-green-orange)',Main.Predictor.name)
    if (length(levels(factor(Main.Component))) >= 4) plot.subtitle = sprintf('color indicates %s (blue-green-orange gradient)',Main.Predictor.name)
  } else {
    Main.Component = rep(1, dim(d$counts)[2])
  }
  
  col.variety = colorRampPalette(c("#0b00c4","#02babc","#25ad00","#cfcd00","#ff8400"))(length(levels(factor(Main.Component))))  
  
  .pardefault <- par(no.readonly = TRUE)
  par(mar=c(1,1,1,1))
  if(dim(d$counts)[2] > 2 & dim(d$counts)[1] > 5){
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
  
  
  if(dim(d$counts)[2] > 2 & dim(d$counts)[1] > 5){
    pdf(file = sprintf('%s/%sMDS.dim.1-2.pdf', out.dir, prefix),
        width=14,height=12,pointsize=12)
    MDS.info = plotMDS(d,dim.plot = c(1,2), col = col.variety[factor(Main.Component)])
    two.dim.MDS.info = cbind(MDS.info$x,MDS.info$y)
    title(main = plot.title, sub = plot.subtitle)
    dev.off()
    
    
    if(dim(d$counts)[2] > 3 & dim(d$counts)[1] > 5){
      pdf(file = sprintf('%s/%sMDS.dim.1-3.pdf', out.dir, prefix),
          width=14,height=11,pointsize=12)
      MDS.info = plotMDS(d,dim.plot = c(1,3), col = col.variety[factor(Main.Component)])
      title(main = plot.title, sub = plot.subtitle)
      dev.off()
      three.dim.MDS.info = cbind(two.dim.MDS.info, MDS.info$y)
      colnames(three.dim.MDS.info) = c('MDS dim1','MDS dim2','MDS dim3')
      write.table.mod(x = three.dim.MDS.info,file = 'MDS info.tsv',sep = '\t')
    }
    
    png.mod(filename = sprintf('%s/%sMDS.dim.1-2.png', out.dir, prefix), units="in",res=600,
            width=8.23,height=6.3,pointsize=12)
    plotMDS(d,dim.plot = c(1,2), col = col.variety[factor(Main.Component)])
    title(main = plot.title, sub = plot.subtitle)
    dev.off()
    
    if(dim(d$counts)[2] > 3 & dim(d$counts)[1] > 5){
      png.mod(filename = sprintf('%s/%sMDS.dim.1-3.png', out.dir, prefix), units="in",res=600,
              width=8.23,height=6.3,pointsize=12)
      a = plotMDS(d,dim.plot = c(1,3), col = col.variety[factor(Main.Component)])
      title(main = plot.title, sub = plot.subtitle)
      dev.off()
    }
  }
  par(.pardefault)
}


Create.Distance.matrices.and.Dendrograms = function(Pars, cpms, schema, Main.Component, GLM.model, Current.Predictors.list, out.dir = NULL, prefix = '', verbose = TRUE){
  if(is.null(out.dir)) out.dir = getwd()
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
        if (Pars$Distance.method == 'canberra.weighted') dist.matrix[F,S] = Weighted.Canberra.dist(cpms[,F], cpms[,S], Pars$Min.CPM.to.include.in.dist, Pars$Canberra.weight.power)
        if (Pars$Distance.method == '1-cor') dist.matrix[F,S] = 1 - suppressWarnings(cor.test(cpms[,F],cpms[,S],method = 'spearman')$estimate)
        ready = ready+1
        if(verbose)
          cat(sprintf('\rCalculating dist matrices. Completed %.2f%%',ready/N/N*100))
        
        #Weighted.Canberra.dist(c(1,2,3,4,5,6,7,8,9,10),c(1,20,30,4,5,6,7,8,9,10),min.value.to.calc.dist = 0.99)
      }
    }
  }
  cat(sprintf('\rCalculating dist matrices completed           \n'))

  tmp=as.dist(dist.matrix)
  dend=hclust(tmp)
  dend=as.dendrogram(dend)
  
  if(!is.null(Main.Component)){
    col.variety = colorRampPalette(c("#0b00c4","#02babc","#25ad00","#cfcd00","#ff8400"))(length(levels(factor(Main.Component))))
    dend.colors = col.variety[factor(Main.Component)]
    names(dend.colors) = schema$'Sample names'
    labels_colors(dend) = dend.colors[labels(dend)]
  }

  png.mod(filename = sprintf('%s/%sClustering, dist as %s.png', out.dir, prefix, Pars$Distance.method),units="in",res=600,
      width=8.23,height=6.3,pointsize=12)
  plot(dend,main = sprintf('%s%s, clustering dendrogram [%s]', prefix, GLM.model, Pars$Distance.method),
       sub = sprintf('distances are calculated with "%s" method', Pars$Distance.method))
  dev.off()
  ##hang=-1
  
  
  #Current.Predictors.list = c('Age','Genotype','Gender')
  info.block = data.frame(as.data.frame(schema)[colnames(cpms),Current.Predictors.list])
  colnames(info.block) = Current.Predictors.list
  info.block.mod = cbind(array(dim = c(dim(info.block)[2],dim(info.block)[2])),t(info.block))
  colnames(info.block.mod) = colnames(cbind(info.block,dist.matrix))
  class(info.block.mod) <- "numeric"

  tsv.file.name = sprintf('%s/%sdistance matrix.tsv', out.dir, prefix)
  simple.tsv.file.name = sprintf('%s/%sdistance matrix simple.tsv', out.dir, prefix)
  simple.xlsx.file.name = sprintf('%s/%sdistance matrix.xlsx', out.dir, prefix)
  write.table.mod(x=rbind(info.block.mod, cbind(info.block,dist.matrix)), file = tsv.file.name, sep="\t", na = "", quote = FALSE)
  write.table.mod(x = dist.matrix, file = simple.tsv.file.name, sep="\t", na = "", quote = FALSE)
  if(Pars$Create.Excel.results){
    CL = sprintf('%s "%s/Any.2.Excel.py" --in "%s" --out "%s" --dist-matrix-mode yes --predictor-rows-count 0',
      Pars$python.bin, Pars$suppl.data.dir, simple.tsv.file.name, simple.xlsx.file.name )
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    tmp = file.remove(simple.tsv.file.name)
  }
}


#Ind.Par.table

Analyze.GLM = function(Startup.Data, GLM.model, forced.parameters = NULL,
                       script.dir = NULL, base.dir = NULL, current.dir = NULL, Create.Excel.results = NULL,
                       include.GLM.model.name.in.Result.names = FALSE,
                       include.GLM.model.name.in.Result.names..as.suffix = TRUE,
                       randomization.test.permut.N = NULL,
                       verbose = TRUE){
  
  if(!is.null(base.dir))     if(endsWith(base.dir, '/'))     base.dir = substr(base.dir, 1, nchar(base.dir) - 1)
  if(!is.null(script.dir))   if(endsWith(script.dir, '/'))   script.dir = substr(script.dir, 1, nchar(script.dir) - 1)
  if(!is.null(current.dir))  if(endsWith(current.dir, '/'))  current.dir = substr(current.dir, 1, nchar(current.dir) - 1)

  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  

  # if(!adjust.bias.in.read.counts)        bias.factors__to.adjust__read.counts = NULL
  # if(!adjust.bias.in.LogFC)              bias.factors__to.adjust__LogFC = NULL
  # if(!include.BF.Quantile.statistics)    bias.factors.for.quantile.statistics = NULL

  # if(adjust.bias.in.read.counts | adjust.bias.in.LogFC)  Evaluate.bias.factors.Associations = TRUE

  # if(is.null(bias.factors.to.analyze))  bias.factors.to.analyze = c()
  # if(!is.null(bias.factors__to.adjust__LogFC))   bias.factors.to.analyze = c(bias.factors.to.analyze, bias.factors__to.adjust__LogFC)
  # if(!is.null(bias.factors__to.adjust__read.counts))   bias.factors.to.analyze = c(bias.factors.to.analyze, bias.factors__to.adjust__read.counts)
  # if(!is.null(bias.factors.for.quantile.statistics))   bias.factors.to.analyze = c(bias.factors.to.analyze, bias.factors.for.quantile.statistics)
  # if(length(bias.factors.to.analyze) > 0)   bias.factors.to.analyze = bias.factors.to.analyze[!duplicated(bias.factors.to.analyze)]

  # if(!Evaluate.bias.factors.Associations)  bias.factors.to.analyze = c()

  # if(!(adjust.bias.in.LogFC.mode %in% c('before GLM', 'after GLM', 'none'))){
  #   stop('incorrect adjust.bias.in.LogFC.mode. SHould be either "before GLM", "after GLM"')
  # }

  # if(!adjust.bias.in.LogFC)  adjust.bias.in.LogFC.mode = 'none'

  
  ####################################################################
  #### Analysing GLM/differential expression
  
  if(is.null(current.dir)) current.dir = getwd()

  GLM.model = gsub("~","",GLM.model,fixed = TRUE)
  
  cat(sprintf("Processing model  \"~ %s\"...\n",GLM.model))
  setwd(Pars$results.dir)
  
  Analysis.Data = new("RTransAnalysisData")
  
  
  ## extracting names of all used individual predictors
  conv.GLM.model = Replace.Preds.with.Abbrebiations(GLM.model,Pars)
  Current.Predictors.list = Extract.predictors(conv.GLM.model)
  Current.Predictors.list = sapply(Current.Predictors.list, function(x){ r = Revert.Preds.from.Abbrebiations(x,Pars) })
  
  #Analysis.name = sprintf("Associations with %s",GLM.model)
  Analysis.Data$GLM.model = GLM.model

  Individual.Pars = NULL
  if(GLM.model %in% colnames(Startup.Data$Ind.Par.table)){
    tmp = subset(Startup.Data$Ind.Par.table, select = GLM.model, subset = !is.na(Startup.Data$Ind.Par.table[,GLM.model]))
    par.table = cbind(as.data.frame(rownames(tmp)), tmp)
    colnames(par.table) = c('variable','value')
    cat(sprintf('Found %d individual parameters:\n', nrow(par.table)))
    for(x in 1:nrow(par.table)){
      cat(sprintf('      %s:    %s\n', par.table[x,1], par.table[x,2]))
    }
    Pars = Process.Parameters.table(Pars, par.table, show.warnings = FALSE, base.dir = base.dir, script.dir = script.dir, current.dir = current.dir)
  }
  
  schema = tryCatch(expr = { suppressMessages(suppressWarnings(readxl::read_excel(Startup.Data$Parameters.xlsx.file.name ,sheet='Sample setup',col_names = TRUE)))
  }, error = function (err){
    stop(sprintf('Sheet "Sample setup" is not found in workbook %s',Startup.Data$Parameters.xlsx.file.name))
  })
  
  schema %<>% as.data.frame
  commented.lines = sapply(X = schema[,1], FUN = function (x) substring(x,first=1,last=1)) %in% '#'
  commented.cols = sapply(X = colnames(schema), FUN = function (x) substring(x,first=1,last=1)) %in% '#'
  schema = schema[!commented.lines, !commented.cols]
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

  all.src..sample.names = schema$'Sample names'
  all.split..sample.names = all.src..sample.names
  if(sum(duplicated(all.src..sample.names)) > 0){
    stop(sprintf('Duplicated sample names are not allowed: %s', all.src..sample.names[duplicated(all.src..sample.names)]))
  }

  paired.tests.are.available = FALSE
  perform.FC.associations.paired.test = FALSE
  perform.classic.paired.test = FALSE
  
  if(Pars$classic.paired.test == '{auto}' & Pars$FC.associations.paired.test == '{auto}'){
    paired.tests.are.available = Check.availability.of.paired.test(schema$'Sample names', verbose = FALSE, verbose.OK = TRUE)
    cat(sprintf('Paired test(s) will %sbe performed\n', ifelse(paired.tests.are.available, '', 'NOT ')))

  } else if(Pars$classic.paired.test == TRUE | Pars$FC.associations.paired.test == TRUE){
    paired.tests.are.available = Check.availability.of.paired.test(schema$'Sample names', verbose = TRUE, verbose.OK = TRUE)
    if(!paired.tests.are.available) stop(sprintf('Paired test is not available for model "%s". However, it is turned on. Check classic.paired.test and FC.associations.paired.test arguments', GLM.model))
  }

  if(paired.tests.are.available){
    new.order = Reorder.sample.names..for.paired.test(schema$'Sample names', return.order = TRUE)
    schema = schema[new.order, ]
    samples..base.names = schema$'Sample names' %>% substr(., 1, nchar(.) - 1)
  }


  Components = tryCatch(expr = {
    data.matrix(schema[,Current.Predictors.list]) 
  }, error = function (err){
    data.matrix(as.numeric(schema[,Current.Predictors.list]))
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
  if(!is.null(randomization.test.permut.N)){
    cat(sprintf('\nRandomization test, permutation %d.\nMain component before shuffling:\n%s\n', randomization.test.permut.N, toString(Components[,Predictors.count])))
    # Components[,Predictors.count] = sample(Components[,Predictors.count])
    cat(sprintf('Main component after shuffling:\n%s\n\n', toString(Components[,Predictors.count])))
  }
  Main.Component = as.numeric(Components[,Predictors.count])
  
  if(paired.tests.are.available & (Pars$classic.paired.test == '{auto}' | Pars$classic.paired.test == TRUE)){
    ## testing whether classic paired test (e.g. Tumor vs Norm paired) is available
    Main.Component..first.part = Main.Component[1:(length(Main.Component)/2)]
    Main.Component..second.part = Main.Component[(length(Main.Component)/2 + 1):length(Main.Component)]
    if((Main.Component..first.part %>% DeDup %>% length) == 1 | (Main.Component..second.part %>% DeDup %>% length) == 1){
      perform.classic.paired.test = TRUE
      cat('Sample setup is compartible with classic paired test. It will be performed\n')
    } else {
      if(Pars$classic.paired.test == TRUE){
        stop(sprintf('Sample setup (Main.Component) seems to be incompartible with classic paired test: each epigroup N/T contains more than 1 corresponding values of Main.Component.\nMain.Component: %s\n',toString(Main.Component)))
      }
      cat(sprintf('Sample setup (Main.Component) seems to be incompartible with classic paired test: each epigroup N/T contains more than 1 corresponding values of Main.Component. Classic paired test will not be performed\n'))
    }
  }

  if(paired.tests.are.available & (Pars$FC.associations.paired.test == '{auto}' | Pars$FC.associations.paired.test == TRUE)){
    ## testing whether associations paired test (e.g. find genes, for which Tumor-to-Norm expression ratios are associatied with some tumor characteristics) is available
    Main.Component..first.part = Main.Component[1:(length(Main.Component)/2)]
    Main.Component..second.part = Main.Component[(length(Main.Component)/2 + 1):length(Main.Component)]
    if(all(Main.Component..first.part == Main.Component..second.part)){
      perform.FC.associations.paired.test = TRUE
      cat('Sample setup is compartible with associations paired test. It will be performed instead of standart tests\n')
    } else {
      if(Pars$FC.associations.paired.test == TRUE){
        stop(sprintf('Sample setup (Main.Component) seems to be incompartible with associations paired test: two epigroups N/T contain non-consistent values of Main.Component.\nMain.Component: %s\n',toString(Main.Component)))
      }
      cat(sprintf('Sample setup (Main.Component) seems to be incompartible with associations paired test: two epigroups N/T contain non-consistent values of Main.Component. Associations paired test will not be performed\n'))
    }
  }

  cat(sprintf('\n\n   Classic paired test:        %s\n', ifelse(perform.classic.paired.test, "YES", "NO")))
  cat(sprintf('   Associations paired test:   %s\n\n', ifelse(perform.FC.associations.paired.test, "YES", "NO")))

  
  if(perform.classic.paired.test){
    if(Pars$Main.test == 'QL'){
      Pars$Main.test = 'QL paired'
    } else if(Pars$Main.test == 'LR'){
      Pars$Main.test = 'LR paired'
    }

    if(!Pars$use.likelihood.ratio.test..paired & !Pars$use.quasi.likelihood.test..paired){
      cat('\nBoth Pars$use.likelihood.ratio.test..paired and Pars$use.quasi.likelihood.test..paired are turned OFF. However, you switched on global paired test setting. Pars$use.quasi.likelihood.test..paired will be enabled\n')
      Pars$use.quasi.likelihood.test..paired = TRUE
    }
  } else {
    Pars$use.quasi.likelihood.test..paired = FALSE
    Pars$use.likelihood.ratio.test..paired = FALSE
  }
  
  if(perform.FC.associations.paired.test) {
    FC.DE.mode = 'FC associations'
    Pars$use.exact.test.in.binary.predictors = FALSE
    Pars$use.quasi.likelihood.test = FALSE
    Pars$use.likelihood.ratio.test = FALSE
    Pars$use.quasi.likelihood.test..paired = FALSE
    Pars$use.likelihood.ratio.test..paired = FALSE
    Pars$Main.test = Pars$Preferred.FC.association.test
    
    if(Pars$Preferred.FC.association.test == 't-test (FC)')
      Pars$use.ttest.for.FC.association.paired.test = TRUE
    if(Pars$Preferred.FC.association.test == 'Mann-Wh. (FC)')
      Pars$use.Mann.Whitney.criterium.for.FC.association.paired.test = TRUE
    if(Pars$Preferred.FC.association.test == 'Lin.Mod. (FC)')
      Pars$use.Linear.Models.for.FC.association.paired.test = TRUE
    
  } else {
    FC.DE.mode = 'standard'
    Pars$use.Mann.Whitney.criterium.for.FC.association.paired.test = FALSE
    Pars$use.ttest.for.FC.association.paired.test = FALSE
    Pars$use.Linear.Models.for.FC.association.paired.test = FALSE
  }
  


  
  #### testing if this step has been preiously completed
  # step.hashmd5 = digest(cbind(schema$'Sample names',Components)) ## + make test if there are ready files
  # Analysis.Data$step.hashmd5 = step.hashmd5
  
  # if (bypass.if.completed & Read.Completed.steps.status(sprintf("%s.%s.%s.complete",GLM.model,Startup.Data$Startup.hashmd5,step.hashmd5),Startup.Data$Completed.steps.file)){
  #   message(sprintf('GLM model "%s" processing was previously completed. Bypassing...\n',GLM.model))
  #   return(Analysis.Data)
  # }
  ####

  if(is.null(script.dir))  script.dir = Get.script.path()


  # if(Pars$gene.bias.factors.file == '')

  if(Pars$Evaluate.bias.factors.Associations){
    if(Pars$gene.bias.factors.file == ''){
      gene.bias.factors.file = sprintf('%s/Gene.Bias.Factors.tsv', script.dir)
      if(file.exists(gene.bias.factors.file)){
        cat('Gene bias factors file is found.\n')
      } else {
        cat('Cannot find file with gene bias factors information (e.g. transcript length or CG content). Bias correction will not be applied. To enable bias correction place Gene.Bias.Factors.tsv (generated with PPLine) in the script directory\n')
        gene.bias.factors.file = NULL
      }
    } else if(!file.exists(gene.bias.factors.file)){
      stop(sprintf('Gene bias factors file %s is NOT found. Check parameter gene.bias.factors.file', gene.bias.factors.file))
    }
  } else {
    gene.bias.factors.file = NULL
  }


  cat(sprintf("%s:   loading and processing RNA-Seq data...\n",GLM.model))
  
  # setting design matrix
  GLM.formula = GLM.model
  #Current.Predictors.list = sprintf('&&&%s&&&',Current.Predictors.list)
  #Current.Predictors.list = mgsub("&&&","",Current.Predictors.list)
  GLM.formula = mgsub(Current.Predictors.list,
                      sprintf('as.numeric(Components[,"%s"])', Pars$available.preds_to_abbreviations[Current.Predictors.list]),
                      GLM.formula,fixed = TRUE)
  
  GLM.formula = Revert.Preds.from.Abbrebiations(GLM.formula, Pars)

  if(perform.classic.paired.test){
    grp = factor(samples..base.names)
    GLM.formula.paired = sprintf("~ grp + %s", GLM.formula)
  }
  GLM.formula = sprintf("~ %s", GLM.formula)
  #for (Pred in Current.Predictors.list){
  #  GLM.formula = gsub(Pred,sprintf('as.numeric(Components[,"%s"])',Pred),GLM.formula,fixed = T)
  #}
  
  Main.Component.levels <- Main.Component %>% .[!duplicated(.)] %>% .[order(.)]
  lev1.cols = Main.Component %in% Main.Component.levels[1]
  lev2.cols = Main.Component %in% Main.Component.levels[2]

  if(perform.FC.associations.paired.test & (length(Current.Predictors.list) != 1))
    stop(sprintf('Now only one predictor is supported in associations paired tests. There are %d predictors in formula: %s',
      length(Current.Predictors.list), toString(GLM.model)))

  if(perform.FC.associations.paired.test & ((Main.Component.levels %>% length) != 2))
    stop(sprintf('Now only two conditions are supported in associaions paired test. These are %d conditions in main predictor: %s',
      length(Main.Component.levels), toString(Main.Component.levels)))


  Single.binary.predictor = FALSE
  Single.binary.predictor__minimal.group.size = NULL
  if (length(Current.Predictors.list) == 1){
    if(length(Main.Component.levels) == 2){
      cat(sprintf('RTrans has detected that GLM model has only one binary predictor\n'))
      Single.binary.predictor__minimal.group.size = min(table(as.numeric(Main.Component)))
      Single.binary.predictor = TRUE
    }
  }

  if(!perform.FC.associations.paired.test){
    if(Pars$Include.DeltaFreq.data == '{auto}' | Pars$Include.DeltaFreq.data == 'auto'){
      if(length(Main.Component.levels) != 2){
        Pars$Include.DeltaFreq.data = FALSE
      } else {
        Pars$Include.DeltaFreq.data = (sum.mod(lev1.cols) >= 7 | sum.mod(lev2.cols) >= 7)
      }
      cat(sprintf('Pars$Include.DeltaFreq.data is set as "%s"\n', toString(Pars$Include.DeltaFreq.data)))

    } else if(Pars$Include.DeltaFreq.data & (length(Main.Component.levels) != 2)){
      msg = sprintf('\nCalculating DeltaFreq.data is not compartible with sample setup for the analysis "%s". There should be only two conditions (really you have %d). DeltaFreq calculating will be disabled\n',
        Cleanup.Model.Name(GLM.model), length(Main.Component.levels))
      Pars$Include.DeltaFreq.data = FALSE
      cat(msg)
      warning(msg)
    }


    if(Pars$Score.calculating.method == '{auto}' | Pars$Score.calculating.method == 'auto'){
      if(length(Main.Component.levels) != 2){
        Pars$Score.calculating.method = 'standard'
      } else {
        if(sum.mod(lev1.cols) >= 7 | sum.mod(lev2.cols) >= 7){
          Pars$Score.calculating.method = 'complex'
        } else Pars$Score.calculating.method = 'standard'
      }
      cat(sprintf('Pars$Score.calculating.method is set as "%s"\n', toString(Pars$Score.calculating.method)))

    } else if(Pars$Score.calculating.method == 'complex' & (length(Main.Component.levels) != 2)){
      msg = sprintf('\nScore.calculating.method is set to "complex" but this is acceptable only for two groups comparison (analysis "%s"). There should be only two conditions (really you have %d). Score.calculating.method is set to standard\n',
        Cleanup.Model.Name(GLM.model), length(Main.Component.levels))
      Pars$Score.calculating.method = 'standard'
      cat(msg)
      warning(msg)
    }


  } else{
    Pars$Include.DeltaFreq.data = FALSE
    if(Pars$Score.calculating.method == 'complex'){
      msg = 'Score calculating method is switched from "complex" to "standard" since paired tests are performed'
      cat(msg)
      warning(msg)
      Pars$Score.calculating.method = 'standard'
    }
  }

  # loading files with counts
  
  some.samples..should.be.merged = FALSE
  if(Pars$Merge.samples.in.sample.setup){
    some.samples..should.be.merged = (length(grep(Pars$Merge.samples..separator, all.src..sample.names)) > 0)
    all.split..sample.names = all.src..sample.names
    if(some.samples..should.be.merged){
      for(x in 1:10)
        all.split..sample.names = gsub(sprintf('%s ', Pars$Merge.samples..separator), Pars$Merge.samples..separator, all.split..sample.names)
      all.split..sample.names = unlist(strsplit(all.split..sample.names, split = ',')) %>% DeDup.na.rm
    } 
  }

  if(Pars$counts.dir != ''){
    if(!dir.exists(Pars$counts.dir))
      stop(sprintf('Directory "%s" does not exist', Pars$counts.dir))
    
    file.names = sprintf("%s/%s%s", Pars$counts.dir, all.split..sample.names, Pars$counts.suffix)
    #file.names = sprintf("%s/%s%s", Pars$counts.dir, schema$'Sample names', Pars$counts.suffix)
    # print()
    if (any(!file.exists(file.names))){
      #warning(sprintf('Warning; The following files do not exist: %s',toString(file.names[!file.exists(file.names)])))
      msg = sprintf('\nTotal %d of %d *.counts files do not exist in the directory "%s": %s\n',
                    sum.mod(!file.exists(file.names)), length(file.names), Pars$counts.dir, toString(
                      sapply(strsplit(file.names[!file.exists(file.names)], split='/'), function(x) { tail(x,1)})
                    ))
      if(some.samples..should.be.merged)
        stop(sprintf('Sample merging is activated. All the samples must be present. %s', msg))
      
      if(perform.classic.paired.test)
        stop(sprintf('Paired test is turned on. All the samples must be present. %s', msg))
      
      cat(msg)
      if (sum.mod(file.exists(file.names)) == 0) stop(sprintf('No *%s files are found in the directory "%s". Please fix "Sample setup" sheet', Pars$counts.suffix, Pars$counts.dir))
      
      if(!some.samples..should.be.merged){
        schema=schema[file.exists(file.names), ]
        all.split..sample.names = all.split..sample.names[file.exists(file.names)]
        all.src..sample.names = all.src..sample.names[file.exists(file.names)]
        rownames(schema) = schema$'Sample names'
        tmp = colnames(Components)
        Components = data.matrix(Components[file.exists(file.names),])
        Main.Component = Main.Component[file.exists(file.names)]
        colnames(Components) = tmp
        file.names = file.names[file.exists(file.names)]
      }
    }

  } else if (Pars$read.counts.table != ''){
    sep = Choose.Separator(Pars$read.counts.table)
    rt = read.table(Pars$read.counts.table, stringsAsFactors = FALSE, sep = sep, header = TRUE, check.names = FALSE)
    src.available.samples = colnames(rt)
    needed.samples = colnames(rt)[colnames(rt) %in% all.split..sample.names]
    rt = rt[,needed.samples, drop = FALSE]
    # print(rownames(rt))
    # stop('')
    absent.samples = all.split..sample.names[!(all.split..sample.names %in% colnames(rt))]
    if(length(absent.samples) > 0){
      msg = sprintf('The following samples are absent in file %s (total %d of %d): %s\n Available samples %s\n', Pars$read.counts.table, length(absent.samples),
        length(all.split..sample.names), toString(absent.samples), toString(src.available.samples))
      if(some.samples..should.be.merged)
        stop(sprintf('Sample merging is activated. %s', msg))
      if(perform.classic.paired.test)
        stop(sprintf('Paired test is turned on. All the samples must be present. There are absent ones. table file %s', Pars$read.counts.table))
      if(length(absent.samples) == length(all.split..sample.names))  stop(sprintf('All samples are absent in the expression data file %s', Pars$read.counts.table))
    }

    sample.is.present = all.split..sample.names %in% colnames(rt)
    if(!some.samples..should.be.merged){
      all.split..sample.names = all.split..sample.names[sample.is.present]
      all.src..sample.names = all.src..sample.names[sample.is.present]
    
      schema = schema[sample.is.present, ]
      rownames(schema) = schema$'Sample names'
      rt = rt[, all.split..sample.names]
      tmp = colnames(Components)
      Components = data.matrix(Components[sample.is.present, ])
      Main.Component = Main.Component[sample.is.present]
      colnames(Components) = tmp
    }

  } else stop('Both Pars$read.counts.table and Pars$counts.dir are not specified')

  if(!perform.FC.associations.paired.test){
    eval(parse(text = sprintf("design <- model.matrix(%s)", GLM.formula)))  
    if(perform.classic.paired.test){
      eval(parse(text = sprintf("design.paired <- model.matrix(%s)", GLM.formula.paired)))
    }
  }


  if(is.null(randomization.test.permut.N)){
    Analysis.name = Cleanup.Model.Name(GLM.model)
    if(!is.null(randomization.test.permut.N))  Analysis.name = sprintf('%s - permutation %d', Analysis.name, randomization.test.permut.N)
    Analysis.Data$Analysis.name = Analysis.name
    current.GLM.results.dir = sprintf("%s/~ %s, results",Pars$results.dir, Analysis.name)
    Analysis.Data$current.GLM.results.dir = current.GLM.results.dir
    dir.create(current.GLM.results.dir, showWarnings = FALSE)
    setwd(current.GLM.results.dir)
    
  } else {
    Analysis.name.root = Cleanup.Model.Name(GLM.model)
    Analysis.name = sprintf('%s - permutation %d', Analysis.name.root, randomization.test.permut.N)
    Analysis.Data$Analysis.name = Analysis.name.root
    current.GLM.results.dir.root = sprintf("%s/~ %s, results",Pars$results.dir, Analysis.name.root)
    dir.create(current.GLM.results.dir.root, showWarnings = FALSE)
    current.GLM.results.dir = sprintf("%s/~ %s, results/perm. %d",Pars$results.dir, Analysis.name.root, randomization.test.permut.N)
    dir.create(current.GLM.results.dir, showWarnings = FALSE)
    Analysis.Data$current.GLM.results.dir = current.GLM.results.dir
    setwd(current.GLM.results.dir)
  }

  Main.Predictor.name = colnames(Components)[Predictors.count]
  Analysis.Data$Main.Predictor.name = Main.Predictor.name
  Analysis.Data$Main.Component = Main.Component
  
  #use.exact.test.in.binary.predictors = FALSE
  
  if(!(Pars$DE.package %in% c('edger')) & !is.null(gene.bias.factors.file) & (Pars$Adjust.Transcript.length.bias__in.Read.counts | Pars$Adjust.Expression.level.bias__in.Read.counts)){
    stop('Bias factor adjustemnt is available only for edgeR')
  }

  if(Pars$DE.package == 'edger'){
    if(Pars$counts.dir != ''){
      samples_counts = suppressMessages(readDGE(file.names))
      counts = samples_counts$counts
      # naming 'counts'
      colnames(counts) = all.split..sample.names  #schema$'Sample names'

    } else if (Pars$read.counts.table != ''){
      sep = Choose.Separator(Pars$read.counts.table)
      rt = read.table(Pars$read.counts.table, stringsAsFactors = FALSE, sep = sep, header = TRUE, check.names = FALSE, comment.char = '#')
      # stop
      # if(all(!is.na(as.numeric(as.character(rownames(rt))))))
      #   stop(sprintf('Incorrect rownames in file "%s". Check the length of header (must be equal to the number of samples - 1)', Pars$read.counts.table))
      # needed.samples = colnames(rt)[colnames(rt) %in% schema$'Sample names']
      rt = rt[ , all.split..sample.names]
      # absent.samples = schema$'Sample names'[!(schema$'Sample names' %in% colnames(rt))]
      # if(length(absent.samples) > 0){
      #   cat(sprintf('The following samples are absent in file %s (total %d of %d): %s\n', Pars$read.counts.table, length(absent.samples),
      #     length(schema$'Sample names'), toString(absent.samples)))
      # }
      # if(length(absent.samples) == length(schema$'Sample names'))  stop(sprintf('All samples are ABSENT in %s', Pars$read.counts.table))
      # schema = schema[schema$'Sample names' %in% colnames(rt), ]
      # rownames(schema) = schema$'Sample names'
      # rt = rt[, schema$'Sample names']
      samples_counts = DGEList(rt)
      counts = samples_counts$counts
      # naming 'counts'
      colnames(counts) = all.split..sample.names  #schema$'Sample names'

    } else stop('Both Pars$read.counts.table and Pars$counts.dir are not specified')

    
    # calculating stats and excluding meta-tags
    meta.tags = rownames(counts)[which(startsWith(rownames(counts), '__'))]
    # ("__no_feature","__ambiguous",'__too_low_aQual','__not_aligned','__alignment_not_unique')
    # meta.tags = meta.tags[meta.tags %in% rownames(counts)]
    if (length(meta.tags) > 0){
      mapping.stats.table = array(dim = c(dim(counts)[2],length(meta.tags) + 2))
      rownames(mapping.stats.table) = all.split..sample.names  #schema$'Sample names'
      colnames(mapping.stats.table) = c(meta.tags,'total hits','useful reads')
      for (mt in meta.tags){
        mapping.stats.table[all.split..sample.names, mt] = counts[mt, all.split..sample.names]
      }
      mapping.stats.table[all.split..sample.names, 'total hits'] = apply(counts,2,sum)
      
      noint = rownames(counts) %in% meta.tags
      counts = counts[!noint,]
      mapping.stats.table[all.split..sample.names, 'useful reads'] = apply(counts,2,sum)
      write.table.mod(mapping.stats.table, 'mapping stats.tsv',sep='\t')
    }
    
    # excluding low expression genes
    cpms = cpm(counts)
    if (Pars$min.samples.with.sufficient.CPM == '{auto}') {
      if (Single.binary.predictor){
        min.samples.with.sufficient.CPM = Calc.min.high.CPM.samples.from.total.samples.count(dim(schema)[1],Single.binary.predictor__minimal.group.size)
      } else {
        min.samples.with.sufficient.CPM = Calc.min.high.CPM.samples.from.total.samples.count(dim(schema)[1])
      }
    } else {
      min.samples.with.sufficient.CPM = Pars$min.samples.with.sufficient.CPM
      if(min.samples.with.sufficient.CPM > ncol(counts))
        stop(sprintf('min.samples.with.sufficient.CPM argument (%d) is greater than the number of samples (%d)', min.samples.with.sufficient.CPM, ncol(counts)))
    }
    
    
    if(Pars$sufficient.CPM.to.analyze %in% c('{auto}', 'auto')){
      min.CPM.single.sample = max(1, 3e+7 / mean(colSums(counts)) / 1.25 ) * Pars$sufficient.CPM..autoadjust..multiplier
    } else {
      min.CPM.single.sample = Pars$sufficient.CPM.to.analyze
    }
    cat(sprintf('\nSetting CPM limit %.1f for at least %d samples', min.CPM.single.sample, min.samples.with.sufficient.CPM))
    
    keep = rowSums(cpms > min.CPM.single.sample) >= min.samples.with.sufficient.CPM
    cat(sprintf('\n%d of %d genes passed CPM threshold\n',sum.mod(keep),length(keep)))
    counts = counts[keep,]
    
    # plotting histograms - before normalization
      
    .pardefault <- par(no.readonly = TRUE)
    
    plot.subtitle = ''
    if (length(levels(factor(Main.Component))) == 2) plot.subtitle = sprintf('color indicates %s (quantity; blue-orange)',Main.Predictor.name)
    if (length(levels(factor(Main.Component))) == 3) plot.subtitle = sprintf('color indicates %s (quantity; blue-green-orange)',Main.Predictor.name)
    if (length(levels(factor(Main.Component))) >= 4) plot.subtitle = sprintf('color indicates %s (quantity; blue-green-orange gradient)',Main.Predictor.name)
    
    col.variety <- colorRampPalette(c("#0b00c4","#02babc","#25ad00","#cfcd00","#ff8400"))(512)
    Main.Component.mod = Main.Component - min(Main.Component)
    cl = col.variety[1 + as.integer(Main.Component.mod/max(Main.Component.mod)*511)]

    if(Pars$minimal.Transcript.length > 0){
      if(is.null(gene.bias.factors.file))  stop('parameter minimal.Transcript.length is set > 0, but gene.bias.factors.file is NULL')
      counts = Filter.counts.by.transcript.length(counts.table = counts, gene.bias.factors.file = gene.bias.factors.file,
                    min.length = Pars$minimal.Transcript.length, max.length = NULL)
    }
    
    adj.counts.table = NULL
    # if(!exclude.bins__in.Read.counts.adjustemnt___BF)  max.bins.to.exclude.num___BF = 0
    results.prefix = ''
    if(!is.null(gene.bias.factors.file) & length(Pars$Bias.factors.to.analyze) > 0){
      if(!Pars$Adjust.Transcript.length.bias__in.Read.counts){
        if(is.null(randomization.test.permut.N))
          tmp = analyze.LogCPM.each.Sample.vs.Bias.Factor.associations(counts.table = counts, out.dir = current.GLM.results.dir,
            forced.parameters = Pars, gene.bias.factors.file = gene.bias.factors.file, bias.factors = Pars$Bias.factors.to.analyze,
            adjustable.bias.factor = NULL, min.CPM.for.association.analysis = 1.0,
            render.not.adjusted.accos.plots = Pars$Draw.bias.associations.plots,
            render.adjustment.CPM.density.plots = Pars$Draw.CPM..to..bias.factor.density.plots,
            render.adjustment.CPM.density.plots.all.genes = Pars$Draw.CPM..to..bias.factor.density.plots,
            script.dir = script.dir, results.prefix = '',
            forced.Analysis.name = Analysis.name, include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
      
      } else {
        adj.counts.table = counts

        if(Pars$Adjust.Transcript.length.bias__in.Read.counts){
          aBF = 'Avg. Transcript Length'
          if(Pars$Adjust.Transcript.length.bias__in.Read.counts__exclude.bins){
            max.bins.to.exclude.num =                              Pars$Adjust.Transcript.length.bias__in.Read.counts__max.bins.to.exclude__count
            max.bins.to.exclude.percent =                          Pars$Adjust.Transcript.length.bias__in.Read.counts__max.bins.to.exclude__percentage
            min.cpms.stdev_abs_value.in.bin.to.exclude =           Pars$Adjust.Transcript.length.bias__in.Read.counts__StDev.abs.value__in.bin.to.exclude
            min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = Pars$Adjust.Transcript.length.bias__in.Read.counts__StDev.ratio__in.bin.to.exclude
          } else {
            max.bins.to.exclude.num = 0
            max.bins.to.exclude.percent = 0
            min.cpms.stdev_abs_value.in.bin.to.exclude  = 100
            min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = 100
          }

          adj.counts.table = adj.counts.table[!(rownames(adj.counts.table) %in% 'dummy'),]
          adj.counts.table = 
            analyze.LogCPM.each.Sample.vs.Bias.Factor.associations(counts.table = adj.counts.table, out.dir = current.GLM.results.dir,
              forced.parameters = Pars, gene.bias.factors.file = gene.bias.factors.file, bias.factors = Pars$Bias.factors.to.analyze,
              adjustable.bias.factor = aBF, min.CPM.for.association.analysis = 1.0,
              render.not.adjusted.accos.plots = (is.null(randomization.test.permut.N) & Pars$Draw.bias.associations.plots),
              render.adjustment.CPM.density.plots.all.genes = (is.null(randomization.test.permut.N) & Pars$Draw.CPM..to..bias.factor.density.plots),
              render.adjustment.CPM.density.plots = (is.null(randomization.test.permut.N) & Pars$Draw.CPM..to..bias.factor.density.plots),
              script.dir = script.dir, results.prefix = results.prefix,
              forced.Analysis.name = Analysis.name, include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
              min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude,
              min.cpms.stdev_abs_value.in.bin.to.exclude = min.cpms.stdev_abs_value.in.bin.to.exclude,
              max.bins.to.exclude.num = max.bins.to.exclude.num,
              max.bins.to.exclude.percent = max.bins.to.exclude.percent)
          results.prefix = sprintf('Adjusted to %s - %s', aBF, results.prefix)
        }

        tmp = analyze.LogCPM.each.Sample.vs.Bias.Factor.associations(counts.table = adj.counts.table, out.dir = current.GLM.results.dir,
          forced.parameters = Pars, gene.bias.factors.file = gene.bias.factors.file, bias.factors = Pars$Bias.factors.to.analyze,
          adjustable.bias.factor = NULL, min.CPM.for.association.analysis = 1.0,
          render.not.adjusted.accos.plots = (is.null(randomization.test.permut.N) & Pars$Draw.bias.associations.plots),
          render.adjustment.CPM.density.plots.all.genes = (is.null(randomization.test.permut.N) & Pars$Draw.CPM..to..bias.factor.density.plots),
          render.adjustment.CPM.density.plots = (is.null(randomization.test.permut.N) & Pars$Draw.CPM..to..bias.factor.density.plots),
          script.dir = script.dir, results.prefix = results.prefix,
          forced.Analysis.name = Analysis.name, include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
      }
    }

    if(Pars$Adjust.Expression.level.bias__in.Read.counts){
      if(!is.null(adj.counts.table)){
        src.pseudo.counts.table = adj.counts.table[!(rownames(adj.counts.table) %in% 'dummy'),]
      } else   src.pseudo.counts.table = counts
      
      if(Pars$Adjust.Expression.level.bias__in.Read.counts__exclude.bins){
        max.bins.to.exclude.num =                              Pars$Adjust.Expression.level.bias__in.Read.counts__max.bins.to.exclude__count
        max.bins.to.exclude.percent =                          Pars$Adjust.Expression.level.bias__in.Read.counts__max.bins.to.exclude__percentage
        min.cpms.stdev_abs_value.in.bin.to.exclude =           Pars$Adjust.Expression.level.bias__in.Read.counts__StDev.abs.value__in.bin.to.exclude
        min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = Pars$Adjust.Expression.level.bias__in.Read.counts__StDev.ratio__in.bin.to.exclude
      } else {
        max.bins.to.exclude.num = 0
        max.bins.to.exclude.percent = 0
        min.cpms.stdev_abs_value.in.bin.to.exclude  = 100
        min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = 100
      }

      adj.counts.table = Normalize.read.counts.by.CPM.bins(src.pseudo.counts.table, Startup.Data,
                                 out.dir = current.GLM.results.dir, forced.parameters = Pars,
                                 bins.count = 8, high.expression.bins.count = 3, std.to.HE.bins.size.ratio = 2.5,
                                 render.adjustment.CPM.density.plots = (is.null(randomization.test.permut.N) & Pars$Draw.CPM..to..bias.factor.density.plots),
                                 render.adjustment.CPM.density.plots.all.genes = (is.null(randomization.test.permut.N) & Pars$Draw.CPM..to..bias.factor.density.plots),
                                 script.dir = script.dir, results.prefix = results.prefix, forced.Analysis.name = Analysis.name,
                                 include.GLM.model.name.in.Result.names = FALSE,
                                 min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude,
                                 min.cpms.stdev_abs_value.in.bin.to.exclude = min.cpms.stdev_abs_value.in.bin.to.exclude,
                                 max.bins.to.exclude.num = max.bins.to.exclude.num,
                                 max.bins.to.exclude.percent = max.bins.to.exclude.percent)
    }


    if(is.null(randomization.test.permut.N)){
      cpms = cpm(counts)
      Write.Density.table(cpms, tsv.file.name = 'CPM density, non-normalized.tsv')
      Render.Density.Plots__color.by.condition(cpms, cl, out.dir = NULL, plot.title = 'Log2(CPM) density, non-normalized', plot.subtitle = plot.subtitle,
        png.file.name = 'CPM density, non-normalized, color by condition.png')
      Render.Density.Plots__color.by.sample(cpms, out.dir = NULL, plot.title = 'Log2(CPM) density, non-normalized', plot.subtitle = '',
        png.file.name = 'CPM density, non-normalized, color by sample.png')
    }    

    # creating DGEList object
    # Calculating Norm Factors
    
    if(some.samples..should.be.merged){
      d = DGEList(counts=counts)
    } else if(Single.binary.predictor & Pars$use.exact.test.in.binary.predictors){
      d = DGEList(counts=counts, group = Main.Component)
    } else {
      d = DGEList(counts=counts)
    }
    
    #RNA.Seq.norm.method = 'RLE'  #,"TMM","upperquartile","none","RLE"
    if(Pars$RNA.Seq.norm.method == 'RG'){
      d = Normalize..Using.Reference.Genes(d, Startup.Data$Reference.genes)
    } else {
      d = calcNormFactors(d, method=Pars$RNA.Seq.norm.method) #,"TMM","upperquartile","none","RLE"
    }
    
    
    Analysis.Data$edgeR_d = d
    
    # plotting histograms - after normalization
    cpms = cpm(d, normalized.lib.sizes = TRUE)
    if(is.null(randomization.test.permut.N)){
      Write.Density.table(cpms, tsv.file.name = sprintf('CPM density, normalized (%s).tsv', Pars$RNA.Seq.norm.method))
      Render.Density.Plots__color.by.sample(cpms, out.dir = NULL, plot.title = sprintf('Log2(CPM) density, normalized (%s)', Pars$RNA.Seq.norm.method), plot.subtitle = '',
        png.file.name = sprintf('CPM density, normalized (%s), color by sample.png', Pars$RNA.Seq.norm.method))
    }    
    Render.Density.Plots__color.by.condition(cpms, cl, out.dir = NULL, plot.title = sprintf('Log2(CPM) density, normalized (%s)', Pars$RNA.Seq.norm.method), plot.subtitle = plot.subtitle,
      png.file.name = sprintf('CPM density, normalized (%s), color by condition.png', Pars$RNA.Seq.norm.method))
    

    
    # plotting histograms - after adjustment
    if(!is.null(adj.counts.table)){
      adj.counts.table.cl = adj.counts.table[!(rownames(adj.counts.table) %in% 'dummy'),]
      adj.cpms = adj.counts.table.cl * 1e+6 / mean(colSums(adj.counts.table.cl))
      if(is.null(randomization.test.permut.N)){
        Write.Density.table(adj.cpms, tsv.file.name = sprintf('CPM density, normalized (%s).tsv', Pars$RNA.Seq.norm.method))
        Render.Density.Plots__color.by.sample(adj.cpms, out.dir = NULL, plot.title = sprintf('Log2(CPM) density,\nadjusted and normalized (%s)', Pars$RNA.Seq.norm.method), plot.subtitle = '',
          png.file.name = sprintf('CPM density, adjusted and normalized (%s), color by sample.png', Pars$RNA.Seq.norm.method))
      }      
      Render.Density.Plots__color.by.condition(adj.cpms, cl, out.dir = NULL, plot.title = sprintf('Log2(CPM) density,\nadjusted and normalized (%s)', Pars$RNA.Seq.norm.method), plot.subtitle = plot.subtitle,
        png.file.name = sprintf('CPM density, adjusted and normalized (%s), color by condition.png', Pars$RNA.Seq.norm.method))
      # cat(sprintf('\n\nUsing counts table that is adjusted to the following bias factors: %s%s\n\n', toString(bias.factors__to.adjust__read.counts),
      #   ifelse(adjust.to.Expression.Level, ", 'absolute' expression level (LogCPM)", "")))
      rm(adj.counts.table.cl)
    }
    par(.pardefault)
  }

  if(!perform.FC.associations.paired.test & !is.null(gene.bias.factors.file) & Pars$Adjust.Transcript.length.bias__in.LogFC & Pars$Adjust.Transcript.length.bias__in.LogFC__mode == 'before GLM'){
    if(some.samples..should.be.merged)
      stop(sprintf('Merging samples is not compartible with Adjust.Transcript.length.bias__in.LogFC__mode mode "before GLM"'))
    
    if(!is.null(adj.counts.table)){
      src.pseudo.counts.table = adj.counts.table[!(rownames(adj.counts.table) %in% 'dummy'),]
    } else {
      if(Pars$DE.package != 'edger')     stop('LogFC bias adjustment (before GLM) is available only for edgeR package')
      src.pseudo.counts.table = t(t(d$counts) / d$samples$lib.size * mean(d$samples$lib.size) / d$samples$norm.factors )
      src.pseudo.counts.table = src.pseudo.counts.table
    }

    if('Expression level' %in% Pars$Bias.factors.to.analyze) {      bias.factors = Pars$Bias.factors.to.analyze
    } else bias.factors = c(Pars$Bias.factors.to.analyze, 'Expression level')

    res.counts.table = Associations.LogFC.with.Bias.Factor__before.GLM(Startup.Data, src.pseudo.counts.table, Main.Component, Analysis.name, gene.bias.factors.file = gene.bias.factors.file,
                                       bias.factors = bias.factors, out.dir = NULL, LogCPM.limits = c(0, 3, 5, 7), main.LogCPM.limit = 3,
                                       adjustable.bias.factor = 'Avg. Transcript Length', additional.plots = is.null(randomization.test.permut.N), script.dir = script.dir,
                                       db.suffix = '', plots.suffix = '', include.GLM.model.name.in.Result.names = FALSE,
                                       LogFC.trimming = 0.03, const.add.read.counts = 0.5, trim = (Pars$lowest.expression.trimming.percents + Pars$highest.expression.trimming.percents) / 100 / 2)
    adj.counts.table = round(res.counts.table)
    cs = colSums(adj.counts.table)
    dummy.row = t(as.data.frame(max(cs) - cs))
    rownames(dummy.row) = 'dummy'
    adj.counts.table = rbind(adj.counts.table, dummy.row)
  }
  
  simple.DE.mode = FALSE
  if(Pars$DE.package !='edger' & perform.FC.associations.paired.test)
    stop('perform.FC.associations.paired.test works only with edgeR package')

  if(Pars$DE.package =='edger'){

    # Creating MDS/PCA plots, Calculating distance matrices
    C.Create.Distance.matrices = Pars$Create.Distance.matrices
    if (Pars$Create.Distance.matrices == '{auto}'){
      C.Create.Distance.matrices = (dim(cpms)[2] <= 50)
    }
    
    Create.MDS.plots(d, Main.Component, Main.Predictor.name, out.dir = NULL, prefix = '')
    if (C.Create.Distance.matrices & is.null(randomization.test.permut.N))      Create.Distance.matrices.and.Dendrograms(Pars, cpms, schema, Main.Component, GLM.model, Current.Predictors.list, prefix = '', verbose = verbose)


    if(!is.null(adj.counts.table)){
      if(some.samples..should.be.merged){
        d = DGEList(counts=adj.counts.table)
      } else if(Single.binary.predictor & Pars$use.exact.test.in.binary.predictors){
        d = DGEList(counts=adj.counts.table, group = Main.Component)
      } else {
        d = DGEList(counts=adj.counts.table)
      }
      
      if(Pars$RNA.Seq.norm.method == 'RG'){
        d = Normalize..Using.Reference.Genes(d, Startup.Data$Reference.genes)
      } else {
        d = calcNormFactors(d, method = 'none') #,"TMM","upperquartile","none","RLE"
      }
      Analysis.Data$edgeR_d = d
      
      ##  for adjusted read counts
      Create.MDS.plots(d, Main.Component, Main.Predictor.name, out.dir = NULL, prefix = 'Bias-adjusted ')
      if (C.Create.Distance.matrices & is.null(randomization.test.permut.N))      Create.Distance.matrices.and.Dendrograms(Pars, adj.cpms, schema, Main.Component, GLM.model, Current.Predictors.list, prefix = 'Bias adjusted ', verbose = verbose)
    }    

    
    if(some.samples..should.be.merged){
      comma.separated..sample.names = all.src..sample.names
      for(x in 1:10)
        comma.separated..sample.names = gsub(sprintf('%s ', Pars$Merge.samples..separator), Pars$Merge.samples..separator, comma.separated..sample.names)
      
      merged.counts..list = vector(mode = 'list')
      sample.names..list = strsplit(comma.separated..sample.names, split = ',')
      n = 1
      for(n in 1:length(comma.separated..sample.names)){
        if(length(sample.names..list[[n]]) == 1){
          merged.counts..list[[all.src..sample.names[n]]] = d$counts[, sample.names..list[[n]], drop = TRUE]
          
        } else {
          cat(sprintf('merging samples %s into "%s" ....  \n', paste(sample.names..list[[n]], collapse = ' + '), all.src..sample.names[n]))
          lib.sizes = d$samples[sample.names..list[[n]], 'lib.size', drop = TRUE]
          norm.factors = d$samples[sample.names..list[[n]], 'norm.factors', drop = TRUE]
          merged.counts..list[[all.src..sample.names[n]]] = 
            apply(d$counts[, sample.names..list[[n]], drop = FALSE], 1, function(x) round(sum(x / norm.factors / lib.sizes * sum(lib.sizes)) / length(x)))
        }
      }
      
      # d = DGEList(counts = as.data.frame(merged.counts..list, check.names = FALSE))
      
      if(Single.binary.predictor & Pars$use.exact.test.in.binary.predictors){    d = DGEList(as.data.frame(merged.counts..list, check.names = FALSE), group = Main.Component)
      } else   d = DGEList(counts = as.data.frame(merged.counts..list, check.names = FALSE))
      
      if(Pars$RNA.Seq.norm.method == 'RG'){
        d = Normalize..Using.Reference.Genes(d, Reference.genes)
      } else {
        d = calcNormFactors(d, method = Pars$RNA.Seq.norm.method) #,"TMM","upperquartile","none"
      }
      Analysis.Data$edgeR_d = d
      
      cat('Merging completed.\n')
      
    }
    
    
    
    if(perform.FC.associations.paired.test){
      norm.CPMs = cpm(d, normalized.lib.sizes = TRUE)
      if(Pars$report.norm.read.counts.instead.of.CPM){
        lib.sizes = d$samples[,'lib.size']
        norm.factors = d$samples[,'norm.factors']
        norm.counts = t(t(d$counts) / lib.sizes / norm.factors * mean(lib.sizes))
      }

      LogFCs.Table..init = (norm.CPMs[, (ncol(norm.CPMs)/2 + 1):ncol(norm.CPMs), drop = FALSE] + 
        Pars$paired.LogFCs.calculating..CPM.const.add) / 
        (norm.CPMs[, 1:(ncol(norm.CPMs)/2), drop = FALSE] + Pars$paired.LogFCs.calculating..CPM.const.add)

      LogFCs.Table..init = log2(LogFCs.Table..init)
      colnames(LogFCs.Table..init) = samples..base.names[1:(length(samples..base.names)/2)]

      Main.Component..LogFC = Main.Component[1:(length(Main.Component)/2)]
      Main.Component..LogFC..levels = Main.Component..LogFC %>% DeDup %>% .[order(.)]

      lev1..LogFC..cols = Main.Component..LogFC %in% Main.Component..LogFC..levels[1]
      lev2..LogFC..cols = Main.Component..LogFC %in% Main.Component..LogFC..levels[2]

      gr.1..mean.LogFC = apply(LogFCs.Table..init[, lev1..LogFC..cols, drop = FALSE], 1, mean)
      gr.2..mean.LogFC = apply(LogFCs.Table..init[, lev2..LogFC..cols, drop = FALSE], 1, mean)
      
      Delta.LogFC = gr.2..mean.LogFC - gr.1..mean.LogFC
      
      LogCPMs = log2((apply(norm.CPMs, 1, sum) / dim(norm.CPMs)[2]) + 0.01)
      
      if(sum.mod(lev1..LogFC..cols) >= Pars$add.Trimmed.LogFC.if.Samples.N.is.greater.than |
         sum.mod(lev2..LogFC..cols) >= Pars$add.Trimmed.LogFC.if.Samples.N.is.greater.than){
        trim = (Pars$lowest.expression.trimming.percents + Pars$highest.expression.trimming.percents) / 100 / 2
        Trimmed.Delta.LogFC = apply(LogFCs.Table..init[, lev2..LogFC..cols, drop = FALSE], 1, function(x) mean(x, trim = trim)) - 
          apply(LogFCs.Table..init[, lev1..LogFC..cols, drop = FALSE], 1, function(x) mean(x, trim = trim))
        Stats.table = cbind(as.data.frame(Delta.LogFC), as.data.frame(Trimmed.Delta.LogFC), as.data.frame(gr.1..mean.LogFC), as.data.frame(gr.2..mean.LogFC), as.data.frame(LogCPMs))
        colnames(Stats.table) = c('delta LogFC','trimmed delta LogFC', 'gr.1 mean LogFC', 'gr.2 mean LogFC', 'LogCPM')
      } else {
        Trimmed.Delta.LogFC = NULL
        Stats.table = cbind(as.data.frame(Delta.LogFC), as.data.frame(gr.1..mean.LogFC), as.data.frame(gr.2..mean.LogFC), as.data.frame(LogCPMs))
        colnames(Stats.table) = c('delta LogFC', 'gr.1 mean LogFC', 'gr.2 mean LogFC', 'LogCPM')
      }
      
      Stats.table = Stats.table[rev(order(abs(Stats.table$'delta LogFC'))), ]
      #rownames(Stats.table)
      
    
    } else if(dim(design)[1] < 3){
      simple.DE.mode = TRUE
      msg = sprintf('Sample setup with only %d samples is found. Dispersion, p-value, FDR and LR values cannot be calculated. Score-based P-value "evaluation" will be performed instead.',dim(design)[1])
      #cat(msg)
      warning(msg)
      if(!(all(Main.Component == c(0,1)) || all(Main.Component == c(1,0)))){
        stop('Only comparison of two samples is supported in "simple" mode')
      }
      cpm.add = 0.4
      norm.CPMs = cpm(d, normalized.lib.sizes = TRUE)
      if(Pars$report.norm.read.counts.instead.of.CPM){
        lib.sizes = d$samples[,'lib.size']
        norm.factors = d$samples[,'norm.factors']
        norm.counts = t(t(d$counts) / lib.sizes / norm.factors * mean(lib.sizes))
      }
      LogFCs = log2((norm.CPMs[,2] + cpm.add)/(norm.CPMs[,1] + cpm.add))
      LogCPMs = log2((apply(norm.CPMs,1,sum)/dim(norm.CPMs)[2]) + 0.01)
      
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
        d.ET = estimateDisp(d)
      }
      if(Pars$use.likelihood.ratio.test | Pars$use.quasi.likelihood.test){
        d.GLM = estimateDisp(d, design)
        #d.GLM = estimateGLMCommonDisp(d, design)
        #d.GLM = estimateGLMTrendedDisp(d.GLM, design)
        #d.GLM = estimateGLMTagwiseDisp(d.GLM, design)
      }
      if(Pars$use.likelihood.ratio.test..paired | Pars$use.quasi.likelihood.test..paired){
        d.GLM.paired = estimateDisp(d, design.paired)
      }
      
      if(Pars$Main.test == 'ET'){
        d = d.ET
      } else {
        d = d.GLM
      }
      
      
      # Plotting Biological coefficient of variation plot
      if(is.null(randomization.test.permut.N)){
        png.mod(filename = 'MeanVar.png',units="in",res=600,
            width=8.23,height=6.3,pointsize=12)
        plotMeanVar(d, show.tagwise.vars = TRUE, NBline = TRUE)
        dev.off()
        png.mod(filename = 'BCV.png',units="in",res=600,
            width=8.23,height=6.3,pointsize=12)
        plotBCV(d)
        dev.off()
      }

      Stats.table = NULL
      Stats.table.ET = NULL
      Stats.table.LR = NULL
      Stats.table.QL = NULL
      Stats.table.LR..paired = NULL
      Stats.table.QL..paired = NULL
      
      if(Single.binary.predictor & Pars$use.exact.test.in.binary.predictors){
        ## exact test
        cat('\nPerforming Exact Tests For Differences Between Two Groups Of Negative-Binomial Counts\n')
        et <- exactTest(d.ET)
        tt <- topTags(et, n=nrow(d.ET))
        Stats.table.ET = tt$table
        rownames(Stats.table.ET) = rownames(tt$table)
        colnames(Stats.table.ET)[which(colnames(Stats.table.ET) %in% 'PValue')] = 'p (ex. test)'
        colnames(Stats.table.ET)[which(colnames(Stats.table.ET) %in% 'FDR')] = 'FDR (ex. test)'
      }

      if (Pars$use.quasi.likelihood.test){
        ## GLM approximation with LR test
        cat('\nPerforming Genewise Quasi-Likelihood Test (Negative Binomial distribution; Generalized Linear Models)\n')
        fit <- glmQLFit (d.GLM, design)
        QLFT <- glmQLFTest(fit)
        tt <- topTags(QLFT, n=nrow(d.GLM))
        Stats.table.QL = tt$table
        colnames(Stats.table.QL)[which(colnames(Stats.table.QL) %in% 'PValue')] = 'p (QLF test)'
        colnames(Stats.table.QL)[which(colnames(Stats.table.QL) %in% 'FDR')] = 'FDR (QLF test)'
      }
      
      if (Pars$use.likelihood.ratio.test){
        ## GLM approximation
        cat('\nPerforming Genewise Likelihood ratio (LR) test (Negative Binomial distribution; Generalized Linear Models)\n')
        fit <- glmFit(d.GLM, design)
        LRT <- glmLRT(fit)
        #lrt <- glmLRT(fit,coef=1)
        tt <- topTags(LRT, n=nrow(d.GLM))
        Stats.table.LR = tt$table
        colnames(Stats.table.LR)[which(colnames(Stats.table.LR) %in% 'PValue')] = 'p (LR test)'
        colnames(Stats.table.LR)[which(colnames(Stats.table.LR) %in% 'FDR')] = 'FDR (LR test)'
      }
      
      if (Pars$use.quasi.likelihood.test..paired){
        ## GLM approximation with LR test
        cat('\nPerforming Genewise Negative Binomial Generalized Linear Models With Quasi-Likelihood Tests (for paired samples)\n')
        fit <- glmQLFit (d.GLM.paired, design.paired)
        QLFT <- glmQLFTest(fit)
        tt <- topTags(QLFT, n=nrow(d.GLM.paired))
        Stats.table.QL..paired = tt$table
        colnames(Stats.table.QL..paired)[which(colnames(Stats.table.QL..paired) %in% 'PValue')] = 'p (QLF test paired)'
        colnames(Stats.table.QL..paired)[which(colnames(Stats.table.QL..paired) %in% 'FDR')] = 'FDR (QLF test paired)'
      }
      
      if (Pars$use.likelihood.ratio.test..paired){
        ## GLM approximation
        cat('\nPerforming Genewise Negative Binomial Generalized Linear Models\n')
        fit <- glmFit(d.GLM.paired, design.paired)
        LRT <- glmLRT(fit)
        #lrt <- glmLRT(fit,coef=1)
        tt <- topTags(LRT, n=nrow(d.GLM.paired))
        Stats.table.LR..paired = tt$table
        colnames(Stats.table.LR..paired)[which(colnames(Stats.table.LR..paired) %in% 'PValue')] = 'p (LR test paired)'
        colnames(Stats.table.LR..paired)[which(colnames(Stats.table.LR..paired) %in% 'FDR')] = 'FDR (LR test paired)'
      }
      
      
      if(Pars$Main.test == 'ET'){
        Stats.table = Stats.table.ET
        if(!is.null(Stats.table.QL)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.QL[rownames(Stats.table), 'p (QLF test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (QLF test)'
        }
        if(!is.null(Stats.table.LR)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.LR[rownames(Stats.table), 'p (LR test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (LR test)'
        }
        if(!is.null(Stats.table.QL..paired)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.QL..paired[rownames(Stats.table), 'p (QLF test paired)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (QLF test paired)'
        }
        if(!is.null(Stats.table.LR..paired)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.LR..paired[rownames(Stats.table), 'p (LR test paired)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (LR test paired)'
        }

      } else if (Pars$Main.test == 'QL'){
        Stats.table = Stats.table.QL
        if(!is.null(Stats.table.LR)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.LR[rownames(Stats.table), 'p (LR test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (LR test)'
        }
        if(!is.null(Stats.table.ET)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.ET[rownames(Stats.table), 'p (ex. test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (ex. test)'
        }
        if(!is.null(Stats.table.QL..paired)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.QL..paired[rownames(Stats.table), 'p (QLF test paired)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (QLF test paired)'
        }
        if(!is.null(Stats.table.LR..paired)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.LR..paired[rownames(Stats.table), 'p (LR test paired)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (LR test paired)'
        }

      } else if (Pars$Main.test == 'LR'){
        Stats.table = Stats.table.LR
        if(!is.null(Stats.table.QL)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.QL[rownames(Stats.table), 'p (QLF test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (QLF test)'
        }
        if(!is.null(Stats.table.ET)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.ET[rownames(Stats.table), 'p (ex. test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (ex. test)'
        }
        if(!is.null(Stats.table.QL..paired)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.QL..paired[rownames(Stats.table), 'p (QLF test paired)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (QLF test paired)'
        }
        if(!is.null(Stats.table.LR..paired)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.LR..paired[rownames(Stats.table), 'p (LR test paired)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (LR test paired)'
        }

      } else if (Pars$Main.test == 'QL paired'){
        Stats.table = Stats.table.QL..paired
        if(!is.null(Stats.table.LR..paired)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.LR..paired[rownames(Stats.table), 'p (LR test paired)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (LR test paired)'
        }
        if(!is.null(Stats.table.QL)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.QL[rownames(Stats.table), 'p (QLF test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (QLF test)'
        }
        if(!is.null(Stats.table.LR)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.LR[rownames(Stats.table), 'p (LR test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (LR test)'
        }
        if(!is.null(Stats.table.ET)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.ET[rownames(Stats.table), 'p (ex. test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (ex. test)'
        }

      } else if (Pars$Main.test == 'LR paired'){
        Stats.table = Stats.table.LR..paired
        if(!is.null(Stats.table.QL..paired)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.QL..paired[rownames(Stats.table), 'p (QLF test paired)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (QLF test paired)'
        }
        if(!is.null(Stats.table.LR)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.LR[rownames(Stats.table), 'p (LR test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (LR test)'
        }
        if(!is.null(Stats.table.QL)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.QL[rownames(Stats.table), 'p (QLF test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (QLF test)'
        }
        if(!is.null(Stats.table.ET)){
          Stats.table = cbind(Stats.table, as.data.frame(Stats.table.ET[rownames(Stats.table), 'p (ex. test)']))
          colnames(Stats.table)[length(colnames(Stats.table))] = 'p (ex. test)'
        }
      }


      norm.CPMs = cpm(d, normalized.lib.sizes = TRUE)
      if(Pars$report.norm.read.counts.instead.of.CPM){
        lib.sizes = d$samples[,'lib.size']
        norm.factors = d$samples[,'norm.factors']
        norm.counts = t(t(d$counts) / lib.sizes / norm.factors * mean(lib.sizes))
      }
      max.samples.in.group = Main.Component %>% .[!duplicated(.)] %>% sapply(., function(x) { sum.mod(Main.Component %in% x)}) %>% max
      if(max.samples.in.group > Pars$add.Trimmed.LogFC.if.Samples.N.is.greater.than && Single.binary.predictor){
        Main.Component.levels <- Main.Component %>% .[!duplicated(.)] %>% .[order(.)]
        lev1.cols = Main.Component %in% Main.Component.levels[1]
        lev2.cols = Main.Component %in% Main.Component.levels[2]
        # apply(norm.CPMs, 1, function(x) {
        #   log2(mean(x[lev2.cols] + Pars$trimmed.LogFC.const.add)/mean(x[lev1.cols] + Pars$trimmed.LogFC.const.add))
        # }) %>% head

        trim = (Pars$lowest.expression.trimming.percents + Pars$highest.expression.trimming.percents) / 100 / 2
        trimmed.LogFC = apply(norm.CPMs, 1, function(x) {
          log2(mean(x[lev2.cols] + Pars$trimmed.LogFC.const.add, trim = trim) / mean(x[lev1.cols] + Pars$trimmed.LogFC.const.add, trim = trim))
        })


        # apply(norm.CPMs, 1, function(x) {
        #   Smoothed_Trimmed_LogFC(x[lev2.cols], x[lev1.cols],
        #                          trim.low = Pars$lowest.expression.trimming.percents / 100,
        #                          trim.high = Pars$highest.expression.trimming.percents / 100,
        #                          const.add = Pars$trimmed.LogFC.const.add)
        # }) %>% head

        Stats.table = cbind(subset(Stats.table, select = c(1)),
              as.data.frame(trimmed.LogFC[rownames(Stats.table)]),
              subset(Stats.table, select = c(2:dim(Stats.table)[2]))
              )
        colnames(Stats.table)[2] = 'trimmed LogFC'
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
    if(is.null(randomization.test.permut.N))  plotMA(res, main="DESeq2", ylim=c(-2,2))
    res = res[order(res$pvalue),]
    
    norm.CPMs = cpm(counts(dds[rownames(res),], normalized=TRUE))
    if(Pars$report.norm.read.counts.instead.of.CPM)
      norm.counts = counts(dds[rownames(res),], normalized=TRUE)
    
    Stats.table = cbind(
      data.frame(res[,'log2FoldChange']),
      data.frame(log2(apply(norm.CPMs,1,mean))),
      data.frame(res[,'stat']*10),
      res[,c('pvalue','padj')])
    colnames(Stats.table) = c('logFC','logCPM','LR','PValue','FDR')
    tt = Stats.table
    
    max.samples.in.group = Main.Component %>% .[!duplicated(.)] %>% sapply(., function(x) { sum.mod(Main.Component %in% x)}) %>% max
    if(max.samples.in.group >= Pars$add.Trimmed.LogFC.if.Samples.N.is.greater.than && Single.binary.predictor){
      Main.Component.levels <- Main.Component %>% .[!duplicated(.)] %>% .[order(.)]
      lev1.cols = Main.Component %in% Main.Component.levels[1]
      lev2.cols = Main.Component %in% Main.Component.levels[2]
      # apply(norm.CPMs, 1, function(x) {
      #   log2(mean(x[lev2.cols] + Pars$trimmed.LogFC.const.add)/mean(x[lev1.cols] + Pars$trimmed.LogFC.const.add))
      # }) %>% head
      
      trimmed.LogFC = apply(norm.CPMs, 1, function(x) {
        log2(mean(x[lev2.cols] + Pars$trimmed.LogFC.const.add, trim = 0.1) / mean(x[lev1.cols] + Pars$trimmed.LogFC.const.add, trim = 0.1))
      })
      
      
      # apply(norm.CPMs, 1, function(x) {
      #   Smoothed_Trimmed_LogFC(x[lev2.cols], x[lev1.cols],
      #                          trim.low = Pars$lowest.expression.trimming.percents / 100,
      #                          trim.high = Pars$highest.expression.trimming.percents / 100,
      #                          const.add = Pars$trimmed.LogFC.const.add)
      # }) %>% head
      
      Stats.table = cbind(subset(Stats.table, select = c(1)),
                          as.data.frame(trimmed.LogFC[rownames(Stats.table)]),
                          subset(Stats.table, select = c(2:dim(Stats.table)[2]))
      )
      colnames(Stats.table)[2] = 'trimmed LogFC'
    }
    
    
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
    #norm.CPMs = cpm(counts(dds,normalized=T))
    
    norm.CPMs = cpm(d, normalized.lib.sizes = TRUE)
    if(Pars$report.norm.read.counts.instead.of.CPM){
      lib.sizes = d$samples[,'lib.size']
      norm.factors = d$samples[,'norm.factors']
      norm.counts = t(t(d$counts) / lib.sizes / norm.factors * mean(lib.sizes))
    }
    
    trimming = 0.03
    
    Cor.table = data.matrix(array(dim=c(dim(norm.CPMs)[1],9)))
    colnames(Cor.table) = c('Spearman r','p (Spearman)','Pearson r','p (Pearson)','logFC','logCPM','LR','PValue','FDR')
    #Pars$DE.package = 'combined_cor'
    for (j in 1:dim(norm.CPMs)[1]){
      if(Pars$DE.package %in% c('combined_cor', 'spearman')){
        ct = suppressWarnings(cor.test(Main.Component, norm.CPMs[j,],method = 'spearman',alternative = 'two.sided'))
        Cor.table[j,'Spearman r'] = ct$estimate
        Cor.table[j,'p (Spearman)'] = ct$p.value
      }
      if(Pars$DE.package %in% c('combined_cor', 'pearson')){
        ct = suppressWarnings(cor.test(Main.Component, norm.CPMs[j,],method = 'pearson',alternative = 'two.sided'))
        Cor.table[j,'Pearson r'] = ct$estimate
        Cor.table[j,'p (Pearson)'] = ct$p.value
      }
    }
    
    #library(chemometrics)
    if(Pars$DE.package == 'spearman'){
      Cor.table[,'PValue'] = Cor.table[,'p (Spearman)']
      Cor.table[,'logFC'] = 2.5*Cor.table[,'Spearman r'] * log(1 + 4*apply(norm.CPMs, 1, rel_sd_trim))
      #Cor.table[,'logFC'] = log((((1 + 1.5*abs(Cor.table[,'Spearman r']))**2
        #*apply(norm.CPMs,1,rel_sd_trim))**0.33)*3,2)*sign(Cor.table[,'Spearman r'])
    } else if(Pars$DE.package == 'pearson'){
      Cor.table[,'PValue'] = Cor.table[,'p (Pearson)']
      Cor.table[,'logFC'] = 2.5*Cor.table[,'Pearson r'] * log(1 + 4*apply(norm.CPMs, 1, rel_sd_trim))
      #Cor.table[,'logFC'] = log((((1 + 1.5*abs(Cor.table[,'Pearson r']))**2
      #                            *apply(norm.CPMs,1,rel_sd_trim))**0.33)*3,2)*sign(Cor.table[,'Pearson r'])
    } else {  ## combined_cor
      integrate.corr.rs <- function(x){
        if(length(levels(factor(sign(x)))) > 1) return(100)
        return (((abs(x[1])*abs(x[1])*abs(x[2]))^0.33)*sign(x[1]))
      }
      Rs = apply(cbind(Cor.table[,'Spearman r'],Cor.table[,'Pearson r']),1,integrate.corr.rs)  ## 'joint' Rs
      Cor.table[,'PValue'] = apply(cbind(Cor.table[,'p (Spearman)'],Cor.table[,'p (Pearson)']),1,max)
      Cor.table[,'logFC'] = 2.5*Rs * log(1 + 4*apply(norm.CPMs,1,rel_sd_trim))
      #Cor.table[,'logFC'] = log((((1 + 1.5*abs(Rs))**2
      #                            *apply(norm.CPMs,1,rel_sd_trim))**0.33)*3,2)*sign(Rs)
    }
    Cor.table[,'LR'] = abs(Cor.table[,'logFC']) * abs(log(Cor.table[,'PValue'],2))*2
    Cor.table[,'logCPM'] = log2(apply(norm.CPMs,1,mean))
    Cor.table[,'FDR'] = p.adjust(Cor.table[,'PValue'],method = 'BH')

    #res <- results(dds)
    #plotMA(res, main="DESeq2", ylim=c(-2,2))
    #res = res[order(res$pvalue),]

    rownames(Cor.table) = rownames(norm.CPMs)
    Stats.table = Cor.table[order(Cor.table[,'LR'],decreasing = T),c('logFC','logCPM','PValue','FDR')]
    #Stats.table = Stats.table[,!(colnames(Stats.table) %in% 'LR')]
    tt = Stats.table
  }

  if('logFC' %in% colnames(Stats.table))
    colnames(Stats.table)[which(colnames(Stats.table) %in%'logFC')] = 'LogFC'
  
  #################################################################
  ########  Adjusting trimmed LogFC according to bias

  if(!perform.FC.associations.paired.test & !is.null(gene.bias.factors.file) & length(Pars$Bias.factors.to.analyze) > 0 & Pars$Evaluate.bias.factors.Associations){

    # if(logFC.type..to..asociate.with.bias %in% c('auto', '{auto}')){
    preferred_columns = c('trimmed LogFC', 'LogFC', 'logFC', 'raw trimmed LogFC', 'raw LogFC')
    # } else {
    #   preferred_columns = logFC.type..to..asociate.with.bias
    # }

    forced.DE.table = Stats.table[ , c(
      Get.LogFC.col.name(colnames(Stats.table), preferred_columns = preferred_columns),
      Get.PValue.col.name(Pars$Main.test, colnames(Stats.table)),
      Get.LogCPM.col.name(colnames(Stats.table))
      )]
    #colnames(forced.DE.table) = c('LogFC', 'PValue')

    if('Expression level' %in% Pars$Bias.factors.to.analyze) {      bias.factors = Pars$Bias.factors.to.analyze
    } else bias.factors = c(Pars$Bias.factors.to.analyze, 'Expression level')

    if(!Pars$Adjust.Transcript.length.bias__in.LogFC){
      tmp = Associations.LogFC.with.Bias.Factor__after.GLM(Startup.Data, forced.parameters = Pars, forced.DE.table = forced.DE.table,
                                                      gene.bias.factors.file = gene.bias.factors.file, bias.factors = bias.factors,
                                                      additional.plots = is.null(randomization.test.permut.N),
                                                      adjustable.bias.factor = NULL, script.dir = script.dir,
                                                      forced.Analysis.name = Analysis.name, db.suffix = '',
                                                      include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
      rm(tmp)
    } else if(Pars$Adjust.Transcript.length.bias__in.LogFC__mode == 'before GLM'){
      tmp = Associations.LogFC.with.Bias.Factor__after.GLM(Startup.Data, forced.parameters = Pars, forced.DE.table = forced.DE.table,
                                                        gene.bias.factors.file = gene.bias.factors.file, bias.factors = bias.factors,
                                                        additional.plots = is.null(randomization.test.permut.N),
                                                        adjustable.bias.factor = NULL, script.dir = script.dir,
                                                        forced.Analysis.name = sprintf('%s, final', Analysis.name), db.suffix = '.adjusted',
                                                        plots.suffix = ', final',
                                                        include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
      rm(tmp)
    } else if(Pars$Adjust.Transcript.length.bias__in.LogFC__mode == 'after GLM'){
      synchronize.pvalues__in.LogFC.adjustemnt = TRUE
      adjusted.DE.table = Associations.LogFC.with.Bias.Factor__after.GLM(Startup.Data, forced.parameters = Pars, forced.DE.table = forced.DE.table,
                                                        gene.bias.factors.file = gene.bias.factors.file, bias.factors = bias.factors,
                                                        additional.plots = is.null(randomization.test.permut.N),
                                                        adjustable.bias.factor = 'Avg. Transcript Length', script.dir = script.dir, forced.Analysis.name = Analysis.name,
                                                        synchronize.pvalues = synchronize.pvalues__in.LogFC.adjustemnt,
                                                        include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
      if(synchronize.pvalues__in.LogFC.adjustemnt){
        if('trimmed LogFC' %in% colnames(Stats.table)){
          Stats.table = cbind(subset(Stats.table, select = c(1, 2)),
              subset(adjusted.DE.table[rownames(Stats.table),], select = c('TBA LogFC', 'TBA p-value', 'TBA FDR')),
              subset(Stats.table, select = c(3:dim(Stats.table)[2])))
        } else {
          Stats.table = cbind(subset(Stats.table, select = c(1)),
              subset(adjusted.DE.table[rownames(Stats.table),], select = c('TBA LogFC', 'TBA p-value', 'TBA FDR')),
              subset(Stats.table, select = c(2:dim(Stats.table)[2])))
        }
      } else {
        if('trimmed LogFC' %in% colnames(Stats.table)){
          Stats.table = cbind(subset(Stats.table, select = c(1, 2)),
              subset(adjusted.DE.table[rownames(Stats.table),], select = c('TBA LogFC')),
              subset(Stats.table, select = c(3:dim(Stats.table)[2])))
        } else {
          Stats.table = cbind(subset(Stats.table, select = c(1)),
              subset(adjusted.DE.table[rownames(Stats.table),], select = c('TBA LogFC')),
              subset(Stats.table, select = c(2:dim(Stats.table)[2])))
        }
      }

      forced.DE.table = Stats.table[ , c(
        'TBA LogFC',
        Get.PValue.col.name(Pars$Main.test, colnames(Stats.table)),
        Get.LogCPM.col.name(colnames(Stats.table))
        )]
      
      tmp = Associations.LogFC.with.Bias.Factor__after.GLM(Startup.Data, forced.parameters = Pars, forced.DE.table = forced.DE.table,
                                                        gene.bias.factors.file = gene.bias.factors.file, bias.factors = bias.factors,
                                                        additional.plots = is.null(randomization.test.permut.N),
                                                        adjustable.bias.factor = NULL, script.dir = script.dir,
                                                        forced.Analysis.name = sprintf('%s, final', Analysis.name), db.suffix = '.adjusted',
                                                        plots.suffix = ', final', synchronize.pvalues = synchronize.pvalues__in.LogFC.adjustemnt,
                                                        include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
      rm(tmp)
    } else {
      stop('Incorrect Pars$Adjust.Transcript.length.bias__in.LogFC__mode')
    }
  }

  Quantile.Statistics.table = NULL
  if(Pars$include.Expression.level__quantile.statistics & is.null(randomization.test.permut.N)){
    cat('Adding quantile statistics [LogCPM]\n')
    Quantile.Statistics.table = Add.Quantile.statistics__LogCPM(Stats.table, Quantile.Statistics.table)
  }

  
  if(Pars$include.Transcript.length__quantile.statistics & !is.null(gene.bias.factors.file) & is.null(randomization.test.permut.N)){
    cat(sprintf('Adding quantile statistics [Avg. Transcript length]\n'))
    Quantile.Statistics.table = Add.Quantile.statistics__BF(Stats.table = Stats.table, Quantile.Statistics.table = Quantile.Statistics.table, BF = 'Avg. Transcript Length',
                                       QStats.colname = 'Tr.Length rank (%)', gene.bias.factors.file = gene.bias.factors.file, script.dir = script.dir)
  }
      # if(bf == 'Avg. Transcript CG content'){
      #   Quantile.Statistics.table = Add.Quantile.statistics__BF(Stats.table = Stats.table, Quantile.Statistics.table = Quantile.Statistics.table, BF = 'Avg. Transcript CG content',
      #                                      QStats.colname = 'CG content rank (%)', gene.bias.factors.file = '{auto}', script.dir = script.dir)
      # }
  
  if(!is.null(Quantile.Statistics.table) & is.null(randomization.test.permut.N)){
    Stats.table = cbind(Quantile.Statistics.table, Stats.table)
  }

 

  ##############################################################
  ############ Creating signle output results table
  
  # combining FRD, LR, P-value info and Read Count Info
  #Scores = abs(tt$table['logFC'])**0.7 * log(tt$table['FDR'])*(-1)
  #colnames(Scores) = c('Score')
  
  
  cat(sprintf("%s:   calculating correlations, scores and writing expression tables...\n",Analysis.Data$GLM.model))
  
  #all.preds = c()
  #for (tmp.Pred.set in c(Analysis.Data$GLM.model,Pars$Models.to.Test))  all.preds = append(all.preds,Extract.predictors(tmp.Pred.set))
  
  #all.preds = levels(as.factor(all.preds))
  
  if (Pars$add.MDS.dims.as.predictors) {
    Analysis.Data$predictor.values = suppressWarnings(as.data.frame(schema, row.names = rownames(schema))[colnames(norm.CPMs),])
    Analysis.Data$predictor.values.nr = Eliminate.Redundant.Predictors.Info(Analysis.Data$predictor.values)
    all.preds.nr = colnames(Analysis.Data$predictor.values.nr)
    
    if(Pars$report.norm.read.counts.instead.of.CPM){
      Counts.Table = suppressWarnings(rbind(t(as.data.frame(schema, row.names = rownames(schema))[colnames(norm.counts), all.preds.nr]),norm.counts[rownames(Stats.table),]))
    } else {
      Counts.Table = suppressWarnings(rbind(t(as.data.frame(schema, row.names = rownames(schema))[colnames(norm.CPMs), all.preds.nr]),norm.CPMs[rownames(Stats.table),]))
    }
    # Counts.Table = suppressWarnings(rbind(t(data.matrix(schema)[colnames(norm.CPMs),all.preds]),norm.CPMs[rownames(tt),]))
    #} else Counts.Table = rbind(t(as.data.frame(schema,row.names = rownames(schema))[colnames(norm.CPMs),all.preds]),norm.CPMs[rownames(tt),])
  } else {
    all.preds = colnames(schema)[!(colnames(schema) %in% 'Sample names')]
    Analysis.Data$predictor.values = suppressWarnings(data.matrix(schema)[colnames(norm.CPMs), all.preds])
    if (length(all.preds) == 1){
      Analysis.Data$predictor.values = as.data.frame(Analysis.Data$predictor.values)
      colnames(Analysis.Data$predictor.values) = c(all.preds)
    }
    Analysis.Data$predictor.values.nr = Eliminate.Redundant.Predictors.Info(Analysis.Data$predictor.values)
    
    all.preds.nr = colnames(Analysis.Data$predictor.values.nr)
    if(Pars$report.norm.read.counts.instead.of.CPM){
      Counts.Table = suppressWarnings(rbind(t(data.matrix(schema)[colnames(norm.counts), all.preds.nr]), norm.counts[rownames(Stats.table),]))
    } else {
      Counts.Table = suppressWarnings(rbind(t(data.matrix(schema)[colnames(norm.CPMs), all.preds.nr]), norm.CPMs[rownames(Stats.table),]))
    }
  }

  if(perform.classic.paired.test | perform.FC.associations.paired.test){
    LogFCs.Table = (Counts.Table[(length(all.preds.nr) + 1):nrow(Counts.Table), (ncol(Counts.Table)/2 + 1):ncol(Counts.Table), drop = FALSE] + 
      Pars$paired.LogFCs.calculating..CPM.const.add) / 
      (Counts.Table[(length(all.preds.nr) + 1):nrow(Counts.Table), 1:(ncol(Counts.Table)/2), drop = FALSE] + Pars$paired.LogFCs.calculating..CPM.const.add)
    LogFCs.Table = log2(LogFCs.Table)
    LogFCs.Table = suppressWarnings(rbind(Counts.Table[1:length(all.preds.nr), 1:(ncol(Counts.Table)/2)], LogFCs.Table))
    colnames(LogFCs.Table) = samples..base.names[1:(length(samples..base.names)/2)]
  } else {
    LogFCs.Table = data.frame(array(dim = c(nrow(Counts.Table), 0)))
    rownames(LogFCs.Table) = rownames(Counts.Table)
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
  
  
  # Main.Component.levels = Main.Component %>% .[!duplicated(.)] %>% .[order(.)]
  #### may be an issue here!!!!
  Main.Component.reordered = Main.Component
  if(!perform.classic.paired.test & !perform.FC.associations.paired.test){
    Main.Component.reordered = Main.Component[order(Main.Component)]
  }
  Main.Component.levels = Main.Component.reordered[!duplicated(Main.Component.reordered)]
  Main.Component.values..to..group.names = Get..Predictor.values..to..group.names(Main.Component, GLM.model)  

  CPMs.groups = vector(mode = 'list')  ## CPMs.groups are reordered according to order(Main.Component) !!!
  LogFCs.groups = vector(mode = 'list')
  if(length(Main.Component.levels) <= Pars$Max.Levels.to.split.Sparkline.groups){
    if(perform.FC.associations.paired.test | perform.classic.paired.test){
      l = Main.Component.levels[1]
      for(l in Main.Component.levels){
        CPMs.groups[[toString(l)]] = which(Main.Component %in% l)
      }
    } else {
      l = Main.Component.levels[1]
      for(l in Main.Component.levels){
        CPMs.groups[[toString(l)]] = which(Main.Component.reordered %in% l)
      }
    }

    if(perform.classic.paired.test){
      LogFCs.groups[['paired LogFCs']] = 1:(length(Main.Component.reordered)/2)
      Main.Component.values..to..group.names[['paired LogFCs']] = 'paired LogFCs'
    } else if(perform.FC.associations.paired.test){
      Main.Component..LogFC = Main.Component[1:(length(Main.Component)/2)]
      Main.Component..LogFC..levels = Main.Component..LogFC %>% DeDup %>% .[order(.)]
      for(l in Main.Component..LogFC..levels){
        LogFCs.groups[[toString(l)]] = which(Main.Component..LogFC %in% l)
      }
    }
  } else {
    CPMs.groups[['rel. expression']] = 1:(length(Main.Component.reordered))
    if(perform.classic.paired.test | perform.FC.associations.paired.test)
      LogFCs.groups[['paired LogFCs']] = 1:(length(Main.Component.reordered)/2)
    Main.Component.values..to..group.names[['rel. expression']] = 'rel. expression'
    Main.Component.values..to..group.names[['paired LogFCs']] = 'paired LogFCs'
  }







  DeltaFreq.table = data.matrix(array(dim=c(dim(norm.CPMs)[1], 0)))
  
  if(perform.FC.associations.paired.test){
    
    Main.Component..LogFC = Main.Component[1:(length(Main.Component)/2)]
    Main.Component..LogFC..levels = Main.Component..LogFC %>% DeDup %>% .[order(.)]

    Two.Pools.tests.available = TRUE
    #KruskWall.test.available = FALSE
    if(length(Main.Component..LogFC..levels) != 2)  Two.Pools.tests.available = FALSE
    if(!all(sapply(Main.Component..LogFC..levels, function(x) length(Main.Component..LogFC %>% .[. %in% x]) >= 2))) Two.Pools.tests.available = FALSE
    
    cor.columns = c()
    if(Pars$use.ttest.for.FC.association.paired.test & Two.Pools.tests.available){
      cor.columns = c(cor.columns, 'p (t-test, FC)')
      if(Pars$Include.ttest.FDR.for.FC.association.paired.test)  cor.columns = c(cor.columns, 'FDR (t-test, FC)')
    }
    
    if(Pars$use.Mann.Whitney.criterium.for.FC.association.paired.test & Two.Pools.tests.available){
      cor.columns = c(cor.columns, 'p (Mann-Wh., FC)')
      if(Pars$Include.Wilcoxon.FDR.for.FC.association.paired.test)  cor.columns = c(cor.columns, 'FDR (Mann-Wh., FC)')
    }
    
    if(Pars$Include.Spearman.corr){
      cor.columns = c(cor.columns, 'Spearman r (FC)')
      if(Pars$Include.Corr.PValues)
        cor.columns = c(cor.columns, 'p (Spearman, FC)')
    }
    
    if(Pars$Include.Pearson.corr){
      cor.columns = c(cor.columns, 'Pearson r (FC)')
      if(Pars$Include.Corr.PValues)
        cor.columns = c(cor.columns, 'p (Pearson, FC)')
    }
    
    if(Two.Pools.tests.available & length(cor.columns) > 0){
      lev1..LogFC..cols = Main.Component..LogFC %in% Main.Component..LogFC..levels[1]
      lev2..LogFC..cols = Main.Component..LogFC %in% Main.Component..LogFC..levels[2]
    } else {
      stop('Now only t-test and Mann-Whitney test are available for FoldChange-associations test')
    }
    
    Cor.table = data.matrix(array(dim=c(dim(norm.CPMs)[1], length(cor.columns))))
    colnames(Cor.table) = cor.columns
    gene.names = rownames(Stats.table)
    
    j = 1
    total.genes = dim(norm.CPMs)[1]
    for (j in 1:total.genes){
      current.LogFCs = LogFCs.Table..init[rownames(Stats.table)[j], ]
      if(Two.Pools.tests.available){
        current.LogFCs.group1 = current.LogFCs[lev1..LogFC..cols]
        current.LogFCs.group2 = current.LogFCs[lev2..LogFC..cols]
        if(Pars$use.ttest.for.FC.association.paired.test){
          tryCatch(expr = {
            Cor.table[j, 'p (t-test, FC)'] = suppressWarnings(t.test(current.LogFCs.group1, current.LogFCs.group2))$p.value
          }, error = function(err){
            Cor.table[j, 'p (t-test, FC)'] = 1
          })
        }
        if(Pars$use.Mann.Whitney.criterium.for.FC.association.paired.test){
          Cor.table[j, 'p (Mann-Wh., FC)'] = suppressWarnings(wilcox.test(current.LogFCs.group1, current.LogFCs.group2))$p.value
        }
        
      }
      
      if(Pars$Include.Spearman.corr){
        ct = suppressWarnings(cor.test(Main.Component..LogFC, current.LogFCs, method = 'spearman', alternative = 'two.sided'))
        Cor.table[j, 'Spearman r (FC)'] = ct$estimate
        if(Pars$Include.Corr.PValues)
          Cor.table[j, 'p (Spearman, FC)'] = ct$p.value
      }
      
      if(Pars$Include.Pearson.corr){
        ct = suppressWarnings(cor.test(Main.Component..LogFC, current.LogFCs, method = 'pearson', alternative = 'two.sided'))
        Cor.table[j ,'Pearson r (FC)'] = ct$estimate
        if(Pars$Include.Corr.PValues)
          Cor.table[j, 'p (Pearson, FC)'] = ct$p.value
      }
      
      if (j %% 25 == 0 & verbose)
        cat(sprintf('\rCalculating correlations and add. stat. tests (paired FC assoc. mode). Completed %.2f%%       ',j/total.genes*100))
    }
    cat(sprintf('\rCalculating correlations and add. stat. tests (paired FC assoc. mode) - completed                                \n'))    
    
    
  } else {
    Two.Pools.tests.available = TRUE
    KruskWall.test.available = FALSE
    if(length(Main.Component.levels) != 2)  Two.Pools.tests.available = FALSE
    if(!all(sapply(Main.Component.levels, function(x) length(Main.Component %>% .[. %in% x]) >= 2))) Two.Pools.tests.available = FALSE
  
    if(!Two.Pools.tests.available){
      KruskWall.test.available = TRUE
      if(length(Main.Component.levels) <= 2)  KruskWall.test.available = FALSE
      if(!any(sapply(Main.Component.levels, function(x) length(Main.Component %>% .[. %in% x]) >= 2))) KruskWall.test.available = FALSE
    }
    
    if(perform.classic.paired.test){
      if(Pars$Include.Mann.Whitney.FDR){
        cor.columns = c('p (Wilcoxon paired)', 'FDR (Wilcoxon paired)')
      } else {
        cor.columns = c('p (Wilcoxon paired)')
      }
    } else cor.columns = c()
  
    if(Two.Pools.tests.available){
      if(Pars$Include.Mann.Whitney.FDR){
        cor.columns = c(cor.columns, 'p (Mann-Wh.)', 'FDR (Mann-Wh.)', 'p (t-test)')
      } else {
        cor.columns = c(cor.columns, 'p (Mann-Wh.)', 'p (t-test)')
      }
    } else if (KruskWall.test.available){
      if(Pars$Include.Kruskal.Wallis.FDR){
        cor.columns = c(cor.columns, 'p (Krusk-W.)', 'FDR (Krusk-W.)')
      } else {
        cor.columns = c(cor.columns, 'p (Krusk-W.)')
      }
    }
  
    if(Pars$Include.Spearman.corr)    cor.columns = c(cor.columns, 'Spearman r')
    if(Pars$Include.Pearson.corr)    cor.columns = c(cor.columns, 'Pearson r')
    if(Pars$Include.Corr.PValues){
      if(Pars$Include.Spearman.corr)
        cor.columns = c(cor.columns, 'p (Spearman)')
      if(Pars$Include.Pearson.corr)
        cor.columns = c(cor.columns, 'p (Pearson)')
    }

    Cor.table = data.matrix(array(dim=c(dim(norm.CPMs)[1], length(cor.columns))))
    colnames(Cor.table) = cor.columns
    gene.names = rownames(Stats.table)
    
    if(!simple.DE.mode){
      j = 1
      total.genes = dim(norm.CPMs)[1]
      for (j in 1:total.genes){
        current.CPMs = norm.CPMs[rownames(Stats.table)[j],]
        
        if(perform.classic.paired.test){
          current.CPMs.group1 = current.CPMs[1:(length(current.CPMs)/2)]
          current.CPMs.group2 = current.CPMs[(length(current.CPMs)/2 + 1) : length(current.CPMs)]
          Cor.table[j, 'p (Wilcoxon paired)'] = suppressWarnings(wilcox.test(pmax(current.CPMs.group1, Pars$Mann.Wh..Kr.Wall..zero.level.CPM),
            pmax(current.CPMs.group2, Pars$Mann.Wh..Kr.Wall..zero.level.CPM), paired = TRUE))$p.value
        }
  
        if(Two.Pools.tests.available){
          current.CPMs.group1 = current.CPMs[Main.Component %in% Main.Component.levels[1]]
          current.CPMs.group2 = current.CPMs[Main.Component %in% Main.Component.levels[2]]
          Cor.table[j, 'p (Mann-Wh.)'] = suppressWarnings(wilcox.test(pmax(current.CPMs.group1, Pars$Mann.Wh..Kr.Wall..zero.level.CPM),
            pmax(current.CPMs.group2, Pars$Mann.Wh..Kr.Wall..zero.level.CPM)))$p.value
          tryCatch(expr = {
            Cor.table[j, 'p (t-test)'] = suppressWarnings(t.test(current.CPMs.group1, current.CPMs.group2))$p.value
          }, error = function (err){
            Cor.table[j, 'p (t-test)'] = 1
          })
  
        } else if(KruskWall.test.available){
          Cor.table[j, 'p (Krusk-W.)'] = suppressWarnings(kruskal.test(pmax(current.CPMs, Pars$Mann.Wh..Kr.Wall..zero.level.CPM), Main.Component))$p.value
        }
        
        if(Pars$Include.Spearman.corr){
          ct = suppressWarnings(cor.test(Main.Component, pmax(current.CPMs, Pars$Mann.Wh..Kr.Wall..zero.level.CPM),method = 'spearman',alternative = 'two.sided'))
          Cor.table[j, 'Spearman r'] = ct$estimate
          if(Pars$Include.Corr.PValues)
            Cor.table[j, 'p (Spearman)'] = ct$p.value
        }
        
        if(Pars$Include.Pearson.corr){
          ct = suppressWarnings(cor.test(Main.Component, current.CPMs,method = 'pearson',alternative = 'two.sided'))
          Cor.table[j ,'Pearson r'] = ct$estimate
          if(Pars$Include.Corr.PValues)
            Cor.table[j, 'p (Pearson)'] = ct$p.value
        }
  
        if (j %% 25 == 0 & verbose)
            cat(sprintf('\rCalculating correlations and add. stat. tests. Completed %.2f%%       ',j/total.genes*100))
      }
      cat(sprintf('\rCalculating correlations and add. stat. tests - completed                                \n'))

      if(Pars$Include.DeltaFreq.data & length(Main.Component.levels) == 2){
        # DeltaFreq.table = data.matrix(array(dim=c(dim(norm.CPMs)[1], 3)))

        # ResTable.gene.part..passed = subset(ResTable.gene.part, subset = keep)

        which.cols..0 = CPMs.groups[[1]] # + total.cols.count - CPMs.cols.count
        which.cols..1 = CPMs.groups[[2]] # + total.cols.count - CPMs.cols.count

        DeltaFreq.table = as.data.frame(t(apply(norm.CPMs[rownames(Stats.table), order(Main.Component)], 1, function(c.row) {
          CPMs..0 = as.numeric(as.character(c.row[which.cols..0]))
          CPMs..1 = as.numeric(as.character(c.row[which.cols..1]))

          mean_logCPM = mean(log2(c(CPMs..0, CPMs..1) + Pars$DeltaFreq..CPM.const.add))
          log.rel.CPMs..0 = log2(CPMs..0 + Pars$DeltaFreq..CPM.const.add) - mean_logCPM
          log.rel.CPMs..1 = log2(CPMs..1 + Pars$DeltaFreq..CPM.const.add) - mean_logCPM
          freq.plus..0 = sum.mod(log.rel.CPMs..0 >= Pars$DeltaFreq..delta.log.CPM.level) / length(log.rel.CPMs..0)
          freq.minus..0 = sum.mod(log.rel.CPMs..0 <= -Pars$DeltaFreq..delta.log.CPM.level) / length(log.rel.CPMs..0)
          freq.plus..1 = sum.mod(log.rel.CPMs..1 >= Pars$DeltaFreq..delta.log.CPM.level) / length(log.rel.CPMs..1)
          freq.minus..1 = sum.mod(log.rel.CPMs..1 <= -Pars$DeltaFreq..delta.log.CPM.level) / length(log.rel.CPMs..1)
          delta.Freq.plus = 100 * abs(freq.plus..1 - freq.plus..0)
          delta.Freq.minus = 100 * abs(freq.minus..1 - freq.minus..0)
          max.delta.Freq = max(delta.Freq.plus, delta.Freq.minus)
          c(delta.Freq.plus, delta.Freq.minus, max.delta.Freq)
          # return(delta.Freq.plus >= Pars$Randomization.test__deltaFreq.positive.threshold &
          #   delta.Freq.minus >= Pars$Randomization.test__deltaFreq.negative.threshold &
          #   max.delta.Freq >= Pars$Randomization.test__max.deltaFreq.threshold)
        } )))

        colnames(DeltaFreq.table) = c('d_Freq+', 'd_Freq-', 'd_Freq max')

      }

    }
  
  }

  if(('p (Mann-Wh.)' %in% cor.columns) & ('FDR (Mann-Wh.)' %in% cor.columns)){
    Cor.table[, 'FDR (Mann-Wh.)'] = p.adjust(Cor.table[, 'p (Mann-Wh.)'], method = 'BH')
  }
  if(('p (Krusk-W.)' %in% cor.columns) & ('FDR (Krusk-W.)' %in% cor.columns)){
    Cor.table[, 'FDR (Krusk-W.)'] = p.adjust(Cor.table[, 'p (Krusk-W.)'], method = 'BH')
  }
  if(('p (Wilcoxon paired)' %in% cor.columns) & ('FDR (Wilcoxon paired)' %in% cor.columns)){
    Cor.table[, 'FDR (Wilcoxon paired)'] = p.adjust(Cor.table[, 'p (Wilcoxon paired)'], method = 'BH')
  }
  if(('p (Mann-Wh., FC)' %in% cor.columns) & ('FDR (Mann-Wh., FC)' %in% cor.columns)){
    Cor.table[, 'FDR (Mann-Wh., FC)'] = p.adjust(Cor.table[, 'p (Mann-Wh., FC)'], method = 'BH')
  }
  if(('p (t-test, FC)' %in% cor.columns) & ('FDR (t-test, FC)' %in% cor.columns)){
    Cor.table[, 'FDR (t-test, FC)'] = p.adjust(Cor.table[, 'p (t-test, FC)'], method = 'BH')
  }
  
  Cor.table.spacer.insertion = as.data.frame(rep(NA, dim(Cor.table)[1]))
  colnames(Cor.table.spacer.insertion) = 'ct_spacer'
  if(ncol(DeltaFreq.table) > 0){
    DeltaFreq.spacer.insertion = as.data.frame(rep(NA, dim(DeltaFreq.table)[1]))
    colnames(DeltaFreq.spacer.insertion) = 'df_spacer'
    Cor.table = cbind(DeltaFreq.spacer.insertion, DeltaFreq.table, Cor.table.spacer.insertion, Cor.table)
  } else {
    Cor.table = cbind(Cor.table.spacer.insertion, Cor.table)
  }
  
  #length(norm.CPMs[rownames(tt)[j],])
  a = array(dim = c(length(all.preds.nr), ncol(Cor.table)))
  colnames(a) = colnames(Cor.table)
  Cor.table.ext = rbind(a, Cor.table)
  
  
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
    colnames(Mart.table) = c("Symbol","Biotype","description","RefSeq_Summary")
  }  else{
    colnames(Mart.table) = c("Symbol","Biotype","description")
  }
  
  if((Pars$Include.GO.annotation.in.Excel.results | Pars$Include.KEGG.annotation.in.Excel.results | Pars$Include.Reactome.annotation.in.Excel.results) & is.null(randomization.test.permut.N)){
    cat('Supplying genes with GO/KEGG/Reactome annotation. this may take a while... \n')
    if(Pars$Include.annotation..top.genes.limit > 0) {
      max.genes.count = min(Pars$Include.annotation..top.genes.limit, nrow(Mart.table))
    } else max.genes.count = nrow(Mart.table)
    DB.annotation.table = data.frame(array(dim = c(max.genes.count, 1)))
    DB.annotation.table[,1] = NA
    rownames(DB.annotation.table) = rownames(Mart.table)[1:max.genes.count]
    
    combine.two.rows = function(x){
      if(!is.na(x[1]) & !is.na(x[2])){
        paste(x, collapse = ', ')
      } else if(is.na(x[1]) & !is.na(x[2])){
        x[2]
      } else if(!is.na(x[1]) & is.na(x[2])){
        x[1]
      } else {
        NA
      }}
    
    if(Pars$Include.KEGG.annotation.in.Excel.results){
      tmp.anno = Get.KEGG.IDs.for.Ensembl.genes(genes = rownames(Mart.table), Startup.Data, forced.parameters = Pars, add.pathway.names = Pars$Include.Anno.entry.names.in.Excel.results)
      tmp.anno = sapply(tmp.anno, function(x) paste(x, collapse = ', '))
      tmp.anno = as.data.frame(tmp.anno)
      DB.annotation.table = apply(cbind(DB.annotation.table, tmp.anno), 1, combine.two.rows)
    }
    
    if(Pars$Include.GO.annotation.in.Excel.results){
      tmp.anno = Get.GO.IDs.for.Ensembl.genes(genes = rownames(Mart.table), Startup.Data, forced.parameters = Pars, add.term.names = Pars$Include.Anno.entry.names.in.Excel.results)
      tmp.anno = sapply(tmp.anno, function(x) paste(x, collapse = ', '))
      tmp.anno = as.data.frame(tmp.anno)
      DB.annotation.table = apply(cbind(DB.annotation.table, tmp.anno), 1, combine.two.rows)
    }

    if(Pars$Include.Reactome.annotation.in.Excel.results){
      tmp.anno = Get.Reactome.IDs.for.Ensembl.genes(genes = rownames(Mart.table), Startup.Data, forced.parameters = Pars, add.pathway.names = Pars$Include.Anno.entry.names.in.Excel.results)
      tmp.anno = sapply(tmp.anno, function(x) paste(x, collapse = ', '))
      tmp.anno = as.data.frame(tmp.anno)
      DB.annotation.table = apply(cbind(DB.annotation.table, tmp.anno), 1, combine.two.rows)
    }
    DB.annotation.table = sapply(DB.annotation.table, function(x) gsub(pattern = ', , ', replacement = ', ', x, fixed = TRUE))
    DB.annotation.table = sapply(DB.annotation.table, function(x) gsub(pattern = '^, ', replacement = '', x))
    DB.annotation.table = sapply(DB.annotation.table, function(x) gsub(pattern = ', $', replacement = '', x))
    DB.annotation.table %<>% as.data.frame
    colnames(DB.annotation.table) = 'Annotation'
    #write.table(x = DB.annotation.table, file = 'anno.table.txt', sep = '\t', quote = F)
    Mart.table = cbind(Mart.table, DB.annotation.table)
  }

  a = array(dim = c(length(all.preds.nr) ,dim(Mart.table)[2]))
  rownames(a) = all.preds.nr
  colnames(a) = colnames(Mart.table)
  Mart.table.ext = rbind(a,Mart.table)
  
  #head(Stats.table)
  ## Calculating scores
  LogFC.col.name = Get.LogFC.col.name(colnames(Stats.table))
  if(!perform.FC.associations.paired.test) PValue.col.name = Get.PValue.col.name(Pars$Main.test, colnames(Stats.table))

  if(perform.FC.associations.paired.test){
    if('trimmed delta LogFC' %in% colnames(Stats.table)){
      Scores = 100 * ((abs(Stats.table$'delta LogFC') * abs(Stats.table$'trimmed delta LogFC'))**0.5) * abs(Cor.table[,"Spearman r (FC)"])
    } else
      Scores = 100 * abs(Stats.table$'delta LogFC') * abs(Cor.table[,"Spearman r (FC)"])
    
  } else if(!simple.DE.mode){
    if(Pars$Score.calculating.method == 'standard'){
      if (Predictors.count > 1) LogFC..score..multiplier = Stats.table[,LogFC.col.name] else 
        LogFC..score..multiplier = 1
        
      if('LR' %in% colnames(Stats.table)){
        if(Pars$Main.test == 'LR'){
          if(!("p (LR test)" %in% colnames(Stats.table)))
            stop(sprintf('[1] Incorrect Stats.table colnames %s', toString(colnames(Stats.table))))
          Scores = (-1)*log2(Stats.table[,"p (LR test)"] + 1e-300)*(2*Stats.table[,"LR"]+10)*LogFC..score..multiplier
        } else if(Pars$Main.test == 'LR paired'){
          if(!("p (LR test paired)" %in% colnames(Stats.table)))
            stop(sprintf('[2] Incorrect Stats.table colnames %s', toString(colnames(Stats.table))))
          Scores = (-1)*log2(Stats.table[,"p (LR test paired)"] + 1e-300)*(2*Stats.table[,"LR"]+10)*LogFC..score..multiplier
        } else {
          stop(sprintf('[3] Incorrect Pars$Main.test %s', Pars$Main.test))
        }

      } else if('F' %in% colnames(Stats.table)){
        if(Pars$Main.test == 'QL'){
          if(!("p (QLF test)" %in% colnames(Stats.table)))
            stop(sprintf('[4] Incorrect Stats.table colnames %s', toString(colnames(Stats.table))))
          Scores = (-1)*log2(Stats.table[,"p (QLF test)"] + 1e-300)*(2*Stats.table[,"F"]+10)*LogFC..score..multiplier
        } else if(Pars$Main.test == 'QL paired'){
          if(!("p (QLF test paired)" %in% colnames(Stats.table)))
            stop(sprintf('[5] Incorrect Stats.table colnames %s', toString(colnames(Stats.table))))
          Scores = (-1)*log2(Stats.table[,"p (QLF test paired)"] + 1e-300)*(2*Stats.table[,"F"]+10)*LogFC..score..multiplier
        } else {
          stop(sprintf('[6] Incorrect Pars$Main.test %s', Pars$Main.test))
        }
        
      } else if('p (ex. test)' %in% colnames(Stats.table)) {
        Scores = (-1)*log2(Stats.table[,"p (ex. test)"] + 1e-300)*50*LogFC..score..multiplier
      } else if('PValue' %in% colnames(Stats.table)) {
        Scores = (-1)*log2(Stats.table[,"PValue"] + 1e-300)*50*LogFC..score..multiplier
      } else if('p' %in% colnames(Stats.table)) {
        Scores = (-1)*log2(Stats.table[,"p"] + 1e-300)*50*LogFC..score..multiplier
      } else {
        stop('[Analyze GLM] Incorrect Stats.table (1)')
      }

      if(Pars$Main.test != 'QL paired' & Pars$Main.test == 'LR paired' & Predictors.count == 1 & "Spearman r" %in% colnames(Cor.table))
        Scores = Scores * ((0.2+abs(Cor.table[,"Spearman r"]))^0.8)
      #Scores = (-1)*log2(Stats.table[,"PValue"])*Stats.table[,"logFC"]
      Scores = 10*abs(Scores)**0.5

    } else if(Pars$Score.calculating.method == 'complex'){
      LogFC.col.name = Get.LogFC.col.name(colnames(Stats.table))
      PValue.col.name = Get.PValue.col.name(Pars$Main.test, colnames(Stats.table))
      FDR.col.name = Get.FDR.col.name(Pars$Main.test, colnames(Stats.table), look.for.paired.test.res = perform.classic.paired.test)
      NP.PValue.col.name = Get.NP.PValue.col.name(colnames(Cor.table))
      LogCPM.col.name = Get.LogCPM.col.name(colnames(Stats.table))
      #Scores = (-1)*log2(Stats.table[,"p (QLF test)"] + 1e-300)*(2*Stats.table[,"F"]+10)*LogFC..score..multiplier
      
      PValue.scores = pmax((-1)*log10(Stats.table[, PValue.col.name] + 1e-200) + log10(0.07), 0) + 0.2
      if(!is.na(NP.PValue.col.name)){
        NP.PValue.scores = pmax((-1)*log10(Cor.table[, NP.PValue.col.name] + 1e-100) + log10(0.07), 0) + 0.3
      } else {
        NP.PValue.scores = 1
      }
      LogFC.scores = abs(Stats.table[, LogFC.col.name])
      LogCPM_scoring_start = 1.8
      LogCPM_scoring_end = 4.4
      LogCPM_scoring_coord = (pmin(pmax(Stats.table[, LogCPM.col.name], LogCPM_scoring_start), LogCPM_scoring_end) - LogCPM_scoring_start) / (LogCPM_scoring_end - LogCPM_scoring_start)
      LogCPM.scores = 0.15 + LogCPM_scoring_coord

      if('d_Freq max' %in% colnames(Cor.table)){
        DeltaFreq.scores = Cor.table[,'d_Freq max'] / 100
        Scores = (PValue.scores^Pars$PValue.score.weigth * 
          (NP.PValue.scores^Pars$NP.PValue.score.weigth) * 
          (DeltaFreq.scores^Pars$DeltaFreq.score.weigth) * 
          (LogFC.scores^Pars$LogFC.score.weigth) * 
          (LogCPM.scores^Pars$LogCPM.score.weigth)) ^ (1/
          (Pars$PValue.score.weigth + Pars$NP.PValue.score.weigth + Pars$DeltaFreq.score.weigth + Pars$LogFC.score.weigth + Pars$LogCPM.score.weigth))
      } else {
        Scores = (PValue.scores^Pars$PValue.score.weigth * 
          (NP.PValue.scores^Pars$NP.PValue.score.weigth) * 
          (LogFC.scores^Pars$LogFC.score.weigth) * 
          (LogCPM.scores^Pars$LogCPM.score.weigth)) ^ (1/
          (Pars$PValue.score.weigth + Pars$NP.PValue.score.weigth + Pars$LogFC.score.weigth + Pars$LogCPM.score.weigth))
      }
      Scores = 100*Scores

    } else {
      stop(sprintf('Incorrect Score.calculating.method %s', Pars$Score.calculating.method))
    }
  }
  #med = median(Scores)
  #Scores = array(Scores,dim = c(length(Scores),1))
  #Scores = apply(Scores,1,function(x) max(0,x-med*5))
  Scores = array(Scores, dim = c(length(Scores),1))
  a = array(dim = c(length(all.preds.nr) ,dim(Scores)[2]))
  Scores.table.ext = rbind(a,Scores)
  colnames(Scores.table.ext) = c('Score')
  
  
  #sum(grepl("^NA",rownames(ResTable.gene.part)))
  
  CPMs.spacer.insertion = as.data.frame(rep(NA, dim(Counts.Table)[1]))

  if(Pars$report.norm.read.counts.instead.of.CPM){
    CPMs.spacer.insertion.colname = 'norm. counts'
  } else  CPMs.spacer.insertion.colname = 'CPMs:'

  colnames(CPMs.spacer.insertion) = CPMs.spacer.insertion.colname

  if(perform.classic.paired.test | perform.FC.associations.paired.test){
    LogFCs.spacer.insertion = as.data.frame(rep(NA, dim(LogFCs.Table)[1]))
    colnames(LogFCs.spacer.insertion) = 'LogFCs:'
  } else {
    LogFCs.spacer.insertion = data.frame(array(dim = c(nrow(Counts.Table), 0)))
    rownames(LogFCs.spacer.insertion) = rownames(Counts.Table)
  }



  #colnames(ResTable) = gsub(pattern = "Pearson.",replacement = "Pearson ",x = colnames(ResTable),fixed = T)
  #colnames(ResTable) = gsub(pattern = "Spearman.",replacement = "Spearman ",x = colnames(ResTable),fixed = T)
  #colnames(ResTable) = gsub(pattern = "Gene.name",replacement = "Gene name",x = colnames(ResTable),fixed = T)
  if(perform.classic.paired.test | perform.FC.associations.paired.test){
    ResTable = cbind(Mart.table.ext, Scores.table.ext, Stats.table.ext, Cor.table.ext, LogFCs.spacer.insertion, LogFCs.Table, CPMs.spacer.insertion, Counts.Table)
  } else {
    ResTable = cbind(Mart.table.ext, Scores.table.ext, Stats.table.ext, Cor.table.ext, CPMs.spacer.insertion, Counts.Table[, order(Main.Component)])
  }
  
  colnames(ResTable)[which(colnames(ResTable) %in% 'logFC')] = 'LogFC'
  colnames(ResTable)[which(colnames(ResTable) %in% 'logCPM')] = 'LogCPM'
  colnames(ResTable)[which(colnames(ResTable) %in% 'Gene name')] = 'Symbol'
  colnames(ResTable)[which(colnames(ResTable) %in% 'description')] = 'Name'

  ResTable.gene.part = ResTable[(length(all.preds.nr)+1):dim(ResTable)[1],]
  ResTable.gene.part = ResTable.gene.part[order(-ResTable.gene.part[,"Score"]), ]

  if(Pars$Genes.sorting.criterium != 'default'){
    if(Pars$Genes.sorting.criterium == 'Spearman'){ 
      ResTable.gene.part = ResTable.gene.part[order(ResTable.gene.part[,"p (Spearman)"]), ]
    } else if(Pars$Genes.sorting.criterium == 'Pearson'){
      ResTable.gene.part = ResTable.gene.part[order(ResTable.gene.part[,"p (Pearson)"]), ]
      # ResTable.gene.part = ResTable.gene.part[order(-abs(ResTable.gene.part[,"p (Pearson)"])), ]
    } else {
      sorting.colname = NA
      if(Pars$Genes.sorting.criterium == 'NP'){
        if('p (Mann-Wh., FC)' %in% colnames(ResTable.gene.part)){
          sorting.colname = 'p (Mann-Wh., FC)'
        } else if('p (Mann-Wh.)' %in% colnames(ResTable.gene.part)){
          sorting.colname = 'p (Mann-Wh.)'
        } else if('p (Krusk-W.)' %in% colnames(ResTable.gene.part)){
          sorting.colname = 'p (Krusk-W.)'
        } else {
          print('Cannot find any non-parametrical stat test')
        }
      } else {
        if(sprintf('p (%s)', Pars$Genes.sorting.criterium) %in% colnames(ResTable.gene.part)){
          sorting.colname = sprintf('p (%s)', Pars$Genes.sorting.criterium)
        }
      }

      if(!is.na(sorting.colname)){
        ResTable.gene.part = ResTable.gene.part[order(ResTable.gene.part[, sorting.colname]), ]
      } else {
        msg = sprintf('Cannot find column name for sorting mode "%s"\n', Pars$Genes.sorting.criterium)
        cat(msg)
        warning(msg)
      }
    }
  }

  CPMs.cols.count = ncol(Counts.Table)
  LogFCs.cols.count = ncol(LogFCs.Table)
  total.cols.count = ncol(ResTable.gene.part)
  source.columns.seq = colnames(ResTable.gene.part)
  current.columns.seq = colnames(ResTable.gene.part)

  CPMs.col.start = total.cols.count - CPMs.cols.count + 1
  LogFCs.col.start = total.cols.count - CPMs.cols.count - 1 - ncol(LogFCs.Table) + 1

  if(Pars$Perform.Randomization.test){
    # fileConn = file("genes.filtering.stats.txt")

    keep = ResTable.gene.part[, 'Score'] >= Pars$Randomization.test__Score.threshold
    msg = sprintf('%d genes passed Score threshold\n', sum.mod(keep))
    cat(msg)
    write(msg, file = 'genes.filtering.stats.txt', append=FALSE)

    PValue.col.name = Get.PValue.col.name(Pars$Main.test, colnames(ResTable.gene.part))
    keep = keep & (ResTable.gene.part[, PValue.col.name] <= Pars$Randomization.test__PValue.threshold)
    msg = sprintf('%d genes passed p-value (%s) threshold\n', sum.mod(keep), Pars$Main.test)
    cat(msg)
    write(msg, file = 'genes.filtering.stats.txt', append = TRUE)

    FDR.col.name = Get.FDR.col.name(Pars$Main.test, colnames(ResTable.gene.part))
    keep = keep & (ResTable.gene.part[, FDR.col.name] <= Pars$Randomization.test__FDR.threshold)
    msg = sprintf('%d genes passed FDR (%s) threshold\n', sum.mod(keep), Pars$Main.test)
    cat(msg)
    write(msg, file = 'genes.filtering.stats.txt', append = TRUE)

    keep = keep & (abs(ResTable.gene.part[, LogFC.col.name]) >= Pars$Randomization.test__abs.LogFC.threshold)
    if(perform.FC.associations.paired.test){
      msg = sprintf('%d genes passed abs. delta LogFC threshold\n', sum.mod(keep))
    } else {
      msg = sprintf('%d genes passed abs. LogFC threshold\n', sum.mod(keep))
    }
    cat(msg)
    write(msg, file = 'genes.filtering.stats.txt', append = TRUE)

    keep = keep & (ResTable.gene.part[, 'LogCPM'] >= Pars$Randomization.test__LogCPM.threshold)
    msg = sprintf('%d genes passed LogCPM threshold\n', sum.mod(keep))
    cat(msg)
    write(msg, file = 'genes.filtering.stats.txt', append = TRUE)

    if(Pars$Randomization.test__Mann.Whitney.PValue.threshold < 1){
      test.passed = TRUE
      if('p (Mann-Wh.)' %in% cor.columns){
        keep = keep & (ResTable.gene.part[, 'p (Mann-Wh.)'] <= Pars$Randomization.test__Mann.Whitney.PValue.threshold)
      } else if('p (Krusk-W.)' %in% cor.columns){
        keep = keep & (ResTable.gene.part[, 'p (Krusk-W.)'] <= Pars$Randomization.test__Mann.Whitney.PValue.threshold)
      } else if('p (Mann-Wh., FC)' %in% cor.columns){
        keep = keep & (ResTable.gene.part[, 'p (Mann-Wh., FC)'] <= Pars$Randomization.test__Mann.Whitney.PValue.threshold)
      } else {
        cat(sprintf('\nRandomization.test__Mann.Whitney.PValue.threshold is set as %g. However, there is no Mann-Whitney/Kruskal-Wallis test performed\n', Pars$Randomization.test__Mann.Whitney.PValue.threshold))
        test.passed = FALSE
      }
      if(test.passed){
        msg = sprintf('%d genes passed Mann-Whitney/Kruskal-Wallis p-value threshold\n', sum.mod(keep))
        cat(msg)
        write(msg, file = 'genes.filtering.stats.txt', append = TRUE)
      }
    }

    if(Pars$Randomization.test__Mann.Whitney.FDR.threshold < 1){
      test.passed = TRUE
      if('FDR (Mann-Wh.)' %in% cor.columns){
        keep = keep & (ResTable.gene.part[, 'FDR (Mann-Wh.)'] <= Pars$Randomization.test__Mann.Whitney.FDR.threshold)
      } else if('FDR (Krusk-W.)' %in% cor.columns){
        keep = keep & (ResTable.gene.part[, 'FDR (Krusk-W.)'] <= Pars$Randomization.test__Mann.Whitney.FDR.threshold)
      } else if('FDR (Mann-Wh., FC)' %in% cor.columns){
        keep = keep & (ResTable.gene.part[, 'FDR (Mann-Wh., FC)'] <= Pars$Randomization.test__Mann.Whitney.FDR.threshold)
      } else {
        cat(sprintf('\nRandomization.test__Mann.Whitney.FDR.threshold is set as %g. However, there is no Mann-Whitney/Kruskal-Wallis test performed\n', Pars$Randomization.test__Mann.Whitney.FDR.threshold))
        test.passed = FALSE
      }
      if(test.passed){
        msg = sprintf('%d genes passed Mann-Whitney/Kruskal-Wallis p-value threshold\n', sum.mod(keep))
        cat(msg)
        write(msg, file = 'genes.filtering.stats.txt', append = TRUE)
      }
    }
    
    if(!perform.FC.associations.paired.test & sum.mod(keep) > 0 & (Pars$Randomization.test__deltaFreq.positive.threshold > 0 | Pars$Randomization.test__deltaFreq.negative.threshold > 0 | Pars$Randomization.test__max.deltaFreq.threshold > 0)){
      if(Pars$Randomization.test__deltaFreq.positive.threshold > 1)  Pars$Randomization.test__deltaFreq.positive.threshold = Pars$Randomization.test__deltaFreq.positive.threshold/100
      if(Pars$Randomization.test__deltaFreq.negative.threshold > 1)  Pars$Randomization.test__deltaFreq.negative.threshold = Pars$Randomization.test__deltaFreq.negative.threshold/100
      if(Pars$Randomization.test__max.deltaFreq.threshold > 1)  Pars$Randomization.test__max.deltaFreq.threshold = Pars$Randomization.test__max.deltaFreq.threshold/100
      if(length(Main.Component.levels) != 2){
        cat(sprintf('\nRandomization.test__deltaFreq.positive.threshold or Randomization.test__deltaFreq.negative.threshold or Randomization.test__max.deltaFreq.threshold are turned on. However, the number of groups (Main.Component) is %d (should be 2)\n',
          length(Main.Component.levels)))
        passed.genes = ResTable.gene.part[keep, 'Symbol']
        passed.genes..IDs = rownames(ResTable.gene.part)[keep]
      } else {
        ResTable.gene.part..passed = subset(ResTable.gene.part, subset = keep)

        which.cols..0 = CPMs.groups[[1]] + total.cols.count - CPMs.cols.count
        which.cols..1 = CPMs.groups[[2]] + total.cols.count - CPMs.cols.count

        c.keep = apply(ResTable.gene.part..passed, 1, function(c.row) {
          CPMs..0 = as.numeric(as.character(c.row[which.cols..0]))
          CPMs..1 = as.numeric(as.character(c.row[which.cols..1]))
          # print(CPMs..0)
          # print(CPMs..1)
          # print('-----------------------')
          mean_logCPM = mean(log2(c(CPMs..0, CPMs..1) + 1))
          log.rel.CPMs..0 = log2(CPMs..0 + 1) - mean_logCPM
          log.rel.CPMs..1 = log2(CPMs..1 + 1) - mean_logCPM
          freq.plus..0 = sum.mod(log.rel.CPMs..0 >= Pars$DeltaFreq..delta.log.CPM.level) / length(log.rel.CPMs..0)
          freq.minus..0 = sum.mod(log.rel.CPMs..0 <= -Pars$DeltaFreq..delta.log.CPM.level) / length(log.rel.CPMs..0)
          freq.plus..1 = sum.mod(log.rel.CPMs..1 >= Pars$DeltaFreq..delta.log.CPM.level) / length(log.rel.CPMs..1)
          freq.minus..1 = sum.mod(log.rel.CPMs..1 <= -Pars$DeltaFreq..delta.log.CPM.level) / length(log.rel.CPMs..1)
          delta.Freq.plus = abs(freq.plus..1 - freq.plus..0)
          delta.Freq.minus = abs(freq.minus..1 - freq.minus..0)
          max.delta.Freq = max(delta.Freq.plus, delta.Freq.minus)
          return(delta.Freq.plus >= Pars$Randomization.test__deltaFreq.positive.threshold &
            delta.Freq.minus >= Pars$Randomization.test__deltaFreq.negative.threshold &
            max.delta.Freq >= Pars$Randomization.test__max.deltaFreq.threshold)
        } )

        msg = sprintf('%d genes passed deltaFreq thresholds\n', sum.mod(c.keep))
        cat(msg)
        write(msg, file = 'genes.filtering.stats.txt', append = TRUE)
        passed.genes = ResTable.gene.part..passed[c.keep, 'Symbol']
        passed.genes..IDs = rownames(ResTable.gene.part..passed)[c.keep]
      }
    } else {
      passed.genes = ResTable.gene.part[keep, 'Symbol']
      passed.genes..IDs = rownames(ResTable.gene.part)[keep]
    }

    if(Startup.Data$miR.mode){  msg = sprintf('Total %d microRNAs passed all thresholds: %s\n', length(passed.genes..IDs), toString(passed.genes..IDs))
    } else msg = sprintf('Total %d genes passed all thresholds: %s\n', length(passed.genes), toString(passed.genes))
    cat(msg)
    write(msg, file = 'genes.filtering.stats.txt', append = TRUE)
    # close(fileConn)

    saveRDS(sum.mod(keep), file = 'passed.genes.count.rds')

    if(!is.null(randomization.test.permut.N))  return(invisible(NULL))
  }

  ResTable.sorted.by.Score = rbind(ResTable[1:length(all.preds.nr),], ResTable.gene.part)
  
  sorted.norm.CPMs = norm.CPMs[rownames(ResTable.gene.part),]
  write.table.mod(sorted.norm.CPMs,'sorted.norm.CPMs.tsv',sep='\t')
  if(Pars$report.norm.read.counts.instead.of.CPM){
    sorted.norm.counts = norm.counts[rownames(ResTable.gene.part),]
    write.table.mod(sorted.norm.counts,'sorted.norm.counts.tsv',sep='\t')
  }
  ResTable.gene.part = ResTable.gene.part[!is.na(ResTable.gene.part[,"Score"]),]
  
  output.tsv.file.name = Verify.path(sprintf("%s_%s_combined.txt",Analysis.Data$Analysis.name,Pars$DE.package))
  write.table.mod(ResTable.sorted.by.Score, file = output.tsv.file.name, sep='\t',na = '', col.names = colnames(ResTable.sorted.by.Score))

  Min.group.size.to..cluster.by.CPMs..filter.passed = FALSE
  for(l in names(CPMs.groups)){
    if(length(CPMs.groups[[l]]) >= Pars$Min.group.size.to..cluster)  Min.group.size.to..cluster.by.CPMs..filter.passed = TRUE
  }

  Min.group.size.to..cluster.by.LogFCs..filter.passed = FALSE
  for(l in names(LogFCs.groups)){
    if(length(LogFCs.groups[[l]]) >= Pars$Min.group.size.to..cluster)  Min.group.size.to..cluster.by.LogFCs..filter.passed = TRUE
  }

  clustering.samples.passed = FALSE
  if(Pars$Cluster.samples.in.Excel.results & Min.group.size.to..cluster.by.CPMs..filter.passed & !(perform.classic.paired.test | perform.FC.associations.paired.test)){
    cat('\nClustering and reordering samples by CPMs profiles...\n')
    # prep.Counts.Table = Counts.Table[,order(Main.Component)]

    l = names(CPMs.groups)[1]
    for(l in names(CPMs.groups)){
      if(length(CPMs.groups[[l]]) < 3)  next
      which.cols = CPMs.col.start + (CPMs.groups[[l]] - 1)
      genes.CPMs = as.data.frame(ResTable.gene.part[1:min(dim(ResTable.gene.part)[1], Pars$use.N.top.Genes.to.Cluster.samples), which.cols, drop = FALSE], stringsAsFactors = FALSE)

      if(type(genes.CPMs[1,]) == "character"){
        rel.LogCPM.table = apply(genes.CPMs, 1, function(all_values) {
          all_values = as.numeric(as.character(all_values))
          gmean_current_CPM = exp(mean(log(all_values + 2)))
          return(log2((all_values + 2) / gmean_current_CPM))
        } )
      } else {
        rel.LogCPM.table = apply(genes.CPMs, 1, function(all_values) {
          gmean_current_CPM = exp(mean(log(all_values + 2)))
          return(log2((all_values + 2) / gmean_current_CPM))
        } )
      }


      dend.order = tryCatch(expr = {
        dend=hclust(dist(rel.LogCPM.table, 'euclidean'), method = Pars$Clustering.method)
        dend$order
      }, error = function (err){
        cat(sprintf('Cannot reorder samples (by CPMs) in the group %s (total %d samples; %d genes)\n', l, dim(genes.CPMs)[2], dim(genes.CPMs)[1]))
        1:dim(genes.CPMs)[2]
      })

      current.columns.seq[which.cols] = current.columns.seq[which.cols][dend.order]
      ResTable[, which.cols] = ResTable[, which.cols][, dend.order]
      ResTable.sorted.by.Score[, which.cols] = ResTable.sorted.by.Score[, which.cols][, dend.order]
      ResTable.gene.part[, which.cols] = ResTable.gene.part[, which.cols][, dend.order]
      clustering.samples.passed = TRUE
    }
    colnames(ResTable) = current.columns.seq
    colnames(ResTable.sorted.by.Score) = current.columns.seq
    colnames(ResTable.gene.part) = current.columns.seq

    if(clustering.samples.passed){
      output.tsv.file.name.reordered = Verify.path(sprintf("%s_%s_combined_reordered.txt",Analysis.Data$Analysis.name, Pars$DE.package))
      write.table.mod(ResTable.sorted.by.Score, file = output.tsv.file.name.reordered, sep='\t',na = '', col.names = colnames(ResTable.sorted.by.Score))
    }
  
  } else if(Pars$Cluster.samples.in.Excel.results & Min.group.size.to..cluster.by.LogFCs..filter.passed & (perform.classic.paired.test | perform.FC.associations.paired.test)){
    cat('\nClustering and reordering samples by LogFCs profiles...\n')
    # prep.Counts.Table = Counts.Table[,order(Main.Component)]

    l = names(LogFCs.groups)[1]
    for(l in names(LogFCs.groups)){
      if(length(LogFCs.groups[[l]]) < 3)  next
      which.cols = LogFCs.col.start + (LogFCs.groups[[l]] - 1)
      which.cols..N = CPMs.col.start + (LogFCs.groups[[l]] - 1)
      which.cols..T = CPMs.col.start + CPMs.cols.count/2 + (LogFCs.groups[[l]] - 1)

      genes.LogFCs = subset(ResTable.gene.part, select = which.cols)[1:min(dim(ResTable.gene.part)[1], Pars$use.N.top.Genes.to.Cluster.samples), ]
      # rel.LogCPM.table = apply(genes.LogFCs, 1, function(all_values) {
      #   gmean_current_CPM = exp(mean(log(all_values + 2)))
      #   return(log2((all_values + 2) / gmean_current_CPM))
      # } )

      dend.order = tryCatch(expr = {
        dend=hclust(dist(t(genes.LogFCs), 'euclidean'), method = Pars$Clustering.method)
        dend$order
      }, error = function (err){
        cat(sprintf('Cannot reorder samples (by LogFCs) in the group %s (total %d samples; %d genes)\n', l, dim(genes.LogFCs)[2], dim(genes.LogFCs)[1]))
        1:dim(genes.LogFCs)[2]
      })

      current.columns.seq[which.cols] = current.columns.seq[which.cols][dend.order]
      ResTable[, which.cols] = ResTable[, which.cols][, dend.order]
      ResTable.sorted.by.Score[, which.cols] = ResTable.sorted.by.Score[, which.cols][, dend.order]
      ResTable.gene.part[, which.cols] = ResTable.gene.part[, which.cols][, dend.order]

      current.columns.seq[which.cols..N] = current.columns.seq[which.cols..N][dend.order]
      ResTable[, which.cols..N] = ResTable[, which.cols..N][, dend.order]
      ResTable.sorted.by.Score[, which.cols..N] = ResTable.sorted.by.Score[, which.cols..N][, dend.order]
      ResTable.gene.part[, which.cols..N] = ResTable.gene.part[, which.cols..N][, dend.order]

      current.columns.seq[which.cols..T] = current.columns.seq[which.cols..T][dend.order]
      ResTable[, which.cols..T] = ResTable[, which.cols..T][, dend.order]
      ResTable.sorted.by.Score[, which.cols..T] = ResTable.sorted.by.Score[, which.cols..T][, dend.order]
      ResTable.gene.part[, which.cols..T] = ResTable.gene.part[, which.cols..T][, dend.order]

      clustering.samples.passed = TRUE
    }
    colnames(ResTable) = current.columns.seq
    colnames(ResTable.sorted.by.Score) = current.columns.seq
    colnames(ResTable.gene.part) = current.columns.seq

    if(clustering.samples.passed){
      output.tsv.file.name.reordered = Verify.path(sprintf("%s_%s_combined_reordered.txt",Analysis.Data$Analysis.name, Pars$DE.package))
      write.table.mod(ResTable.sorted.by.Score, file = output.tsv.file.name.reordered, sep='\t',na = '', col.names = colnames(ResTable.sorted.by.Score))
    }
  }


  source.samples.seq = source.columns.seq[(total.cols.count - CPMs.cols.count + 1) : total.cols.count]
  current.samples.seq = current.columns.seq[(total.cols.count - CPMs.cols.count + 1) : total.cols.count]

  if(is.null(Create.Excel.results))  Create.Excel.results = Pars$Create.Excel.results
  
  if (Create.Excel.results){
    cat('Writing output files....\n')

    forced_col_types_by_names..CL.insertion = '--forced-col-types-by-names "LogFC:&LogFC&,&trimmed LogFC&,&TBA LogFC&" "logfc_array2:&gr.1 mean LogFC&,&gr.2 mean LogFC&" "Delta_LogFC:&delta LogFC&,&trimmed delta LogFC&" "LogCPM:&LogCPM&" "p:&PValue&,&TBA p-value&,&p (QLF test)&,&p (LR test)&,&p (QLF test paired)&,&p (LR test paired)&,&p (ex. test)&,&p (Mann-Wh.)&,&p (Wilcoxon paired)&,&p (t-test)&,&p (Mann-Wh., FC)&,&p (t-test, FC)&,&p (Spearman)&,&p (Pearson)&,&p (Spearman, FC)&,&p (Pearson, FC)&" "FDR:&TBA FDR&,&FDR (QLF test)&,&FDR (LR test)&,&FDR (QLF test paired)&,&FDR (LR test paired)&,&FDR (ex. test)&,&FDR (Wilcoxon paired)&,&FDR (Mann-Wh.)&,&FDR (Mann-Wh., FC)&,&FDR (t-test, FC)&,&FDR&" "score:Score" "corr:&Spearman r&,&Pearson r&,&Spearman r (FC)&,&Pearson r (FC)&" "Biotype:Biotype"  "Gene_Symbol:Symbol" "Gene_Name:Name" "Spacer:&ct_spacer&,&df_spacer&,&LogFCs&" "rank_1:&CG content rank (%)&,&Tr.Length rank (%)&,&LogCPM rank (%)&" "warm_gradient:&d_Freq+&" "cold_gradient:&d_Freq-&" "gray_gradient:&d_Freq max&"'
    if(Startup.Data$miR.mode)
      forced_col_types_by_names..CL.insertion = sprintf('%s "Hidden:&Symbol&,&Biotype&,&Name&,&RefSeq_Summary&,&Annotation&"', forced_col_types_by_names..CL.insertion)
    if(perform.FC.associations.paired.test)
      forced_col_types_by_names..CL.insertion = sprintf("%s ", forced_col_types_by_names..CL.insertion)

    if(perform.classic.paired.test | perform.FC.associations.paired.test){
      LogFCs.cols.numbers.text = paste(sprintf('%d', LogFCs.col.start:(CPMs.col.start-1)), collapse = ',')
      paired.test...CL.insertion = sprintf('"logfc_array:%s"', LogFCs.cols.numbers.text)
    } else paired.test...CL.insertion = ''

    CPMs.cols.numbers.text = paste(sprintf('%d', CPMs.col.start:ncol(ResTable)), collapse = ',')
    if(Pars$Create.CPM.profiles.in.Excel.results){
      if(include.GLM.model.name.in.Result.names){
        out.Excel.file = Verify.path(sprintf('%s, DE results, CPM profiles.xlsx',Analysis.Data$Analysis.name))
      } else if(include.GLM.model.name.in.Result.names..as.suffix){
        out.Excel.file = Verify.path(sprintf('DE results, CPM profiles - %s.xlsx',Analysis.Data$Analysis.name))
      } else {
        out.Excel.file = Verify.path(sprintf('DE results, CPM profiles.xlsx'))
      }
      
      CL = sprintf('%s "%s/Any.2.Excel.py" --in "%s" --out-excel "%s" --one-book yes --sheet-names "%s" --cpm-cells-formatting--logfc-max 2.5 --cpm-cells-formatting--cpm-max 500 --cpm-cells-formatting--gradient-color-schema 3 --cpm-cells-formatting--const-add 0.5 --predictor-rows-count %d %s --entry-type "Gene ID" --forced-cols "cpm_array:%s" "cpm:%s" %s --disable-coltype-autoassign yes --first-col-width 19 --first-col-italic yes --cpm-sparklines-groups %s --cpm-sparklines-mode %s --cpm-sparklines-start-col %d --logfc-sparklines-groups %s --logfc-sparklines-start-col %d %s --max-strings %s',
                   Pars$python.bin, Pars$suppl.data.dir, output.tsv.file.name, out.Excel.file, GLM.model, length(all.preds.nr),
                   forced_col_types_by_names..CL.insertion,
                   CPMs.cols.numbers.text, CPMs.cols.numbers.text,
                   paired.test...CL.insertion,
                   sapply(names(CPMs.groups), function(x) sprintf('"%s:%s"', Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(CPMs.col.start - 1 + CPMs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                   'relative', 7,
                   sapply(names(LogFCs.groups), function(x) sprintf('"%s:%s"', Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(LogFCs.col.start - 1 + LogFCs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                   7,
                   ifelse('TBA LogFC' %in% colnames(ResTable.sorted.by.Score), '--forced-col-widths-by-names 1:"&trimmed LogFC&,&Trimmed LogFC&"', ''),
                   as.character(Pars$Max.genes.in.Excel.results))
      
      if(Pars$bidirectional.bars.in.LogFC.cells.in.Excel.reports){
        CL = sprintf("%s --bidirectional-bars-in-logfc-cells yes --simple-logfc-format no", CL)
      } else if (Pars$colored.LogFC.cells.in.Excel.reports | perform.FC.associations.paired.test){
        CL = sprintf("%s --simple-logfc-format yes --bidirectional-bars-in-logfc-cells no", CL)
      } else {
        CL = sprintf("%s --simple-logfc-format no --bidirectional-bars-in-logfc-cells no", CL)
      }

      if(Pars$white.background.in.Excel.reports){  CL = sprintf("%s --whitespace-mode yes", CL)
      } else                                       CL = sprintf("%s --whitespace-mode no",  CL)

      if(Pars$freeze.header.in.Excel.reports)  CL = sprintf("%s --freeze-rows 1", CL)
      if(Pars$freeze.gene.names.in.Excel.reports)  CL = sprintf("%s --freeze-cols 2", CL)
      
      if(Pars$Include.small.heatmaps.in.Excel.reports){
        CL = sprintf("%s --include-cpm-heatmaps yes --include-logfc-heatmaps yes", CL)
      }

      code = system(CL)
      if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
      if (!file.exists(out.Excel.file)){
        message(sprintf('Excel worksheet generation (DE results) for ~%s FAILED', Analysis.Data$Analysis.name))
      }

      if(clustering.samples.passed){
        if(include.GLM.model.name.in.Result.names){
          out.Excel.file = Verify.path(sprintf('%s, DE results, CPM profiles [reordered].xlsx',Analysis.Data$Analysis.name))
        } else if(include.GLM.model.name.in.Result.names..as.suffix){
          out.Excel.file = Verify.path(sprintf('DE results, CPM profiles [reordered] - %s.xlsx',Analysis.Data$Analysis.name))
        } else {
          out.Excel.file = Verify.path(sprintf('DE results, CPM profiles [reordered].xlsx'))
        }
        
        CL = sprintf('%s "%s/Any.2.Excel.py" --in "%s" --out-excel "%s" --one-book yes --sheet-names "%s" --cpm-cells-formatting--logfc-max 2.5 --cpm-cells-formatting--cpm-max 500 --cpm-cells-formatting--gradient-color-schema 3 --cpm-cells-formatting--const-add 0.5 --predictor-rows-count %d %s --entry-type "Gene ID" --forced-cols "cpm_array:%s" "cpm:%s" %s --disable-coltype-autoassign yes --first-col-width 19 --first-col-italic yes --cpm-sparklines-groups %s --cpm-sparklines-mode %s --cpm-sparklines-start-col %d --logfc-sparklines-groups %s --logfc-sparklines-start-col %d %s --max-strings %s',
                     Pars$python.bin, Pars$suppl.data.dir, output.tsv.file.name.reordered, out.Excel.file, GLM.model, length(all.preds.nr),
                     forced_col_types_by_names..CL.insertion,
                     CPMs.cols.numbers.text, CPMs.cols.numbers.text,
                     paired.test...CL.insertion,
                     sapply(names(CPMs.groups), function(x) sprintf('"%s:%s"', Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(CPMs.col.start - 1 + CPMs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     # sapply(names(CPMs.groups), function(x) sprintf('"%s:%s"', Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(CPMs.groups[[x]] + dim(ResTable)[2] - dim(Counts.Table)[2], collapse = ','))) %>% paste(., collapse = ' '),
                     'relative', 7,
                     sapply(names(LogFCs.groups), function(x) sprintf('"%s:%s"', Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(LogFCs.col.start - 1 + LogFCs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     7,
                     ifelse('TBA LogFC' %in% colnames(ResTable.sorted.by.Score), '--forced-col-widths-by-names 1:"&trimmed LogFC&,&Trimmed LogFC&"', ''),
                     as.character(Pars$Max.genes.in.Excel.results))

      if(Pars$bidirectional.bars.in.LogFC.cells.in.Excel.reports){
        CL = sprintf("%s --bidirectional-bars-in-logfc-cells yes --simple-logfc-format no", CL)
      } else if (Pars$colored.LogFC.cells.in.Excel.reports | perform.FC.associations.paired.test){
        CL = sprintf("%s --simple-logfc-format yes --bidirectional-bars-in-logfc-cells no", CL)
      } else {
        CL = sprintf("%s --simple-logfc-format no --bidirectional-bars-in-logfc-cells no", CL)
      }

      if(Pars$white.background.in.Excel.reports){  CL = sprintf("%s --whitespace-mode yes", CL)
      } else                                       CL = sprintf("%s --whitespace-mode no",  CL)

      if(Pars$freeze.header.in.Excel.reports)  CL = sprintf("%s --freeze-rows 1", CL)
      if(Pars$freeze.gene.names.in.Excel.reports)  CL = sprintf("%s --freeze-cols 2", CL)
      
      if(Pars$Include.small.heatmaps.in.Excel.reports){
        CL = sprintf("%s --include-cpm-heatmaps yes --include-logfc-heatmaps yes", CL)
      }

      code = system(CL)
      if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
        if (!file.exists(out.Excel.file)){
          message(sprintf('Excel worksheet generation (DE results) for ~%s FAILED', Analysis.Data$Analysis.name))
          cat(sprintf('\nCommand line FAILED:\n%s\n', CL))
        }
      }
    }

    if(Pars$Create.rel.LogCPM.profiles.in.Excel.results){
      if(include.GLM.model.name.in.Result.names){
        out.Excel.file = Verify.path(sprintf('%s, DE results, Log.rel.CPM profiles.xlsx',Analysis.Data$Analysis.name))
      } else if(include.GLM.model.name.in.Result.names..as.suffix){
        out.Excel.file = Verify.path(sprintf('DE results, Log.rel.CPM profiles - %s.xlsx',Analysis.Data$Analysis.name))
      } else {
        out.Excel.file = Verify.path(sprintf('DE results, Log.rel.CPM profiles.xlsx'))
      }

      CL = sprintf('%s "%s/Any.2.Excel.py" --in "%s" --out-excel "%s" --one-book yes --sheet-names "%s" --cpm-cells-formatting--logfc-max 2.5 --cpm-cells-formatting--cpm-max 500 --cpm-cells-formatting--gradient-color-schema 3 --cpm-cells-formatting--const-add 0.5 --predictor-rows-count %d %s --entry-type "Gene ID" --forced-cols "cpm_array:%s" "cpm:%s" %s --disable-coltype-autoassign yes --first-col-width 19 --first-col-italic yes --cpm-sparklines-groups %s --cpm-sparklines-mode %s --cpm-sparklines-start-col %d --logfc-sparklines-groups %s --logfc-sparklines-start-col %d %s --max-strings %s',
                   Pars$python.bin, Pars$suppl.data.dir, output.tsv.file.name, out.Excel.file, GLM.model, length(all.preds.nr),
                   forced_col_types_by_names..CL.insertion,
                   CPMs.cols.numbers.text, CPMs.cols.numbers.text,
                   paired.test...CL.insertion,
                   sapply(names(CPMs.groups), function(x) sprintf('"%s:%s"', Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(CPMs.col.start - 1 + CPMs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                   'logfc', 7,
                   sapply(names(LogFCs.groups), function(x) sprintf('"%s:%s"', Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(LogFCs.col.start - 1 + LogFCs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                   7,
                   ifelse('TBA LogFC' %in% colnames(ResTable.sorted.by.Score), '--forced-col-widths-by-names 1:"&trimmed LogFC&,&Trimmed LogFC&"', ''),
                   as.character(Pars$Max.genes.in.Excel.results))

      if(Pars$bidirectional.bars.in.LogFC.cells.in.Excel.reports){
        CL = sprintf("%s --bidirectional-bars-in-logfc-cells yes --simple-logfc-format no", CL)
      } else if (Pars$colored.LogFC.cells.in.Excel.reports | perform.FC.associations.paired.test){
        CL = sprintf("%s --simple-logfc-format yes --bidirectional-bars-in-logfc-cells no", CL)
      } else {
        CL = sprintf("%s --simple-logfc-format no --bidirectional-bars-in-logfc-cells no", CL)
      }

      if(Pars$white.background.in.Excel.reports){  CL = sprintf("%s --whitespace-mode yes", CL)
      } else                                       CL = sprintf("%s --whitespace-mode no",  CL)

      if(Pars$freeze.header.in.Excel.reports)  CL = sprintf("%s --freeze-rows 1", CL)
      if(Pars$freeze.gene.names.in.Excel.reports)  CL = sprintf("%s --freeze-cols 2", CL)
      
      if(Pars$Include.small.heatmaps.in.Excel.reports){
        CL = sprintf("%s --include-cpm-heatmaps yes --include-logfc-heatmaps yes", CL)
      }

      code = system(CL)
      if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
      if (!file.exists(out.Excel.file)){
        message(sprintf('Excel worksheet generation (DE results) for ~%s FAILED', Analysis.Data$Analysis.name))
      }

      if(clustering.samples.passed){
        if(include.GLM.model.name.in.Result.names){
          out.Excel.file = Verify.path(sprintf('%s, DE results, Log.rel.CPM profiles [reordered].xlsx',Analysis.Data$Analysis.name))
        } else if(include.GLM.model.name.in.Result.names..as.suffix){
          out.Excel.file = Verify.path(sprintf('DE results, Log.rel.CPM profiles [reordered] - %s.xlsx',Analysis.Data$Analysis.name))
        } else {
          out.Excel.file = Verify.path(sprintf('DE results, Log.rel.CPM profiles [reordered].xlsx'))
        }

        CL = sprintf('%s "%s/Any.2.Excel.py" --in "%s" --out-excel "%s" --one-book yes --sheet-names "%s" --cpm-cells-formatting--logfc-max 2.5 --cpm-cells-formatting--cpm-max 500 --cpm-cells-formatting--gradient-color-schema 3 --cpm-cells-formatting--const-add 0.5 --predictor-rows-count %d %s --entry-type "Gene ID" --forced-cols "cpm_array:%s" "cpm:%s" %s --disable-coltype-autoassign yes --first-col-width 19 --first-col-italic yes --cpm-sparklines-groups %s --cpm-sparklines-mode %s --cpm-sparklines-start-col %d --logfc-sparklines-groups %s --logfc-sparklines-start-col %d %s --max-strings %s',
                     Pars$python.bin, Pars$suppl.data.dir, output.tsv.file.name.reordered, out.Excel.file, GLM.model, length(all.preds.nr),
                     forced_col_types_by_names..CL.insertion,
                     CPMs.cols.numbers.text, CPMs.cols.numbers.text,
                     paired.test...CL.insertion,
                     sapply(names(CPMs.groups), function(x) sprintf('"%s:%s"', Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(CPMs.col.start - 1 + CPMs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     'logfc', 7,
                     sapply(names(LogFCs.groups), function(x) sprintf('"%s:%s"', Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(LogFCs.col.start - 1 + LogFCs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     7,
                     ifelse('TBA LogFC' %in% colnames(ResTable.sorted.by.Score), '--forced-col-widths-by-names 1:"&trimmed LogFC&,&Trimmed LogFC&"', ''),
                     as.character(Pars$Max.genes.in.Excel.results))

      if(Pars$bidirectional.bars.in.LogFC.cells.in.Excel.reports){
        CL = sprintf("%s --bidirectional-bars-in-logfc-cells yes --simple-logfc-format no", CL)
      } else if (Pars$colored.LogFC.cells.in.Excel.reports | perform.FC.associations.paired.test){
        CL = sprintf("%s --simple-logfc-format yes --bidirectional-bars-in-logfc-cells no", CL)
      } else {
        CL = sprintf("%s --simple-logfc-format no --bidirectional-bars-in-logfc-cells no", CL)
      }

      if(Pars$white.background.in.Excel.reports){  CL = sprintf("%s --whitespace-mode yes", CL)
      } else                                       CL = sprintf("%s --whitespace-mode no",  CL)

      if(Pars$freeze.header.in.Excel.reports)  CL = sprintf("%s --freeze-rows 1", CL)
      if(Pars$freeze.gene.names.in.Excel.reports)  CL = sprintf("%s --freeze-cols 2", CL)
      
      if(Pars$Include.small.heatmaps.in.Excel.reports){
        CL = sprintf("%s --include-cpm-heatmaps yes --include-logfc-heatmaps yes", CL)
      }

      code = system(CL)
      if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
        if (!file.exists(out.Excel.file)){
          message(sprintf('Excel worksheet generation (DE results) for ~%s FAILED', Analysis.Data$Analysis.name))
        }
      }
    }    
  }
  
  
  #cor.test(counts['FBgn0004242',],counts['FBgn0013972',],method = 'spearman')

  saveRDS(ResTable.gene.part, file='DE.results.rds')
  Analysis.Data$ResTable.gene.part = ResTable.gene.part
  Analysis.Data$ResTable.pred.part = ResTable[1:length(all.preds.nr),]
  Analysis.Data$sorted.norm.CPMs = sorted.norm.CPMs
  if(Pars$report.norm.read.counts.instead.of.CPM)
    Analysis.Data$sorted.norm.counts = sorted.norm.counts
  Analysis.Data$Main.Component = Main.Component
  Analysis.Data$schema = schema
  Analysis.Data$Current.Predictors.list = Current.Predictors.list
  Analysis.Data$samples..sorted.by.main.component = source.samples.seq
  Analysis.Data$samples..reclustered = current.samples.seq
  Analysis.Data$clustering.samples.passed = clustering.samples.passed
  Analysis.Data$Main.Component.values..to..group.names = Main.Component.values..to..group.names
  Analysis.Data$perform.classic.paired.test = perform.classic.paired.test
  Analysis.Data$perform.FC.associations.paired.test = perform.FC.associations.paired.test
  Analysis.Data$Main.test = Pars$Main.test
  Analysis.Data$CPMs.groups = CPMs.groups
  # Analysis.Data$miR.mode = Startup.Data$miR.mode

  suppressWarnings(rm(counts))
  suppressWarnings(rm(noint))
  suppressWarnings(rm(cpms))
  suppressWarnings(rm(keep))
  suppressWarnings(rm(d))
  suppressWarnings(rm(ResTable))
  rm(ResTable.gene.part)
  rm(sorted.norm.CPMs)
  rm(norm.CPMs)
  if(Pars$report.norm.read.counts.instead.of.CPM){
    rm(sorted.norm.counts)
    rm(norm.counts)
  }

  
  Analysis.Data$GLM.analysis.performed = TRUE
  if(include.GLM.model.name.in.Result.names){
    saveRDS(Analysis.Data, file = sprintf('%s, analysis data.rds',Analysis.Data$Analysis.name))
  } else {
    saveRDS(Analysis.Data, file = sprintf('analysis data.rds'))
  }
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
#   #sorted.norm.CPMs = ResTable.gene.part[,14:dim(ResTable.gene.part)[2]]
#   rn.refseq = Analysis.Data$ResTable.gene.part[rownames(Analysis.Data$sorted.norm.CPMs),"Gene name"]
#   rn.ensembl = rownames(Analysis.Data$sorted.norm.CPMs)
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
#   sorted.norm.CPMs.for.heatmaps = Analysis.Data$sorted.norm.CPMs
#   rownames(sorted.norm.CPMs.for.heatmaps) = rn.refseq
#   
#   
#   
#   #x=50
#   for (x in Pars$Top.genes.to.include.in.heatmaps.list){
#     x = min(x,dim(sorted.norm.CPMs.for.heatmaps)[1])
#     if (x == dim(sorted.norm.CPMs.for.heatmaps)[1]){
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
#     limit = quantile(sorted.norm.CPMs.for.heatmaps[1:x,],0.983)
#     par(cex.main=0.6)
#     max.chars = max(nchar(colnames(sorted.norm.CPMs.for.heatmaps)))
#     heatmap.2(main = Heatmap_Title, data.matrix(pmin(sorted.norm.CPMs.for.heatmaps[1:x,],limit))**0.5,Colv = FALSE,Rowv = T,scale = 'none',
#               col = hmcols, key=TRUE, symkey=FALSE,
#               density.info="none", trace="none",
#               cexRow=0.55/x*50,cexCol=min(1, 0.8/max.chars*10),dendrogram='row')
#     dev.off()
#     
#     
#     if (x == dim(sorted.norm.CPMs.for.heatmaps)[1]){
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
#     intra.norm = sorted.norm.CPMs.for.heatmaps/rowSums(sorted.norm.CPMs.for.heatmaps)
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
                           Analysis.Data.RDS.file = NULL,
                           forced.parameters = NULL,
                           limit.to.genes = NULL, out.dir.root = NULL, max.Genes.counts.list = NULL,
                           add.Heatmaps.with.Custom.samples.order = NULL, Heatmaps.with.Custom.samples.order.dir = '.',
                           add.grid = FALSE,
                           annotation.logFC.range = 4, annotation.log10P.range = 12,
                           include.GLM.model.name.in.Result.names = FALSE,
                           use.provided.Sample.annotation..if.available = TRUE){
  
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model, or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  if(is.null(GLM.model))  GLM.model = Analysis.Data$GLM.model
  
  if (!Analysis.Data$GLM.analysis.performed) stop('Create.Heatmaps can be processed only after GLM/DE analysis. It is not performed.')

  Individual.Pars = NULL
  if(GLM.model %in% colnames(Startup.Data$Ind.Par.table)){
    tmp = subset(Startup.Data$Ind.Par.table, select = GLM.model, subset = !is.na(Startup.Data$Ind.Par.table[,GLM.model]))
    par.table = cbind(as.data.frame(rownames(tmp)), tmp)
    colnames(par.table) = c('variable','value')
    cat(sprintf('Found %d individual parameters:\n', nrow(par.table)))
    for(x in 1:nrow(par.table)){
      cat(sprintf('      %s:    %s\n', par.table[x,1], par.table[x,2]))
    }
    Pars = Process.Parameters.table(Pars, par.table, show.warnings = FALSE, suppress.base.dir.errors = TRUE)
  }

  if(is.null(out.dir.root))  out.dir.root = Analysis.Data$current.GLM.results.dir

  dir.create(out.dir.root, showWarnings = FALSE)

  #!is.null(custom.sample.order) + cluster.all.samples + !is.null(predictor.to.sort) + use.GLM.sorted.samples + use.GLM.clustered.samples

  if(Pars$Create.heatmaps..with.Src.sample.order){
    cat(sprintf("~ %s:   creating heatmaps - intact samples order...\n", GLM.model))
    out.dir = sprintf("%s", out.dir.root)
    Create.Heatmaps__internal(Startup.Data = Startup.Data, Analysis.Data = Analysis.Data, GLM.model = NULL, GLM.working.dir = GLM.working.dir,
                             Analysis.Data.RDS.file = NULL, forced.parameters = Pars,
                             limit.to.genes = limit.to.genes, out.dir = out.dir, max.Genes.counts.list = max.Genes.counts.list,
                             cluster.all.samples = FALSE,
                             predictor.to.sort = NULL,
                             custom.sample.order = NULL,
                             use.GLM.sorted.samples = FALSE,
                             use.GLM.clustered.samples = FALSE,
                             add.grid = add.grid,
                             annotation.logFC.range = annotation.logFC.range, annotation.log10P.range = annotation.log10P.range,
                             include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                             use.provided.Sample.annotation..if.available = use.provided.Sample.annotation..if.available)
  }

  # print(Analysis.Data$samples..sorted.by.main.component)
  if(Pars$Create.heatmaps..with.Sorted.samples & (any(Analysis.Data$samples..sorted.by.main.component != colnames(Analysis.Data$sorted.norm.CPMs)) | (!Pars$Create.heatmaps..with.Src.sample.order & !(Pars$Create.heatmaps..with.Clustered.samples & Analysis.Data$clustering.samples.passed)))){
    cat(sprintf("~ %s:   creating heatmaps - sorted by main predictor...\n", GLM.model))
    out.dir = sprintf("%s/Heatmaps - sorted samples", out.dir.root)
    Create.Heatmaps__internal(Startup.Data = Startup.Data, Analysis.Data = Analysis.Data, GLM.model = NULL, GLM.working.dir = GLM.working.dir,
                             Analysis.Data.RDS.file = NULL, forced.parameters = Pars,
                             limit.to.genes = limit.to.genes, out.dir = out.dir, max.Genes.counts.list = max.Genes.counts.list,
                             cluster.all.samples = FALSE,
                             predictor.to.sort = NULL,
                             custom.sample.order = NULL,
                             use.GLM.sorted.samples = TRUE,
                             use.GLM.clustered.samples = FALSE,
                             add.grid = add.grid,
                             annotation.logFC.range = annotation.logFC.range, annotation.log10P.range = annotation.log10P.range,
                             include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                             use.provided.Sample.annotation..if.available = use.provided.Sample.annotation..if.available)
  }

  if(Pars$Create.heatmaps..with.Clustered.samples & Pars$Cluster.samples.in.Excel.results & (Analysis.Data$clustering.samples.passed | (!Pars$Create.heatmaps..with.Src.sample.order & !Pars$Create.heatmaps..with.Sorted.samples))){
    cat(sprintf("~ %s:   creating heatmaps - clustered samples...\n", GLM.model))
    out.dir = sprintf("%s/Heatmaps - clustered samples", out.dir.root)
    Create.Heatmaps__internal(Startup.Data = Startup.Data, Analysis.Data = Analysis.Data, GLM.model = NULL, GLM.working.dir = GLM.working.dir,
                             Analysis.Data.RDS.file = NULL, forced.parameters = Pars,
                             limit.to.genes = limit.to.genes, out.dir = out.dir, max.Genes.counts.list = max.Genes.counts.list,
                             cluster.all.samples = FALSE,
                             predictor.to.sort = NULL,
                             custom.sample.order = NULL,
                             use.GLM.sorted.samples = FALSE,
                             use.GLM.clustered.samples = TRUE,
                             add.grid = add.grid,
                             annotation.logFC.range = annotation.logFC.range, annotation.log10P.range = annotation.log10P.range,
                             include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                             use.provided.Sample.annotation..if.available = use.provided.Sample.annotation..if.available)
  }

  if(Pars$Create.heatmaps..with.All.samples & nrow(Analysis.Data$schema) < nrow(Startup.Data$schema)){
  # if(Pars$Create.heatmaps..with.All.samples){
    cat(sprintf("~ %s:   creating heatmaps - complete samples list...\n", GLM.model))
    out.dir = sprintf("%s/Heatmaps - all samples", out.dir.root)
    Create.Heatmaps__internal(Startup.Data = Startup.Data, Analysis.Data = Analysis.Data, GLM.model = NULL, GLM.working.dir = GLM.working.dir,
                             Analysis.Data.RDS.file = NULL, forced.parameters = Pars,
                             limit.to.genes = limit.to.genes, out.dir = out.dir, max.Genes.counts.list = max.Genes.counts.list,
                             cluster.all.samples = FALSE,
                             predictor.to.sort = NULL,
                             custom.sample.order = NULL,
                             use.GLM.sorted.samples = FALSE,
                             use.GLM.clustered.samples = FALSE,
                             custom.CPMs.data.file = sprintf('%s/CPMs.rds', Pars$results.dir),
                             custom.sample.annotations = Startup.Data$predictor.values.nr,
                             add.grid = add.grid,
                             annotation.logFC.range = annotation.logFC.range, annotation.log10P.range = annotation.log10P.range,
                             include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                             use.provided.Sample.annotation..if.available = use.provided.Sample.annotation..if.available)
  }

  if(!is.null(add.Heatmaps.with.Custom.samples.order)){
    cat(sprintf("~ %s:   creating heatmaps - custom samples order...\n", GLM.model))
    out.dir = sprintf("%s/%s", out.dir.root, Heatmaps.with.Custom.samples.order.dir)
    Create.Heatmaps__internal(Startup.Data = Startup.Data, Analysis.Data = Analysis.Data, GLM.model = NULL, GLM.working.dir = GLM.working.dir,
                             Analysis.Data.RDS.file = NULL, forced.parameters = Pars,
                             limit.to.genes = limit.to.genes, out.dir = out.dir, max.Genes.counts.list = max.Genes.counts.list,
                             cluster.all.samples = FALSE,
                             predictor.to.sort = NULL,
                             custom.sample.order = add.Heatmaps.with.Custom.samples.order,
                             use.GLM.sorted.samples = FALSE,
                             use.GLM.clustered.samples = FALSE,
                             add.grid = add.grid,
                             annotation.logFC.range = annotation.logFC.range, annotation.log10P.range = annotation.log10P.range,
                             include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                             use.provided.Sample.annotation..if.available = use.provided.Sample.annotation..if.available)
  }

}


Create.Heatmaps__internal = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL,
                           Analysis.Data.RDS.file = NULL, forced.parameters = NULL,
                           limit.to.genes = NULL, out.dir = NULL, max.Genes.counts.list = NULL,
                           cluster.all.samples = FALSE, predictor.to.sort = NULL, custom.sample.order = NULL,
                           use.GLM.sorted.samples = FALSE, use.GLM.clustered.samples = FALSE,
                           custom.CPMs.data.file = NULL, custom.sample.annotations = NULL,
                           add.grid = FALSE, annotation.logFC.range = 4, annotation.log10P.range = 12,
                           include.GLM.model.name.in.Result.names = FALSE,
                           use.provided.Sample.annotation..if.available = TRUE){
  # print(custom.sample.order)
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  if(is.null(GLM.model))  GLM.model = Analysis.Data$GLM.model
  
  if (!Analysis.Data$GLM.analysis.performed) stop('Create.Heatmaps can be processed only after GLM/DE analysis. It is not performed.')

  # if (bypass.if.completed & Read.Completed.steps.status(sprintf("%s.%s.%s.heatmaps",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
  #   message(sprintf('\nCreating heatmaps for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
  #   return()
  # }

  
  if(!is.null(custom.sample.order) + cluster.all.samples + !is.null(predictor.to.sort) + use.GLM.sorted.samples + use.GLM.clustered.samples > 1 ){
    msg = '\nCreate.Heatmaps(...):  Only one of the following parameters should not be provided/turned ON: custom.sample.order, cluster.all.samples, predictor.to.sort, use.GLM.sorted.samples, use.GLM.clustered.samples\n'
    stop(msg)
  }

    
  #############################################################
  ########### Creating Heatmaps
  
  if(is.null(out.dir))  out.dir = Analysis.Data$current.GLM.results.dir
  dir.create(out.dir, showWarnings = FALSE)
  
  # cat(sprintf("\n~ %s:   creating heatmaps...\n", GLM.model))
  #sorted.norm.CPMs = ResTable.gene.part[,14:dim(ResTable.gene.part)[2]]
  
  if(!is.null(limit.to.genes)){
    keep.by.Gene.list = (rownames(Analysis.Data$sorted.norm.CPMs) %in% limit.to.genes)
  } else {
    keep.by.Gene.list = rep(TRUE, nrow(Analysis.Data$sorted.norm.CPMs))
  }

  LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
  # print(Pars$Main.test)
  PValue.col.name = Get.PValue.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part))
  FDR.col.name = Get.FDR.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part), look.for.paired.test.res = Analysis.Data$perform.classic.paired.test)
  NP.PValue.col.name = Get.NP.PValue.col.name(colnames(Analysis.Data$ResTable.gene.part))
  LogCPM.col.name = Get.LogCPM.col.name(colnames(Analysis.Data$ResTable.gene.part))

  keep.by.LogFC =  abs(Analysis.Data$ResTable.gene.part[rownames(Analysis.Data$sorted.norm.CPMs), LogFC.col.name]) >= Pars$Heatmaps.abs.LogFC.threshold
  keep.by.PValue = Analysis.Data$ResTable.gene.part[rownames(Analysis.Data$sorted.norm.CPMs), PValue.col.name] <= Pars$Heatmaps.PValue.threshold
  keep.by.FDR =    Analysis.Data$ResTable.gene.part[rownames(Analysis.Data$sorted.norm.CPMs), FDR.col.name] <= Pars$Heatmaps.FDR.threshold
  keep.by.LogCPM = Analysis.Data$ResTable.gene.part[rownames(Analysis.Data$sorted.norm.CPMs), LogCPM.col.name] >= Pars$Heatmaps.LogCPM.threshold

  keep = keep.by.Gene.list & keep.by.LogFC & keep.by.PValue & keep.by.FDR & keep.by.LogCPM
  if(!is.na(NP.PValue.col.name))
    keep = keep & (Analysis.Data$ResTable.gene.part[rownames(Analysis.Data$sorted.norm.CPMs), NP.PValue.col.name] <= Pars$Heatmaps.NP.PValue.threshold)

  keep[which(is.na(keep))] = FALSE
  
  if(sum(keep) == 0){
    cat('\n[Create.Heatmaps]  No genes to render\n')
    return(invisible())
  } else {
    cat(sprintf('%d of %d genes passed thresholds p < %g, FDR < %g, non-parametric p < %g, LogCPM > %g, abs(LogFC) > %g\n', sum.mod(keep), length(keep),
      Pars$Heatmaps.PValue.threshold, Pars$Heatmaps.FDR.threshold, Pars$Heatmaps.NP.PValue.threshold,
      Pars$Heatmaps.LogCPM.threshold, Pars$Heatmaps.abs.LogFC.threshold))
  }
  
  rn.refseq = Analysis.Data$ResTable.gene.part[rownames(Analysis.Data$sorted.norm.CPMs)[keep],"Symbol"]
  rn.ensembl = rownames(Analysis.Data$sorted.norm.CPMs)[keep]
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

  names(rn.refseq) = rn.ensembl

  if(!is.null(custom.CPMs.data.file)){
    sorted.norm.CPMs.for.heatmaps = readRDS(file = custom.CPMs.data.file)
    desired.genes = rownames(Analysis.Data$sorted.norm.CPMs[keep,])
    common.genes = intersect(rownames(sorted.norm.CPMs.for.heatmaps), desired.genes)
    rn.refseq = rn.refseq[desired.genes %in% common.genes]
    desired.genes = desired.genes[desired.genes %in% common.genes]
    sorted.norm.CPMs.for.heatmaps = sorted.norm.CPMs.for.heatmaps[desired.genes, ]
    rownames(sorted.norm.CPMs.for.heatmaps) = rn.refseq
    corresponding.Ensembl.IDs = desired.genes
  } else {
    sorted.norm.CPMs.for.heatmaps = Analysis.Data$sorted.norm.CPMs[keep,]
    rownames(sorted.norm.CPMs.for.heatmaps) = rn.refseq
    corresponding.Ensembl.IDs = rownames(Analysis.Data$sorted.norm.CPMs[keep,])
  }

  
  LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
  PValue.col.name = Get.PValue.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part))
   
  genes.annotation = as.data.frame(Analysis.Data$ResTable.gene.part[corresponding.Ensembl.IDs, c(LogFC.col.name, 'LogCPM')])
  colnames(genes.annotation)[which(colnames(genes.annotation) %in% LogFC.col.name)] = 'LogFC'

  if(is.null(annotation.logFC.range) | 'auto' %in% annotation.logFC.range){
    annotation.logFC.range = (mean(abs(genes.annotation[1:50, 'LogFC'])) + quantile(abs(genes.annotation[ ,'LogFC']), 0.99)) / 2
    annotation.logFC.range %<>% min(.,6) %>% max(., 1.5) %>% round(.,1)
  }
  if(is.null(annotation.log10P.range) | 'auto' %in% annotation.log10P.range){
    annotation.log10P.range = 15
    #annotation.log10P.range = (mean(abs(genes.annotation[1:50,'logFC'])) + quantile(abs(genes.annotation[,'logFC']), 0.99)) / 2
    #annotation.log10P.range %<>% min(.,6) %>% max(., 1.5) %>% round(.,1)
  }
  
  genes.annotation[,'LogFC'] %<>% pmin(. , annotation.logFC.range) %>% pmax (., -annotation.logFC.range)
  genes.annotation[,'LogCPM'] %<>% pmin(. , 10) %>% pmax (., 0)
  genes.annotation = cbind(genes.annotation, (-1)*log10(Analysis.Data$ResTable.gene.part[corresponding.Ensembl.IDs, PValue.col.name] + 1e-100))
  colnames(genes.annotation)[3] <- '-Log10(p)'
  genes.annotation[,'-Log10(p)'] %<>% pmin(. , annotation.log10P.range) %>% pmax (., 0)

  #head(Analysis.Data$ResTable.gene.part$'Spearman r')
  Spearman_r.is.present = FALSE
  if('Spearman r' %in% colnames(Analysis.Data$ResTable.gene.part)){
    if(any(!is.na(Analysis.Data$ResTable.gene.part[corresponding.Ensembl.IDs, 'Spearman r']))){
      genes.annotation = cbind(genes.annotation, Analysis.Data$ResTable.gene.part[corresponding.Ensembl.IDs, 'Spearman r'])
      colnames(genes.annotation)[4] <- 'Spearman_r'
      Spearman_r.is.present = TRUE
    }
  }
  
  rownames(genes.annotation) = rownames(sorted.norm.CPMs.for.heatmaps)

  #sample.order = 1:ncol(sorted.norm.CPMs.for.heatmaps)

  if(!is.null(custom.sample.order)){
    sorted.norm.CPMs.for.heatmaps = sorted.norm.CPMs.for.heatmaps[, custom.sample.order]
  } else if(!is.null(predictor.to.sort)){
    sample.annotations = as.data.frame(Analysis.Data$predictor.values)
    if(predictor.to.sort %in% colnames(sample.annotations)){
      sample.order = order(sample.annotations[,predictor.to.sort])
      sorted.norm.CPMs.for.heatmaps = sorted.norm.CPMs.for.heatmaps[,sample.order]
    } else if (!msg.sent){
      msg = sprintf('\nPredictor %s is not found among sample annotations. Sorting will not be applied\n', predictor.to.sort)
      stop(msg)
    }
  } else if(use.GLM.clustered.samples){
    sorted.norm.CPMs.for.heatmaps = sorted.norm.CPMs.for.heatmaps[, Analysis.Data$samples..reclustered]
  } else if(use.GLM.sorted.samples){
    sorted.norm.CPMs.for.heatmaps = sorted.norm.CPMs.for.heatmaps[, Analysis.Data$samples..sorted.by.main.component]
  }

  if(is.null(max.Genes.counts.list))  max.Genes.counts.list = Pars$Top.genes.to.include.in.heatmaps.list
  max.Genes.counts.list = sapply(max.Genes.counts.list, function(x) { min(x, nrow(sorted.norm.CPMs.for.heatmaps)) })
  max.Genes.counts.list = max.Genes.counts.list[!duplicated(max.Genes.counts.list)]

  ref.samples = c()
  
  sorted.Log.Rel.CPM.for.heatmaps = NULL

  all.sorted.norm.CPMs.for.heatmaps__ncol = ncol(sorted.norm.CPMs.for.heatmaps)

  if(Pars$Create.heatmaps..Log.rel.CPM){
    if(Pars$Create.heatmaps..Log.rel.CPM..using.ext.samples){
      src.CPMs.filename = sprintf('%s/CPMs.rds', Pars$results.dir)
      if(file.exists(src.CPMs.filename)){
        all.sorted.norm.CPMs.for.heatmaps = readRDS(file = sprintf('%s/CPMs.rds', Pars$results.dir))

        all.sorted.norm.CPMs.for.heatmaps = all.sorted.norm.CPMs.for.heatmaps[rownames(all.sorted.norm.CPMs.for.heatmaps) %in% names(rn.refseq), ]
        rownames(all.sorted.norm.CPMs.for.heatmaps) = rn.refseq[rownames(all.sorted.norm.CPMs.for.heatmaps)]

        common.genes = intersect(rownames(all.sorted.norm.CPMs.for.heatmaps), rownames(sorted.norm.CPMs.for.heatmaps))
        if(length(common.genes) == 0){
          stop('Something strange. Gene lists in CPMs.rds and sorted.norm.CPMs.for.heatmaps are completely different\n')
        } else {
          all.sorted.norm.CPMs.for.heatmaps = all.sorted.norm.CPMs.for.heatmaps[rownames(sorted.norm.CPMs.for.heatmaps)[rownames(sorted.norm.CPMs.for.heatmaps) %in% rownames(all.sorted.norm.CPMs.for.heatmaps)], ]
        }
      } else {
        stop(sprintf('Source file %s is not found. Skipping loading src CPMs info\n', src.CPMs.filename))
      }
    } else {
      all.sorted.norm.CPMs.for.heatmaps = sorted.norm.CPMs.for.heatmaps
    }

    ref.samples = c()
    if(Pars$Create.heatmaps..Log.rel.CPM..using.ext.samples){
      if(Pars$Create.heatmaps..Log.rel.CPM..include.ref.samples.with.suffix != ""){
        Pars$Create.heatmaps..Log.rel.CPM..include.ref.samples.with.suffix = gsub('"', "", Pars$Create.heatmaps..Log.rel.CPM..include.ref.samples.with.suffix, fixed = TRUE)
        Pars$Create.heatmaps..Log.rel.CPM..include.ref.samples.with.suffix = gsub("'", "", Pars$Create.heatmaps..Log.rel.CPM..include.ref.samples.with.suffix, fixed = TRUE)
        keep = endsWith(colnames(all.sorted.norm.CPMs.for.heatmaps), Pars$Create.heatmaps..Log.rel.CPM..include.ref.samples.with.suffix)
        if(sum.mod(keep) > 0){
          ref.samples = colnames(all.sorted.norm.CPMs.for.heatmaps)[keep]
        } else {
          cat(sprintf('No reference samples that end with suffix "%s" are found. Log.rel.CPM heatmaps with custom normal samples will not be created\n', Pars$Create.heatmaps..Log.rel.CPM..include.ref.samples.with.suffix))
          ref.samples = c()
        }
      } else {
        ref.samples = colnames(all.sorted.norm.CPMs.for.heatmaps)
      }

      if(length(ref.samples) > 0){
        c.sorted.norm.CPMs.for.heatmaps = sorted.norm.CPMs.for.heatmaps[rownames(sorted.norm.CPMs.for.heatmaps) %in% rownames(all.sorted.norm.CPMs.for.heatmaps), ]
        if(any(rownames(all.sorted.norm.CPMs.for.heatmaps) != rownames(c.sorted.norm.CPMs.for.heatmaps))){
          stop(sprintf('Create.Heatmaps__internal needs debug. %d of %d rownames are different', sum.mod(any(rownames(all.sorted.norm.CPMs.for.heatmaps) != rownames(c.sorted.norm.CPMs.for.heatmaps))), dim(all.sorted.norm.CPMs.for.heatmaps)[1]))
        }
        ref.samples.norm.CPMs = subset(all.sorted.norm.CPMs.for.heatmaps, select = which(colnames(all.sorted.norm.CPMs.for.heatmaps) %in% ref.samples))
        
        if(length(ref.samples) >= Pars$add.Trimmed.LogFC.if.Samples.N.is.greater.than){
          ref.CPM.values = apply(ref.samples.norm.CPMs, 1, function(x) {  mean(x, trim = 0.1)  })
        } else {
          ref.CPM.values = apply(ref.samples.norm.CPMs, 1, mean)
        }
        
        sorted.Log.Rel.CPM.for.heatmaps = log2((c.sorted.norm.CPMs.for.heatmaps + Pars$Create.heatmaps..Log.rel.CPM..const.CPM.add) / (ref.CPM.values + Pars$Create.heatmaps..Log.rel.CPM..const.CPM.add))
        sorted.Log.Rel.CPM.for.heatmaps = pmin(pmax(sorted.Log.Rel.CPM.for.heatmaps, -Pars$Create.heatmaps..Log.rel.CPM..log.range), Pars$Create.heatmaps..Log.rel.CPM..log.range)
        rm(c.sorted.norm.CPMs.for.heatmaps)
      }
      all.sorted.norm.CPMs.for.heatmaps__ncol = ncol(all.sorted.norm.CPMs.for.heatmaps)
      rm(all.sorted.norm.CPMs.for.heatmaps)
    }
  }
  
  
  #cluster.samples = TRUE
  x=max.Genes.counts.list[length(max.Genes.counts.list)]
  x=120
  for (x.src in max.Genes.counts.list){
    cat(sprintf("\r  processing heatmap for top %d DE genes...", x.src))
    x = min(x.src, dim(sorted.norm.CPMs.for.heatmaps)[1])
    
    #genes.annotation = rownames(sorted.norm.CPMs.for.heatmaps[1:x,])
    
    hmcols<-colorRampPalette(c("white","yellow","orange","red","blue"))(256)
    limit = quantile(sorted.norm.CPMs.for.heatmaps[1:x,],0.983)
    # par(cex.main=0.4)
    max.chars = max(nchar(colnames(sorted.norm.CPMs.for.heatmaps)))
    
    if(use.provided.Sample.annotation..if.available && !is.null(Startup.Data$Heatmaps.annotation)){
      common.samples = colnames(sorted.norm.CPMs.for.heatmaps)[colnames(sorted.norm.CPMs.for.heatmaps) %in% rownames(Startup.Data$Heatmaps.annotation)]
      if(length(common.samples) == 0){  sample.annotations = NULL
      } else sample.annotations = Startup.Data$Heatmaps.annotation[common.samples,, drop = FALSE]
    } else if(!is.null(custom.sample.annotations)){
      sample.annotations = as.data.frame(custom.sample.annotations)
    } else {
      sample.annotations = as.data.frame(Analysis.Data$predictor.values.nr)
    }
    #sample.annotations = sample.annotations[,!endsWith(colnames(sample.annotations),' D')]
    
    # keep = !(apply(sample.annotations, 2, function(x) { length(levels(as.factor(x)))}) %in% 1)
    # sample.annotations.clf = sample.annotations[, keep, drop = FALSE]
    if(!is.null(sample.annotations)){
      sample.annotations.clf = sample.annotations[, apply(sample.annotations, 2, function(x) { !all(is.na(x)) & length(x[!duplicated(x)]) > 1 }), drop = FALSE]
    } else  sample.annotations.clf = NULL
    

    current.genes.annotation = genes.annotation[1:x,]
    logFC.colors <- colorRampPalette(c("#4875e6","white","#f95136"))(512)
    Spearman_r.colors <- colorRampPalette(c("#0e8bec","white","#f13004"))(512)
    P.colors <- colorRampPalette(c("white","#f1ab04","#69cd0f"))(512)
    logCPM.colors <- colorRampPalette(c("white", "#ecab0e"))(512)
    
    logFC.color.coord.start <- (length(logFC.colors) * (min(current.genes.annotation$LogFC) + annotation.logFC.range) / 2 / annotation.logFC.range) %>% round %>% min(., length(logFC.colors) - 1) %>% max(., 0)
    logFC.color.coord.end <- (length(logFC.colors) * (max(current.genes.annotation$LogFC) + annotation.logFC.range) / 2 / annotation.logFC.range) %>% round %>% min(., length(logFC.colors)) %>% max(., 0)
    P.color.coord.start <- (length(P.colors) * (min(current.genes.annotation$`-Log10(p)`) + 0) / annotation.log10P.range) %>% round %>% min(., length(P.colors) - 1) %>% max(., 0)
    P.color.coord.end <- (length(P.colors) * (max(current.genes.annotation$`-Log10(p)`) + 0) / annotation.log10P.range) %>% round %>% min(., length(P.colors)) %>% max(., 0)
    logCPM.color.coord.start <- (length(logCPM.colors) * (min(current.genes.annotation$LogCPM) + 0) / 10) %>% round %>% min(., length(logCPM.colors) - 1) %>% max(., 0)
    logCPM.color.coord.end <- (length(logCPM.colors) * (max(current.genes.annotation$LogCPM) + 0) / 10) %>% round %>% min(., length(logCPM.colors)) %>% max(., 0)
    if(Spearman_r.is.present){
      current.genes.annotation$Spearman_r[is.na(current.genes.annotation$Spearman_r)] = 0
      Speaman_r.color.coord.start <- (length(Spearman_r.colors) * (min(current.genes.annotation$Spearman_r) + 1) / 2) %>% round %>% min(., length(Spearman_r.colors) - 1) %>% max(., 0)
      Speaman_r.color.coord.end <- (length(Spearman_r.colors) * (max(current.genes.annotation$Spearman_r) + 1) / 2) %>% round %>% min(., length(Spearman_r.colors)) %>% max(., 0)
    }
    
    ann.colors = list(
      LogFC = logFC.colors [logFC.color.coord.start : logFC.color.coord.end],
      LogCPM = logCPM.colors[logCPM.color.coord.start : logCPM.color.coord.end]
    )
    ann.colors[['-Log10(p)']] = P.colors[P.color.coord.start : P.color.coord.end]
    if(Spearman_r.is.present){
      ann.colors[['Spearman_r']] = Spearman_r.colors[Speaman_r.color.coord.start : Speaman_r.color.coord.end]
    }
    
    width.min = 8
    width.max = 19
    hm.width = max(width.min, min(width.max, dim(sorted.norm.CPMs.for.heatmaps)[2]*17/29))
    if(is.null(sample.annotations)){  hm.height = min(16, max(8, 8*(ncol(genes.annotation))/9))
    } else hm.height = min(16, max(8, 8*(ncol(sample.annotations) + ncol(genes.annotation))/9))
    
    if(add.grid == TRUE){
      border_color = "grey60"
    } else if(add.grid == FALSE){
      border_color = NA
    } else {
      border_color = if (x <= 70) "grey60" else NA
    }

    
    if(is.null(sample.annotations.clf))  sample.annotations.clf.char = NULL
    if(Pars$Discretize.sample.annotations & !is.null(sample.annotations.clf)){
      #base.colors = c('#4f8cff', '#ff9319', '#5fd852', '#95e84c', '#0dd6d6', '#77a53e', '#ba801d', '#e23d76', '#9553e0')
      base.colors..pallete = c('#f59342', '#f0bf0d', '#a1d825', '#79dcf1', '#81aeea')
      base.colors = colorRampPalette(base.colors..pallete)(ncol(sample.annotations.clf))
      if(ncol(sample.annotations.clf) == 1)
        base.colors = '#90d000'
      sample.annotations.clf.char = sample.annotations.clf

      predN = 1
      for(predN in 1:ncol(sample.annotations.clf)){
        c.lev = DeDup(sample.annotations.clf[, predN])
        zero.color = colorRampPalette(c("white", base.colors[predN]))(8)[2]
        if(length(c.lev) <= Pars$Max.Levels.to.discretize | !suppressWarnings(all(is.na(c.lev) | !is.na(as.numeric(as.character(c.lev)))))){
          c.lev = c.lev[order(c.lev)]
          conv.pred.anno = as.character(sample.annotations.clf[, predN])
          Predictor.values..to..group.names = Get..Predictor.values..to..group.names(conv.pred.anno[!is.na(conv.pred.anno)], colnames(sample.annotations.clf)[predN])

          sample.annotations.clf.char[, predN] = sapply(conv.pred.anno, function(x){ ifelse(is.na(x), NA, Predictor.values..to..group.names[[as.character(x)]][[1]]) })
          conv.pred.anno[which(is.na(conv.pred.anno))] = 'NA'


          if(predN <= length(base.colors)){
            if(NA %in% c.lev){
              c.lev.without.NA = c.lev[!is.na(c.lev)]
              c.col = c(colorRampPalette(c(zero.color, base.colors[predN]))(length(c.lev.without.NA)), '#d8d8d8')
              names(c.col) = c(sapply(c.lev.without.NA, function(x){ ifelse(is.na(x), NA, Predictor.values..to..group.names[[as.character(x)]][[1]]) }), 'NA')
            } else {
              c.col = colorRampPalette(c(zero.color, base.colors[predN]))(length(c.lev))
              names(c.col) = sapply(c.lev, function(x){ ifelse(is.na(x), NA, Predictor.values..to..group.names[[as.character(x)]][[1]]) })
            }
            ann.colors[[colnames(sample.annotations.clf)[predN]]] = c.col
          }
        } else {
          if(predN <= length(base.colors)){
            c.pal = colorRampPalette(c(zero.color, base.colors[predN]))(10)
            if(any(is.na(sample.annotations.clf[, predN]))){
              ann.colors[[colnames(sample.annotations.clf)[predN]]] = c(c.pal[2], c.pal[10])
            } else {
              ann.colors[[colnames(sample.annotations.clf)[predN]]] = c(c.pal[1], c.pal[10])
            }
          }
        }
      }
      #rownames(sample.annotations.clf.char) = rownames(sample.annotations.clf)
      sample.annotations.clf = sample.annotations.clf.char
    }
    
    if(Pars$Create.heatmaps..CPM){
      if (x == dim(sorted.norm.CPMs.for.heatmaps)[1]){
        insert_text = sprintf('all %d genes', x)
      } else {
        insert_text = sprintf('top %d DE genes', x)
      }
      Heatmap_Title = sprintf("~ %s, %s. sqrt(CPM)", GLM.model, insert_text, x)

      if(include.GLM.model.name.in.Result.names){
        png.mod(filename = sprintf('%s/%s, Heatmap, sqrt(CPM), %s.png', out.dir, Analysis.name, insert_text),
            units="in", width=hm.width, height = hm.height, pointsize=12, res=300)
      } else {
        png.mod(filename = sprintf('%s/Heatmap, sqrt(CPM), %s.png', out.dir, insert_text),
            units="in", width=hm.width, height = hm.height, pointsize=12, res=300)
      }

      # sorted.norm.CPMs.for.heatmaps = sorted.norm.CPMs.for.heatmaps[apply(sorted.norm.CPMs.for.heatmaps, 1, function(r) (length(DeDup.na.rm(r)) > 1)),,drop = FALSE]
      # print(sorted.norm.CPMs.for.heatmaps[(nrow(sorted.norm.CPMs.for.heatmaps) - 40):nrow(sorted.norm.CPMs.for.heatmaps), ])
      # print(limit)
      tryCatch(expr = {
        pheatmap(mat = data.matrix(pmin(sorted.norm.CPMs.for.heatmaps[1:x, ],limit))**0.5, scale = 'none', color = hmcols,
                 cluster_cols = cluster.all.samples, border_color = border_color,
                 main = Heatmap_Title, fontsize = 7, fontsize_row = min(9, 9/x*50) * hm.height / 8,
                 fontsize_col = min(22*min(0.5, 0.8/max.chars*10), 23*min(1, 30/dim(sorted.norm.CPMs.for.heatmaps)[2])),
                 annotation_col = sample.annotations.clf, annotation_row = current.genes.annotation,
                 annotation_colors = ann.colors)
      }, error = function (err){
        cat(sprintf('Heatmap creating has FAILED\n'))
      })
      
      dev.off()
    }
    # heatmap.2(main = Heatmap_Title, data.matrix(pmin(sorted.norm.CPMs.for.heatmaps[1:x,],limit))**0.5,Colv = FALSE,Rowv = T,scale = 'none',
    #           col = hmcols, key=TRUE, symkey=FALSE,
    #           density.info="none", trace="none",
    #           cexRow=0.55/x*50, cexCol=min(1, 0.8/max.chars*10),dendrogram='row')
    
    
    
    if(Pars$Create.heatmaps..ZScore){
      if (x == dim(sorted.norm.CPMs.for.heatmaps)[1]){
        insert_text = sprintf('all %d genes', x)
      } else {
        insert_text = sprintf('top %d DE genes', x)
      }
      Heatmap_Title = sprintf("~ %s, %s. CPM z-score", GLM.model, insert_text, x)

      if(include.GLM.model.name.in.Result.names){
        png.mod(filename = sprintf('%s/%s, Heatmap, CPM Z-scores, %s.png', out.dir, Analysis.name, insert_text),
            units="in", width=hm.width, height = hm.height, pointsize=12, res=300)
      } else {
        png.mod(filename = sprintf('%s/Heatmap, CPM Z-scores, %s.png', out.dir, insert_text),
            units="in", width=hm.width, height = hm.height, pointsize=12, res=300)
      }
      
      intra.norm = sorted.norm.CPMs.for.heatmaps[1:x,]
      intra.norm = intra.norm[apply(intra.norm, 1, function(r) (length(DeDup.na.rm(r)) > 1)),,drop = FALSE]
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
      
      # par(cex.main=0.4)
      hmcols<-colorRampPalette(c("#7372b2","#78abf8",'white',"#ffb31f","#ff6521"))(512)
      #print(intra.norm[(nrow(intra.norm) - 40):nrow(intra.norm), ])
      # intra.norm = intra.norm[!(rownames(intra.norm) %in% 'cel-miR-39'),,drop=FALSE]
      tryCatch(expr = {
        pheatmap(mat = intra.norm, scale = 'row', color = hmcols, cluster_cols = cluster.all.samples,
                 border_color = border_color,
                 main = Heatmap_Title, fontsize = 7, fontsize_row = min(9, 9/x*50) * hm.height / 8,
                 fontsize_col = min(22*min(0.5, 0.8/max.chars*10), 23*min(1, 30/dim(sorted.norm.CPMs.for.heatmaps)[2])),
                 annotation_col = sample.annotations.clf, annotation_row = current.genes.annotation,
                 annotation_colors = ann.colors)
      }, error = function (err){
        cat(sprintf('Heatmap creating has FAILED\n'))
      })
      
      # heatmap.2(main = Heatmap_Title, intra.norm[1:x,]**1,Colv = FALSE,Rowv = T,scale = 'row',
      #           col = hmcols, key=TRUE, symkey=FALSE,
      #           density.info="none", trace="none",
      #           cexRow=0.55/x*50,cexCol=min(1, 0.8/max.chars*10),dendrogram='row')
      dev.off()
    }

    if(Pars$Create.heatmaps..Log.rel.CPM){
      x = min(x.src, dim(sorted.norm.CPMs.for.heatmaps)[1])
      #Pars$Create.heatmaps..Log.rel.CPM..log.range

      current.sorted.norm.CPMs.for.heatmaps = sorted.norm.CPMs.for.heatmaps[1:x, ]

      if(ncol(current.sorted.norm.CPMs.for.heatmaps) >= Pars$add.Trimmed.LogFC.if.Samples.N.is.greater.than){
        ref.CPM.values = apply(current.sorted.norm.CPMs.for.heatmaps, 1, function(x) {  mean(x, trim = 0.1)  })
      } else {
        ref.CPM.values = apply(current.sorted.norm.CPMs.for.heatmaps, 1, mean)
      }
      
      current.sorted.Log.Rel.CPM.for.heatmaps = log2((current.sorted.norm.CPMs.for.heatmaps + Pars$Create.heatmaps..Log.rel.CPM..const.CPM.add) / (ref.CPM.values + Pars$Create.heatmaps..Log.rel.CPM..const.CPM.add))
      current.sorted.Log.Rel.CPM.for.heatmaps = pmin(pmax(current.sorted.Log.Rel.CPM.for.heatmaps, -Pars$Create.heatmaps..Log.rel.CPM..log.range), Pars$Create.heatmaps..Log.rel.CPM..log.range)

      if (x == dim(sorted.norm.CPMs.for.heatmaps)[1]){
        insert_text = sprintf('all %d genes', x)
      } else {
        insert_text = sprintf('top %d DE genes', x)
      }

      Heatmap_Title = sprintf("~ %s, %s. Log2(CPM / avg.CPM)", GLM.model, insert_text)

      if(include.GLM.model.name.in.Result.names){
        png.mod(filename = sprintf('%s/%s, Heatmap, Log.rel.CPM, %s.png', out.dir, Analysis.name, insert_text),
            units="in", width=hm.width, height = hm.height, pointsize=12, res=300)
      } else {
        png.mod(filename = sprintf('%s/Heatmap, Log.rel.CPM, %s.png', out.dir, insert_text),
            units="in", width=hm.width, height = hm.height, pointsize=12, res=300)
      }

      # print(dim(sorted.Log.Rel.CPM.for.heatmaps))
      # print(dim(intra.norm))
      # print(x)
      # print(data.matrix(sorted.Log.Rel.CPM.for.heatmaps[1:x, ]))
      #stop('')
      # print(intra.norm)
      ####  encorrect x!!!!

      # hmcols2 = colorRampPalette(c("#7372b2","#78abf8",'white',"#ffb31f","#ff6521"))(512)
      Log.rel.CPM..log.range = Pars$Create.heatmaps..Log.rel.CPM..log.range
      breaksList = seq(-Log.rel.CPM..log.range, Log.rel.CPM..log.range, by = 2*Log.rel.CPM..log.range/(length(hmcols) - 1))
      
      tryCatch(expr = {
        pheatmap(mat = data.matrix(current.sorted.Log.Rel.CPM.for.heatmaps), scale = 'none', color = hmcols,
                 cluster_cols = cluster.all.samples, border_color = border_color,
                 main = Heatmap_Title, fontsize = 7, fontsize_row = min(9, 9/x*50) * hm.height / 8,
                 fontsize_col = min(22*min(0.5, 0.8/max.chars*10), 23*min(1, 30/dim(current.sorted.Log.Rel.CPM.for.heatmaps)[2])),
                 annotation_col = sample.annotations.clf, annotation_row = current.genes.annotation,
                 annotation_colors = ann.colors,
                 breaks = breaksList)
      }, error = function (err){
        cat(sprintf('Heatmap creating has FAILED\n'))
      })
      
      dev.off()
    }

    if(Pars$Create.heatmaps..Log.rel.CPM..using.ext.samples & !is.null(sorted.Log.Rel.CPM.for.heatmaps)){
      if(!Pars$Create.heatmaps..Log.rel.CPM | Pars$Create.heatmaps..Log.rel.CPM..include.ref.samples.with.suffix != "" | all.sorted.norm.CPMs.for.heatmaps__ncol != ncol(sorted.Log.Rel.CPM.for.heatmaps)){
        x = min(x.src, dim(sorted.Log.Rel.CPM.for.heatmaps)[1])
        #Pars$Create.heatmaps..Log.rel.CPM..log.range
        if (x == dim(sorted.norm.CPMs.for.heatmaps)[1]){
          insert_text = sprintf('all %d genes', x)
        } else {
          insert_text = sprintf('top %d DE genes', x)
        }
        if(Pars$Create.heatmaps..Log.rel.CPM..include.ref.samples.with.suffix != ""){
          Heatmap_Title = sprintf("~ %s, %s. Log2(CPM / avg.CPM). %d samples as ref (**%s)", GLM.model, insert_text, length(ref.samples), Pars$Create.heatmaps..Log.rel.CPM..include.ref.samples.with.suffix)
        } else {
          Heatmap_Title = sprintf("~ %s, %s. Log2(CPM / avg.CPM), relatively all samples", GLM.model, insert_text)
        }

        if(include.GLM.model.name.in.Result.names){
          png.mod(filename = sprintf('%s/%s, Heatmap, Log.rel.CPM (+ext.samples), %s.png', out.dir, Analysis.name, insert_text),
              units="in", width=hm.width, height = hm.height, pointsize=12, res=300)
        } else {
          png.mod(filename = sprintf('%s/Heatmap, Log.rel.CPM (+ext.samples), %s.png', out.dir, insert_text),
              units="in", width=hm.width, height = hm.height, pointsize=12, res=300)
        }

        # print(dim(sorted.Log.Rel.CPM.for.heatmaps))
        # print(dim(intra.norm))
        # print(x)
        # print(data.matrix(sorted.Log.Rel.CPM.for.heatmaps[1:x, ]))
        #stop('')
        # print(intra.norm)
        ####  encorrect x!!!!
        tryCatch(expr = {
          pheatmap(mat = data.matrix(sorted.Log.Rel.CPM.for.heatmaps[1:x, ]), scale = 'none', color = hmcols,
                   cluster_cols = cluster.all.samples, border_color = border_color,
                   main = Heatmap_Title, fontsize = 7, fontsize_row = min(9, 9/x*50) * hm.height / 8,
                   fontsize_col = min(22*min(0.5, 0.8/max.chars*10), 23*min(1, 30/dim(sorted.norm.CPMs.for.heatmaps)[2])),
                   annotation_col = sample.annotations.clf, annotation_row = current.genes.annotation,
                   annotation_colors = ann.colors)
        }, error = function (err){
          cat(sprintf('Heatmap creating has FAILED\n'))
        })
        
        dev.off()
      }
    }

    


  }
  # Add.Completed.step(sprintf("%s.%s.%s.heatmaps",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
  Analysis.Data$Heatmaps.created = TRUE
  cat('\nCreating heatmaps completed.\n')
}




Create.Heatmaps.for.WGCNA = function(Startup.Data, WGCNA.group = NULL, WGCNA.working.dir = NULL,
                                     WGCNA.Data.RDS.file = NULL, forced.parameters = NULL,
                                     limit.to.genes = NULL, out.dir = NULL, max.Genes.counts.list = NULL,
                                     cluster.all.samples = FALSE, predictor.to.sort = NULL, custom.sample.order = NULL,
                                     use.GLM.sorted.samples = FALSE, use.GLM.clustered.samples = FALSE,
                                     custom.CPMs.data.file = NULL, custom.sample.annotations = NULL,
                                     add.grid = FALSE, annotation.logFC.range = 4, annotation.log10P.range = 12,
                                     include.GLM.model.name.in.Result.names = FALSE){
  # print(custom.sample.order)
  if(sum(c(!is.null(WGCNA.group), !is.null(WGCNA.working.dir), !is.null(WGCNA.Data.RDS.file))) != 1) {
    stop('Please specify either WGCNA.group, WGCNA.working.dir or WGCNA.Data.RDS.file arguments')
  }

#  WGCNA.Data$edgeR_d
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(!is.null(WGCNA.group)){
    WGCNA.Analysis.name = Cleanup.Model.Name(WGCNA.group)  
    WGCNA.working.dir = sprintf("%s/WGCNA - %s, results", Pars$results.dir, WGCNA.Analysis.name)
    WGCNA.Data.RDS.file = file.path(WGCNA.working.dir, "WGCNA.data.rds")
    if(!file.exists(WGCNA.Data.RDS.file))
        stop(sprintf('WGCNA data file %s is not found', WGCNA.Data.RDS.file))
    WGCNA.Data = readRDS(file = file.path(WGCNA.working.dir, "WGCNA.data.rds"))
    
  } else if (!is.null(WGCNA.working.dir)){
    if(!file.exists(WGCNA.Data.RDS.file))
        stop(sprintf('WGCNA data file %s is not found', WGCNA.Data.RDS.file))
    WGCNA.Data = readRDS(file = file.path(WGCNA.working.dir, "WGCNA.data.rds"))
    WGCNA.Analysis.name = WGCNA.Data$WGCNA.Analysis.name
    
  } else if (!is.null(WGCNA.Data.RDS.file)){
    if(!file.exists(WGCNA.Data.RDS.file))
        stop(sprintf('WGCNA data file %s is not found', WGCNA.Data.RDS.file))
    WGCNA.Data = readRDS(file = file.path(WGCNA.working.dir, "WGCNA.data.rds"))
    WGCNA.Analysis.name = WGCNA.Data$WGCNA.Analysis.name
    WGCNA.working.dir = dirname(WGCNA.Data.RDS.file)

  }


  # WGCNA.Data$WGCNA.group = WGCNA.group
  # WGCNA.Data$MEs = MEs
  # WGCNA.Data$datExpr.ord = datExpr.ord
  # WGCNA.Data$moduleLabels.ord = moduleLabels.ord
  # WGCNA.Data$moduleColors.ord = moduleColors.ord
  # WGCNA.Data$WGCNA.Analysis.name = WGCNA.Analysis.name


  
  if(!is.null(custom.sample.order) + cluster.all.samples + !is.null(predictor.to.sort) + use.GLM.sorted.samples + use.GLM.clustered.samples > 1 ){
    msg = '\nCreate.Heatmaps.for.WGCNA(...):  Only one of the following parameters should not be provided/turned ON: custom.sample.order, cluster.all.samples, predictor.to.sort, use.GLM.sorted.samples, use.GLM.clustered.samples\n'
    stop(msg)
  }






}





Get.LogFCs.for.Ensembl.IDs.vector <- function(Analysis.Data, Pars, genes.vector, PValue.thr = 1, logCPM.thr = 0, sort.logFCs = TRUE){
  if(is.null(Analysis.Data$Main.test)){
    Main.test = Pars$Main.test
  } else Main.test = Analysis.Data$Main.test

  LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
  PValue.col.name = Get.PValue.col.name(Main.test, colnames(Analysis.Data$ResTable.gene.part))
  # LogFC.col.name = 'LogFC'
  # LogFC.col.name = 'p (QLF test)'
  # cat(sprintf('Using "%s" and "%s" indicators\n', LogFC.col.name, PValue.col.name))
  
  #genes.vector = Genes.per.DB.entry
  
  gene.ID.array = genes.vector[[1]]
  return(
    sapply(genes.vector, function(gene.ID.array){
      gene.ID.array = gene.ID.array[!duplicated(gene.ID.array)]
      gene.ID.array = gene.ID.array[gene.ID.array %in% rownames(Analysis.Data$ResTable.gene.part)]
      sub.table = Analysis.Data$ResTable.gene.part[gene.ID.array,]
      sub.table = sub.table[sub.table$LogCPM >= logCPM.thr,]
      sub.table[,LogFC.col.name][which(sub.table[,PValue.col.name] > PValue.thr)] = 0
      
      if(sort.logFCs)  return(sub.table[,LogFC.col.name][rev(order(sub.table[,LogFC.col.name]))])
      return(sub.table[,LogFC.col.name])
    })
  )  
}


# Get.LogFCs.for.Ensembl.IDs.vector <- function(Analysis.Data, Pars, genes.vector, PValue.thr = 1, logCPM.thr = 0, sort.logFCs = TRUE){
#   
#   LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
#   PValue.col.name = Get.PValue.col.name(Pars)
#   
#   #genes.vector=Genes.per.DB.entry
#   #genes.vector[['R-HSA-977441']]
#   
#   return(
#     sapply(genes.vector, function(gene.ID.array){
#       gene.ID.array = gene.ID.array[!duplicated(gene.ID.array)]
#       if(sum.mod(gene.ID.array %in% rownames(Analysis.Data$ResTable.gene.part)) == 0){
#         return(Analysis.Data$ResTable.gene.part[c(),])
#       }
#       gene.ID.array = gene.ID.array[gene.ID.array %in% rownames(Analysis.Data$ResTable.gene.part)]
#       sub.table = subset(Analysis.Data$ResTable.gene.part, subset = (rownames(Analysis.Data$ResTable.gene.part) %in% gene.ID.array))
#       if(dim(sub.table)[1] == 0)  return(NA)
#       sub.table = subset(sub.table, subset = (sub.table[,'LogCPM'] >= logCPM.thr))
#       if(dim(sub.table)[1] == 0)  return(NA)
#       sub.table[which(sub.table[,PValue.col.name] > PValue.thr), LogFC.col.name] = 0
#       
#       LogFCs = sub.table[,LogFC.col.name]
#       if(sort.logFCs)  return(LogFCs[rev(order(LogFCs))])
#       return(LogFCs)
#     })
#   )  
# }
# 

Get.complete.LogFCs.with.Entrez.IDs <- function(Startup.Data, Pars, Analysis.Data){
  #if(Pars$Species == 'dme')  gene_ID_type = 'flybasecgid_gene' else  gene_ID_type = 'entrezgene_id'
  # gene_ID_type = 'entrezgene_id'
  if(is.null(Analysis.Data$Main.test)){
    Main.test = Pars$Main.test
  } else Main.test = Analysis.Data$Main.test

  if('entrezgene_id' %in% colnames(Startup.Data$General.maRt.table)){
    gene_ID_type = 'entrezgene_id'
  } else if('entrezgene' %in% colnames(Startup.Data$General.maRt.table)){
    gene_ID_type = 'entrezgene'
  } else {
    stop(sprintf('Cannot find Entrez gene coulmn in Startup.Data$General.maRt.table. Available columns: %s', toString(colnames(Startup.Data$General.maRt.table))))
  }

  EntrezGene.IDs = Startup.Data$General.maRt.table[!(Startup.Data$General.maRt.table[,gene_ID_type] %in% NA),]
  rownames(EntrezGene.IDs) = EntrezGene.IDs[,"ensembl_gene_id"]
  
  LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
  PValue.col.name = Get.PValue.col.name(Main.test, colnames(Analysis.Data$ResTable.gene.part))

  Entrez.Stats.table = Analysis.Data$ResTable.gene.part[,c(LogFC.col.name, PValue.col.name, 'Score')]
  Entrez.Stats.table = Entrez.Stats.table[(rownames(Entrez.Stats.table) %in% rownames(EntrezGene.IDs)),]
  Entrez.Stats.table = cbind(Entrez.Stats.table,data.frame(EntrezGene.IDs[rownames(Entrez.Stats.table),c(gene_ID_type)]))
  colnames(Entrez.Stats.table)[4] = gene_ID_type
  
  my_geneList = Entrez.Stats.table[,LogFC.col.name]
  names(my_geneList) = Entrez.Stats.table[,gene_ID_type]
  my_geneList = my_geneList[rev(order(my_geneList))]
  return(my_geneList)
}


Do.end.run = FALSE

Get.Entrez.IDs <- function(Startup.Data, Pars, Analysis.Data,
                           max.DE.genes = NULL, DE.type = 'both', min.abs.LogFC = 0, max.PValue = 1, max.NP.PValue = 1, min.LogCPM.Value = NA,
                           min.Score = 0, disable.filters = FALSE, use.Org.Db__instead.of__biomaRt.table = FALSE, getLogFC = FALSE){

  if(use.Org.Db__instead.of__biomaRt.table && Pars$Species == 'dme'){
    cat('\nUsing org.Dm.eg.db is not compartible with enrichKEGG / enrichPathway / enrichGO functions. biomaRt tables will be used instead of org.Dm.eg.db\n')
    use.Org.Db__instead.of__biomaRt.table = FALSE
  }
  
  if(is.null(Analysis.Data$Main.test)){
    Main.test = Pars$Main.test
  } else Main.test = Analysis.Data$Main.test

  LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
  PValue.col.name = Get.PValue.col.name(Main.test, colnames(Analysis.Data$ResTable.gene.part))
  NP.PValue.col.name = Get.NP.PValue.col.name(colnames(Analysis.Data$ResTable.gene.part))
  LogCPM.col.name = Get.LogCPM.col.name(colnames(Analysis.Data$ResTable.gene.part))

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
    
    #if(Pars$Species == 'dme')  gene_ID_type = 'flybasecgid_gene' else  gene_ID_type = 'entrezgene_id'
    
    if(max.NP.PValue < 1 & !is.na(NP.PValue.col.name)){
      Ensembl.Stats.table = Analysis.Data$ResTable.gene.part[,c(LogFC.col.name, PValue.col.name, NP.PValue.col.name, LogCPM.col.name, 'Score')]
    } else {
      Ensembl.Stats.table = Analysis.Data$ResTable.gene.part[,c(LogFC.col.name, PValue.col.name, LogCPM.col.name, 'Score')]
    }

    Entrez.IDs = mapIds(org.DB, keys=rownames(Ensembl.Stats.table), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
    Entrez.IDs = Entrez.IDs[!(Entrez.IDs %in% NA)]
    #cat(sprintf('\n%d of %d genes unmapped\n',sum.mod(Entrez.IDs %in% NA),length(Entrez.IDs)))
    Ensembl.Stats.table = Ensembl.Stats.table[(rownames(Ensembl.Stats.table) %in% names(Entrez.IDs)),]
    Ensembl.Stats.table = cbind(Ensembl.Stats.table,data.frame(Entrez.IDs[rownames(Ensembl.Stats.table)]))
    # head(Ensembl.Stats.table)

    gene_ID_type = 'entrezgene_id'
    colnames(Ensembl.Stats.table)[ncol(Ensembl.Stats.table)] = gene_ID_type
    
  } else {
    #if(Pars$Species == 'dme')  gene_ID_type = 'flybasecgid_gene' else  gene_ID_type = 'entrezgene_id'
    if('entrezgene_id' %in% colnames(Startup.Data$General.maRt.table)){
      gene_ID_type = 'entrezgene_id'
    } else if('entrezgene' %in% colnames(Startup.Data$General.maRt.table)){
      gene_ID_type = 'entrezgene'
    } else {
      stop(sprintf('Cannot find Entrez gene coulmn in Startup.Data$General.maRt.table. Available columns: %s', toString(colnames(Startup.Data$General.maRt.table))))
    }

    EntrezGene.IDs = Startup.Data$General.maRt.table[!(Startup.Data$General.maRt.table[,gene_ID_type] %in% NA),]
    rownames(EntrezGene.IDs) = EntrezGene.IDs[,"ensembl_gene_id"]
    if(max.NP.PValue < 1 & !is.na(NP.PValue.col.name)){
      Ensembl.Stats.table = Analysis.Data$ResTable.gene.part[,c(LogFC.col.name, PValue.col.name, NP.PValue.col.name, LogCPM.col.name, 'Score')]
    } else {
      Ensembl.Stats.table = Analysis.Data$ResTable.gene.part[,c(LogFC.col.name, PValue.col.name, LogCPM.col.name, 'Score')]
    }
    Ensembl.Stats.table = Ensembl.Stats.table[(rownames(Ensembl.Stats.table) %in% rownames(EntrezGene.IDs)),]
    Ensembl.Stats.table = cbind(Ensembl.Stats.table,data.frame(EntrezGene.IDs[rownames(Ensembl.Stats.table),c(gene_ID_type)]))
    #head(data.frame(EntrezGene.IDs[rownames(Ensembl.Stats.table),c("entrezgene_id")]))
    
    colnames(Ensembl.Stats.table)[ncol(Ensembl.Stats.table)] = gene_ID_type
  }
  

  
  if(!disable.filters){
    if (DE.type == 'de' | DE.type == 'up-down' | DE.type == 'upreg-downreg' | DE.type == 'both'){
      Passed = (abs(Ensembl.Stats.table[,LogFC.col.name]) >= min.abs.LogFC) & 
        (Ensembl.Stats.table[,PValue.col.name] <= max.PValue) & 
        (Ensembl.Stats.table[,"Score"] >= min.Score)
    }
    if (DE.type == 'upreg' | DE.type == 'up'){
      Passed = (Ensembl.Stats.table[,LogFC.col.name] > 0) & 
        (abs(Ensembl.Stats.table[,LogFC.col.name]) >= min.abs.LogFC) & 
        (Ensembl.Stats.table[,PValue.col.name] <= max.PValue) & 
        (Ensembl.Stats.table[,"Score"] >= min.Score)
    }
    if (DE.type == 'downreg' | DE.type == 'down'){
      Passed = (Ensembl.Stats.table[,LogFC.col.name] < 0) & 
        (abs(Ensembl.Stats.table[,LogFC.col.name]) >= min.abs.LogFC) & 
        (Ensembl.Stats.table[,PValue.col.name] <= max.PValue) & 
        (Ensembl.Stats.table[,"Score"] >= min.Score)
    }

    if(max.NP.PValue < 1 & !is.na(NP.PValue.col.name))
      Passed = Passed & (Ensembl.Stats.table[, NP.PValue.col.name] <= max.NP.PValue)

    if(!is.na(min.LogCPM.Value))
      Passed = Passed & (Ensembl.Stats.table[, LogCPM.col.name] >= min.LogCPM.Value)

    if (max.DE.genes >= length(Passed[Passed])){
      Do.end.run <<- TRUE
    }
    
    max.DE.genes = min(length(Passed[Passed])-1,max.DE.genes)
    Ensembl.Stats.table = Ensembl.Stats.table[Passed,]
    min.Score.threshold.effective = max(min.Score, sort(Ensembl.Stats.table[,"Score"],decreasing = T)[max.DE.genes+1])
    Ensembl.Stats.table = Ensembl.Stats.table[Ensembl.Stats.table[,"Score"] > min.Score.threshold.effective,]
  }
  ## if(Pars$Species == 'dme')  Entrez.IDs = sprintf('Dmel_CG%s',Entrez.IDs)
  
  if(getLogFC){
    LogFCs.Ensembl = Ensembl.Stats.table[,LogFC.col.name]
    names(LogFCs.Ensembl) = rownames(Ensembl.Stats.table)
    Entrez.names = Startup.Data$General.maRt.table[rownames(Ensembl.Stats.table), gene_ID_type]
    LogFCs.Ensembl = LogFCs.Ensembl[!duplicated(Entrez.names)]
    Entrez.names = Entrez.names[!duplicated(Entrez.names)]
    names(LogFCs.Ensembl) = Entrez.names
    return(LogFCs.Ensembl)

  } else{
    Entrez.IDs = levels(as.factor(Startup.Data$General.maRt.table[rownames(Ensembl.Stats.table), gene_ID_type]))
    return(Entrez.IDs)
  }
}


Get.Ensembl.IDs <- function(Startup.Data, Pars, Analysis.Data,
                           max.DE.genes = NULL, DE.type = 'both', min.abs.LogFC = 0, max.PValue = 1, max.NP.PValue = 1, min.LogCPM.Value = NA,
                           min.Score = 0, disable.filters = FALSE, getLogFC = FALSE){
  
  if(is.null(Analysis.Data$Main.test)){
    Main.test = Pars$Main.test
  } else Main.test = Analysis.Data$Main.test

  LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
  PValue.col.name = Get.PValue.col.name(Main.test, colnames(Analysis.Data$ResTable.gene.part))
  NP.PValue.col.name = Get.NP.PValue.col.name(colnames(Analysis.Data$ResTable.gene.part))
  LogCPM.col.name = Get.LogCPM.col.name(colnames(Analysis.Data$ResTable.gene.part))

  if(max.NP.PValue < 1 & !is.na(NP.PValue.col.name)){
    Ensembl.Stats.table = Analysis.Data$ResTable.gene.part[,c(LogFC.col.name, PValue.col.name, NP.PValue.col.name, LogCPM.col.name, 'Score')]
  } else {
    Ensembl.Stats.table = Analysis.Data$ResTable.gene.part[,c(LogFC.col.name, PValue.col.name, LogCPM.col.name, 'Score')]
  }
  rownames(Ensembl.Stats.table) = rownames(Analysis.Data$ResTable.gene.part)

  if(!disable.filters){
    if (DE.type == 'de' | DE.type == 'up-down' | DE.type == 'upreg-downreg' | DE.type == 'both'){
      Passed = (abs(Ensembl.Stats.table[,LogFC.col.name]) >= min.abs.LogFC) & 
        (Ensembl.Stats.table[,PValue.col.name] <= max.PValue) & 
        (Ensembl.Stats.table[,"Score"] >= min.Score)
    }
    if (DE.type == 'upreg' | DE.type == 'up'){
      Passed = (Ensembl.Stats.table[,LogFC.col.name] > 0) & 
        (abs(Ensembl.Stats.table[,LogFC.col.name]) >= min.abs.LogFC) & 
        (Ensembl.Stats.table[,PValue.col.name] <= max.PValue) & 
        (Ensembl.Stats.table[,"Score"] >= min.Score)
    }
    if (DE.type == 'downreg' | DE.type == 'down'){
      Passed = (Ensembl.Stats.table[,LogFC.col.name] < 0) & 
        (abs(Ensembl.Stats.table[,LogFC.col.name]) >= min.abs.LogFC) & 
        (Ensembl.Stats.table[,PValue.col.name] <= max.PValue) & 
        (Ensembl.Stats.table[,"Score"] >= min.Score)
    }

    if(max.NP.PValue < 1 & !is.na(NP.PValue.col.name))
      Passed = Passed & (Ensembl.Stats.table[, NP.PValue.col.name] <= max.NP.PValue)

    if(!is.na(min.LogCPM.Value))
      Passed = Passed & (Ensembl.Stats.table[, LogCPM.col.name] >= min.LogCPM.Value)
    

    if (max.DE.genes >= length(Passed[Passed])){
      Do.end.run <<- TRUE
    }
    
    max.DE.genes = min(length(Passed[Passed])-1,max.DE.genes)
    Ensembl.Stats.table = Ensembl.Stats.table[Passed,]
    min.Score.threshold.effective = max(min.Score, sort(Ensembl.Stats.table[,"Score"],decreasing = TRUE)[max.DE.genes+1])
    Ensembl.Stats.table = Ensembl.Stats.table[Ensembl.Stats.table[,"Score"] > min.Score.threshold.effective,]
  }
  
  if(getLogFC){
    res = Ensembl.Stats.table[, LogFC.col.name]
    names(res) = rownames(Ensembl.Stats.table)
    return(res)
  } else return(rownames(Ensembl.Stats.table))
}


Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table = function(GSEA.summary.table, Startup.Data, Pars, gene.col, Ensembl.mode = FALSE){
  if(nrow(GSEA.summary.table)==0) return(GSEA.summary.table)
  # if(!Ensembl.mode){
  #   #if(Pars$Species == 'dme')  gene_ID_type = 'flybasecgid_gene' else  gene_ID_type = 'entrezgene_id'
  #   gene_ID_type = 'entrezgene_id'
  # } else {
  #   gene_ID_type = 'ensembl_gene_id'
  # }
  
  if('entrezgene_id' %in% colnames(Startup.Data$General.maRt.table)){
    gene_ID_type = 'entrezgene_id'
  } else if('entrezgene' %in% colnames(Startup.Data$General.maRt.table)){
    gene_ID_type = 'entrezgene'
  } else {
    stop(sprintf('Cannot find Entrez gene coulmn in Startup.Data$General.maRt.table. Available columns: %s', toString(colnames(Startup.Data$General.maRt.table))))
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

Convert.Gene.IDs.in.values.Vector = function(values.by.gene.IDs, Startup.Data, Pars, Ensembl.mode = FALSE){
  if(length(values.by.gene.IDs)==0) return(c())
  # if(!Ensembl.mode){
  #   #if(Pars$Species == 'dme')  gene_ID_type = 'flybasecgid_gene' else  gene_ID_type = 'entrezgene_id'
  #   gene_ID_type = 'entrezgene_id'
  # } else {
  #   gene_ID_type = 'ensembl_gene_id'
  # }
  
  if('entrezgene_id' %in% colnames(Startup.Data$General.maRt.table)){
    gene_ID_type = 'entrezgene_id'
  } else if('entrezgene' %in% colnames(Startup.Data$General.maRt.table)){
    gene_ID_type = 'entrezgene'
  } else {
    stop(sprintf('Cannot find Entrez gene coulmn in Startup.Data$General.maRt.table. Available columns: %s', toString(colnames(Startup.Data$General.maRt.table))))
  }

  tmp.maRt.table = Startup.Data$General.maRt.table[!(Startup.Data$General.maRt.table[,gene_ID_type] %in% NA),]
  tmp.maRt.table = tmp.maRt.table[!duplicated(tmp.maRt.table[,gene_ID_type]),]
  rownames(tmp.maRt.table) = sapply(tmp.maRt.table[,gene_ID_type],function(y) { as.character(y) })
  
  values.by.gene.IDs = values.by.gene.IDs[names(values.by.gene.IDs) %in% rownames(tmp.maRt.table)]
  gene.Names = tmp.maRt.table[names(values.by.gene.IDs), 'external_gene_name']
  
  values.by.gene.IDs = values.by.gene.IDs[!duplicated(gene.Names)]
  gene.Names = gene.Names[!duplicated(gene.Names)]
  names(values.by.gene.IDs) = gene.Names
  return(values.by.gene.IDs)
}


Convert.Gene.IDs.in.values.Vector__ENSEMBL = function(values.by.gene.IDs, Startup.Data, Pars){
  if(length(values.by.gene.IDs)==0) return(c())
  Ensembl.IDs__to__gene.symbols = Startup.Data$Ensembl.IDs__to__gene.symbols

  values.by.gene.IDs = values.by.gene.IDs[names(values.by.gene.IDs) %in% names(Ensembl.IDs__to__gene.symbols)]
  gene.Names = Ensembl.IDs__to__gene.symbols[names(values.by.gene.IDs)]
  
  values.by.gene.IDs = values.by.gene.IDs[!duplicated(gene.Names)]
  gene.Names = gene.Names[!duplicated(gene.Names)]
  names(values.by.gene.IDs) = gene.Names
  return(values.by.gene.IDs)
}

Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table__ENSEMBL = function(GSEA.summary.table, Startup.Data, Pars, gene.col){
  if(nrow(GSEA.summary.table)==0) return(GSEA.summary.table)
  Ensembl.IDs__to__gene.symbols = Startup.Data$Ensembl.IDs__to__gene.symbols

  for(x in 1:dim(GSEA.summary.table)[1]){
    gene.IDs = strsplit(GSEA.summary.table[x,gene.col],'/',fixed = TRUE)[[1]]
    gene.IDs = gene.IDs[gene.IDs %in% names(Ensembl.IDs__to__gene.symbols)]
    gene.Names = Ensembl.IDs__to__gene.symbols[gene.IDs]
    GSEA.summary.table[x,gene.col] = paste(gene.Names, collapse = '/')
  }
  return(GSEA.summary.table)
}

Get.Gene.Ensembl.IDs__for__GO.term.IDs = function(terms, Startup.Data, forced.parameters = NULL, remove.notfound.terms = TRUE){
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
  
  # terms = c('GO:0042416','GO:0051971','GO:0022008')
  # use Startup.Data$GO.full.descriptions to get to know GO terms IDs

  eval(parse(text = sprintf("GO.ID__to__Entrez = as.list(%s::%sGO2ALLEGS)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  if(remove.notfound.terms)  terms = terms[terms %in% names(GO.ID__to__Entrez)]
  eval(parse(text = sprintf("Entrez__to__Ensembl = as.list(%s::%sENSEMBL)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  #gene.list = sapply(Reactome.ID__to__Entrez[terms],function(x) { unlist(mapIds(org.Mm.eg.db, keys = x, column="ENSEMBL", keytype="ENTREZID", multiVals="list")) })
  return(sapply(GO.ID__to__Entrez[terms], function(x) { unlist(Entrez__to__Ensembl[x]) }))
}


Get.GO.IDs.for.Ensembl.genes = function(genes, Startup.Data, forced.parameters = NULL, add.term.names = FALSE){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters

  
  
  if (Pars$Species %in% names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
    org.DB = Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes[Pars$Species]
  } else {
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s). Get.GO.IDs.for.Ensembl.genes cannot run\n',
                  Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                  paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                  paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
    cat(msg)
    warning(msg)
    org.DB = NULL
    return(vector(mode = 'list'))
  }
  eval(parse(text = sprintf("Ensembl__to__Entrez = as.list(%s::%sENSEMBL2EG)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  eval(parse(text = sprintf("Entrez__to__GO = as.list(%s::%sGO)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  Entrez.gene.names = Ensembl__to__Entrez[genes]
  not.found.genes.count = sum.mod(is.na(names(Entrez.gene.names)))
  cat(sprintf('\n%d of %d genes cannot be assigned with Entrez ID\n', not.found.genes.count, length(genes)))
  #tmp.sel = Entrez.gene.names[which(sapply(Entrez.gene.names, length) == 0)]
  #x = tmp.sel[[2]]
  res = sapply(Entrez.gene.names, function(x){
    DeDup(unlist(sapply(Entrez__to__GO[x], names)))
  })

  if(add.term.names){
    GO.ID__to__GO.name = Startup.Data$GO.full.descriptions[, 'name']
    names(GO.ID__to__GO.name) = rownames(Startup.Data$GO.full.descriptions)
    res = sapply(res, function(x){
      sprintf('%s (%s)', x, GO.ID__to__GO.name[x])
    })
  }
  return(res)
}


Get.KEGG.IDs.for.Ensembl.genes = function(genes, Startup.Data, forced.parameters = NULL, add.pathway.names = FALSE){
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
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s). Get.KEGG.IDs.for.Ensembl.genes will not be run\n',
                  Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                  paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                  paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
    cat(msg)
    warning(msg)
    org.DB = NULL
    return(vector(mode = 'list'))
  }

  # KEGG.ID__to__KEGG.Name = KEGG_DATA$PATHID2NAME

  source(sprintf('%s/kegg.clusterProfiler.R',Pars$suppl.data.dir))
  if(Pars$Species == 'dme'){
    Startup.Data$KEGG_DATA <- prepare_KEGG(species = Pars$Species, keyType = 'ncbi-geneid')
  } else {
    Startup.Data$KEGG_DATA <- prepare_KEGG(species = Pars$Species, keyType = 'kegg')
  }

  eval(parse(text = sprintf("Ensembl__to__Entrez = as.list(%s::%sENSEMBL2EG)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  Entrez.gene.names = Ensembl__to__Entrez[genes]
  not.found.genes.count = sum.mod(is.na(names(Entrez.gene.names)))
  cat(sprintf('\n%d of %d genes cannot be assigned with Entrez ID\n', not.found.genes.count, length(genes)))
  #tmp.sel = Entrez.gene.names[which(sapply(Entrez.gene.names, length) == 2)]
  Entres__to__KEGG = as.list(Startup.Data$KEGG_DATA$EXTID2PATHID)

  res = sapply(Entrez.gene.names, function(x){
    DeDup(unlist(Entres__to__KEGG[x]))
  })
  
  if(add.pathway.names){
    KEGG.ID__to__KEGG.name = Startup.Data$KEGG_DATA$PATHID2NAME
    res = sapply(res, function(x){
      sprintf('%s (%s)', x, KEGG.ID__to__KEGG.name[x])
    })
  }
  return(res)
}


Get.Reactome.IDs.for.Ensembl.genes = function(genes, Startup.Data, forced.parameters = NULL, add.pathway.names = FALSE){
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
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s). Get.Reactome.IDs.for.Ensembl.genes will not be run\n',
                  Startup.Data$DB.data$Taxons.by.KEGG.codes[Pars$Species], Pars$Species,
                  paste(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                  paste(Startup.Data$DB.data$Taxons.by.KEGG.codes[names(Startup.Data$DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
    cat(msg)
    warning(msg)
    org.DB = NULL
    return(vector(mode = 'list'))
  }


  eval(parse(text = sprintf("Ensembl__to__Entrez = as.list(%s::%sENSEMBL2EG)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  Entrez.gene.names = Ensembl__to__Entrez[genes]
  not.found.genes.count = sum.mod(is.na(names(Entrez.gene.names)))
  cat(sprintf('\n%d of %d genes cannot be assigned with Entrez ID\n', not.found.genes.count, length(genes)))
  #tmp.sel = Entrez.gene.names[which(sapply(Entrez.gene.names, length) == 2)]
  Entrez__to__Reactome.ID = as.list(reactomeEXTID2PATHID)

  res = sapply(Entrez.gene.names, function(x){
    DeDup(unlist(Entrez__to__Reactome.ID[x]))
  })
  
  if(add.pathway.names){
    Reactome.ID__to__Reactome.Name = as.list(reactomePATHID2NAME)
    x = res$ENSG00000172936
    res = sapply(res, function(x){
      x.ok = x[x %in% names(Reactome.ID__to__Reactome.Name)]
      x.not.ok = x[!(x %in% names(Reactome.ID__to__Reactome.Name))]
      c(sprintf('%s (%s)', x.ok, Reactome.ID__to__Reactome.Name[x.ok]), x.not.ok)
    })
  }
  return(res)
}


Get.Gene.Ensembl.IDs__for__GO.term.IDs.pack = function(terms.vector, species, DB.data){
  if (species %in% names(DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
    org.DB = DB.data$Bioconductor.Org.DBs.by.KEGG.codes[species]
  } else {
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Available only: %s (%s). Get.Gene.Ensembl.IDs__for__GO.term.IDs.pack cannot run\n',
                  DB.data$Taxons.by.KEGG.codes[species], species,
                  paste(DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                  paste(DB.data$Taxons.by.KEGG.codes[names(DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
    cat(msg)
    warning(msg)
    org.DB = NULL
    return(vector(mode = 'list'))
  }
  
  # terms = c('GO:0042416','GO:0051971','GO:0022008')
  # use Startup.Data$GO.full.descriptions to get to know GO terms IDs

  eval(parse(text = sprintf("GO.ID__to__Entrez = as.list(%s::%sGO2ALLEGS)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  eval(parse(text = sprintf("Entrez__to__Ensembl = as.list(%s::%sENSEMBL)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  not.found.terms.vector = vector(mode = 'list')
  genes.vector = vector(mode = 'list')
  for(group.name in names(terms.vector)){
    terms = terms.vector[[group.name]]
    not.found.terms = terms[!(terms %in% names(GO.ID__to__Entrez))]
    not.found.terms.vector[[group.name]] = not.found.terms
    terms = terms[terms %in% names(GO.ID__to__Entrez)]
    genes.by.terms = sapply(GO.ID__to__Entrez[terms], function(x) {
      genes = unlist(Entrez__to__Ensembl[x])
      names(genes) = NULL
      genes = genes[!is.na(genes)]
    })
    all.genes = unlist(genes.by.terms)
    all.genes = all.genes[!duplicated(all.genes)]
    genes.vector[[group.name]] = all.genes
  }
  res = vector(mode = 'list')
  res[['genes.vector']] = genes.vector
  res[['not.found.terms.vector']] = not.found.terms.vector
  return(res)
  
  #gene.list = sapply(Reactome.ID__to__Entrez[term.list],function(x) { unlist(mapIds(org.Mm.eg.db, keys = x, column="ENSEMBL", keytype="ENTREZID", multiVals="list")) })
  # return(sapply(GO.ID__to__Entrez[term.list], function(x) { unlist(Entrez__to__Ensembl[x]) }))

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
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s). Get.Gene.Ensembl.IDs__for__KEGG.pathway.IDs___not.using.clusterProfiler will not run\n',
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
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s). Get.Gene.Ensembl.IDs__for__KEGG.pathway.IDs will not be run\n',
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

  source(sprintf('%s/kegg.clusterProfiler.R',Pars$suppl.data.dir))
  if(Pars$Species == 'dme'){
    Startup.Data$KEGG_DATA <- prepare_KEGG(species = Pars$Species, keyType = 'ncbi-geneid')
  } else {
    Startup.Data$KEGG_DATA <- prepare_KEGG(species = Pars$Species, keyType = 'kegg')
  }
  if(remove.notfound.pathways)  pathway.list = pathway.list[pathway.list %in% names(Startup.Data$KEGG_DATA$PATHID2EXTID)]

  
  if(!is.null(org.DB)){
    eval(parse(text = sprintf("Entrez__to__Ensembl = as.list(%s::%sENSEMBL)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
    #gene.list = sapply(Reactome.ID__to__Entrez[pathway.list],function(x) { unlist(mapIds(org.Mm.eg.db, keys = x, column="ENSEMBL", keytype="ENTREZID", multiVals="list")) })
    return(sapply(Startup.Data$KEGG_DATA$PATHID2EXTID[pathway.list], function(x) { unlist(Entrez__to__Ensembl[x]) %>% .[!is.na(.)] }))
  } else {
    return(Startup.Data$KEGG_DATA$PATHID2EXTID[pathway.list])
  }
}



Get.Gene.Ensembl.IDs__for__KEGG.pathway.IDs.pack = function(pathway.list, species, DB.data, suppl.data.dir){
  if (!(species %in% rownames(DB.data$KEGG.code.table))){
    msg = sprintf('\nSpecies "%s" is not found among available KEGG species. Ensure that you provided adequate KEGG code (e.g. hsa, mmu, dme, ...)\n', species)
    cat(msg)
    warning(msg)
  }

  if (species %in% names(DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
    org.DB = DB.data$Bioconductor.Org.DBs.by.KEGG.codes[species]
  } else {
    msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s). Get.Gene.Ensembl.IDs__for__KEGG.pathway.IDs will not be run\n',
                  DB.data$Taxons.by.KEGG.codes[species], species,
                  paste(DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
                  paste(DB.data$Taxons.by.KEGG.codes[names(DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
    cat(msg)
    warning(msg)
    org.DB = NULL
    return(vector(mode = 'list'))
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

  source(sprintf('%s/kegg.clusterProfiler.R', suppl.data.dir))
  if(species == 'dme'){
    KEGG_DATA <- prepare_KEGG(species = species, keyType = 'ncbi-geneid')
  } else {
    KEGG_DATA <- prepare_KEGG(species = species, keyType = 'kegg')
  }
  # if(remove.notfound.pathways)  pathway.list = pathway.list[pathway.list %in% names(KEGG_DATA$PATHID2EXTID)]

  
  # if(!is.null(org.DB)){
  #   eval(parse(text = sprintf("Entrez__to__Ensembl = as.list(%s::%sENSEMBL)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
  #   #gene.list = sapply(Reactome.ID__to__Entrez[pathway.list],function(x) { unlist(mapIds(org.Mm.eg.db, keys = x, column="ENSEMBL", keytype="ENTREZID", multiVals="list")) })
  #   return(sapply(KEGG_DATA$PATHID2EXTID[pathway.list], function(x) { unlist(Entrez__to__Ensembl[x]) %>% .[!is.na(.)] }))
  # } else {
  #   return(KEGG_DATA$PATHID2EXTID[pathway.list])
  # }

  not.found.pathway.list = vector(mode = 'list')
  genes.vector = vector(mode = 'list')
  for(group.name in names(pathway.list)){
    pathways = pathway.list[[group.name]]
    not.found.pathways = pathways[!(pathways %in% names(KEGG_DATA$PATHID2EXTID))]
    not.found.pathway.list[[group.name]] = not.found.pathways
    pathways = pathways[pathways %in% names(KEGG_DATA$PATHID2EXTID)]
    genes.by.pathways = sapply(KEGG_DATA$PATHID2EXTID[pathways], function(x) {
      genes = unlist(Entrez__to__Ensembl[x])
      names(genes) = NULL
      genes = genes[!is.na(genes)]
    })
    all.genes = unlist(genes.by.pathways)
    all.genes = all.genes[!duplicated(all.genes)]
    genes.vector[[group.name]] = all.genes
  }
  res = vector(mode = 'list')
  res[['genes.vector']] = genes.vector
  res[['not.found.pathway.list']] = not.found.pathway.list
  return(res)

}




# Get.Gene.Ensembl.IDs__for__KEGG.pathway.IDs.pack = function(pathway.list, DB.data, KEGG_DATA, species, remove.notfound.pathways = TRUE){
#   if (!(species %in% rownames(DB.data$KEGG.code.table))){
#     msg = sprintf('\nSpecies "%s" is not found among available KEGG species. Ensure that you provided adequate KEGG code (e.g. hsa, mmu, dme, ...)\n', species)
#     cat(msg)
#     warning(msg)
#   }

#   if (species %in% names(DB.data$Bioconductor.Org.DBs.by.KEGG.codes)){
#     org.DB = DB.data$Bioconductor.Org.DBs.by.KEGG.codes[species]
#   } else {
#     msg = sprintf('\nSpecies "%s" (%s) is not found among available Bioconductor org.DB. Allowed only: %s (%s). Gene ID conversion will not be performed\n',
#                   DB.data$Taxons.by.KEGG.codes[species], species,
#                   paste(DB.data$Bioconductor.Org.DBs.by.KEGG.codes, collapse = ', '),
#                   paste(DB.data$Taxons.by.KEGG.codes[names(DB.data$Bioconductor.Org.DBs.by.KEGG.codes)], collapse = ', '))
#     cat(msg)
#     warning(msg)
#     org.DB = NULL
#     #return(pathway.list)
#   }
  
#   # Startup.Data$KEGG_DATA = KEGG_DATA
#   # pathway.list = c('mmu04142',
#   #   'mmu00511',
#   #   'mmu04915',
#   #   'mmu04960',
#   #   'mmu04150',
#   #   'mmu04213',
#   #   'mmu05221',
#   #   'mmu04211',
#   #   'left')
  

#   # KEGG.ID__to__KEGG.Name = KEGG_DATA$PATHID2NAME

#   if(remove.notfound.pathways)  pathway.list = pathway.list[pathway.list %in% names(KEGG_DATA$PATHID2EXTID)]

  
#   if(!is.null(org.DB)){
#     eval(parse(text = sprintf("Entrez__to__Ensembl = as.list(%s::%sENSEMBL)",org.DB,substr(org.DB,1,nchar(org.DB)-3))))
#     #gene.list = sapply(Reactome.ID__to__Entrez[pathway.list],function(x) { unlist(mapIds(org.Mm.eg.db, keys = x, column="ENSEMBL", keytype="ENTREZID", multiVals="list")) })
#     return(sapply(KEGG_DATA$PATHID2EXTID[pathway.list], function(x) { unlist(Entrez__to__Ensembl[x]) %>% .[!is.na(.)] }))
#   } else {
#     return(KEGG_DATA$PATHID2EXTID[pathway.list])
#   }



#   not.found.terms.vector = vector(mode = 'list')
#   genes.vector = vector(mode = 'list')
#   for(group.name in names(terms.vector)){
#     terms = terms.vector[[group.name]]
#     not.found.terms = terms[!(terms %in% names(GO.ID__to__Entrez))]
#     not.found.terms.vector[[group.name]] = not.found.terms
#     terms = terms[terms %in% names(GO.ID__to__Entrez)]
#     genes.by.terms = sapply(GO.ID__to__Entrez[terms], function(x) {
#       genes = unlist(Entrez__to__Ensembl[x])
#       names(genes) = NULL
#       genes = genes[!is.na(genes)]
#     })
#     all.genes = unlist(genes.by.terms)
#     all.genes = all.genes[!duplicated(all.genes)]
#     genes.vector[[group.name]] = all.genes
#   }
#   res = vector(mode = 'list')
#   res[['genes.vector']] = genes.vector
#   res[['not.found.terms.vector']] = not.found.terms.vector
#   return(res)



# }




# Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table.Ensembl = function(GSEA.summary.table, Startup.Data, Pars, gene.col){
#   if(Pars$Species == 'dme')  gene_ID_type = 'flybasecgid_gene' else  gene_ID_type = 'entrezgene_id'
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
#   save(KEGG.classic.Enrichment__results__by.DE.type__by.max.DE.genes, file='KEGG classic enrichment.rds')
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
                              forced.parameters = NULL, additional.plots = TRUE,
                              make.simplify = TRUE, simplify.cutoff = 0.55, make.simplify.plots = TRUE,
                              include.GLM.model.name.in.Result.names = FALSE, old.compartibility = FALSE){
  
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
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
  
  # if (bypass.if.completed & 
  #     Read.Completed.steps.status(sprintf("%s.%s.%s.%s.classic.Enrichment", database.plus, Analysis.Data$GLM.model, Startup.Data$Startup.hashmd5, Analysis.Data$step.hashmd5), Startup.Data$Completed.steps.file) &
  #     Read.Completed.steps.status(sprintf("%s.%s.%s.pathways.visualization", Analysis.Data$GLM.model, Startup.Data$Startup.hashmd5, Analysis.Data$step.hashmd5), Startup.Data$Completed.steps.file)){
  #   message(sprintf('\n %s.classic.Enrichment pathway enrichment and visualization for GLM model "%s" was previously completed. Bypassing...', database.plus, Analysis.Data$GLM.model))
  #   return(invisible(Analysis.Data))
  # }
  
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
  eval(parse(text = sprintf('gene.max.NP.PValue.threshold = Pars$%s.classic.Enrich__with.clusterProfiler___gene.max.NP.PValue.threshold', database)))
  eval(parse(text = sprintf('gene.min.LogCPM.threshold = Pars$%s.classic.Enrich__with.clusterProfiler___gene.min.LogCPM.threshold', database)))  
  eval(parse(text = sprintf('gene.min.Score.threshold = Pars$%s.classic.Enrich__with.clusterProfiler___gene.min.Score.threshold', database)))
  
  eval(parse(text = sprintf('Render.Summary.Plot = Pars$%s.classic.Enrich__with.clusterProfiler___Render.Summary.Plot', database)))

  eval(parse(text = sprintf('EEP___limit.PValues.to = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___limit.PValues.to', database)))
  eval(parse(text = sprintf('EEP___Max.PValue.threshold = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___Max.PValue.threshold__classic.test', database)))
  eval(parse(text = sprintf('EEP___Max.NP.PValue.threshold = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___Max.NP.PValue.threshold__classic.test', database)))
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
        
        my_geneList..LogFC = Get.Entrez.IDs(Startup.Data, Pars, Analysis.Data,
                                    max.DE.genes = max.DE.genes,  DE.type = DE.mode,
                                    min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                    max.PValue = gene.max.PValue.threshold,
                                    max.NP.PValue = gene.max.NP.PValue.threshold,
                                    min.LogCPM.Value = gene.min.LogCPM.threshold,
                                    min.Score = gene.min.Score.threshold, getLogFC = TRUE)
        if(is.null(my_geneList..LogFC)) return(invisible(Analysis.Data))
        my_geneList = names(my_geneList..LogFC)
        cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
        #### fun = "gse.res = gsePathway(my_geneList, organism = org, pAdjustMethod = 'BH', pvalueCutoff = FDR.Cutoff)"
        enrich.res = enrichPathway(my_geneList, organism=org, pvalueCutoff=1.00, pAdjustMethod="BH", qvalueCutoff=1)
      
      } else if (database == 'KEGG'){
        # stop(my_geneList)
        my_geneList..LogFC = Get.Entrez.IDs(Startup.Data, Pars, Analysis.Data,
                                     max.DE.genes = max.DE.genes,  DE.type = DE.mode,
                                     min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                     max.PValue = gene.max.PValue.threshold,
                                     max.NP.PValue = gene.max.NP.PValue.threshold,
                                     min.LogCPM.Value = gene.min.LogCPM.threshold,
                                     min.Score = gene.min.Score.threshold, getLogFC = TRUE)
        if(is.null(my_geneList..LogFC)) return(invisible(Analysis.Data))
        my_geneList = names(my_geneList..LogFC)
        cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
        #### fun = "gse.res = gseKEGG(my_geneList, organism = Pars$Species, pAdjustMethod = 'BH', pvalueCutoff = FDR.Cutoff)"
        if(Pars$Species == 'dme'){
          enrich.res = enrichKEGG(my_geneList, organism=Pars$Species, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0, keyType = 'ncbi-geneid')
        } else {
          enrich.res = enrichKEGG(my_geneList, organism=Pars$Species, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)
        }
        # print(my_geneList)
        # stop()
      
      } else if(database == 'Disease_Ont'){
        my_geneList..LogFC = Get.Entrez.IDs(Startup.Data, Pars, Analysis.Data,
                                     max.DE.genes = max.DE.genes,  DE.type = DE.mode,
                                     min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                     max.PValue = gene.max.PValue.threshold,
                                     max.NP.PValue = gene.max.NP.PValue.threshold,
                                     min.LogCPM.Value = gene.min.LogCPM.threshold,
                                     min.Score = gene.min.Score.threshold, getLogFC = TRUE)
        if(is.null(my_geneList..LogFC)) return(invisible(Analysis.Data))
        my_geneList = names(my_geneList..LogFC)
        cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
        # if(Pars$Species == 'dme'){
        #   enrich.res = enrichDO(my_geneList, organism=Pars$Species, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0, keyType = 'ncbi-geneid')
        # } else {
        enrich.res = enrichDO(my_geneList, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)
        # }

      } else if(database == 'NCG'){
        my_geneList..LogFC = Get.Entrez.IDs(Startup.Data, Pars, Analysis.Data,
                                     max.DE.genes = max.DE.genes,  DE.type = DE.mode,
                                     min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                     max.PValue = gene.max.PValue.threshold,
                                     max.NP.PValue = gene.max.NP.PValue.threshold,
                                     min.LogCPM.Value = gene.min.LogCPM.threshold,
                                     min.Score = gene.min.Score.threshold, getLogFC = TRUE)
        if(is.null(my_geneList..LogFC)) return(invisible(Analysis.Data))
        my_geneList = names(my_geneList..LogFC)
        cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
        # if(Pars$Species == 'dme'){
        #   enrich.res = enrichDO(my_geneList, organism=Pars$Species, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0, keyType = 'ncbi-geneid')
        # } else {
        enrich.res = enrichNCG(my_geneList, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)
        # }

      } else if(database == 'DisGeNET'){
        my_geneList..LogFC = Get.Entrez.IDs(Startup.Data, Pars, Analysis.Data,
                                     max.DE.genes = max.DE.genes,  DE.type = DE.mode,
                                     min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                     max.PValue = gene.max.PValue.threshold,
                                     max.NP.PValue = gene.max.NP.PValue.threshold,
                                     min.LogCPM.Value = gene.min.LogCPM.threshold,
                                     min.Score = gene.min.Score.threshold, getLogFC = TRUE)
        if(is.null(my_geneList..LogFC)) return(invisible(Analysis.Data))
        my_geneList = names(my_geneList..LogFC)
        cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
        # if(Pars$Species == 'dme'){
        #   enrich.res = enrichDO(my_geneList, organism=Pars$Species, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0, keyType = 'ncbi-geneid')
        # } else {
        enrich.res = enrichDGN(my_geneList, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)
        # }

      } else if(database == 'WikiPathways'){
        my_geneList..LogFC = Get.Entrez.IDs(Startup.Data, Pars, Analysis.Data,
                                     max.DE.genes = max.DE.genes,  DE.type = DE.mode,
                                     min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                     max.PValue = gene.max.PValue.threshold,
                                     max.NP.PValue = gene.max.NP.PValue.threshold,
                                     min.LogCPM.Value = gene.min.LogCPM.threshold,
                                     min.Score = gene.min.Score.threshold, getLogFC = TRUE)
        if(is.null(my_geneList..LogFC)) return(invisible(Analysis.Data))
        my_geneList = names(my_geneList..LogFC)
        cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
        # if(Pars$Species == 'dme'){
        #   enrich.res = enrichDO(my_geneList, organism=Pars$Species, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0, keyType = 'ncbi-geneid')
        # } else {
        # enrich.res = enrichDGN(my_geneList, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)
        # }
        #wpgmtfile <- system.file(Pars$WikiPathways.GMT.file, package="clusterProfiler")
        if(is.null(Pars$WikiPathways.GMT.file))  return(invisible(Analysis.Data))
        if(!file.exists(Pars$WikiPathways.GMT.file)){
          msg = sprintf('WikiPathways GMT file "%s" does not exist. Please download it from http://data.wikipathways.org/current/gmt/\n',
            Pars$WikiPathways.GMT.file)
          cat(msg)
          warning(msg)
          return(invisible(Analysis.Data))
        }
        wp2gene <- read.gmt(Pars$WikiPathways.GMT.file)
        wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
        wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
        wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

        enrich.res <- enricher(my_geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)

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
          my_geneList..LogFC = Get.Entrez.IDs(Startup.Data, Pars, Analysis.Data,
                                       max.DE.genes = max.DE.genes,  DE.type = DE.mode,
                                       min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                       max.PValue = gene.max.PValue.threshold,
                                       max.NP.PValue = gene.max.NP.PValue.threshold,
                                       min.LogCPM.Value = gene.min.LogCPM.threshold,
                                       min.Score = gene.min.Score.threshold, getLogFC = TRUE)
          if(is.null(my_geneList..LogFC)) return(invisible(Analysis.Data))
          my_geneList = names(my_geneList..LogFC)
          cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
          #### fun = "gse.res = gseGO(geneList = my_geneList, keyType = 'ENTREZID', OrgDb = org.DB, ont = GO.type, pAdjustMethod='BH', pvalueCutoff = FDR.Cutoff)"
          enrich.res = enrichGO(geneList = my_geneList, OrgDb = org.DB, keyType = 'ENTREZID', ont = GO.type, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)
        
        } else {
          my_geneList..LogFC = Get.Ensembl.IDs(Startup.Data, Pars, Analysis.Data,
                                      max.DE.genes = max.DE.genes, DE.type = DE.mode,
                                      min.abs.LogFC = gene.min.abs.LogFC.threshold,
                                      max.PValue = gene.max.PValue.threshold,
                                      max.NP.PValue = gene.max.NP.PValue.threshold,
                                      min.LogCPM.Value = gene.min.LogCPM.threshold,
                                      min.Score = gene.min.Score.threshold, getLogFC = TRUE)
          if(is.null(my_geneList..LogFC)) return(invisible(Analysis.Data))
          my_geneList = names(my_geneList..LogFC)
          cat(sprintf('\n"~ %s":   Performing %s classic enrichment for %d %s genes', Analysis.Data$GLM.model, printed.name, length(my_geneList), DE.mode))
          #### fun = "gse.res = gseGO(geneList = my_geneList, keyType = 'ENSEMBL', OrgDb = org.DB, ont = GO.type, pAdjustMethod='BH', pvalueCutoff = FDR.Cutoff)"
          if(old.compartibility){
            enrich.res = enrichGO(gene = my_geneList, OrgDb = org.DB, keytype = 'ENSEMBL', ont = GO.type, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)
          } else {
            enrich.res = enrichGO(gene = my_geneList, OrgDb = org.DB, keyType = 'ENSEMBL', ont = GO.type, pvalueCutoff=1.0, pAdjustMethod="BH", qvalueCutoff=1.0)
          }
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
            current.p = max(current.p, EEP___limit.PValues.to)
            current.score = 0
            if (current.p > EEP___Max.PValue.threshold){
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
            current.p = max(current.p, EEP___limit.PValues.to)
            current.score = 0
            if (current.p > EEP___Max.PValue.threshold){
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
      
      eval(parse(text = sprintf('pvalue.lim = Pars$%s.classic.Enrich__with.clusterProfiler___term.max.PValue.threshold', database)))
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
        my_geneList..LogFC..native.gene.names = Convert.Gene.IDs.in.values.Vector__ENSEMBL(my_geneList..LogFC, Startup.Data, Pars)
      } else {
        enrich.res.converted@result = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table(enrich.res@result, Startup.Data, Pars, gene.col = 'geneID')
        my_geneList..LogFC..native.gene.names = Convert.Gene.IDs.in.values.Vector(my_geneList..LogFC, Startup.Data, Pars)
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
          if(Pars$classic.Enrich__with.clusterProfiler___Render.CNE.plots){
            png.mod(filename = sprintf("%s[cneplot] %s classic enrich. for top-%d %s genes.png", dirs.list[t], printed.name, max.DE.genes, DE.mode),
                units="in", width=12, height=12, pointsize=10, res=300)
            print({
              cnetplot(res.list[[t]], showCategory=20, foldChange = my_geneList..LogFC..native.gene.names)
            })
            dev.off()
          }
          
          if(Pars$classic.Enrich__with.clusterProfiler___Render.Dot.plots){
            png.mod(filename = sprintf("%s[dotplot] %s classic enrich. for top-%d %s genes.png", dirs.list[t], printed.name, max.DE.genes, DE.mode),
                units="in", width=12, height=7, pointsize=8, res=300)
            print({
              dotplot(res.list[[t]], showCategory=30) + ggplot2::xlab(sprintf("Gene ratio; %s, top-%d %s genes", printed.name, max.DE.genes, DE.mode)) #+ scale_colour_gradient(limits=c(0, 0.05), high = 'blue', low="red")
            })
            dev.off()
          }
          
          if(Pars$classic.Enrich__with.clusterProfiler___Render.Heat.plots){

            png.mod(filename = sprintf("%s[heatplot] %s classic enrich. for top-%d %s genes.png", dirs.list[t], printed.name, max.DE.genes, DE.mode),
                units="in", width=14, height=7, pointsize=10, res=300)
            print({
              heatplot(res.list[[t]], showCategory=20, foldChange = my_geneList..LogFC..native.gene.names)
            })
            dev.off()
          }
          
          ### optionally, one can use:
          ## p = dotplot(enrich.res.converted,showCategory=30)
          ## ggplot2::ggsave(file="plot.png", plot=p)
          
          #par(mar=c(1,1,1,1))
          if(Pars$classic.Enrich__with.clusterProfiler___Render.Bar.plots){
            png.mod(filename = sprintf("%s[barplot] %s classic enrich. for top-%d %s genes.png", dirs.list[t], printed.name, max.DE.genes, DE.mode),
                units="in", width=12, height=7, pointsize=8, res=300)
            print({
              barplot(res.list[[t]], showCategory=30) + ggplot2::ylab(sprintf("Gene count (per %s) among top-%d %s ones", substr(printed.name, 1, nchar(printed.name)-1), max.DE.genes, DE.mode))
            })
            dev.off()
          }
        }
      }

      # head(as.data.frame(enrich.res.converted))
      # if(grepl('GO',database) & !convert.ENSEMBL.ids.to.ENTREZID){
      #   classic.Enrich.summary.table = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table__ENSEMBL(as.data.frame(enrich.res), Startup.Data, Pars, gene.col = 'geneID')
      # } else {
      #   classic.Enrich.summary.table = Convert.Gene.IDs.in.CusterProfiler.Enrichment.results.table(as.data.frame(enrich.res), Startup.Data, Pars, gene.col = 'geneID')
      # }
      
      if(include.GLM.model.name.in.Result.names){
        file.name = Verify.path(sprintf("%s, %s classic enrichment for top-%d %s genes.tsv", Analysis.Data$Analysis.name, printed.name, max.DE.genes, DE.mode))
      } else {
        file.name = Verify.path(sprintf("%s classic enrichment for top-%d %s genes.tsv", printed.name, max.DE.genes, DE.mode))
      }

      classic.Enrichment.output.file.names = c(classic.Enrichment.output.file.names,file.name)
      write.table.mod(x = as.data.frame(enrich.res.converted), file = file.name, sep='\t')
      
      if(make.simplify){
        if(include.GLM.model.name.in.Result.names){
          file.name = Verify.path(sprintf("%s, %s classic enrichment for top-%d %s genes (RR).tsv", Analysis.Data$Analysis.name, printed.name, max.DE.genes, DE.mode))
        } else {
          file.name = Verify.path(sprintf("%s classic enrichment for top-%d %s genes (RR).tsv", printed.name, max.DE.genes, DE.mode))
        }
        simplified.classic.Enrichment.output.file.names = c(simplified.classic.Enrichment.output.file.names, file.name)
        write.table.mod(x = as.data.frame(simplified.enrich.res.converted), file = file.name, sep='\t')
      }
      
      # save(enrich.res.complete,sprintf("%s, %s classic enrichment for top-%d %s genes.rds", Analysis.Data$Analysis.name, printed.name, max.DE.genes, DE.mode))
      
      Enriched.Pathways.Classic = c(Enriched.Pathways.Classic, rownames(as.data.frame(enrich.res))) # This is a cumulative list of all found enriched pathways
      
    }
  }

  if(Render.Summary.Plot){
    ck = tryCatch(expr = {
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

      } else if (database == 'Disease_Ont') {
        ## enrichKEGG(my_geneList, organism=Pars$Species, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=1)
        # if (Pars$Species == 'dme'){
        #   ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
        #                        fun = "enrichKEGG", organism=Pars$Species,  pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1,keyType = 'ncbi-geneid')
        # } else {
          ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
                               fun = "enrichDO", pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1)
        # }

      } else if (database == 'NCG') {
        ## enrichKEGG(my_geneList, organism=Pars$Species, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=1)
        # if (Pars$Species == 'dme'){
        #   ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
        #                        fun = "enrichKEGG", organism=Pars$Species,  pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1,keyType = 'ncbi-geneid')
        # } else {
          ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
                               fun = "enrichNCG", pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1)
        # }

      } else if (database == 'DisGeNET') {
        ## enrichKEGG(my_geneList, organism=Pars$Species, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=1)
        # if (Pars$Species == 'dme'){
        #   ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
        #                        fun = "enrichKEGG", organism=Pars$Species,  pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1,keyType = 'ncbi-geneid')
        # } else {
          ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
                               fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1)
        # }

      } else if (database == 'WikiPathways') {
        ## enrichKEGG(my_geneList, organism=Pars$Species, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=1)
        # if (Pars$Species == 'dme'){
        #   ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
        #                        fun = "enrichKEGG", organism=Pars$Species,  pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1,keyType = 'ncbi-geneid')
        # } else {
          ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
                               fun = "enricher", TERM2GENE = wpid2gene, TERM2NAME = wpid2name,
                               pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1)
        # }

      } else if (grepl('GO', database)){
        if(old.compartibility){
          if(convert.ENSEMBL.ids.to.ENTREZID){
            ## enrichGO(geneList = my_geneList, OrgDb = org.DB, keytype = 'ENTREZID', ont = GO.type, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=1)
            ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
                               fun = "enrichGO", OrgDb = org.DB, keytype = 'ENTREZID', ont = GO.type, pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1)
          } else {
            ## enrichGO(geneList = my_geneList, OrgDb = org.DB, keytype = 'ENSEMBL', ont = GO.type, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=1)
            ck <- compareCluster(geneCluster = Gene.Lists__by.DE.type.and.max.DE.genes,
                                 fun = "enrichGO", OrgDb = org.DB, keytype = 'ENSEMBL', ont = GO.type, pvalueCutoff=0.05, pAdjustMethod="none", qvalueCutoff=1)
          }
        } else {

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
      }
      ck
    }, error = function (err){
      warning('Creating summary plot has failed')
      return(NULL)
    })
        
    if(!is.null(ck)){
      ck.mod = ck
      ck.mod@compareClusterResult = ck.mod@compareClusterResult[!duplicated(ck.mod@compareClusterResult$ID),]
      ids.of.duplicated.descriptions = ck.mod@compareClusterResult$ID[duplicated(tolower(ck.mod@compareClusterResult$Description))]
      ck@compareClusterResult = ck@compareClusterResult[!(ck@compareClusterResult$ID %in% ids.of.duplicated.descriptions),]
      
      max.term.name.length = 55
      ck@compareClusterResult$Description = sapply(ck@compareClusterResult$Description, function(x){ if (nchar(x) > max.term.name.length) sprintf("%s...", substr(x, 1, max.term.name.length)) else x })
    
      ### ck@compareClusterResult$Description
      
      font.size = min(14, max(6, 1500/dim(ck@compareClusterResult)[1]))
      font.size.for.min = min(14, max(4, 1000/dim(ck@compareClusterResult)[1]))
    
      if(include.GLM.model.name.in.Result.names){
        png.filename = sprintf("[MS dotplot] %s, Enriched %s multi.png",Analysis.Data$Analysis.name, printed.name)
      } else {
        png.filename = sprintf("[MS dotplot] Enriched %s multi.png", printed.name)
      }
      png.mod(filename = png.filename,
          units="in", width=14, height=8, pointsize=8, res=300)
      print({
        suppressMessages(dotplot(ck, showCategory = 9, color = 'pvalue') + theme_dose(font.size = font.size) + scale_colour_gradient(limits=c(0, 0.05), high = 'blue', low="red"))
      })
      dev.off()

      file.copy(png.filename, sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, png.filename))
      
      if(include.GLM.model.name.in.Result.names){
        png.filename = sprintf("[MS dotplot] %s, Enriched %s multi-set, max. terms count.png",Analysis.Data$Analysis.name, printed.name)
      } else {
        png.filename = sprintf("[MS dotplot] Enriched %s multi-set, max. terms count.png", printed.name)
      }

      png.mod(filename = png.filename,
          units="in", width=14, height=8, pointsize=8, res=300)
      print({
        suppressMessages(dotplot(ck, showCategory = 70, color = 'pvalue') + theme_dose(font.size = font.size.for.min) + scale_colour_gradient(limits=c(0, 0.05), high = 'blue', low="red"))
      })
      dev.off()
    }
  }  
  eval(parse(text = sprintf('%s.classic.Enrichment__results__by.DE.type__by.max.DE.genes = Enrichment__results__by.DE.type__by.max.DE.genes', database.plus)))
  if(include.GLM.model.name.in.Result.names){
    eval(parse(text = sprintf('saveRDS(%s.classic.Enrichment__results__by.DE.type__by.max.DE.genes, file="%s, %s classic enrichment.rds")', database.plus, Analysis.Data$Analysis.name, database.plus)))
  } else {
    eval(parse(text = sprintf('saveRDS(%s.classic.Enrichment__results__by.DE.type__by.max.DE.genes, file="%s classic enrichment.rds")', database.plus, database.plus)))
  }
  
  
  if (Pars$Create.Excel.results){
    CL = sprintf('%s "%s/Classic.enrichment.results.to.Excel.py" %s', Pars$python.bin, Pars$suppl.data.dir, database)
    # CL = gsub('/','\\',CL,fixed = TRUE)
    if(include.GLM.model.name.in.Result.names){
      res.file.name = sprintf('%s, %s classic enrichment.xlsx', Analysis.Data$Analysis.name, printed.name)
    } else {
      res.file.name = sprintf('%s classic enrichment.xlsx', printed.name)
    }
    CL = paste(c(CL, sprintf('"%s"', Verify.path(res.file.name)), sprintf('"%s"', classic.Enrichment.output.file.names)),collapse = ' ')
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    if (!file.exists(Verify.path(res.file.name))){
      message(sprintf('Excel worksheet generation (%s classic enrichment results) for "~ %s" FAILED', printed.name, Analysis.Data$Analysis.name))
    } else {
      file.remove(classic.Enrichment.output.file.names)
      file.copy(Verify.path(res.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, res.file.name)), overwrite = TRUE)
    }
  }
  if (Pars$Create.Excel.results && make.simplify){
    CL = sprintf('%s "%s/Classic.enrichment.results.to.Excel.py" %s', Pars$python.bin, Pars$suppl.data.dir, database)
    # CL = gsub('/','\\',CL,fixed = TRUE)
    if(include.GLM.model.name.in.Result.names){
      res.file.name = sprintf('%s, %s classic enrichment (RR).xlsx', Analysis.Data$Analysis.name, printed.name)
    } else {
      res.file.name = sprintf('%s classic enrichment (RR).xlsx', printed.name)
    }
    CL = paste(c(CL, sprintf('"%s"',Verify.path(res.file.name)), sprintf('"%s"', simplified.classic.Enrichment.output.file.names)),collapse = ' ')
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    if (!file.exists(Verify.path(res.file.name))){
      message(sprintf('Excel worksheet generation (%s classic enrichment results, RR) for "~ %s" FAILED', printed.name, Analysis.Data$Analysis.name))
    } else {
      file.remove(simplified.classic.Enrichment.output.file.names)
      file.copy(Verify.path(res.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, res.file.name)), overwrite = TRUE)
    }
  }
  
  all.DB.entries = c(names(enrichment.scores_by.term.id.UP),names(enrichment.scores_by.term.id.DOWN))
  all.DB.entries = all.DB.entries[!duplicated(all.DB.entries)]
  if(Create.DE.plots.from.clusterProfile.res & length(all.DB.entries) > 0){
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
    
    if(include.GLM.model.name.in.Result.names){
      write.table.mod(x=Enrich.info.split, file=sprintf('%s, %s classic Enrich. stats.tsv', Analysis.Data$Analysis.name, database.plus), sep='\t', quote = FALSE, na = "")
      saveRDS(Enrich.info.split, file = sprintf('%s, %s classic Enrich. stats.rds', Analysis.Data$Analysis.name, database.plus))
    } else {
      write.table.mod(x=Enrich.info.split, file=sprintf('%s classic Enrich. stats.tsv', database.plus), sep='\t', quote = FALSE, na = "")
      saveRDS(Enrich.info.split, file = sprintf('%s classic Enrich. stats.rds', database.plus))
    }
    
  }

  if(include.GLM.model.name.in.Result.names){
    saveRDS(Enriched.Pathways.Classic, file = sprintf('%s, %s enriched pathways (classic).rds', Analysis.Data$Analysis.name, database.plus))
  } else {
    saveRDS(Enriched.Pathways.Classic, file = sprintf('%s enriched pathways (classic).rds', database.plus))
  }
  
  # Add.Completed.step(sprintf("%s.%s.%s.%s.classic.Enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5,database.plus),Startup.Data$Completed.steps.file)
  eval(parse(text = sprintf("Analysis.Data$%s.classic.Enriched.Pathways = levels(factor(Enriched.Pathways.Classic))", database.plus)))
  eval(parse(text = sprintf('Analysis.Data$%s.classic.Enrichment.performed = TRUE', database.plus)))
  
  
  setwd(Analysis.Data$current.GLM.results.dir)
  cat(sprintf('\n%s enrichment (classic) completed.\n', printed.name))
  return(invisible(Analysis.Data))
}



Reactome.classic.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                       forced.parameters = NULL, additional.plots = TRUE,
                                       include.GLM.model.name.in.Result.names = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = classic.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'Reactome',
                                    printed.name = 'Reactome Pathways',
                                    forced.parameters = forced.parameters, additional.plots = additional.plots, make.simplify = FALSE,
                                    include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
  
  return(invisible(Analysis.Data))
}

KEGG.classic.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                   forced.parameters = NULL, additional.plots = TRUE,
                                   include.GLM.model.name.in.Result.names = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = classic.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'KEGG',
                                    printed.name = 'KEGG Pathways',
                                    forced.parameters = forced.parameters, additional.plots = additional.plots, make.simplify = FALSE,
                                    include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
  
  return(invisible(Analysis.Data))
}

Disease_Ont.classic.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                   forced.parameters = NULL, additional.plots = TRUE,
                                   include.GLM.model.name.in.Result.names = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(!(Pars$Species %in% c('hsa', 'mmu', 'rno'))){
    msg = sprintf('Disease_Ont.classic.Enrichment works properly only for human and (maybe!) for mammals\n')
    cat(msg)
    warning(msg)
    return()
  }

  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = classic.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'Disease_Ont',
                                    printed.name = 'Disease Ontology DB',
                                    forced.parameters = forced.parameters, additional.plots = additional.plots, make.simplify = FALSE,
                                    include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
  
  return(invisible(Analysis.Data))
}

NCG.classic.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                   forced.parameters = NULL, additional.plots = TRUE,
                                   include.GLM.model.name.in.Result.names = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(!(Pars$Species %in% c('hsa', 'mmu', 'rno'))){
    msg = sprintf('NCG.classic.Enrichment works properly only for human and (maybe!) for mammals\n')
    cat(msg)
    warning(msg)
    return()
  }

  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = classic.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'NCG',
                                    printed.name = 'Network of Cancer Gene DB',
                                    forced.parameters = forced.parameters, additional.plots = additional.plots, make.simplify = FALSE,
                                    include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
  
  return(invisible(Analysis.Data))
}

DisGeNET.classic.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                   forced.parameters = NULL, additional.plots = TRUE,
                                   include.GLM.model.name.in.Result.names = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(!(Pars$Species %in% c('hsa', 'mmu', 'rno'))){
    msg = sprintf('DisGeNET.classic.Enrichment works properly only for human and (maybe!) for mammals\n')
    cat(msg)
    warning(msg)
    return()
  }

  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = classic.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'DisGeNET',
                                    printed.name = 'DisGeNET',
                                    forced.parameters = forced.parameters, additional.plots = additional.plots, make.simplify = FALSE,
                                    include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
  
  return(invisible(Analysis.Data))
}

WikiPathways.classic.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                   forced.parameters = NULL, additional.plots = TRUE,
                                   include.GLM.model.name.in.Result.names = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  # if(!(Pars$Species %in% c('hsa', 'mmu', 'rno'))){
  #   msg = sprintf('Disease_Ont.classic.Enrichment works properly only for human and (maybe!) for mammals\n')
  #   cat(msg)
  #   warning(msg)
  #   return()
  # }

  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = classic.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'WikiPathways',
                                    printed.name = 'WikiPathways',
                                    forced.parameters = forced.parameters, additional.plots = additional.plots, make.simplify = FALSE,
                                    include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
  
  return(invisible(Analysis.Data))
}

GO.classic.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                 GO.types = '{from.parameters}', forced.parameters = NULL, additional.plots = TRUE,
                                 include.GLM.model.name.in.Result.names = FALSE, old.compartibility = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
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
                                      forced.parameters = forced.parameters, additional.plots = additional.plots, make.simplify = TRUE,
                                      include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names, old.compartibility = old.compartibility)
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
#   save(gse.res, file='KEGG trends enrichment.rds')
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
#   save(gse.res, file='Reactome trends enrichment.rds')
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
                             forced.parameters = NULL, additional.plots = TRUE, make.simplify = FALSE, simplify.cutoff = 0.55, make.simplify.plots = FALSE,
                             include.GLM.model.name.in.Result.names = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
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
  
  # if (bypass.if.completed & 
  #     Read.Completed.steps.status(sprintf("%s.%s.%s.%s.trends.Enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5,database.plus),Startup.Data$Completed.steps.file) &
  #     Read.Completed.steps.status(sprintf("%s.%s.%s.pathways.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
  #   message(sprintf('\n %s.trends.Enrichment pathway enrichment and visualization for GLM model "~ %s" was previously completed. Bypassing...',database.plus,Analysis.Data$GLM.model))
  #   return(invisible(Analysis.Data))
  # }
  
  ######################################################
  ### trends enrichment with clusterProfiler
  
  wd = sprintf("%s/%s - trends Enrichment Analysis", Analysis.Data$current.GLM.results.dir, printed.name)
  dir.create(wd,showWarnings = F)
  setwd(wd)
  eval(parse(text = sprintf('Analysis.Data$%s.trends.Enrich.working.dir = wd', database.plus)))
  
  LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
  PValue.col.name = Get.PValue.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part))
  
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

  } else if (database == 'Disease_Ont'){
    my_geneList = Get.complete.LogFCs.with.Entrez.IDs(Startup.Data, Pars, Analysis.Data)
    if(is.null(my_geneList)) return(invisible(Analysis.Data))
    fun = "gse.res = gseDO(my_geneList, pAdjustMethod = 'BH', pvalueCutoff = 1.0)"

  } else if (database == 'NCG'){
    my_geneList = Get.complete.LogFCs.with.Entrez.IDs(Startup.Data, Pars, Analysis.Data)
    if(is.null(my_geneList)) return(invisible(Analysis.Data))
    fun = "gse.res = gseNCG(my_geneList, pAdjustMethod = 'BH', pvalueCutoff = 1.0)"

  } else if (database == 'DisGeNET'){
    my_geneList = Get.complete.LogFCs.with.Entrez.IDs(Startup.Data, Pars, Analysis.Data)
    if(is.null(my_geneList)) return(invisible(Analysis.Data))
    fun = "gse.res = gseDGN(my_geneList, pAdjustMethod = 'BH', pvalueCutoff = 1.0)"

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
      my_geneList = Analysis.Data$ResTable.gene.part[,LogFC.col.name]
      names(my_geneList) = rownames(Analysis.Data$ResTable.gene.part)
      my_geneList = my_geneList[rev(order(my_geneList))]
      fun = "gse.res = gseGO(geneList = my_geneList, keyType = 'ENSEMBL', OrgDb = org.DB, ont = GO.type, pAdjustMethod='BH', pvalueCutoff = 1.0)"
    }
  } else {
    stop('Unknown database')
  }
  
  #gse.res = gseGO(geneList = my_geneList, keyType = 'ENSEMBL', OrgDb = org.DB, ont = GO.type, pAdjustMethod='BH', pvalueCutoff = FDR.Cutoff)
  my_geneList = my_geneList[!duplicated(names(my_geneList))]

  failed = FALSE
  tryCatch(expr = {
    eval(parse(text = fun))
  }, error = function (err){
    msg = sprintf('Trends enrichment (%s) failed\n', database)
    cat(msg)
    warning(msg)
    failed <<- TRUE    
  })
  if(failed)  return(invisible(NULL))

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
  
  
  eval(parse(text = sprintf('EEP___limit.PValues.to = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___limit.PValues.to', database)))
  eval(parse(text = sprintf('EEP___Max.PValue.threshold = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___Max.PValue.threshold__trends.test', database)))
  eval(parse(text = sprintf('EEP___Max.FDR.threshold = Pars$%s.clusterProfiler.Enriched.Expression.Profiles___Max.FDR.threshold__trends.test', database)))

  if(dim(gse.res@result)[1] == 0){
    msg = sprintf('trends.Enrichment was unable to find out any enriched terms/pathways/entries\n')
    cat(msg)
    warning(msg)
    return(invisible(NULL))
  }

  x = 1
  for(x in 1:dim(gse.res@result)[1]){
    current.p = max(gse.res@result$pvalue[x], EEP___limit.PValues.to)
    current.FDR = gse.res@result$p.adjust[x]
    if (current.p > EEP___Max.PValue.threshold){
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
  saveRDS(Enrich.info.split, file = sprintf('%s, %s trends Enrich. stats.rds', Analysis.Data$Analysis.name, database.plus))
  
  eval(parse(text = sprintf('pvalue.lim = Pars$%s.trends.Enrich__with.clusterProfiler___term.max.PValue.threshold', database)))
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
  saveRDS(gse.res, file=sprintf('%s, %s trends enrichment.rds', Analysis.Data$Analysis.name, database.plus))
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
      #edox <- setReadable(res.list[[t]], 'org.Hs.eg.db', 'ENTREZID')
      if(Pars$trends.Enrich__with.clusterProfiler___Render.EnrichMap){
        png.mod(filename = sprintf("%s[enrichMap] %s, %s trends enrich.png", prefixes[[t]], Analysis.Data$Analysis.name, printed.name),
            units="in", width=12, height=7, pointsize=8, res=300)
        #print({
        emapplot(res.list[[t]], color = 'pvalue')
        #})
        dev.off()
      }
      
      if(Pars$trends.Enrich__with.clusterProfiler___Render.Ridge.plot){
        png.mod(filename = sprintf("%s[ridgeplot] %s, %s trends enrich.png", prefixes[[t]], Analysis.Data$Analysis.name, printed.name),
            units="in", width=12, height=7, pointsize=3, res=300)
        ridgeplot(res.list[[t]], fill = 'pvalue')
        dev.off()
      }

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
  
  if(include.GLM.model.name.in.Result.names){
    enrich.file.name = Verify.path(sprintf("%s, %s trends enrichment.tsv",Analysis.Data$Analysis.name, printed.name))
    write.table.mod(x = as.data.frame(gse.res.converted), file = enrich.file.name,sep='\t')
    if(make.simplify){
      simplified.enrich.file.name = Verify.path(sprintf("%s, %s trends enrichment (RR).tsv",Analysis.Data$Analysis.name, printed.name))
      write.table.mod(x = as.data.frame(simplified.gse.res.converted), file = simplified.enrich.file.name, sep='\t')
    }
  } else {
    enrich.file.name = Verify.path(sprintf("%s trends enrichment.tsv", printed.name))
    write.table.mod(x = as.data.frame(gse.res.converted), file = enrich.file.name,sep='\t')
    if(make.simplify){
      simplified.enrich.file.name = Verify.path(sprintf("%s trends enrichment (RR).tsv", printed.name))
      write.table.mod(x = as.data.frame(simplified.gse.res.converted), file = simplified.enrich.file.name, sep='\t')
    }
  }
  
  Enriched.Pathways.Trended = rownames(as.data.frame(gse.res.converted))
  
  if(additional.plots & (nrow(gse.res@result) > 0)){
    gse.res.table = as.data.frame(gse.res)
    dir.create('GSEA plots', showWarnings = FALSE)
    max_GSEA_plots_per_direction = 10
    pathway.ID.list.UP = gse.res.table$ID[gse.res.table$NES > 0]
    pathway.ID.list.UP = pathway.ID.list.UP[1:min(max_GSEA_plots_per_direction, length(pathway.ID.list.UP))]
    pathway.ID.list.DOWN = gse.res.table$ID[gse.res.table$NES < 0]
    pathway.ID.list.DOWN = pathway.ID.list.DOWN[1:min(max_GSEA_plots_per_direction, length(pathway.ID.list.DOWN))]
    pathway.ID.list = c(pathway.ID.list.UP, pathway.ID.list.DOWN)
    pathway.ID.list = pathway.ID.list[!is.na(pathway.ID.list)]
    
    n = 1
    for (pathway.ID in pathway.ID.list){
      cat(sprintf('\r Rendering GSEA plot %d of %d...',n,min(max_GSEA_plots_per_direction*2,dim(gse.res.table)[1])))
      pathway.name = gsub("[^[:alnum:] ]", "", gse.res.table[pathway.ID,'Description'])
      pathway.ID.str = gsub("[^[:alnum:] ]", "_", pathway.ID)
      if(nchar(pathway.name) > 50) pathway.name = sprintf('%s...', substr(pathway.name,1,45))
      pathway.pvalue = gse.res.table[pathway.ID,'pvalue']
      
      png.mod(filename = sprintf("GSEA plots/%s - %s, DE genes distrib.png", pathway.ID.str, pathway.name),
          units="in", width=12, height=7, pointsize=8, res=150)
      print({gseaplot(gse.res, geneSetID = pathway.ID)})
      dev.off()
      n = n+1
    }
  }
  
  
  if (Pars$Create.Excel.results){
    CL = sprintf('%s "%s/Trends.enrichment.results.to.Excel.py" %s',Pars$python.bin, Pars$suppl.data.dir, gsub('Disease_Ont', 'DO', database, fixed = TRUE))
    # CL = gsub('/','\\',CL,fixed = TRUE)
    if(include.GLM.model.name.in.Result.names){
      res.file.name = sprintf('%s, %s trends enrichment.xlsx',Analysis.Data$Analysis.name, printed.name)
    } else {
      res.file.name = sprintf('%s trends enrichment.xlsx', printed.name)
    }
    CL = paste(c(CL, sprintf('"%s"',Verify.path(res.file.name)), sprintf('"%s"', enrich.file.name)), collapse = ' ')
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    if (!file.exists(Verify.path(res.file.name))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (%s trends enrichment results) for "~ %s" FAILED', printed.name, Analysis.Data$Analysis.name))
    } else {
      file.remove(enrich.file.name)
      file.copy(Verify.path(res.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, res.file.name)), overwrite = TRUE)
    }
  }
  
  if (Pars$Create.Excel.results && make.simplify){
    CL = sprintf('%s "%s/Trends.enrichment.results.to.Excel.py" %s',Pars$python.bin, Pars$suppl.data.dir, gsub('Disease_Ont', 'DO', database, fixed = TRUE))
    # CL = gsub('/','\\',CL,fixed = T)
    if(include.GLM.model.name.in.Result.names){
      res.file.name = sprintf('%s, %s trends enrichment (RR).xlsx',Analysis.Data$Analysis.name, printed.name)
    } else {
      res.file.name = sprintf('%s trends enrichment (RR).xlsx', printed.name)
    }
    CL = paste(c(CL, sprintf('"%s"',Verify.path(res.file.name)), sprintf('"%s"', simplified.enrich.file.name)), collapse = ' ')
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    if (!file.exists(Verify.path(res.file.name))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (%s trends enrichment results) for "~ %s" FAILED', printed.name, Analysis.Data$Analysis.name))
    } else {
      file.remove(simplified.enrich.file.name)
      file.copy(Verify.path(res.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, res.file.name)), overwrite = TRUE)
    }
  }
  
  # sprintf('%s, %s trends Enrich. stats.rds', Analysis.Data$Analysis.name, printed.name)
  if(include.GLM.model.name.in.Result.names){
    saveRDS(Enriched.Pathways.Trended, file = sprintf('%s, %s enriched pathways (trends).rds', Analysis.Data$Analysis.name, database.plus))
  } else {
    saveRDS(Enriched.Pathways.Trended, file = sprintf('%s enriched pathways (trends).rds', database.plus))
  }

  
  # Add.Completed.step(sprintf("%s.%s.%s.%s.trends.Enrichment", 
  #                            Analysis.Data$GLM.model, Startup.Data$Startup.hashmd5, Analysis.Data$step.hashmd5, database.plus), Startup.Data$Completed.steps.file)
  tmp = sprintf("Analysis.Data$%s.trends.Enriched.Pathways = levels(factor(Enriched.Pathways.Trended))", database.plus)
  eval(parse(text = tmp))
  
  tmp = sprintf("Analysis.Data$%s.trends.Enrichment.performed = TRUE", database.plus)
  eval(parse(text = tmp))
  
  setwd(Analysis.Data$current.GLM.results.dir)
  
  cat(sprintf('\r%s enrichment (trends) completed.       \n',printed.name))
  return(invisible(Analysis.Data))
}



Reactome.trends.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                      forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = trends.Enrichment(Startup.Data, Analysis.Data,
                             database = 'Reactome',
                             printed.name = 'Reactome Pathways',
                             forced.parameters = forced.parameters, additional.plots = additional.plots,
                             make.simplify = FALSE)
  
  return(invisible(Analysis.Data))
}

KEGG.trends.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                  forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = trends.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'KEGG',
                                    printed.name = 'KEGG Pathways',
                                    forced.parameters = forced.parameters, additional.plots = additional.plots,
                                    make.simplify = FALSE)
  
  return(invisible(Analysis.Data))
}

GO.trends.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                GO.types = '{from.parameters}', forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
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
                                      forced.parameters = forced.parameters, additional.plots = additional.plots,
                                      make.simplify = FALSE)
  }
  return(invisible(Analysis.Data))
}


DisGeNET.trends.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                  forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(!(Pars$Species %in% c('hsa', 'mmu', 'rno'))){
    msg = sprintf('DisGeNET.trends.Enrichment works properly only for human and (maybe!) for mammals\n')
    cat(msg)
    warning(msg)
    return()
  }

  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = trends.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'DisGeNET',
                                    printed.name = 'DisGeNET',
                                    forced.parameters = forced.parameters, additional.plots = additional.plots,
                                    make.simplify = FALSE)
  
  return(invisible(Analysis.Data))
}

Disease_Ont.trends.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                  forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(!(Pars$Species %in% c('hsa', 'mmu', 'rno'))){
    msg = sprintf('Disease_Ont.trends.Enrichment works properly only for human and (maybe!) for mammals\n')
    cat(msg)
    warning(msg)
    return()
  }

  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = trends.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'Disease_Ont',
                                    printed.name = 'Disease Ontology DB',
                                    forced.parameters = forced.parameters, additional.plots = additional.plots,
                                    make.simplify = FALSE)
  
  return(invisible(Analysis.Data))
}

NCG.trends.Enrichment = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                  forced.parameters = NULL, additional.plots = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
  }
  
  if(typeof(Analysis.Data) == 'character'){
    GLM.model = Analysis.Data
    Analysis.Data = NULL
  }
  
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(!(Pars$Species %in% c('hsa', 'mmu', 'rno'))){
    msg = sprintf('NCG.trends.Enrichment works properly only for human and (maybe!) for mammals\n')
    cat(msg)
    warning(msg)
    return()
  }

  if(is.null(Analysis.Data)){
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  Analysis.Data = trends.Enrichment(Startup.Data, Analysis.Data,
                                    database = 'NCG',
                                    printed.name = 'Network of Cancer Gene DB',
                                    forced.parameters = forced.parameters, additional.plots = additional.plots,
                                    make.simplify = FALSE)
  
  return(invisible(Analysis.Data))
}




KEGG.Pathways.Visualization = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                       Add.enriched.pathways = TRUE,
                                       forced.parameters = NULL, forced.pathways = NULL,
                                       same.layer = 'adaptive', kegg.native = TRUE, write.de.info = TRUE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  if (!Analysis.Data$GLM.analysis.performed) stop('Pathway.Visualization can be processed only after GLM/DE analysis. It is not performed.')
  # if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  # } else Pars = forced.parameters
  
  # if (bypass.if.completed & Read.Completed.steps.status(sprintf("%s.%s.%s.pathways.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
  #   message(sprintf('\nKEGG pathway visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
  #   return()
  # }
  
  ######################################################
  ### Pathway visualization
  
  if(is.null(forced.pathways)){
    Pathways.to.Visualize = Startup.Data$Custom.Pathways.to.Visualize
    if(Add.enriched.pathways){
      if(!is.null(Analysis.Data$KEGG.trends.Enriched.Pathways))      Pathways.to.Visualize = c(Pathways.to.Visualize, Analysis.Data$KEGG.trends.Enriched.Pathways)
      if(!is.null(Analysis.Data$KEGG.classic.Enriched.Pathways))      Pathways.to.Visualize = c(Pathways.to.Visualize, Analysis.Data$KEGG.classic.Enriched.Pathways)
      
      file.name = sprintf("%s/KEGG Pathways - classic Enrichment Analysis/%s, KEGG enriched pathways (classic).rds", Analysis.Data$current.GLM.results.dir, Analysis.Data$Analysis.name)
      if(file.exists(file.name)){
        ep = readRDS(file.name)
        Pathways.to.Visualize = c(Pathways.to.Visualize, ep)
      }
      
      file.name = sprintf("%s/KEGG Pathways - trends Enrichment Analysis/%s, KEGG enriched pathways (trends).rds", Analysis.Data$current.GLM.results.dir, Analysis.Data$Analysis.name)
      if(file.exists(file.name)){
        ep = readRDS(file.name)
        Pathways.to.Visualize = c(Pathways.to.Visualize, ep)
      }
    }
    
    Pathways.to.Visualize = levels(as.factor(Pathways.to.Visualize))
    #Pathways.to.Visualize = c("hsa01100","mmu01100","dme01100","hsa00511","hsa00514","hsa00533","hsa04215")
    
    forbidden.pathways = c(Startup.Data$Omit.pathways,"01100","00511","00514","00533","04215")
    # forbidden.pathways = c("01100","00511","00514","00533","04215")
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
    
    KEGG.Pathways.Visualization(Startup.Data, Analysis.Data = Analysis.Data,
                                     forced.parameters = forced.parameters, forced.pathways = Pathways.to.Visualize[pathway.IDs.num < 2000],
                                     same.layer = FALSE, kegg.native = TRUE, write.de.info = write.de.info)
    KEGG.Pathways.Visualization(Startup.Data, Analysis.Data = Analysis.Data,
                                     forced.parameters = forced.parameters, forced.pathways = Pathways.to.Visualize[pathway.IDs.num >= 2000],
                                     same.layer = TRUE, kegg.native = TRUE, write.de.info = write.de.info)
    # Add.Completed.step(sprintf("%s.%s.%s.pathways.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
    Analysis.Data$Pathway.Visualization.performed = TRUE
    cat('\nPathway visualization completed.\n')
    
  } else {
    cat(sprintf("\n%s:   Performing pathway visualization...",Analysis.Data$GLM.model))
    Analysis.Data$PA.working.dir = sprintf("%s/KEGG Pathways - Visualization",Analysis.Data$current.GLM.results.dir)
    dir.create(Analysis.Data$PA.working.dir, showWarnings = FALSE)
    setwd(Analysis.Data$PA.working.dir)
    
    DE.info.file.names.by.pathway.id = vector(mode="list")
    
    LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
    PValue.col.name = Get.PValue.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part))

    if(Pars$Pathway.Visualization___CPM.aware){
      logFCs = Analysis.Data$ResTable.gene.part[,LogFC.col.name]
      CPMs = 2^Analysis.Data$ResTable.gene.part$LogCPM
      PValues = Analysis.Data$ResTable.gene.part[,PValue.col.name]
      logFCs [PValues > Pars$Pathway.Visualization___max.gene.PValue] = 0
      names(logFCs) = rownames(Analysis.Data$ResTable.gene.part)
      names(CPMs) = rownames(Analysis.Data$ResTable.gene.part)
      #logFCs = logFCs[tt$table$PValue < Pathway.Visualization___max.gene.PValue]
      #CPMs = CPMs[tt$table$PValue < Pathway.Visualization___max.gene.PValue]
      
      cat(sprintf("\nTotal %d genes has been selected for pathway analysis.\nNow rendering pathway images (modified pathview)....\n",sum.mod(abs(logFCs) > 0)))
      
      #Pathways.to.Visualize=c('hsa04110')
      pv.out <- suppressWarnings(suppressMessages(pathviewmod(gene.data = logFCs,
                                             gene.idtype = "ENSEMBL",
                                             pathway.id = Pathways.to.Visualize, species = Pars$Species, same.layer = same.layer,
                                             out.suffix = Analysis.Data$Analysis.name, keys.align = "y", kegg.native = kegg.native,
                                             limit = list(gene = Pars$Pathway.Visualization___LogFC.limits), bins = 20, cpm.data = CPMs,
                                             write.de.info = write.de.info, cached.dir=sprintf('%s/KEGG.cache', Pars$suppl.data.dir))))
    
    } else {
      logFCs = Analysis.Data$ResTable.gene.part[,LogFC.col.name]
      names(logFCs) = rownames(Analysis.Data$ResTable.gene.part)
      logFCs = logFCs[Analysis.Data$ResTable.gene.part[,PValue.col.name] < Pars$Pathway.Visualization___max.gene.PValue]
      
      cat(sprintf("\nTotal %d genes has been selected for pathway analysis.\nNow rendering pathway images (standard pathview)....\n",sum.mod(abs(logFCs) > 0)))

      pv.out <- suppressWarnings(suppressMessages(pathview(gene.data = logFCs,
                                          gene.idtype = "ENSEMBL",
                                          pathway.id = Pathways.to.Visualize, species = Pars$Species, same.layer = same.layer,
                                          out.suffix = Analysis.Data$Analysis.name, keys.align = "y", kegg.native = kegg.native,
                                          limit = list(gene = Pars$Pathway.Visualization___LogFC.limits), bins = list(gene = 12))))
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
                            forced.parameters = NULL, include.GLM.model.name.in.Result.names = FALSE){
  for(GO.type in GO.types){
    topGO.Enrichment.single.ontology(Startup.Data, Analysis.Data = Analysis.Data, GLM.model = GLM.model, GLM.working.dir = GLM.working.dir,
                                     Analysis.Data.RDS.file = Analysis.Data.RDS.file, GO.type = GO.type,
                                     forced.parameters = forced.parameters,
                                     include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
  }
}

topGO.Enrichment.single.ontology = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                            GO.type = 'BP',
                            forced.parameters = NULL, include.GLM.model.name.in.Result.names = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }
  
  LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
  PValue.col.name = Get.PValue.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part))
  NP.PValue.col.name = Get.NP.PValue.col.name(colnames(Analysis.Data$ResTable.gene.part))
  LogCPM.col.name = Get.LogCPM.col.name(colnames(Analysis.Data$ResTable.gene.part))
  
  # if (!Analysis.Data$GLM.analysis.performed) stop('topGO.Enrichment can be processed only after GLM/DE analysis. It is not performed.')
  # if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  # } else Pars = forced.parameters
  
  
  # if (bypass.if.completed & 
  #     Read.Completed.steps.status(sprintf("%s.%s.%s.GO.expression.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file) &
  #     Read.Completed.steps.status(sprintf("%s.%s.%s.GO.enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
  #   message(sprintf('\nGO enrichment and GO-centric expression profiles visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
  #   return(Analysis.Data)
  # }
  
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
        Passed = (abs(Analysis.Data$ResTable.gene.part[,LogFC.col.name]) >= Pars$GO.Enrich__with.topGO___gene.min.abs.LogFC.threshold) & 
          (Analysis.Data$ResTable.gene.part[,PValue.col.name] <= Pars$GO.Enrich__with.topGO___gene.max.PValue.threshold) & 
          (Analysis.Data$ResTable.gene.part[,"Score"] >= Pars$GO.Enrich__with.topGO___gene.min.Score.threshold)
      }
      if (GO.Enrich__with.topGO___mode == 'upreg'){
        Passed = (Analysis.Data$ResTable.gene.part[,LogFC.col.name] > 0) & 
          (abs(Analysis.Data$ResTable.gene.part[,LogFC.col.name]) >= Pars$GO.Enrich__with.topGO___gene.min.abs.LogFC.threshold) & 
          (Analysis.Data$ResTable.gene.part[,PValue.col.name] <= Pars$GO.Enrich__with.topGO___gene.max.PValue.threshold) & 
          (Analysis.Data$ResTable.gene.part[,"Score"] >= Pars$GO.Enrich__with.topGO___gene.min.Score.threshold)
      }
      if (GO.Enrich__with.topGO___mode == 'downreg'){
        Passed = (Analysis.Data$ResTable.gene.part[,LogFC.col.name] < 0) & 
          (abs(Analysis.Data$ResTable.gene.part[,LogFC.col.name]) >= Pars$GO.Enrich__with.topGO___gene.min.abs.LogFC.threshold) & 
          (Analysis.Data$ResTable.gene.part[,PValue.col.name] <= Pars$GO.Enrich__with.topGO___gene.max.PValue.threshold) & 
          (Analysis.Data$ResTable.gene.part[,"Score"] >= Pars$GO.Enrich__with.topGO___gene.min.Score.threshold)
      }

      if(!is.na(NP.PValue.col.name) & Pars$GO.Enrich__with.topGO___gene.max.NP.PValue.threshold < 1)
        Passed = Passed & (Analysis.Data$ResTable.gene.part[,NP.PValue.col.name] <= Pars$GO.Enrich__with.topGO___gene.max.NP.PValue.threshold)
      
      if(!is.null(Pars$GO.Enrich__with.topGO___gene.min.LogCPM.threshold))
        Passed = Passed & (Analysis.Data$ResTable.gene.part[,LogCPM.col.name] >= Pars$GO.Enrich__with.topGO___gene.min.LogCPM.threshold)

      if (sum.mod(Passed) <= GO.Enrich__with.topGO___max.DE.genes)  Do.end.run = TRUE
      
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
                } else {  current.p.value = max(allRes[i,'P.Fisher'],Pars$topGO.Enriched.Expression.Profiles___limit.PValues.to)  }
                
                current.score = 0
                if (current.p.value > Pars$topGO.Enriched.Expression.Profiles___Max.PValue.threshold){
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
            current.p.value = max(allRes[i,'P.Fisher'],Pars$topGO.Enriched.Expression.Profiles___limit.PValues.to)  }
          
          current.score = 0
          if (current.p.value > Pars$topGO.Enriched.Expression.Profiles___Max.PValue.threshold){
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
        if(include.GLM.model.name.in.Result.names){
          filename = sprintf('%s, GO GSEA for top-%s %s genes, top %d GO.terms.png',Analysis.Data$Analysis.name,GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode,N)
        } else {
          filename = sprintf('GO GSEA for top-%s %s genes, top %d GO.terms.png', GO.Enrich__with.topGO___max.DE.genes, GO.Enrich__with.topGO___mode, N)
        }
        png.mod(filename = filename, units="in",width=8.23,height=6.3,pointsize=12,
            res=600)
        suppressWarnings(showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = N, useInfo = 'all'))
        dev.off()
      }
      
      
      if(include.GLM.model.name.in.Result.names){
        filename = sprintf('%s, GO GSEA for top-%d %s genes.pdf',Analysis.Data$Analysis.name,GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode)
      } else {
        filename = sprintf('GO GSEA for top-%d %s genes.pdf', GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode)
      }
      pdf(file = filename, width=8.23, height=6.3,onefile = T)
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
        all.GT.logFCs = Analysis.Data$ResTable.gene.part[names(gene.scores),LogFC.col.name]
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
      
      if(include.GLM.model.name.in.Result.names){
        filename = Verify.path(sprintf("%s, GO GSEA for top-%s %s genes.txt",Analysis.Data$Analysis.name,GO.Enrich__with.topGO___max.DE.genes,GO.Enrich__with.topGO___mode))
      } else {
        filename = Verify.path(sprintf("GO GSEA for top-%s %s genes.txt", GO.Enrich__with.topGO___max.DE.genes, GO.Enrich__with.topGO___mode))
      }
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
    # CL = gsub('/','\\',CL,fixed = T)
    if(include.GLM.model.name.in.Result.names){
      file.base.name = sprintf('%s, GO Terms (%s) topGO enrichment.xlsx', Analysis.Data$Analysis.name, GO.type)
    } else {
      file.base.name = sprintf('GO Terms (%s) topGO enrichment.xlsx', GO.type)
    }
    CL = paste(c(CL, Verify.path(sprintf('"%s"', file.base.name)),sprintf('"%s"',output.file.names)),collapse = ' ')
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
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

  if(include.GLM.model.name.in.Result.names){
    write.table.mod(x=GO.terms.GSEA.info,file=sprintf('%s, GO terms GSEA stats.tsv',Analysis.Data$Analysis.name),sep='\t',quote = F,na = "")
  } else {
    write.table.mod(x=GO.terms.GSEA.info,file=sprintf('GO terms GSEA stats.tsv'),sep='\t',quote = F,na = "")
  }

  setwd(Analysis.Data$current.GLM.results.dir)
  # Add.Completed.step(sprintf("%s.%s.%s.GO.enrichment",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
  
  cat('\nGene Ontology enrichment completed.\n')
  return(invisible(Analysis.Data))
}

topGO.Expression.Profiles.Custom = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                            GO.type = 'BP', profile.type = 'custom',
                                            forced.parameters = NULL, forced.terms = NULL,
                                            out.dir = NULL, forced.Analysis.name = NULL,GO.stats.file.name = NULL, sparklines.axis.limits = 2, min.term.size = NULL,
                                            include.GLM.model.name.in.Result.names = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
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
  
  # if (bypass.if.completed & 
  #     Read.Completed.steps.status(sprintf("%s.%s.%s.GO.expression.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
  #   message(sprintf('\nGO-centric expression profiles visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
  #   return()
  # }
  
  LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
  PValue.col.name = Get.PValue.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part))
  
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
    cf1 = sprintf('%s/GO Terms (%s) - topGO - Enrichment Analysis/GO terms GSEA stats.tsv', Analysis.Data$current.GLM.results.dir, GO.type)
    cf2 = sprintf('%s/GO Terms (%s) - topGO - Enrichment Analysis/%s, GO terms GSEA stats.tsv', Analysis.Data$current.GLM.results.dir, GO.type, Analysis.name)
    if(file.exists(cf1)){
      GO.stats.file.name = cf1
    } else if(file.exists(cf2)){
      GO.stats.file.name = cf2
    }
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
              if (length(unavailable.terms) < 5){
                cat(sprintf('\n Term %s is not found in DB',gt))
              } else if (length(unavailable.terms) == 5){
                cat('\n .... ')
              }
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
        all.GT.logFCs = Analysis.Data$ResTable.gene.part[gene.names,LogFC.col.name]
        names(all.GT.logFCs) = gene.names
        all.GT.logFCs = all.GT.logFCs[Analysis.Data$ResTable.gene.part[gene.names,"LogCPM"] >= topGO.Expression.Profiles___Min.gene.logCPM]
        
        gene.names = names(all.GT.logFCs)
        all.GT.logFCs = all.GT.logFCs * (Analysis.Data$ResTable.gene.part[gene.names, PValue.col.name] <= topGO.Expression.Profiles___Max.P)
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
      if(include.GLM.model.name.in.Result.names){
        file.name = Verify.path(sprintf('%s, GO-DE info, logCPM.gt.%g, P.lt.%g.tsv', Analysis.name, topGO.Expression.Profiles___Min.gene.logCPM, topGO.Expression.Profiles___Max.P))
      } else {
        file.name = Verify.path(sprintf('GO-DE info, logCPM.gt.%g, P.lt.%g.tsv', topGO.Expression.Profiles___Min.gene.logCPM, topGO.Expression.Profiles___Max.P))
      }
      
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
    # CL = gsub('/','\\',CL,fixed = T)
    if(include.GLM.model.name.in.Result.names){
      res.file.name.base = sprintf('%s, GO Terms (%s) - topGO - %s DE profiles.xlsx', Analysis.name, GO.type, profile.type)
    } else {
      res.file.name.base = sprintf('GO Terms (%s) - topGO - %s DE profiles.xlsx', GO.type, profile.type)
    }

    CL = paste(c(CL,sprintf('--out-excel "%s"', Verify.path(res.file.name.base)), '--in', sprintf('"%s"',output.file.names),
                 '--models', sprintf('"%s"', rep(Analysis.name, length(output.file.names))),
                 '--maximal-genes-count',as.character(Pars$topGO.Expression.Profiles___Maximal.genes.in.terms.to.visualize),
                 sprintf('--sparklines-axis-limit %g',sparklines.axis.limits)),
               collapse = ' ')
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    if (!file.exists(Verify.path(res.file.name.base))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (GO-centric DE profiles) for ~%s FAILED',Analysis.name))
      cat(CL)
    } else {
      file.copy(Verify.path(res.file.name.base), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, res.file.name.base)), overwrite = TRUE)
    }
    
    CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces yes',Pars$python.bin,Pars$suppl.data.dir)
    # CL = gsub('/','\\',CL,fixed = T)
    if(include.GLM.model.name.in.Result.names){
      res.file.name.base = sprintf('%s, GO Terms (%s) - topGO - %s DE profiles (with spaces).xlsx', Analysis.name, GO.type, profile.type)
    } else {
      res.file.name.base = sprintf('GO Terms (%s) - topGO - %s DE profiles (with spaces).xlsx', GO.type, profile.type)
    }
    CL = paste(c(CL,sprintf('--out-excel "%s"', Verify.path(res.file.name.base)),'--in',sprintf('"%s"',output.file.names),
                 '--models', sprintf('"%s"', rep(Analysis.name, length(output.file.names))),
                 '--maximal-genes-count',as.character(Pars$topGO.Expression.Profiles___Maximal.genes.in.terms.to.visualize),
                 sprintf('--sparklines-axis-limit %g',sparklines.axis.limits)),collapse = ' ')
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    if (!file.exists(Verify.path(res.file.name.base))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (GO-centric DE profiles) for ~%s FAILED',Analysis.name))
      cat(CL)
    }
  }
  else {
    file.copy(Verify.path(res.file.name.base), sprintf('%s/%s', Verify.path(Analysis.Data$current.GLM.results.dir, res.file.name.base)), overwrite = TRUE)
  }
  
  # Add.Completed.step(sprintf("%s.%s.%s.GO.expression.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)
  
  
  if (is.null(out.dir)) setwd(Analysis.Data$current.GLM.results.dir)
  Analysis.Data$GO.Expression.Profiles.Created = T
  cat('\nCreating GO-centric gene expression profiles completed.\n')
}

topGO.Expression.Profiles.Enriched = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                              GO.type = 'BP',
                                              forced.parameters = NULL, out.dir = NULL, forced.Analysis.name = NULL,
                                              GO.stats.file.name = NULL,sparklines.axis.limits = 2,
                                              include.GLM.model.name.in.Result.names = FALSE){
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
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
  
  # if (bypass.if.completed & 
  #     Read.Completed.steps.status(sprintf("%s.%s.%s.GO.expression.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
  #   message(sprintf('\nGO-centric expression profiles visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
  #   return()
  # }
  
  
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
    cf1 = sprintf('%s/GO Terms (%s) - topGO - Enrichment Analysis/GO terms GSEA stats.tsv', Analysis.Data$current.GLM.results.dir, GO.type)
    cf2 = sprintf('%s/GO Terms (%s) - topGO - Enrichment Analysis/%s, GO terms GSEA stats.tsv', Analysis.Data$current.GLM.results.dir, GO.type, Analysis.name)
    if(file.exists(cf1)){
      GO.stats.file.name = cf1
    } else if(file.exists(cf2)){
      GO.stats.file.name = cf2
    }
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
                                   sparklines.axis.limits = sparklines.axis.limits, min.term.size = Pars$topGO.Enriched.Expression.Profiles___minimal.term.size,
                                   include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
}


Collect.GLM.results = function(Startup.Data, GLM.models = NULL, forced.parameters = NULL, out.dir = NULL,
                               include.GLM.model.name.in.Result.names = FALSE,
                               include.GLM.model.name.in.Result.names..as.suffix = TRUE){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if (!Pars$Create.Excel.results){
    warning('Create.Excel.results parameter is turned OFF. "Collect.GLM.results" will not run')
    return()
  }
  
  if(is.null(GLM.models)) GLM.models = Pars$Models.to.Test

  if(is.null(out.dir))  out.dir = Pars$results.dir
  
  for(DE.out.file.type in c('Log.rel.CPM profiles', 'CPM profiles')){
    file.names = c()
    for (GLM.model in GLM.models){
      GLM.model = Delete.prefix.in.model(GLM.model)
      Analysis.name  = Cleanup.Model.Name(GLM.model)
      if(include.GLM.model.name.in.Result.names){
        f = sprintf("%s/~ %s, results/%s, DE results, %s.xlsx", Pars$results.dir, Analysis.name, Analysis.name, DE.out.file.type)
      } else if(include.GLM.model.name.in.Result.names..as.suffix){
        f = sprintf("%s/~ %s, results/DE results, %s - %s.xlsx", Pars$results.dir, Analysis.name, DE.out.file.type, Analysis.name)
      } else {
        f = sprintf("%s/~ %s, results/DE results, %s.xlsx", Pars$results.dir, Analysis.name, DE.out.file.type)
      }
      # f = gsub("/","\\",f,fixed = T)
      f = Verify.path(f)

      if(file.exists(f))  temp = file.copy(f, sprintf("%s/[DE results, %s]    %s.xlsx", out.dir, DE.out.file.type, Analysis.name))

      file.names = append(file.names, f)
    }
    if (sum.mod(!file.exists(file.names))>0){
      warning(sprintf('The following files with DE/GLM info do not exist (%s): %s', DE.out.file.type, toString(file.names[!file.exists(file.names)])))
    }
  }
}



Summarize.GLM.results = function(Startup.Data, GLM.models = NULL, forced.parameters = NULL, out.file.name.base = '@ Summary expression info'){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if (!Pars$Create.Excel.results){
    warning('Create.Excel.results parameter is turned OFF. "Summarize.GLM.results" will not run')
    return()
  }
  
  if(is.null(GLM.models)) GLM.models = Pars$Models.to.Test
  file.names = c()
  Analysis.names = c()
  for (GLM.model in GLM.models){
    GLM.model = Delete.prefix.in.model(GLM.model)
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    f = sprintf("%s/~ %s, results/%s_%s_combined.txt",Pars$results.dir,Analysis.name,Analysis.name,Pars$DE.package)
    # f = gsub("/","\\",f,fixed = T)
    f = Verify.path(f)
    file.names = c(file.names, f)
    Analysis.names = c(Analysis.names, Analysis.name)
  }

  if (sum.mod(!file.exists(file.names))>0){
    warning(sprintf('The following files with DE/GLM info do not exist: %s',toString(file.names[!file.exists(file.names)])))
    Analysis.names = Analysis.names[file.exists(file.names)]
    file.names = file.names[file.exists(file.names)]
  }

  
  
  
  ##### 
  # CL = sprintf('%s "%s/DE.results.to.Excel.py"',Pars$python.bin, Pars$suppl.data.dir)
  # CL = gsub('/','\\',CL,fixed = T)
  #if(is.null(out.dir)) out.dir = Pars$results.dir
  if(is.null(out.file.name.base)){
    out.file.name = sprintf('%s/Summary expression info.xlsx', Pars$results.dir)
    out.file.name.ff = sprintf('%s/Summary expression info, joint.xlsx', Pars$results.dir)
    out.file.name.lite = sprintf('%s/Summary expression info, joint, lite.xlsx', Pars$results.dir)
    out.file.name.lite.logfconly = sprintf('%s/Summary expression info, joint, lite, LogFC table.xlsx', Pars$results.dir)
  } else {
    out.file.name.ff = sprintf("%s, joint.xlsx", out.file.name.base)
    out.file.name.lite = sprintf("%s, joint lite.xlsx", out.file.name.base)
    out.file.name.lite.logfconly = sprintf("%s, joint lite, LogFC table.xlsx", out.file.name.base)
  }
  # CL = paste(c(CL,sprintf('"%s"',Verify.path(out.file.name)),sprintf('"%s"',file.names)),collapse = ' ')
  # system(CL)
  # if (!file.exists(Verify.path(out.file.name))){
  #   #file.remove(output.file.names)
  #   stop(sprintf('Excel worksheet generation FAILED [non-joint], %s',out.file.name))
  # }
  
  
  ### running merging txt > Excel conversion
  CL_base = sprintf('%s "%s/DE.results.to.Excel.joint.py"', Pars$python.bin, Pars$suppl.data.dir)
  # CL_base = gsub('/','\\', CL_base, fixed = TRUE)
  
  CL.ff = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs.tsv" --lr no --score no --logcpm yes --logcpm-common yes --spearmanp no --pearsonr no --pearsonp no',
                  CL_base, paste(sprintf('"%s"',file.names),collapse=' '), Verify.path(out.file.name.ff), Pars$results.dir)
  CL.lite = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs.tsv" --lr no --score no --logcpm yes --logcpm-common yes --spearmanp no --spearmanr no --fdr no --pearsonr no --pearsonp no --simple-logfc-format-mode yes',
                  CL_base, paste(sprintf('"%s"',file.names),collapse=' '), Verify.path(out.file.name.lite), Pars$results.dir)
  CL.lite.logfconly = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs.tsv" --p no --fdr no --lr no --score no --logcpm no --logcpm-common yes --spearmanr no --spearmanp no --pearsonr no --pearsonp no --simple-logfc-format-mode yes',
                    CL_base, paste(sprintf('"%s"',file.names),collapse=' '), Verify.path(out.file.name.lite.logfconly), Pars$results.dir)
  
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

Create.Euler.Venn.diagrams = function(Startup.Data, GLM.models = NULL, forced.parameters = NULL,
  output.dir = NULL, bypass..too.large.comparisons = TRUE, reported.analysis.group.name = '', write.intersections = TRUE, img.resolution = NULL){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(GLM.models)) GLM.models = Pars$Models.to.Test
  RDS.file.names = c()
  Analysis.names = c()
  for (GLM.model in GLM.models){
    GLM.model = Delete.prefix.in.model(GLM.model)
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    f = sprintf("%s/~ %s, results/analysis data.rds", Pars$results.dir, Analysis.name)
    f = Verify.path(f)
    RDS.file.names = c(RDS.file.names, f)
    Analysis.names = c(Analysis.names, Analysis.name)
  }

  if(is.null(img.resolution))  img.resolution = Pars$Venn.diagrams___image.resolution

  if (sum.mod(!file.exists(RDS.file.names))>0){
    warning(sprintf('The following files with DE/GLM info do not exist: %s', toString(RDS.file.names[!file.exists(RDS.file.names)])))
    Analysis.names = Analysis.names[file.exists(RDS.file.names)]
    RDS.file.names = RDS.file.names[file.exists(RDS.file.names)]
  }

  if(is.null(output.dir)){
    output.dir = getwd()
  } else dir.create(output.dir, showWarnings = FALSE)

  all.tsv.file.names = c()
  all.Excel.sheet.names = c()

  for(threshold.n in 1:length(Pars$Venn.diagrams___threshold.names.list)){
    threshold.name = Pars$Venn.diagrams___threshold.names.list[threshold.n]

    t1 = Pars$Venn.diagrams___PValue.higher.bounds.list[threshold.n]
    t2 = Pars$Venn.diagrams___PValue.lower.bounds.list[threshold.n]
    Venn.diagrams___PValue.higher.bound = max(t1, t2)
    Venn.diagrams___PValue.lower.bound = min(t1, t2)

    t1 = Pars$Venn.diagrams___FDR.higher.bounds.list[threshold.n]
    t2 = Pars$Venn.diagrams___FDR.lower.bounds.list[threshold.n]
    Venn.diagrams___FDR.higher.bound = max(t1, t2)
    Venn.diagrams___FDR.lower.bound = min(t1, t2)

    t1 = Pars$Venn.diagrams___NP.PValue.higher.bounds.list[threshold.n]
    t2 = Pars$Venn.diagrams___NP.PValue.lower.bounds.list[threshold.n]
    Venn.diagrams___NP.PValue.higher.bound = max(t1, t2)
    Venn.diagrams___NP.PValue.lower.bound = min(t1, t2)
    
    t1 = Pars$Venn.diagrams___abs.LogFC.higher.bounds.list[threshold.n]
    t2 = Pars$Venn.diagrams___abs.LogFC.lower.bounds.list[threshold.n]
    Venn.diagrams___abs.LogFC.higher.bound = max(t1, t2)
    Venn.diagrams___abs.LogFC.lower.bound = min(t1, t2)

    t1 = Pars$Venn.diagrams___LogCPM.higher.bounds.list[threshold.n]
    t2 = Pars$Venn.diagrams___LogCPM.lower.bounds.list[threshold.n]
    Venn.diagrams___LogCPM.higher.bound = max(t1, t2)
    Venn.diagrams___LogCPM.lower.bound = max(t1, t2)

    cat(sprintf('[Create.Euler.Venn.diagrams]  %s - %s - Loading RDS data and deriving the complete gene list....\n', threshold.name, reported.analysis.group.name))
    all.DE.genes = c()
    all.DE.genes..NP = c()
    include.NP.diagrams = FALSE
    analyses.n.to.include = c()
    max.avg.group.size = 0
    for(n in 1:length(Analysis.names)){
      Analysis.Data = readRDS(RDS.file.names[n])
      avg.group.size = mean(sapply(Analysis.Data$CPMs.groups, length))
      NP.PValue.col.name = Get.NP.PValue.col.name(colnames(Analysis.Data$ResTable.gene.part))
      if(avg.group.size >= Pars$Venn.diagrams___Include.NP.diagrams.if.avg.samples.N.per.group.exceeds && !is.na(NP.PValue.col.name))
        include.NP.diagrams = TRUE
      max.avg.group.size = max(avg.group.size, max.avg.group.size)

      PValue.col.name = Get.PValue.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part))
      FDR.col.name = Get.FDR.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part), look.for.paired.test.res = Analysis.Data$perform.classic.paired.test)
      LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
      LogCPM.col.name = Get.LogCPM.col.name(colnames(Analysis.Data$ResTable.gene.part))
      if(is.na(PValue.col.name) | is.na(FDR.col.name) | is.na(LogFC.col.name) | is.na(LogCPM.col.name)){
        msg = sprintf('[Create.Euler.Venn.diagrams]  %s - Analysis "%s": cannot determine column names for P-value, FDR, LogFC or LogCPM. This analysis will be bypassed', reported.analysis.group.name, Analysis.names[n])
        cat(msg)
        warning(msg)
        next
      }
      analyses.n.to.include = c(analyses.n.to.include, n)
      keep.by.PValue = Analysis.Data$ResTable.gene.part[,PValue.col.name] <= Venn.diagrams___PValue.lower.bound
      keep.by.FDR = Analysis.Data$ResTable.gene.part[,FDR.col.name] <= Venn.diagrams___FDR.lower.bound
      keep.by.LogFC = abs(Analysis.Data$ResTable.gene.part[,LogFC.col.name]) >= Venn.diagrams___abs.LogFC.higher.bound
      keep.by.LogCPM = Analysis.Data$ResTable.gene.part[,LogCPM.col.name] >= Venn.diagrams___LogCPM.higher.bound
      keep = keep.by.PValue & keep.by.FDR & keep.by.LogFC & keep.by.LogCPM
      all.DE.genes = c(all.DE.genes, rownames(Analysis.Data$ResTable.gene.part)[keep])

      if(!is.na(NP.PValue.col.name)){
        keep.by.NP.PValue = Analysis.Data$ResTable.gene.part[, NP.PValue.col.name] <= Venn.diagrams___NP.PValue.lower.bound
        all.DE.genes..NP = c(all.DE.genes..NP, rownames(Analysis.Data$ResTable.gene.part)[keep & keep.by.NP.PValue])
      }
    }
    all.DE.genes = sort(DeDup.na.rm(all.DE.genes))
    all.DE.genes..NP = sort(DeDup.na.rm(all.DE.genes..NP))

    Analysis.names = Analysis.names[analyses.n.to.include]
    RDS.file.names = RDS.file.names[analyses.n.to.include]

    if(length(all.DE.genes) == 0){
      cat('No genes have passed differential expression thresholds. Venn diagrams will not be created\n')
      next
    }

    cat(sprintf('[Create.Euler.Venn.diagrams]  %s - %s - Euler-Venn diagrams for genes with non-parametric test filtering will%s be created (max. avg. sample group size is %.1f)\n',
      reported.analysis.group.name, threshold.name, ifelse(include.NP.diagrams, '', ' NOT'), max.avg.group.size))

    cat(sprintf('[Create.Euler.Venn.diagrams]  %s - %s - Total %d genes are considered as differentially expressed among all %d analyses\n',
      reported.analysis.group.name, threshold.name, length(all.DE.genes), length(Analysis.names)))
    if(include.NP.diagrams){
      cat(sprintf('[Create.Euler.Venn.diagrams]  %s - %s - Total %d genes are considered as differentially expressed (plus non-parametric tests) among all %d analyses\n',
        reported.analysis.group.name, threshold.name, length(all.DE.genes..NP), length(Analysis.names)))
      if(length(all.DE.genes..NP) == 0){
        cat('NP diagrams will be disabled\n')
        include.NP.diagrams = FALSE
      }
    }

    gene.vs.analysis.DE.table = data.frame(array(dim = c(length(all.DE.genes), length(Analysis.names))))
    rownames(gene.vs.analysis.DE.table) = all.DE.genes
    colnames(gene.vs.analysis.DE.table) = Analysis.names
    gene.vs.analysis.DE.table[] = FALSE
    if(include.NP.diagrams){
      gene.vs.analysis.DE.table..NP = data.frame(array(dim = c(length(all.DE.genes..NP), length(Analysis.names))))
      rownames(gene.vs.analysis.DE.table..NP) = all.DE.genes..NP
      colnames(gene.vs.analysis.DE.table..NP) = Analysis.names
      gene.vs.analysis.DE.table..NP[] = FALSE
    }
    print(all.DE.genes..NP)

    cat(sprintf('[Create.Euler.Venn.diagrams]  %s - %s - Loading RDS data and constructing gene-centric table....\n', reported.analysis.group.name, threshold.name))
    for(n in 1:length(Analysis.names)){
      Analysis.Data = readRDS(RDS.file.names[n])
      PValue.col.name = Get.PValue.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part))
      FDR.col.name = Get.FDR.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part), look.for.paired.test.res = Analysis.Data$perform.classic.paired.test)
      LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
      LogCPM.col.name = Get.LogCPM.col.name(colnames(Analysis.Data$ResTable.gene.part))
      keep.by.PValue = Analysis.Data$ResTable.gene.part[,PValue.col.name] <= Venn.diagrams___PValue.higher.bound
      keep.by.FDR = Analysis.Data$ResTable.gene.part[,FDR.col.name] <= Venn.diagrams___FDR.higher.bound
      keep.by.LogFC = abs(Analysis.Data$ResTable.gene.part[,LogFC.col.name]) >= Venn.diagrams___abs.LogFC.lower.bound
      keep.by.LogCPM = Analysis.Data$ResTable.gene.part[,LogCPM.col.name] >= Venn.diagrams___LogCPM.lower.bound
      keep = keep.by.PValue & keep.by.FDR & keep.by.LogFC & keep.by.LogCPM
      keep.genes = rownames(Analysis.Data$ResTable.gene.part)[keep]
      gene.vs.analysis.DE.table[which(all.DE.genes %in% keep.genes), Analysis.names[n]] = TRUE

      if(include.NP.diagrams){
        NP.PValue.col.name = Get.NP.PValue.col.name(colnames(Analysis.Data$ResTable.gene.part))
        if(!is.na(NP.PValue.col.name)){
          keep.by.NP.PValue = Analysis.Data$ResTable.gene.part[,NP.PValue.col.name] <= Venn.diagrams___NP.PValue.higher.bound
          keep.genes..NP = rownames(Analysis.Data$ResTable.gene.part)[keep & keep.by.NP.PValue]
          gene.vs.analysis.DE.table..NP[which(all.DE.genes %in% keep.genes..NP), Analysis.names[n]] = TRUE
        }
      }
    }


    for(Venn.diagrams___type in Pars$Venn.diagrams___types){
      if((Venn.diagrams___type == 'Euler_NP' | Venn.diagrams___type == 'Venn_NP') & !include.NP.diagrams)  next
      cat(sprintf('[Create.Euler.Venn.diagrams]  %s - %s - Rendering %s diagrams....\n', reported.analysis.group.name, threshold.name, Venn.diagrams___type))
      cat(sprintf('%d genes passed "%s" threshold\n', nrow(gene.vs.analysis.DE.table), threshold.name))
      if(include.NP.diagrams)
        cat(sprintf('%d genes passed "%s" NP-threshold\n', nrow(gene.vs.analysis.DE.table..NP), threshold.name))
      if(Venn.diagrams___type == 'Euler' | Venn.diagrams___type == 'Venn'){
        GvAt = gene.vs.analysis.DE.table
      } else if(Venn.diagrams___type == 'Euler_NP' | Venn.diagrams___type == 'Venn_NP'){
        GvAt = gene.vs.analysis.DE.table..NP
      } else {
        stop(sprintf('Incorrect Venn.diagrams___type "%s"', Venn.diagrams___type))
      }
      Venn.diagrams___type


      venn.setup..list = vector(mode = 'list')
      intersected.taxons..list = vector(mode = 'list')

      # GvAt = as.data.frame(array(sample(c(TRUE,FALSE), nrow(tax.data)*ncol(tax.data), TRUE), dim = c(nrow(tax.data), ncol(tax.data))))
      # GvAt = GvAt[, 1:(ncol(GvAt)-1)]
      # print(dim(GvAt))
      
      if(ncol(GvAt) == 2){
        venn.setup..list$x1 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE))) %>% sum
        venn.setup..list$x2 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE))) %>% sum
        venn.setup..list$x1.x2 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE))) %>% sum


        intersected.taxons..list$x1 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE))) %>% rownames(GvAt)[.]
        
      } else if(ncol(GvAt) == 3){
        venn.setup..list$x1 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE))) %>% sum
        venn.setup..list$x2 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE))) %>% sum
        venn.setup..list$x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE))) %>% sum
        venn.setup..list$x1.x2 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE))) %>% sum
        venn.setup..list$x1.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE))) %>% sum
        venn.setup..list$x2.x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE))) %>% sum
        venn.setup..list$x1.x2.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE))) %>% sum


        intersected.taxons..list$x1 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]

      } else if(ncol(GvAt) == 4){
        venn.setup..list$x1 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE))) %>% sum
        venn.setup..list$x2 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE))) %>% sum
        venn.setup..list$x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE))) %>% sum
        venn.setup..list$x1.x2 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE))) %>% sum
        venn.setup..list$x1.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE))) %>% sum
        venn.setup..list$x2.x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE))) %>% sum
        venn.setup..list$x1.x2.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE))) %>% sum

        venn.setup..list$x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE))) %>% sum
        venn.setup..list$x1.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE))) %>% sum
        venn.setup..list$x2.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE))) %>% sum
        venn.setup..list$x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE))) %>% sum
        venn.setup..list$x1.x2.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE))) %>% sum
        venn.setup..list$x1.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE))) %>% sum
        venn.setup..list$x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE))) %>% sum
        venn.setup..list$x1.x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE))) %>% sum


        intersected.taxons..list$x1 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]

        intersected.taxons..list$x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]


      } else if(ncol(GvAt) == 5){

        venn.setup..list$x1 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, FALSE))) %>% sum
        venn.setup..list$x2 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, FALSE))) %>% sum
        venn.setup..list$x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, FALSE))) %>% sum
        venn.setup..list$x1.x2 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, FALSE))) %>% sum
        venn.setup..list$x1.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, FALSE))) %>% sum
        venn.setup..list$x2.x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, FALSE))) %>% sum
        venn.setup..list$x1.x2.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, FALSE))) %>% sum

        venn.setup..list$x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, FALSE))) %>% sum
        venn.setup..list$x1.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, FALSE))) %>% sum
        venn.setup..list$x2.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, FALSE))) %>% sum
        venn.setup..list$x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, FALSE))) %>% sum
        venn.setup..list$x1.x2.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, FALSE))) %>% sum
        venn.setup..list$x1.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, FALSE))) %>% sum
        venn.setup..list$x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, FALSE))) %>% sum
        venn.setup..list$x1.x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, FALSE))) %>% sum


        venn.setup..list$x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, FALSE, TRUE))) %>% sum
        venn.setup..list$x1.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, TRUE))) %>% sum
        venn.setup..list$x2.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, TRUE))) %>% sum
        venn.setup..list$x3.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, TRUE))) %>% sum
        venn.setup..list$x1.x2.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, TRUE))) %>% sum
        venn.setup..list$x1.x3.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, TRUE))) %>% sum
        venn.setup..list$x2.x3.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, TRUE))) %>% sum
        venn.setup..list$x1.x2.x3.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, TRUE))) %>% sum

        venn.setup..list$x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, TRUE))) %>% sum
        venn.setup..list$x1.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, TRUE))) %>% sum
        venn.setup..list$x2.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, TRUE))) %>% sum
        venn.setup..list$x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, TRUE))) %>% sum
        venn.setup..list$x1.x2.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, TRUE))) %>% sum
        venn.setup..list$x1.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, TRUE))) %>% sum
        venn.setup..list$x2.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, TRUE))) %>% sum
        venn.setup..list$x1.x2.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, TRUE))) %>% sum



        intersected.taxons..list$x1 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]

        intersected.taxons..list$x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]


        intersected.taxons..list$x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]

        intersected.taxons..list$x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]


      } else if(ncol(GvAt) == 6 | ncol(GvAt) == 7){

        for(nc in 1:ncol(GvAt)){
          venn.setup..list[[colnames(GvAt)[nc]]] = rownames(GvAt)[GvAt[,nc]]
        }

        # venn.setup..list$x1 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x2 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x1.x2 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x1.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x2.x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x1.x2.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))) %>% sum

        # venn.setup..list$x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x1.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x2.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x1.x2.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x1.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE))) %>% sum
        # venn.setup..list$x1.x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))) %>% sum


        # venn.setup..list$x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x1.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x2.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x3.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x1.x2.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x1.x3.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x2.x3.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x1.x2.x3.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE))) %>% sum

        # venn.setup..list$x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x1.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x2.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x1.x2.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x1.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x2.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE))) %>% sum
        # venn.setup..list$x1.x2.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))) %>% sum



        # venn.setup..list$x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x1.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x2.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x3.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x1.x2.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x1.x3.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x2.x3.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x1.x2.x3.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE))) %>% sum

        # venn.setup..list$x4.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x1.x4.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x2.x4.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x3.x4.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x1.x2.x4.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x1.x3.x4.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x2.x3.x4.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE))) %>% sum
        # venn.setup..list$x1.x2.x3.x4.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE))) %>% sum


        # venn.setup..list$x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x1.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x2.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x3.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x1.x2.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x1.x3.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x2.x3.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x1.x2.x3x.5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE))) %>% sum

        # venn.setup..list$x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x1.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x2.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x3.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x1.x2.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x1.x3.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x2.x3.x4.x5x.6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE))) %>% sum
        # venn.setup..list$x1.x2.x3.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE))) %>% sum





        intersected.taxons..list$x1 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))) %>% rownames(GvAt)[.]

        intersected.taxons..list$x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3.x4 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))) %>% rownames(GvAt)[.]


        intersected.taxons..list$x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE))) %>% rownames(GvAt)[.]

        intersected.taxons..list$x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3.x4.x5 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))) %>% rownames(GvAt)[.]



        intersected.taxons..list$x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE))) %>% rownames(GvAt)[.]

        intersected.taxons..list$x4.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x4.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x4.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3.x4.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x4.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3.x4.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3.x4.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3.x4.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE))) %>% rownames(GvAt)[.]


        intersected.taxons..list$x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3x.5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE))) %>% rownames(GvAt)[.]

        intersected.taxons..list$x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x3.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x3.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x2.x3.x4.x5x.6 = apply(GvAt, 1, function(x) all(x == c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]
        intersected.taxons..list$x1.x2.x3.x4.x5.x6 = apply(GvAt, 1, function(x) all(x == c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE))) %>% rownames(GvAt)[.]

      } else {
        msg = sprintf('[Create.Euler.Venn.diagrams]  %s - Too large number of groups in Euler/Venn diagram - %d groups. Only 2-7 groups are suppoted now\n', reported.analysis.group.name, ncol(GvAt))
        if(!bypass..too.large.comparisons)  stop(msg)
        cat(msg)
        next
      }

      print(GvAt)
      print(venn.setup..list)
      venn.setup = unlist(venn.setup..list)
      nn = gsub(".", "&", names(venn.setup), fixed = TRUE)
      nn = gsub("x", "JKsid5bs", nn, fixed = TRUE)
      for(xn in 1:ncol(GvAt)){
        nn = gsub(sprintf("JKsid5bs%d", xn), colnames(GvAt)[xn], nn, fixed = TRUE)
      }
      names(venn.setup) = nn


      nn = gsub(".", " + ", names(intersected.taxons..list), fixed = TRUE)
      nn = gsub("x", "JKsid5bs", nn, fixed = TRUE)
      for(xn in 1:ncol(GvAt)){
        nn = gsub(sprintf("JKsid5bs%d", xn), colnames(GvAt)[xn], nn, fixed = TRUE)
      }
      names(intersected.taxons..list) = nn        

      
      if(ncol(GvAt) <= 5){
        if(Venn.diagrams___type == 'Venn' | Venn.diagrams___type == 'Venn_NP'){
          g = tryCatch(expr = {
            fit.venn <- eulerr::venn(venn.setup)
            plot(fit.venn, legend = TRUE, labels = TRUE, quantities = TRUE, fills = suppressWarnings(brewer.pal(ncol(GvAt), 'Spectral')))
          }, error = function (err){
            cat(sprintf('[Create.Euler.Venn.diagrams]  %s - %s - creating %s diagrams has failed (1) with error: "%s"\n', reported.analysis.group.name, threshold.name, Venn.diagrams___type, trimws(err)))
            # stop('ddd')
            NULL
          })
          if(!is.null(g))
            if(Venn.diagrams___type %in% c('Euler', 'Venn')){
              ggsave(filename = sprintf('%s/Venn, %s.png', output.dir, threshold.name), plot = g, dpi = img.resolution)
            } else {
              ggsave(filename = sprintf('%s/Venn (NP), %s.png', output.dir, threshold.name), plot = g, dpi = img.resolution)
            }
          # dev.off()
        }

        if(Venn.diagrams___type == 'Euler' | Venn.diagrams___type == 'Euler_NP'){
          g = tryCatch(expr = {
            fit.euler <- euler(venn.setup, shape = "ellipse")
            g = plot(fit.euler, legend = TRUE, labels = TRUE, quantities = TRUE, fills = suppressWarnings(brewer.pal(ncol(GvAt), 'Spectral')))
          }, error = function (err){
            cat(sprintf('[Create.Euler.Venn.diagrams]  %s - %s - creating %s diagrams has failed (2) with error: "%s"\n', reported.analysis.group.name, threshold.name, Venn.diagrams___type, trimws(err)))
            # stop('ddd')
            NULL
          })
          if(!is.null(g))
            if(Venn.diagrams___type %in% c('Euler', 'Venn')){
              ggsave(filename = sprintf('%s/Euler, %s.png', output.dir, threshold.name), plot = g, dpi = img.resolution)
            } else {
              ggsave(filename = sprintf('%s/Euler (NP), %s.png', output.dir, threshold.name), plot = g, dpi = img.resolution)
            }
        }
      } else if(Venn.diagrams___type == 'Venn' | Venn.diagrams___type == 'Venn_NP'){
        tryCatch(expr = {
          if(Venn.diagrams___type %in% c('Euler', 'Venn')){
            png(filename = sprintf('%s/Venn, %s.png', output.dir, threshold.name), units="in", width=7, height=7, pointsize=16, res = img.resolution)
          } else {
            png(filename = sprintf('%s/Venn (NP), %s.png', output.dir, threshold.name), units="in", width=7, height=7, pointsize=16, res = img.resolution)
          }
          venn::venn(venn.setup..list, ilabels = TRUE, zcolor = "style")
          dev.off()          
        }, error = function (err){
          cat(sprintf('[Create.Euler.Venn.diagrams]  %s - %s - creating %s diagrams has failed (3) with error: "%s"\n', reported.analysis.group.name, threshold.name, Venn.diagrams___type, trimws(err)))
          # stop('ddd')
          NULL
        })
      } else {
        cat(sprintf('[Create.Euler.Venn.diagrams]  %s - Only 5 groups are supported for Euler diagrams. However, you specified %d groups. This diagram will be bypassed\n', reported.analysis.group.name, ncol(GvAt)))
      }
    
      if(write.intersections){
        max.len = max(sapply(intersected.taxons..list, length))
        intersected.taxons..df = as.data.frame(lapply(intersected.taxons..list, function(x) c(x, rep("", max.len - length(x)))))

        if(Venn.diagrams___type %in% c('Euler', 'Venn')){
          tsv.file.name = sprintf('%s/Venn Intersections, %s.tsv', output.dir, threshold.name)
          Excel.sheet.name = threshold.name
        } else {
          tsv.file.name = sprintf('%s/Venn Intersections NP, %s.tsv', output.dir, threshold.name)
          Excel.sheet.name = sprintf("%s (NP)", threshold.name)
        }
        if(!(Excel.sheet.name %in% all.Excel.sheet.names)){
          write.table(x = intersected.taxons..df, file = tsv.file.name, sep = '\t', quote = FALSE)
          all.tsv.file.names %<>% c(., tsv.file.name)
          all.Excel.sheet.names %<>% c(., Excel.sheet.name)
        }
      }  
    }
  }
  if(write.intersections){
    Excel.file.name = sprintf('%s/Venn intersections.xlsx', output.dir)
    CL = sprintf('%s %s/Any.2.Excel.py --in %s --one-book yes --sheet-names %s --out-excel "%s" --disable-coltype-autoassign yes --entry-type N',
               Pars$python.bin, Pars$suppl.data.dir,
               sprintf('"%s"', all.tsv.file.names) %>% paste(., collapse = ' '),
               sprintf('"%s"', all.Excel.sheet.names) %>% paste(., collapse = ' '),
               Excel.file.name
    )
    exit.code = system(CL)
    if(exit.code > 0){
      cat(sprintf('\n\nCommand:\n      %s\nFAILED with exit code %d\n\n', CL, exit.code))
      stop(sprintf('Excel generation for %s returned exit code %d', Excel.file.name, exit.code))
    }
    if(length(all.tsv.file.names) > 0)
      invisible(suppressWarnings(file.remove(all.tsv.file.names)))
  }
}

Collect.inDetails.GLM.results = function(Startup.Data, GLM.models = NULL, forced.parameters = NULL,
                                           out.dir = NULL,
                                           inDetails.entries = NULL, add.biotypes = NULL){ # c('lncRNA + lincRNA + antisense_RNA + macro_lincRNA')
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(GLM.models)) GLM.models = Pars$Models.to.Test
  if(is.null(inDetails.entries))  inDetails.entries = Startup.Data$inDetails.entries
  
  if(length(inDetails.entries) == 0){
    cat('\nNo inDetails groups are declared - see "GO terms (inDetails)" sheet\n')
    return(invisible())
  }
  
  if(is.null(add.biotypes))  add.biotypes = Pars$inDetails.additional.biotypes
  all.group.names = c(names(inDetails.entries), names(add.biotypes))
  all.group.Name_s = gsub("[^[:alnum:] ]", "", all.group.names)
  if(is.null(out.dir))  out.dir = sprintf('%s/inDetails DE summary', Pars$results.dir)
  dir.create(out.dir, showWarnings = FALSE)
  
  present.src.files = 0
  absent.src.files = 0
  for(group.Name_s in all.group.Name_s){
    for(GLM.model in GLM.models){
      Analysis.name  = Cleanup.Model.Name(GLM.model)
      GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
      inDetails.main.wd = sprintf("%s/inDetails", GLM.working.dir)
      
      if(Pars$Create.CPM.profiles.in.Excel.results){
        inDetails.Excel.DE.file.name = sprintf("%s/%s/[%s], DE results, CPM profiles.xlsx",
                                        inDetails.main.wd, group.Name_s, group.Name_s)
        inDetails.Excel.DE.file.name = Verify.path(inDetails.Excel.DE.file.name)
        dest.inDetails.Excel.DE.file.name = sprintf("%s/[%s]  [DE results, CPM profiles]    %s.xlsx", out.dir, group.Name_s, Analysis.name)

        if(file.exists(inDetails.Excel.DE.file.name)){
          tmp = file.copy(inDetails.Excel.DE.file.name, dest.inDetails.Excel.DE.file.name)
          present.src.files = c(present.src.files, inDetails.Excel.DE.file.name)
        } else absent.src.files = c(absent.src.files, inDetails.Excel.DE.file.name)
      }
      
      if(Pars$Create.rel.LogCPM.profiles.in.Excel.results){
        inDetails.Excel.DE.file.name = sprintf("%s/%s/[%s], DE results, Log.rel.CPM profiles.xlsx",
                                        inDetails.main.wd, group.Name_s, group.Name_s)
        inDetails.Excel.DE.file.name = Verify.path(inDetails.Excel.DE.file.name)
        dest.inDetails.Excel.DE.file.name = sprintf("%s/[%s]  [DE results, Log.rel.CPM profiles]    %s.xlsx", out.dir, group.Name_s, Analysis.name)

        if(file.exists(inDetails.Excel.DE.file.name)){
          tmp = file.copy(inDetails.Excel.DE.file.name, dest.inDetails.Excel.DE.file.name)
          present.src.files = c(present.src.files, inDetails.Excel.DE.file.name)
        } else absent.src.files = c(absent.src.files, inDetails.Excel.DE.file.name)
      }
    }
  }
  if(length(absent.src.files) > 0)   cat(sprintf('The following inDetails Excel results do not exist (total %d of %d): %s\n', length(absent.src.files), length(absent.src.files) + length(present.src.files), toString(absent.src.files)))
}




Summarize.inDetails.GLM.results = function(Startup.Data, GLM.models = NULL, forced.parameters = NULL,
                                           out.dir = NULL,
                                           inDetails.entries = NULL, add.biotypes = NULL){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  
  if(is.null(GLM.models)) GLM.models = Pars$Models.to.Test
  if(is.null(inDetails.entries))  inDetails.entries = Startup.Data$inDetails.entries
  
  if(length(inDetails.entries) == 0){
    cat('\nNo inDetails groups are declared - see "GO terms (inDetails)" sheet\n')
    return(invisible())
  }
  
  if(is.null(add.biotypes))  add.biotypes = Pars$inDetails.additional.biotypes
  all.group.names = c(names(inDetails.entries), names(add.biotypes))
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
      cat(sprintf('[%s]: %d of %d DE analysis results (*.tsv) are not found\n', group.Name_s, length(keep) - sum(keep), length(keep)))
      next
    } else {
      cat(sprintf('[%s]: all %d DE analysis results are present\n', group.Name_s, length(keep)))
    }
    
    all.inDetails.DE.file.names = all.inDetails.DE.file.names[keep]
    # all.inDetails.DE.file.names = gsub("/", "\\", all.inDetails.DE.file.names, fixed = TRUE)
    
    out.file.name = Verify.path(sprintf("%s/%s - Summary expression info.xlsx", out.dir, group.Name_s))
    
    ### running non-merging txt > Excel conversion
    # CL = sprintf('%s "%s/DE.results.to.Excel.py"',Pars$python.bin,Pars$suppl.data.dir)
    # # CL = gsub('/','\\',CL,fixed = T)
    # CL = paste(c(CL,sprintf('"%s"',out.file.name),sprintf('"%s"',all.inDetails.DE.file.names)),collapse = ' ')
    # system(CL)
    # if (!file.exists(out.file.name)){
    #   stop(sprintf('Excel worksheet generation FAILED [non-joint], %s',out.file.name)) }

        
    ### running merging txt > Excel conversion
    out.file.name.ff = Verify.path(sprintf("%s, joint.xlsx",gsub('.{5}$', '', out.file.name)))
    out.file.name.lite = Verify.path(sprintf("%s, joint lite.xlsx",gsub('.{5}$', '', out.file.name)))
    out.file.name.lite.logfconly = Verify.path(sprintf("%s, joint lite, LogFC table.xlsx",gsub('.{5}$', '', out.file.name)))
    CL_base = sprintf('%s "%s/DE.results.to.Excel.joint.py"',Pars$python.bin,Pars$suppl.data.dir)
    # CL_base = gsub('/','\\',CL_base,fixed = T)
    CL.ff = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs.tsv" --lr no --score no --logcpm yes --logcpm-common yes --spearmanp no --pearsonr no --pearsonp no --exclude-genes-without-de-info yes',
                    CL_base, paste(sprintf('"%s"',all.inDetails.DE.file.names),collapse=' '), out.file.name.ff, Pars$results.dir)
    CL.lite = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs.tsv" --lr no --score no --logcpm yes --logcpm-common yes --spearmanp no --spearmanr no --fdr no --pearsonr no --pearsonp no --simple-logfc-format-mode yes --exclude-genes-without-de-info yes',
                      CL_base, paste(sprintf('"%s"',all.inDetails.DE.file.names),collapse=' '), out.file.name.lite, Pars$results.dir)
    CL.lite.logfconly = sprintf('%s -i %s -o "%s" --cpm-table "%s/CPMs.tsv" --p no --fdr no --lr no --score no --logcpm no --logcpm-common yes --spearmanr no --spearmanp no --pearsonr no --pearsonp no --simple-logfc-format-mode yes --exclude-genes-without-de-info yes',
                                CL_base, paste(sprintf('"%s"',all.inDetails.DE.file.names),collapse=' '), out.file.name.lite.logfconly, Pars$results.dir)
    
    cat(CL.ff)
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
                                                      forced.parameters = NULL,
                                                      out.dir = NULL, out.file.prefix = '',
                                                      include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                                                      include.with.spaces = FALSE){
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
    if(include.GLM.model.name.in.Result.names){
      f = sprintf("%s/~ %s, results/GO Terms (%s) - topGO - custom terms DE/%s, GO-DE info, logCPM.gt.%g, P.lt.%g.tsv",
                  Pars$results.dir, Analysis.name, GO.type, Analysis.name,
                  Pars$topGO.Expression.Profiles___Min.gene.logCPM.for.summary.across.models, Pars$topGO.Expression.Profiles___Max.gene.PValue.for.summary.across.models)
    } else {
      f = sprintf("%s/~ %s, results/GO Terms (%s) - topGO - custom terms DE/GO-DE info, logCPM.gt.%g, P.lt.%g.tsv",
                  Pars$results.dir, Analysis.name, GO.type,
                  Pars$topGO.Expression.Profiles___Min.gene.logCPM.for.summary.across.models, Pars$topGO.Expression.Profiles___Max.gene.PValue.for.summary.across.models)
    }

    # f = gsub("/","\\", f, fixed = TRUE)
    f = Verify.path(f)
    file.names = append(file.names, f)
  }

  if (sum.mod(!file.exists(file.names))>0){
    cat(sprintf('[Summarize.topGO.Expression.Profiles.Custom] Total %d files with custom GO-terms expression profiles do not exist: %s\n',sum.mod(!file.exists(file.names)),toString(file.names[!file.exists(file.names)])))
    GLM.models = GLM.models[file.exists(file.names)]
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
               '--models', sprintf('"%s"', GLM.models),
               '--maximal-genes-count',as.character(Pars$topGO.Expression.Profiles___Maximal.genes.in.terms.to.visualize)),collapse = ' ')
  # CL = gsub('/','\\',CL,fixed = TRUE)
  code = system(CL)
  if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
  if (!file.exists(Verify.path(out.file.name)))    stop(sprintf('[Summarize.topGO.Expression.Profiles.Custom] Excel worksheet generation FAILED'))
  
  if(include.with.spaces){
    CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces yes ',Pars$python.bin,Pars$suppl.data.dir)
    out.file.name = sprintf('%s/%sSummary topGO (%s) DE profiles (with spaces).xlsx', out.dir, out.file.prefix, GO.type)
    CL = paste(c(CL,sprintf('--out-excel "%s"',Verify.path(out.file.name)),'--in',sprintf('"%s"',file.names),
                 '--models', sprintf('"%s"', GLM.models),
                 '--maximal-genes-count',as.character(Pars$topGO.Expression.Profiles___Maximal.genes.in.terms.to.visualize)),collapse = ' ')
    # CL = gsub('/','\\',CL,fixed = TRUE)
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    if (!file.exists(Verify.path(out.file.name)))    stop(sprintf('[Summarize.topGO.Expression.Profiles.Custom] Excel worksheet generation FAILED'))
  }
}


Summarize.clusterProfiler.Expression.Profiles.Custom = function(Startup.Data, GLM.models = NULL, database = 'GO', printed.name = 'GO Terms',
                                                      forced.parameters = NULL, out.dir = NULL, out.file.prefix = "",
                                                      include.GLM.model.name.in.Result.names = FALSE, include.with.spaces = FALSE){
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
    if(include.GLM.model.name.in.Result.names){
      f = sprintf("%s/~ %s, results/%s - custom terms DE/%s, %s-DE info, logCPM.gt.%g, P.lt.%g.tsv", Pars$results.dir, Analysis.name, printed.name, Analysis.name, database, logCPM.thr, PValue.thr)
    } else {
      f = sprintf("%s/~ %s, results/%s - custom terms DE/%s-DE info, logCPM.gt.%g, P.lt.%g.tsv", Pars$results.dir, Analysis.name, printed.name, database, logCPM.thr, PValue.thr)
    }

    # f = gsub("/","\\", f, fixed = TRUE)
    f = Verify.path(f)
    file.names = append(file.names, f)
  }
  if (sum.mod(!file.exists(file.names))>0){
    cat(sprintf('\n[Summarize.clusterProfiler.Expression.Profiles.Custom]  Total %d files with custom GO-terms expression profiles do not exist: %s\n',sum.mod(!file.exists(file.names)),toString(file.names[!file.exists(file.names)])))
  }
  
  if(is.null(out.dir)) out.dir = Pars$results.dir
  
  GLM.models = GLM.models[file.exists(file.names)]
  file.names = file.names[file.exists(file.names)]
  if(length(file.names) == 0){
    cat('No files found\n')
    return(invisible())
  }
  
  CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces no --database %s ',Pars$python.bin,Pars$suppl.data.dir,database)
  # CL = gsub('/','\\',CL,fixed = TRUE)
  out.file.name = sprintf('%s/%sSummary %s DE profiles.xlsx', out.dir, out.file.prefix, database)
  CL = paste(c(CL,sprintf('--out-excel "%s"',Verify.path(out.file.name)),'--in',sprintf('"%s"',file.names),
               '--models', sprintf('"%s"', GLM.models),
               '--maximal-genes-count', as.character(max.Genes)),collapse = ' ')
  code = system(CL)
  if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
  if (!file.exists(Verify.path(out.file.name)))    stop(sprintf('Excel worksheet generation FAILED'))
  
  if(include.with.spaces){
    CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces yes --database %s ',Pars$python.bin,Pars$suppl.data.dir,database)
    # CL = gsub('/','\\',CL,fixed = TRUE)
    out.file.name = sprintf('%s/%sSummary %s DE profiles (with spaces).xlsx', out.dir, out.file.prefix, database)
    CL = paste(c(CL,sprintf('--out-excel "%s"',Verify.path(out.file.name)),'--in',sprintf('"%s"',file.names),
                 '--models', sprintf('"%s"', GLM.models),
                 '--maximal-genes-count', as.character(max.Genes)),collapse = ' ')
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    if (!file.exists(Verify.path(out.file.name)))    stop(sprintf('Excel worksheet generation FAILED'))
  }
}


Summarize.GO.Expression.Profiles.Custom = function(Startup.Data, GLM.models = NULL, forced.parameters = NULL,
                                                   out.dir = NULL, out.file.prefix = '', include.GLM.model.name.in.Result.names = FALSE,
                                                   include.with.spaces = FALSE){
  Summarize.clusterProfiler.Expression.Profiles.Custom(Startup.Data = Startup.Data, GLM.models = GLM.models,
                                                       database = 'GO', printed.name = 'GO Terms', out.dir = out.dir, out.file.prefix = out.file.prefix,
                                                       forced.parameters = forced.parameters,
                                                       include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                                                       include.with.spaces = include.with.spaces)
}

Summarize.KEGG.Expression.Profiles.Custom = function(Startup.Data, GLM.models = NULL, forced.parameters = NULL,
                                                     out.dir = NULL, out.file.prefix = "", include.GLM.model.name.in.Result.names = FALSE,
                                                     include.with.spaces = FALSE){
  Summarize.clusterProfiler.Expression.Profiles.Custom(Startup.Data = Startup.Data, GLM.models = GLM.models,
                                                       database = 'KEGG', printed.name = 'KEGG Pathways', out.dir = out.dir, out.file.prefix = out.file.prefix,
                                                       forced.parameters = forced.parameters,
                                                       include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                                                       include.with.spaces = include.with.spaces)
}

Summarize.Reactome.Expression.Profiles.Custom = function(Startup.Data, GLM.models = NULL, forced.parameters = NULL,
                                                         out.dir = NULL, out.file.prefix = '', include.GLM.model.name.in.Result.names = FALSE,
                                                         include.with.spaces = FALSE){
  Summarize.clusterProfiler.Expression.Profiles.Custom(Startup.Data = Startup.Data, GLM.models = GLM.models,
                                                       database = 'Reactome', printed.name = 'Reactome Pathways', out.dir = out.dir, out.file.prefix = out.file.prefix,
                                                       forced.parameters = forced.parameters,
                                                       include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                                                       include.with.spaces = include.with.spaces)
}



Summarize.KEGG.DE.info = function(Startup.Data, GLM.models = NULL, forced.parameters = NULL,
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
  
  if(length(existing.all.src.files) == 0) return(invisible(NULL))
  
  
  #setwd(KEGG.de.summary.dir)
  CL=sprintf('%s "%s/Integrate.KEGG.pathway.summaries.py" --input-file-names "%s" --out-excel "%s/%sKEGG DE summary.xlsx" --bar-max-value-override 100',
             Pars$python.bin, Pars$suppl.data.dir, paste(existing.all.src.files, collapse = ';'), out.dir, out.file.prefix)
  code = system(CL)
  if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
  
  #file.copy(from = 'KEGG.DE.summary.xlsx', to = sprintf('%s/KEGG.DE.summary.xlsx',Pars$results.dir),overwrite = TRUE)
  CL=sprintf('%s "%s/Integrate.KEGG.pathway.summaries.py" --input-file-names "%s" --out-excel "%s/%sKEGG DE summary combined.xlsx" --bar-max-value-override 100 --combine-node-count-and-perc yes',
             Pars$python.bin, Pars$suppl.data.dir, paste(existing.all.src.files, collapse = ';'), out.dir, out.file.prefix)
  code = system(CL)
  if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
  #file.copy(from = 'KEGG.DE.summary.combined.xlsx', to = sprintf('%s/KEGG.DE.summary.combined.xlsx',Pars$results.dir),overwrite = TRUE)
  #setwd(Pars$results.dir)
  
  #"KEGG.DE.summary.to.Excel.py"
}



Process.ext.DE.data = function(Startup.Data, ext.DE.data.file, custom.analysis.name = NULL, forced.parameters = NULL){
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
  
  # Analysis.Data$step.hashmd5 = digest(ext.DE.data.file)
  
  ResTable.gene.part = Startup.Data$ext.DE.data[[ext.DE.data.file]][,c('gene','score','logfc','pvalue','fdr','logcpm')]
  colnames(ResTable.gene.part) = c('Symbol','Score','logFC','PValue','FDR','logCPM')
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
    colnames(Mart.table) = c("Symbol","Biotype","description","RefSeq_Summary")
  }  else{
    colnames(Mart.table) = c("Symbol","Biotype","description")
  }

  ResTable.gene.part = cbind(Mart.table,ResTable.gene.part[,-1])

  output.file.name = sprintf("%s_external_data.txt",Analysis.Data$Analysis.name)
  write.table.mod(ResTable.gene.part, file = output.file.name,sep='\t',na = '',col.names = colnames(ResTable.gene.part))
  
  if (Pars$Create.Excel.results){
    CL = sprintf('%s "%s/DE.results.to.Excel.py"',Pars$python.bin,Pars$suppl.data.dir)
    # CL = gsub('/','\\',CL,fixed = T)
    CL = paste(c(CL,sprintf('"%s, DE results.xlsx"',Analysis.Data$Analysis.name),sprintf('"%s"',c(output.file.name))),collapse = ' ')
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    if (!file.exists(sprintf('%s, DE results.xlsx',Analysis.Data$Analysis.name))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (DE results) for ~%s FAILED',Analysis.Data$Analysis.name)) }
    
  }
  
  Analysis.Data$ResTable.gene.part = ResTable.gene.part
  saveRDS(ResTable.gene.part,file='DE.results.rds')

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
                              forced.parameters = NULL, out.dir = NULL,
                              forced.Analysis.name = NULL, stats.file.name = NULL, sparklines.axis.limits = 2,
                              include.GLM.model.name.in.Result.names = FALSE){
  
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
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
      all.stats = c('min.p.value_by.term.id.UP','min.FDR_by.term.id.UP','min.p.value_by.term.id.DOWN','min.FDR_by.term.id.DOWN',
                    'enrichment.score_by.term.id.final',
                    'enrichment.score_by.term.id.final.UP','enrichment.score_by.term.id.final.DOWN')
      if(database == 'GO'){
        Enrich.info.split =  data.frame(array(dim=c(0,length(all.stats))))
        colnames(Enrich.info.split) = all.stats
        for (GO.type in c('BP', 'MF', 'CC')){
          file.name = sprintf("%s/GO Terms (%s) - %s Enrichment Analysis/GO.%s %s Enrich. stats.rds",
                              Analysis.Data$current.GLM.results.dir, GO.type, enrich.test.type, GO.type, enrich.test.type)
          # print(file.name)
          # print(file.exists(file.name))
          # stop('jjjjj')
          if(file.exists(file.name)){
            Enrich.info.split = rbind(Enrich.info.split, readRDS(file = file.name))
          }
        }
        if(dim(Enrich.info.split)[1] > 0)  Enrich.info.split__is.found = TRUE
      } else {
        if (database == 'KEGG'){
          file.name = sprintf('%s/KEGG Pathways - %s Enrichment Analysis/KEGG %s Enrich. stats.rds',
                              Analysis.Data$current.GLM.results.dir, enrich.test.type, enrich.test.type)
        } else if (database == 'Reactome'){
          file.name = sprintf('%s/Reactome Pathways - %s Enrichment Analysis/Reactome %s Enrich. stats.rds',
                              Analysis.Data$current.GLM.results.dir, enrich.test.type, enrich.test.type)
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
    Enrich.info = read.table(file = stats.file.name,header = TRUE,sep = '\t')
    
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
      all.LogFCs.per.DB.entry = Get.LogFCs.for.Ensembl.IDs.vector(Analysis.Data, Pars, Genes.per.DB.entry, PValue.thr = PValue.thr, logCPM.thr = logCPM.thr, sort.logFCs = TRUE)
      #all.LogFCs.per.DB.entry[['R-HSA-977441']]
      #all.LogFCs.per.DB.entry[['R-HSA-5368287']]
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
      
      DE.profile.table = DE.profile.table[!is.na(all.LogFCs.per.DB.entry),]
      DE.profile.table = DE.profile.table[which(sapply(all.LogFCs.per.DB.entry, length) > 0),]
      
      if(include.GLM.model.name.in.Result.names){
        file.name = Verify.path(sprintf('%s, %s-DE info, logCPM.gt.%g, P.lt.%g.tsv', Analysis.name, database, logCPM.thr, PValue.thr))
      } else {
        file.name = Verify.path(sprintf('%s-DE info, logCPM.gt.%g, P.lt.%g.tsv', database, logCPM.thr, PValue.thr))
      }

      write.table.mod(DE.profile.table, file = file.name ,sep='\t',quote = FALSE)
      output.file.names = c(output.file.names, file.name)
    }
  }
  
  if (Pars$Create.Excel.results){
    CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces no --database %s',Pars$python.bin,Pars$suppl.data.dir,database)
    # CL = gsub('/','\\',CL,fixed = T)
    if(include.GLM.model.name.in.Result.names){
      base.res.file.name = sprintf('%s, %s %s DE profiles.xlsx', Analysis.name, printed.name, profile.type)
    } else {
      base.res.file.name = sprintf('%s %s DE profiles.xlsx', printed.name, profile.type)
    }
    CL = paste(c(CL,sprintf('--out-excel "%s"', Verify.path(base.res.file.name)),'--in',sprintf('"%s"',output.file.names),
                 '--models', sprintf('"%s"', rep(GLM.model, length(output.file.names))),
                 '--maximal-genes-count',as.character(Pars$clusterProfiler.Expression.Profiles___Maximal.genes.in.terms.to.visualize),
                 sprintf('--sparklines-axis-limit %g',sparklines.axis.limits)),
               collapse = ' ')
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    if (!file.exists(Verify.path(base.res.file.name))){
      #file.remove(output.file.names)
      message(sprintf('Excel worksheet generation (%s-centric DE profiles) for ~%s FAILED', database, Analysis.name))
    } else {
      file.copy(Verify.path(base.res.file.name), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, base.res.file.name)))
    }
    
    CL = sprintf('%s "%s/GO.expression.to.Excel.py" --insert-spaces yes --database %s',Pars$python.bin,Pars$suppl.data.dir,database)
    # CL = gsub('/','\\',CL,fixed = T)

    if(include.GLM.model.name.in.Result.names){
      base.res.file.name = sprintf('%s, %s %s DE profiles (with spaces).xlsx', Analysis.name, printed.name, profile.type)
    } else {
      base.res.file.name = sprintf('%s %s DE profiles (with spaces).xlsx', printed.name, profile.type)
    }
    CL = paste(c(CL,sprintf('--out-excel "%s"', Verify.path(base.res.file.name)),'--in',sprintf('"%s"',output.file.names),
                 '--models', sprintf('"%s"', rep(GLM.model, length(output.file.names))),
                 '--maximal-genes-count',as.character(Pars$clusterProfiler.Expression.Profiles___Maximal.genes.in.terms.to.visualize),
                 sprintf('--sparklines-axis-limit %g',sparklines.axis.limits)),collapse = ' ')
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
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
  forced.parameters = NULL, out.dir = NULL,
  forced.Analysis.name = NULL, stats.file.name = NULL, sparklines.axis.limits = 2, include.GLM.model.name.in.Result.names = FALSE){
    
    DB.entries = Startup.Data$GO.clusterProfiler.Expression.Profiles___Custom.DB.entries
    clusterProfiler.Expression.Profiles(Startup.Data, DB.entries, Analysis.Data = Analysis.Data, GLM.model = GLM.model,
                                        GLM.working.dir = GLM.working.dir, Analysis.Data.RDS.file = Analysis.Data.RDS.file,
                                        database = 'GO', printed.name = 'GO Terms', profile.type = 'custom',
                                        forced.parameters = forced.parameters, out.dir = out.dir,
                                        forced.Analysis.name = forced.Analysis.name, stats.file.name = stats.file.name,
                                        sparklines.axis.limits = sparklines.axis.limits,
                                        include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
}


KEGG.Expression.Profiles.Custom =  function(
  Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
  forced.parameters = NULL, out.dir = NULL,
  forced.Analysis.name = NULL, stats.file.name = NULL, sparklines.axis.limits = 2, include.GLM.model.name.in.Result.names = FALSE){
  
  DB.entries = Startup.Data$KEGG.clusterProfiler.Expression.Profiles___Custom.DB.entries
  clusterProfiler.Expression.Profiles(Startup.Data, DB.entries, Analysis.Data = Analysis.Data, GLM.model = GLM.model,
                                      GLM.working.dir = GLM.working.dir, Analysis.Data.RDS.file = Analysis.Data.RDS.file,
                                      database = 'KEGG', printed.name = 'KEGG Pathways', profile.type = 'custom',
                                      forced.parameters = forced.parameters, out.dir = out.dir,
                                      forced.Analysis.name = forced.Analysis.name, stats.file.name = stats.file.name,
                                      sparklines.axis.limits = sparklines.axis.limits,
                                      include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
}


Reactome.Expression.Profiles.Custom =  function(
  Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
  forced.parameters = NULL, out.dir = NULL,
  forced.Analysis.name = NULL, stats.file.name = NULL, sparklines.axis.limits = 2, include.GLM.model.name.in.Result.names = FALSE){
  
  DB.entries = Startup.Data$Reactome.clusterProfiler.Expression.Profiles___Custom.DB.entries
  clusterProfiler.Expression.Profiles(Startup.Data, DB.entries, Analysis.Data = Analysis.Data, GLM.model = GLM.model,
                                      GLM.working.dir = GLM.working.dir, Analysis.Data.RDS.file = Analysis.Data.RDS.file,
                                      database = 'Reactome', printed.name = 'Reactome Pathways', profile.type = 'custom',
                                      forced.parameters = forced.parameters, out.dir = out.dir,
                                      forced.Analysis.name = forced.Analysis.name, stats.file.name = stats.file.name,
                                      sparklines.axis.limits = sparklines.axis.limits,
                                      include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
}


inDetails = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                     forced.parameters = NULL, out.dir = NULL,
                     add.biotypes = NULL,
                     include.GLM.model.name.in.Result.names = FALSE, verbose = TRUE, do.create.Heatmaps = TRUE){
  
  if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
    stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
    if(is.null(Analysis.Data.RDS.file)){
      cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
      cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(file.exists(cf1)){
        Analysis.Data.RDS.file = cf1
      } else if(file.exists(cf2)){
        Analysis.Data.RDS.file = cf2
      } else {
        stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
      }
    }
    Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
    Analysis.name = Analysis.Data$Analysis.name
    Analysis.Data$current.GLM.results.dir = GLM.working.dir
  } else {
    Analysis.name = Analysis.Data$Analysis.name
  }

  if(is.null(add.biotypes)) add.biotypes = Pars$inDetails.additional.biotypes
  
  if(length(Startup.Data$inDetails.entries) == 0){
    cat('\nNo inDetails groups are declared - please check "inDetails" sheet in parameters Excel file\n')
    return(invisible())
  }
  
  if(is.null(out.dir)){
    main.wd = sprintf("%s/inDetails", Analysis.Data$current.GLM.results.dir)
  } else {
    main.wd = out.dir
  }
  dir.create(main.wd, showWarnings = FALSE)
  setwd(main.wd)
  
  if(Pars$re.Cluster.samples.in.Details.Excel.results < Pars$Create.heatmaps..with.re.Clustered.samples.inDetails){
    msg = sprintf('You turned %s re.Cluster.samples.in.Details.Excel.results but turned %s Create.heatmaps..with.re.Clustered.samples.inDetails. Clustered heatmaps cannot be created if the sample clustering is turned FALSE. It will be turned TRUE\n',
                  toString(Pars$re.Cluster.samples.in.Details.Excel.results), toString(Pars$Create.heatmaps..with.re.Clustered.samples.inDetails))
    cat(msg)
    Pars$re.Cluster.samples.in.Details.Excel.results = TRUE
  }


  all__group.Types = vector(mode = 'list')
  all__genes.of.interest = vector(mode = 'list')
  for(group.Name in names(Startup.Data$inDetails.entries)){
    genes.of.interest = c()
    if(!is.null(Startup.Data$inDetails.GO.genes[[group.Name]]))    genes.of.interest = c(genes.of.interest, Startup.Data$inDetails.GO.genes[[group.Name]])
    if(!is.null(Startup.Data$inDetails.Standalone.genes[[group.Name]]))    genes.of.interest = c(genes.of.interest, Startup.Data$inDetails.Standalone.genes[[group.Name]])
    genes.of.interest = genes.of.interest[!duplicated(genes.of.interest)]
    if(length(genes.of.interest) == 0){
      cat(sprintf('No genes were found for group "%s"\n',group.Name))
      next
    } else if(length(genes.of.interest) < 2){
      cat(sprintf('Group "%s" is bypassed since only 1 gene is found within\n',group.Name))
      next
    }
    all__genes.of.interest[[group.Name]] = genes.of.interest
    all__group.Types[[group.Name]] = 'user_declared'
  }
  
  if(length(add.biotypes) > 0){
    #biotypes = add.biotypes[1]
    for(group.Name in names(add.biotypes)){
      biotypes = add.biotypes[[group.Name]]
      #biotypes= 'lincRNA + antisense_RNA + macro_lincRNA'
      #biotypes %<>%  gsub(" ", "", .) %>% strsplit(., split = '+', fixed = TRUE) %>% unlist
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
  }
  
  all.output.file.names = c()
  all.output.file.names.reord = c()

  
  group.Name = names(all__group.Types)[1]
  for(group.Name in names(all__group.Types)){
    genes.of.interest = all__genes.of.interest[[group.Name]]
    cat(sprintf('inDetails analysis of group "%s" (%d genes)...\n', group.Name, length(genes.of.interest)))
    if(length(genes.of.interest) < 2){
      cat('The analysis is bypassed since the genes count is less than 2\n')
      next
    }
    
    group.Name_s = gsub("[^[:alnum:] ]", "", group.Name)
    dir.create(group.Name_s, showWarnings = FALSE)
    setwd(group.Name_s)
    
    if('edgeR_d' %in% names(Analysis.Data) & FALSE){
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
      } else if(Pars$Keep.original..in.Details..gene.order) {
        d$counts = d$counts[match(genes.of.interest, rownames(d$counts)) %>% .[!is.na(.)], ]
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
      
      
      if(dim(d$counts)[2] > 2 & dim(d$counts)[1] > 5){
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
        
        if(dim(d$counts)[2] > 3 & dim(d$counts)[1] > 5){
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
        
        if(dim(d$counts)[2] > 3 & dim(d$counts)[1] > 5){
          png.mod(filename = 'MDS.dim.2-3.png',units="in",res=600,
              width=8.23,height=6.3,pointsize=12)
          a = plotMDS(d,dim.plot = c(2,3), col = col.variety[factor(Main.Component)])
          title(main = plot.title, sub = plot.subtitle)
          dev.off()
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
                if(verbose)
                  cat(sprintf('\rCalculating dist matrices. Completed %.2f%%',ready/N/N*100))
                
                #Weighted.Canberra.dist(c(1,2,3,4,5,6,7,8,9,10),c(1,20,30,4,5,6,7,8,9,10),min.value.to.calc.dist = 0.99)
              }
            }
          }
          
          tmp=as.dist(dist.matrix)
          tryCatch(expr = {
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
          }, error = function (err){
            print('Clustering samples FAILED')
          })
        
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
    }
    
    # if (Pars$Create.Excel.results){
    #   CL = sprintf('%s "%s/DE.results.to.Excel.py"', Pars$python.bin, Pars$suppl.data.dir)
    #   CL = gsub('/','\\',CL,fixed = T)
    #   CL = paste(c(CL, Verify.path(sprintf('"%s [%s], DE results.xlsx"', Analysis.name, group.Name_s)), sprintf('"%s"',Verify.path(output.file.name))),collapse = ' ')
    #   system(CL)
    #   if (!file.exists(Verify.path(sprintf('%s [%s], DE results.xlsx', Analysis.name, group.Name_s)))){
    #     #file.remove(output.file.names)
    #     message(sprintf('Excel worksheet generation (DE results) for ~%s FAILED', Analysis.name))
    #   }
    # }
    
    Main.Component = Analysis.Data$Main.Component
    Main.Component.reordered = Main.Component
    if(!Analysis.Data$perform.classic.paired.test & !Analysis.Data$perform.FC.associations.paired.test){
      Main.Component.reordered = Main.Component[order(Main.Component)]
    }
    Main.Component.levels = Main.Component.reordered[!duplicated(Main.Component.reordered)]
    

    CPMs.groups = vector(mode = 'list')  ## CPMs.groups are reordered according to order(Main.Component) !!!
    LogFCs.groups = vector(mode = 'list')
    if(length(Main.Component.levels) <= Pars$Max.Levels.to.split.Sparkline.groups){
      if(Analysis.Data$perform.FC.associations.paired.test | Analysis.Data$perform.classic.paired.test){
        for(l in Main.Component.levels){
          CPMs.groups[[toString(l)]] = which(Main.Component %in% l)
        }
      } else {
        for(l in Main.Component.levels){
          CPMs.groups[[toString(l)]] = which(Main.Component.reordered %in% l)
        }
      }

      if(Analysis.Data$perform.classic.paired.test){
        LogFCs.groups[['paired LogFCs']] = 1:(length(Main.Component.reordered)/2)
      } else if(Analysis.Data$perform.FC.associations.paired.test){
        Main.Component..LogFC = Main.Component[1:(length(Main.Component)/2)]
        Main.Component..LogFC..levels = Main.Component..LogFC %>% DeDup %>% .[order(.)]
        for(l in Main.Component..LogFC..levels){
          LogFCs.groups[[toString(l)]] = which(Main.Component..LogFC %in% l)
        }
      }
    } else {
      CPMs.groups[['rel. expression']] = 1:(length(Main.Component.reordered))
      if(Analysis.Data$perform.classic.paired.test | Analysis.Data$perform.FC.associations.paired.test)
        LogFCs.groups[['paired LogFCs']] = 1:(length(Main.Component.reordered)/2)
    }

    if(Pars$Keep.original..in.Details..gene.order) {
      gene.numbers = match(genes.of.interest, rownames(Analysis.Data$ResTable.gene.part)) %>% .[!is.na(.)]
      if(length(gene.numbers) > 1){
        ResTable.gene.part = Analysis.Data$ResTable.gene.part[gene.numbers, ]
      } else if(length(gene.numbers) == 1) {
        ResTable.gene.part = subset(Analysis.Data$ResTable.gene.part, subset = (rownames(Analysis.Data$ResTable.gene.part) %in% genes.of.interest))
      }
    } else {
      # ResTable.gene.part = Analysis.Data$ResTable.gene.part[rownames(Analysis.Data$ResTable.gene.part) %in% genes.of.interest, ]
      ResTable.gene.part = subset(Analysis.Data$ResTable.gene.part, subset = (rownames(Analysis.Data$ResTable.gene.part) %in% genes.of.interest))
    }
    
    ResTable.pred.part = Analysis.Data$ResTable.pred.part

    if(Pars$report.norm.read.counts.instead.of.CPM){
      CPMs.spacer.insertion.colname = 'norm. counts'
    } else  CPMs.spacer.insertion.colname = 'CPMs:'

    total.cols.count = ncol(ResTable.pred.part)
    if(Analysis.Data$perform.classic.paired.test | Analysis.Data$perform.FC.associations.paired.test){
      stats.cols.count = which(colnames(ResTable.pred.part) %in% 'LogFCs:')[1]
      LogFCs.col.start = which(colnames(ResTable.pred.part) %in% 'LogFCs:')[1] + 1
      CPMs.col.start = which(colnames(ResTable.pred.part) %in% CPMs.spacer.insertion.colname)[1] + 1
      CPMs.cols.count = total.cols.count - CPMs.col.start + 1
      LogFCs.cols.count = CPMs.col.start - LogFCs.col.start - 1
      if(LogFCs.cols.count < 0)  stop('LogFCs.cols.count cannot be less than 0')

    } else {
      print('---------------------------------------------------')
      print(CPMs.spacer.insertion.colname)
      print('--------------------')
      print(colnames(ResTable.pred.part))
      print('--------------------')
      print(colnames(Analysis.Data$ResTable.gene.part))
      print('--------------------')
      # stop('')
      stats.cols.count = which(colnames(ResTable.pred.part) %in% CPMs.spacer.insertion.colname)[1]
      CPMs.col.start = which(colnames(ResTable.pred.part) %in% CPMs.spacer.insertion.colname)[1] + 1
      CPMs.cols.count = total.cols.count - CPMs.col.start + 1
      LogFCs.col.start = NA
      LogFCs.cols.count = NA
    }

    ResTable = rbind(ResTable.pred.part, ResTable.gene.part)

    output.tsv.file.name = sprintf("%s - %s_%s_combined.txt", group.Name_s, Analysis.name, Pars$DE.package)
    output.tsv.file.name = Verify.path(output.tsv.file.name)
    write.table.mod(ResTable, file = output.tsv.file.name,sep='\t',na = '',col.names = colnames(ResTable))
    all.output.file.names = c(all.output.file.names, output.tsv.file.name)

    inDetails..reordered.samples.seq = NULL
    inDetails..reclustering.samples.passed = FALSE

    Min.group.size.to..cluster.by.CPMs..filter.passed = FALSE
    for(l in names(CPMs.groups)){
        if(length(CPMs.groups[[l]]) >= Pars$Min.group.size.to..cluster)  Min.group.size.to..cluster.by.CPMs..filter.passed = TRUE
    }

    Min.group.size.to..cluster.by.LogFCs..filter.passed = FALSE
    for(l in names(LogFCs.groups)){
      if(length(LogFCs.groups[[l]]) >= Pars$Min.group.size.to..cluster)  Min.group.size.to..cluster.by.LogFCs..filter.passed = TRUE
    }

    if(Pars$re.Cluster.samples.in.Details.Excel.results & Min.group.size.to..cluster.by.CPMs..filter.passed & !(Analysis.Data$perform.classic.paired.test | Analysis.Data$perform.FC.associations.paired.test)){
      cat('re-Clustering samples by CPMs profiles...\n')
      # prep.Counts.Table = Counts.Table[,order(Main.Component)]
      # CPMs.cols.count = ncol(Counts.Table)
      
      l = names(CPMs.groups)[1]
      for(l in names(CPMs.groups)){
        if(length(CPMs.groups[[l]]) < 3)  next
        # which.cols = CPMs.groups[[l]] + stats.cols.count
        which.cols = CPMs.col.start + (CPMs.groups[[l]] - 1)
        genes.CPMs = subset(ResTable.gene.part, select = which.cols)[1:min(dim(ResTable.gene.part)[1], Pars$use.N.top.Genes.to.Cluster.samples), ]
        

        if(type(genes.CPMs[1,]) == "character"){
          rel.LogCPM.table = apply(genes.CPMs, 1, function(all_values) {
            all_values = as.numeric(as.character(all_values))
            gmean_current_CPM = exp(mean(log(all_values + 2)))
            return(log2((all_values + 2) / gmean_current_CPM))
          } )
        } else {
          rel.LogCPM.table = apply(genes.CPMs, 1, function(all_values) {
            gmean_current_CPM = exp(mean(log(all_values + 2)))
            return(log2((all_values + 2) / gmean_current_CPM))
          } )
        }

        # rel.LogCPM.table = apply(genes.CPMs, 1, function(all_values) {
        #   gmean_current_CPM = exp(mean(log(all_values + 2)))
        #   return(log2((all_values + 2) / gmean_current_CPM))
        # } )
        # rel.LogCPM.table = t(scale(t(genes.CPMs)))

        dend.order = tryCatch(expr = {
          dend=hclust(dist(rel.LogCPM.table, 'euclidean'), method = Pars$Clustering.method)
          dend$order
        }, error = function (err){
          cat(sprintf('Cannot reorder samples (by CPMs) in the group %s (total %d samples; %d genes)\n', l, dim(genes.CPMs)[2], dim(genes.CPMs)[1]))
          1:dim(genes.CPMs)[2]
        })
        # print(dend.order)
        # print(which.cols)
        # print(genes.CPMs)
        ResTable[, which.cols] = ResTable[, which.cols][, dend.order]
        colnames(ResTable)[which.cols] = colnames(ResTable)[which.cols][dend.order]
        inDetails..reclustering.samples.passed = TRUE
      }

      inDetails..reordered.samples.seq = colnames(ResTable)[(total.cols.count - CPMs.cols.count + 1) : total.cols.count]

      output.tsv.file.name.reord = sprintf("%s - %s_%s_combined, reord.txt", group.Name_s, Analysis.name, Pars$DE.package)
      output.tsv.file.name.reord = Verify.path(output.tsv.file.name.reord)
      write.table.mod(ResTable, file = output.tsv.file.name.reord, sep='\t',na = '',col.names = colnames(ResTable))
      all.output.file.names.reord = c(all.output.file.names.reord, output.tsv.file.name.reord)
    
    } else if(Pars$re.Cluster.samples.in.Details.Excel.results & Min.group.size.to..cluster.by.LogFCs..filter.passed & (Analysis.Data$perform.classic.paired.test | Analysis.Data$perform.FC.associations.paired.test)){
      cat('re-Clustering samples by LogFCs profiles...\n')
      # prep.Counts.Table = Counts.Table[,order(Main.Component)]
      # CPMs.cols.count = ncol(Counts.Table)
      
      l = names(LogFCs.groups)[1]
      for(l in names(LogFCs.groups)){
        if(length(LogFCs.groups[[l]]) < 3)  next
        # which.cols = CPMs.groups[[l]] + stats.cols.count
        which.cols = LogFCs.col.start + (LogFCs.groups[[l]] - 1)
        which.cols..N = CPMs.col.start + (LogFCs.groups[[l]] - 1)
        which.cols..T = CPMs.col.start + CPMs.cols.count/2 + (LogFCs.groups[[l]] - 1)

        genes.LogFCs = subset(ResTable.gene.part, select = which.cols)[1:min(dim(ResTable.gene.part)[1], Pars$use.N.top.Genes.to.Cluster.samples), ]
        
        # rel.LogCPM.table = apply(genes.CPMs, 1, function(all_values) {
        #   gmean_current_CPM = exp(mean(log(all_values + 2)))
        #   return(log2((all_values + 2) / gmean_current_CPM))
        # } )
        # # rel.LogCPM.table = t(scale(t(genes.CPMs)))

        dend.order = tryCatch(expr = {
          dend=hclust(dist(t(genes.LogFCs), 'euclidean'), method = Pars$Clustering.method)
          dend$order
        }, error = function (err){
          cat(sprintf('Cannot reorder samples (by LogFCs) in the group %s (total %d samples; %d genes)\n', l, dim(genes.LogFCs)[2], dim(genes.LogFCs)[1]))
          1:dim(genes.LogFCs)[2]
        })
        # print(dend.order)
        # print(which.cols)
        # print(genes.CPMs)
        ResTable[, which.cols] = ResTable[, which.cols][, dend.order]
        colnames(ResTable)[which.cols] = colnames(ResTable)[which.cols][dend.order]

        ResTable[, which.cols..N] = ResTable[, which.cols..N][, dend.order]
        colnames(ResTable)[which.cols..N] = colnames(ResTable)[which.cols..N][dend.order]

        ResTable[, which.cols..T] = ResTable[, which.cols..T][, dend.order]
        colnames(ResTable)[which.cols..T] = colnames(ResTable)[which.cols..T][dend.order]

        inDetails..reclustering.samples.passed = TRUE
      }

      inDetails..reordered.samples.seq = colnames(ResTable)[(total.cols.count - CPMs.cols.count + 1) : total.cols.count]

      output.tsv.file.name.reord = sprintf("%s - %s_%s_combined, reord.txt", group.Name_s, Analysis.name, Pars$DE.package)
      output.tsv.file.name.reord = Verify.path(output.tsv.file.name.reord)
      write.table.mod(ResTable, file = output.tsv.file.name.reord, sep='\t',na = '',col.names = colnames(ResTable))
      all.output.file.names.reord = c(all.output.file.names.reord, output.tsv.file.name.reord)
    }

    if (Pars$Create.Excel.results){
      # CPMs.cols.numbers.text = paste(sprintf('%d', (dim(ResTable)[2] - CPMs.cols.count + 1):dim(ResTable)[2]), collapse = ',')
      CPMs.cols.numbers.text = paste(sprintf('%d', CPMs.col.start:ncol(ResTable)), collapse = ',')

      forced_col_types_by_names..CL.insertion = '--forced-col-types-by-names "LogFC:&LogFC&,&trimmed LogFC&,&TBA LogFC&,&gr.1 mean LogFC&,&gr.2 mean LogFC&" "Delta_LogFC:&delta LogFC&,&trimmed delta LogFC&" "LogCPM:&LogCPM&" "p:&PValue&,&TBA p-value&,&p (QLF test)&,&p (LR test)&,&p (QLF test paired)&,&p (LR test paired)&,&p (ex. test)&,&p (Mann-Wh.)&,&p (Wilcoxon paired)&,&p (t-test)&,&p (Mann-Wh., FC)&,&p (t-test, FC)&,&p (Spearman)&,&p (Pearson)&,&p (Spearman, FC)&,&p (Pearson, FC)&" "FDR:&TBA FDR&,&FDR (QLF test)&,&FDR (LR test)&,&FDR (QLF test paired)&,&FDR (LR test paired)&,&FDR (ex. test)&,&FDR (Wilcoxon paired)&,&FDR (Mann-Wh.)&,&FDR (Mann-Wh., FC)&,&FDR (t-test, FC)&,&FDR&" "score:Score" "corr:&Spearman r&,&Pearson r&,&Spearman r (FC)&,&Pearson r (FC)&" "Biotype:Biotype"  "Gene_Symbol:Symbol" "Gene_Name":"Name" "Spacer:&ct_spacer&,&df_spacer&,&LogFCs&" "rank_1:&CG content rank (%)&,&Tr.Length rank (%)&,&LogCPM rank (%)&" "warm_gradient:&d_Freq+&" "cold_gradient:&d_Freq-&" "gray_gradient:&d_Freq max&"'
      if(Startup.Data$miR.mode)
        forced_col_types_by_names..CL.insertion = sprintf('%s "Hidden:&Symbol&,&Biotype&,&Name&,&RefSeq_Summary&,&Annotation&"', forced_col_types_by_names..CL.insertion)
      if(Analysis.Data$perform.FC.associations.paired.test)
        forced_col_types_by_names..CL.insertion = sprintf("%s ", forced_col_types_by_names..CL.insertion)

      if(Analysis.Data$perform.classic.paired.test | Analysis.Data$perform.FC.associations.paired.test){
        LogFCs.cols.numbers.text = paste(sprintf('%d', LogFCs.col.start:(CPMs.col.start-1)), collapse = ',')
        paired.test...CL.insertion = sprintf('"logfc_array:%s"', LogFCs.cols.numbers.text)
      } else paired.test...CL.insertion = ''

      if(Pars$Create.CPM.profiles.in.Excel.results){
        out.Excel.file = Verify.path(sprintf('[%s], DE results, CPM profiles.xlsx', group.Name_s))
        CL = sprintf('%s "%s/Any.2.Excel.py" --in "%s" --out-excel "%s" --one-book yes --sheet-names "%s" --cpm-cells-formatting--logfc-max 2.5 --cpm-cells-formatting--cpm-max 500 --cpm-cells-formatting--gradient-color-schema 3 --cpm-cells-formatting--const-add 0.5 --predictor-rows-count %d %s --entry-type "Gene ID" --forced-cols "cpm_array:%s" "cpm:%s" %s --disable-coltype-autoassign yes --first-col-width 19 --first-col-italic yes --cpm-sparklines-groups %s --cpm-sparklines-mode %s --cpm-sparklines-start-col %d --logfc-sparklines-groups %s --logfc-sparklines-start-col %d %s --max-strings %s',
                     Pars$python.bin, Pars$suppl.data.dir, output.tsv.file.name, out.Excel.file, GLM.model, dim(Analysis.Data$predictor.values.nr)[2],
                     forced_col_types_by_names..CL.insertion,
                     CPMs.cols.numbers.text, CPMs.cols.numbers.text,
                     paired.test...CL.insertion,
                     sapply(names(CPMs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(CPMs.col.start - 1 + CPMs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     'relative', 7,
                     sapply(names(LogFCs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(LogFCs.col.start - 1 + LogFCs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     7,
                     ifelse('TBA LogFC' %in% colnames(ResTable), '--forced-col-widths-by-names 1:"&trimmed LogFC&,&Trimmed LogFC&"', ''),
                     as.character(Pars$Max.genes.in.Excel.results))

        if(Pars$bidirectional.bars.in.LogFC.cells.in.Excel.reports){
          CL = sprintf("%s --bidirectional-bars-in-logfc-cells yes --simple-logfc-format no", CL)
        } else if (Pars$colored.LogFC.cells.in.Excel.reports | Analysis.Data$perform.FC.associations.paired.test){
          CL = sprintf("%s --simple-logfc-format yes --bidirectional-bars-in-logfc-cells no", CL)
        } else {
          CL = sprintf("%s --simple-logfc-format no --bidirectional-bars-in-logfc-cells no", CL)
        }

        if(Pars$white.background.in.Excel.reports){  CL = sprintf("%s --whitespace-mode yes", CL)
        } else                                       CL = sprintf("%s --whitespace-mode no",  CL)

        if(Pars$freeze.header.in.Excel.reports)  CL = sprintf("%s --freeze-rows 1", CL)
        if(Pars$freeze.gene.names.in.Excel.reports)  CL = sprintf("%s --freeze-cols 2", CL)
        
        if(Pars$Include.small.heatmaps.in.Excel.reports){
          CL = sprintf("%s --include-cpm-heatmaps yes --include-logfc-heatmaps yes", CL)
        }

        code = system(CL)
        if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
        if (!file.exists(out.Excel.file)){
          message(sprintf('Excel worksheet generation (DE results) for ~%s [%s] FAILED', Analysis.name, group.Name_s))
        }

        if(Pars$re.Cluster.samples.in.Details.Excel.results & Min.group.size.to..cluster.by.CPMs..filter.passed){
          out.Excel.file = Verify.path(sprintf('[%s], DE results, CPM profiles [reordered].xlsx', group.Name_s))
          CL = sprintf('%s "%s/Any.2.Excel.py" --in "%s" --out-excel "%s" --one-book yes --sheet-names "%s" --cpm-cells-formatting--logfc-max 2.5 --cpm-cells-formatting--cpm-max 500 --cpm-cells-formatting--gradient-color-schema 3 --cpm-cells-formatting--const-add 0.5 --predictor-rows-count %d %s --entry-type "Gene ID" --forced-cols "cpm_array:%s" "cpm:%s" %s --disable-coltype-autoassign yes --first-col-width 19 --first-col-italic yes --cpm-sparklines-groups %s --cpm-sparklines-mode %s --cpm-sparklines-start-col %d --logfc-sparklines-groups %s --logfc-sparklines-start-col %d %s --max-strings %s',
                       Pars$python.bin, Pars$suppl.data.dir, output.tsv.file.name.reord, out.Excel.file, GLM.model, dim(Analysis.Data$predictor.values.nr)[2],
                       forced_col_types_by_names..CL.insertion,
                       CPMs.cols.numbers.text, CPMs.cols.numbers.text,
                       paired.test...CL.insertion,
                       sapply(names(CPMs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(CPMs.col.start - 1 + CPMs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                       'relative', 7,
                       sapply(names(LogFCs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(LogFCs.col.start - 1 + LogFCs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                       7,
                       ifelse('TBA LogFC' %in% colnames(ResTable), '--forced-col-widths-by-names 1:"&trimmed LogFC&,&Trimmed LogFC&"', ''),
                       as.character(Pars$Max.genes.in.Excel.results))

          if(Pars$bidirectional.bars.in.LogFC.cells.in.Excel.reports){
            CL = sprintf("%s --bidirectional-bars-in-logfc-cells yes --simple-logfc-format no", CL)
          } else if (Pars$colored.LogFC.cells.in.Excel.reports | Analysis.Data$perform.FC.associations.paired.test){
            CL = sprintf("%s --simple-logfc-format yes --bidirectional-bars-in-logfc-cells no", CL)
          } else {
            CL = sprintf("%s --simple-logfc-format no --bidirectional-bars-in-logfc-cells no", CL)
          }

          if(Pars$white.background.in.Excel.reports){  CL = sprintf("%s --whitespace-mode yes", CL)
          } else                                       CL = sprintf("%s --whitespace-mode no",  CL)

          if(Pars$freeze.header.in.Excel.reports)  CL = sprintf("%s --freeze-rows 1", CL)
          if(Pars$freeze.gene.names.in.Excel.reports)  CL = sprintf("%s --freeze-cols 2", CL)
          
          if(Pars$Include.small.heatmaps.in.Excel.reports){
            CL = sprintf("%s --include-cpm-heatmaps yes --include-logfc-heatmaps yes", CL)
          }

          code = system(CL)
          if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
          if (!file.exists(out.Excel.file)){
            message(sprintf('Excel worksheet generation (DE results) for ~%s [%s] FAILED', Analysis.name, group.Name_s))
          }
        }
      }
      
      if(Pars$Create.rel.LogCPM.profiles.in.Excel.results){
        out.Excel.file = Verify.path(sprintf('[%s], DE results, Log.rel.CPM profiles.xlsx', group.Name_s))
        CL = sprintf('%s "%s/Any.2.Excel.py" --in "%s" --out-excel "%s" --one-book yes --sheet-names "%s" --cpm-cells-formatting--logfc-max 2.5 --cpm-cells-formatting--cpm-max 500 --cpm-cells-formatting--gradient-color-schema 3 --cpm-cells-formatting--const-add 0.5 --predictor-rows-count %d %s --entry-type "Gene ID" --forced-cols "cpm_array:%s" "cpm:%s" %s --disable-coltype-autoassign yes --first-col-width 19 --first-col-italic yes --cpm-sparklines-groups %s --cpm-sparklines-mode %s --cpm-sparklines-start-col %d --logfc-sparklines-groups %s --logfc-sparklines-start-col %d %s --max-strings %s',
                     Pars$python.bin, Pars$suppl.data.dir, output.tsv.file.name, out.Excel.file, GLM.model, dim(Analysis.Data$predictor.values.nr)[2],
                     forced_col_types_by_names..CL.insertion,
                     CPMs.cols.numbers.text, CPMs.cols.numbers.text,
                     paired.test...CL.insertion,
                     sapply(names(CPMs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(CPMs.col.start - 1 + CPMs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     'logfc', 7,
                     sapply(names(LogFCs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(LogFCs.col.start - 1 + LogFCs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     7,
                     ifelse('TBA LogFC' %in% colnames(ResTable), '--forced-col-widths-by-names 1:"&trimmed LogFC&,&Trimmed LogFC&"', ''),
                     as.character(Pars$Max.genes.in.Excel.results))

        if(Pars$bidirectional.bars.in.LogFC.cells.in.Excel.reports){
          CL = sprintf("%s --bidirectional-bars-in-logfc-cells yes --simple-logfc-format no", CL)
        } else if (Pars$colored.LogFC.cells.in.Excel.reports | Analysis.Data$perform.FC.associations.paired.test){
          CL = sprintf("%s --simple-logfc-format yes --bidirectional-bars-in-logfc-cells no", CL)
        } else {
          CL = sprintf("%s --simple-logfc-format no --bidirectional-bars-in-logfc-cells no", CL)
        }

        if(Pars$white.background.in.Excel.reports){  CL = sprintf("%s --whitespace-mode yes", CL)
        } else                                       CL = sprintf("%s --whitespace-mode no",  CL)

        if(Pars$freeze.header.in.Excel.reports)  CL = sprintf("%s --freeze-rows 1", CL)
        if(Pars$freeze.gene.names.in.Excel.reports)  CL = sprintf("%s --freeze-cols 2", CL)
        
        if(Pars$Include.small.heatmaps.in.Excel.reports){
          CL = sprintf("%s --include-cpm-heatmaps yes --include-logfc-heatmaps yes", CL)
        }

        code = system(CL)
        if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
        if (!file.exists(out.Excel.file)){
          message(sprintf('Excel worksheet generation (DE results) for ~%s [%s] FAILED', Analysis.name, group.Name_s))
        }

        if(Pars$re.Cluster.samples.in.Details.Excel.results & Min.group.size.to..cluster.by.CPMs..filter.passed){
          out.Excel.file = Verify.path(sprintf('[%s], DE results, Log.rel.CPM profiles [reordered].xlsx', group.Name_s))
          CL = sprintf('%s "%s/Any.2.Excel.py" --in "%s" --out-excel "%s" --one-book yes --sheet-names "%s" --cpm-cells-formatting--logfc-max 2.5 --cpm-cells-formatting--cpm-max 500 --cpm-cells-formatting--gradient-color-schema 3 --cpm-cells-formatting--const-add 0.5 --predictor-rows-count %d %s --entry-type "Gene ID" --forced-cols "cpm_array:%s" "cpm:%s" %s --disable-coltype-autoassign yes --first-col-width 19 --first-col-italic yes --cpm-sparklines-groups %s --cpm-sparklines-mode %s --cpm-sparklines-start-col %d --logfc-sparklines-groups %s --logfc-sparklines-start-col %d %s --max-strings %s',
                       Pars$python.bin, Pars$suppl.data.dir, output.tsv.file.name.reord, out.Excel.file, GLM.model, dim(Analysis.Data$predictor.values.nr)[2],
                       forced_col_types_by_names..CL.insertion,
                       CPMs.cols.numbers.text, CPMs.cols.numbers.text,
                       paired.test...CL.insertion,
                       sapply(names(CPMs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(CPMs.col.start - 1 + CPMs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                       'logfc', 7,
                       sapply(names(LogFCs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(LogFCs.col.start - 1 + LogFCs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                       7,
                       ifelse('TBA LogFC' %in% colnames(ResTable), '--forced-col-widths-by-names 1:"&trimmed LogFC&,&Trimmed LogFC&"', ''),
                       as.character(Pars$Max.genes.in.Excel.results))

          if(Pars$bidirectional.bars.in.LogFC.cells.in.Excel.reports){
            CL = sprintf("%s --bidirectional-bars-in-logfc-cells yes --simple-logfc-format no", CL)
          } else if (Pars$colored.LogFC.cells.in.Excel.reports | Analysis.Data$perform.FC.associations.paired.test){
            CL = sprintf("%s --simple-logfc-format yes --bidirectional-bars-in-logfc-cells no", CL)
          } else {
            CL = sprintf("%s --simple-logfc-format no --bidirectional-bars-in-logfc-cells no", CL)
          }

          if(Pars$white.background.in.Excel.reports){  CL = sprintf("%s --whitespace-mode yes", CL)
          } else                                       CL = sprintf("%s --whitespace-mode no",  CL)

          if(Pars$freeze.header.in.Excel.reports)  CL = sprintf("%s --freeze-rows 1", CL)
          if(Pars$freeze.gene.names.in.Excel.reports)  CL = sprintf("%s --freeze-cols 2", CL)
          
          if(Pars$Include.small.heatmaps.in.Excel.reports){
            CL = sprintf("%s --include-cpm-heatmaps yes --include-logfc-heatmaps yes", CL)
          }

          code = system(CL)
          if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
          if (!file.exists(out.Excel.file)){
            message(sprintf('Excel worksheet generation (DE results) for ~%s [%s] FAILED', Analysis.name, group.Name_s))
          }
        }
      }      
    }
    
    
    if(do.create.Heatmaps){
      if(!inDetails..reclustering.samples.passed)  inDetails..reordered.samples.seq = NULL
      Create.Heatmaps(Startup.Data, Analysis.Data = Analysis.Data,
                    limit.to.genes = genes.of.interest, out.dir.root = '.', max.Genes.counts.list = Pars$Top.genes.to.include.inDetails.heatmaps.list,
                    add.Heatmaps.with.Custom.samples.order = inDetails..reordered.samples.seq,
                    Heatmaps.with.Custom.samples.order.dir = 'Heatmaps - re-clustered samples')
    }
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
  
  # stats.cols.count = which(colnames(ResTable) %in% 'CPMs:')[1]
  # CPMs.cols.count = dim(ResTable)[2] - stats.cols.count
  # CPMs.col.start
  
  if(Pars$report.norm.read.counts.instead.of.CPM){
    CPMs.spacer.insertion.colname = 'norm. counts'
  } else  CPMs.spacer.insertion.colname = 'CPMs:'

  total.cols.count = ncol(ResTable)
  if(Analysis.Data$perform.classic.paired.test | Analysis.Data$perform.FC.associations.paired.test){
    stats.cols.count = which(colnames(ResTable) %in% 'LogFCs:')[1]
    LogFCs.col.start = which(colnames(ResTable) %in% 'LogFCs:')[1] + 1
    CPMs.col.start = which(colnames(ResTable) %in% 'CPMs.spacer.insertion.colname')[1] + 1
    CPMs.cols.count = total.cols.count - CPMs.col.start + 1
    LogFCs.cols.count = CPMs.col.start - LogFCs.col.start - 1
    if(LogFCs.cols.count < 0)  stop('LogFCs.cols.count cannot be less than 0')

  } else {
    stats.cols.count = which(colnames(ResTable) %in% 'CPMs.spacer.insertion.colname')[1]
    CPMs.col.start = which(colnames(ResTable) %in% 'CPMs.spacer.insertion.colname')[1] + 1
    CPMs.cols.count = total.cols.count - CPMs.col.start + 1
    LogFCs.col.start = NA
    LogFCs.cols.count = NA
  }

  Main.Component = Analysis.Data$Main.Component
  Main.Component.reordered = Main.Component
  if(!Analysis.Data$perform.classic.paired.test & !Analysis.Data$perform.FC.associations.paired.test){
    Main.Component.reordered = Main.Component[order(Main.Component)]
  }
  Main.Component.levels = Main.Component.reordered[!duplicated(Main.Component.reordered)]

  CPMs.groups = vector(mode = 'list')  ## CPMs.groups are reordered according to order(Main.Component) !!!
  LogFCs.groups = vector(mode = 'list')

  if(length(Main.Component.levels) <= Pars$Max.Levels.to.split.Sparkline.groups){
    if(Analysis.Data$perform.FC.associations.paired.test | Analysis.Data$perform.classic.paired.test){
      for(l in Main.Component.levels){
        CPMs.groups[[toString(l)]] = which(Main.Component %in% l)
      }
    } else {
      for(l in Main.Component.levels){
        CPMs.groups[[toString(l)]] = which(Main.Component.reordered %in% l)
      }
    }

    if(Analysis.Data$perform.classic.paired.test){
      LogFCs.groups[['paired LogFCs']] = 1:(length(Main.Component.reordered)/2)
    } else if(Analysis.Data$perform.FC.associations.paired.test){
      Main.Component..LogFC = Main.Component[1:(length(Main.Component)/2)]
      Main.Component..LogFC..levels = Main.Component..LogFC %>% DeDup %>% .[order(.)]
      for(l in Main.Component..LogFC..levels){
        LogFCs.groups[[toString(l)]] = which(Main.Component..LogFC %in% l)
      }
    }
  } else {
    CPMs.groups[['rel. expression']] = 1:(length(Main.Component.reordered))
    if(Analysis.Data$perform.classic.paired.test | Analysis.Data$perform.FC.associations.paired.test)
      LogFCs.groups[['paired LogFCs']] = 1:(length(Main.Component.reordered)/2)
  }

  # CPMs.cols.numbers.text = paste(sprintf('%d', (dim(ResTable)[2] - CPMs.cols.count + 1):dim(ResTable)[2]), collapse = ',')
  CPMs.cols.numbers.text = paste(sprintf('%d', CPMs.col.start:ncol(ResTable)), collapse = ',')

  if (Pars$Create.Excel.results){
    forced_col_types_by_names..CL.insertion = '--forced-col-types-by-names "LogFC:&LogFC&,&trimmed LogFC&,&TBA LogFC&,&gr.1 mean LogFC&,&gr.2 mean LogFC&" "Delta_LogFC:&delta LogFC&,&trimmed delta LogFC&" "LogCPM:&LogCPM&" "p:&PValue&,&TBA p-value&,&p (QLF test)&,&p (LR test)&,&p (QLF test paired)&,&p (LR test paired)&,&p (ex. test)&,&p (Mann-Wh.)&,&p (Wilcoxon paired)&,&p (t-test)&,&p (Mann-Wh., FC)&,&p (t-test, FC)&,&p (Spearman)&,&p (Pearson)&,&p (Spearman, FC)&,&p (Pearson, FC)&" "FDR:&TBA FDR&,&FDR (QLF test)&,&FDR (LR test)&,&FDR (QLF test paired)&,&FDR (LR test paired)&,&FDR (ex. test)&,&FDR (Wilcoxon paired)&,&FDR (Mann-Wh.)&,&FDR (Mann-Wh., FC)&,&FDR (t-test, FC)&,&FDR&" "score:Score" "corr:&Spearman r&,&Pearson r&,&Spearman r (FC)&,&Pearson r (FC)&" "Biotype:Biotype"  "Gene_Symbol:Symbol" "Gene_Name":"Name" "Spacer:&ct_spacer&,&df_spacer&,&LogFCs&" "rank_1:&CG content rank (%)&,&Tr.Length rank (%)&,&LogCPM rank (%)&" "warm_gradient:&d_Freq+&" "cold_gradient:&d_Freq-&" "gray_gradient:&d_Freq max&"'

    if(Startup.Data$miR.mode)
      forced_col_types_by_names..CL.insertion = sprintf('%s "Hidden:&Symbol&,&Biotype&,&Name&,&RefSeq_Summary&,&Annotation&"', forced_col_types_by_names..CL.insertion)
    if(Analysis.Data$perform.FC.associations.paired.test)
      forced_col_types_by_names..CL.insertion = sprintf("%s ", forced_col_types_by_names..CL.insertion)

    if(Analysis.Data$perform.classic.paired.test | Analysis.Data$perform.FC.associations.paired.test){
      LogFCs.cols.numbers.text = paste(sprintf('%d', LogFCs.col.start:(CPMs.col.start-1)), collapse = ',')
      paired.test...CL.insertion = sprintf('"logfc_array:%s"', LogFCs.cols.numbers.text)
    } else paired.test...CL.insertion = ''

    if(Pars$Create.CPM.profiles.in.Excel.results){
      out.Excel.file = Verify.path(sprintf('inDetails - DE results, CPM profiles.xlsx'))

      CL = sprintf('%s "%s/Any.2.Excel.py" --in %s --out-excel "%s" --one-book yes --sheet-names %s --cpm-cells-formatting--logfc-max 2.5 --cpm-cells-formatting--cpm-max 500 --cpm-cells-formatting--gradient-color-schema 3 --cpm-cells-formatting--const-add 0.5 --predictor-rows-count %d %s --entry-type "Gene ID" --forced-cols "cpm_array:%s" "cpm:%s" %s --disable-coltype-autoassign yes --first-col-width 19 --first-col-italic yes --cpm-sparklines-groups %s --cpm-sparklines-mode %s --cpm-sparklines-start-col %d --logfc-sparklines-groups %s --logfc-sparklines-start-col %d %s --max-strings %s',
                     Pars$python.bin, Pars$suppl.data.dir,  paste(sprintf('"%s"', all.output.file.names), collapse = ' '), out.Excel.file, paste(sprintf('"%s"', all.group.Name_s), collapse = ' '),
                     dim(Analysis.Data$predictor.values.nr)[2],
                     forced_col_types_by_names..CL.insertion,
                     CPMs.cols.numbers.text, CPMs.cols.numbers.text,
                     paired.test...CL.insertion,
                     sapply(names(CPMs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(CPMs.col.start - 1 + CPMs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     'relative', 7,
                     sapply(names(LogFCs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(LogFCs.col.start - 1 + LogFCs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     7,
                     ifelse('TBA LogFC' %in% colnames(ResTable), '--forced-col-widths-by-names 1:"&trimmed LogFC&,&Trimmed LogFC&"', ''),
                     as.character(Pars$Max.genes.in.Excel.results))

      if(Pars$bidirectional.bars.in.LogFC.cells.in.Excel.reports){
        CL = sprintf("%s --bidirectional-bars-in-logfc-cells yes --simple-logfc-format no", CL)
      } else if (Pars$colored.LogFC.cells.in.Excel.reports | Analysis.Data$perform.FC.associations.paired.test){
        CL = sprintf("%s --simple-logfc-format yes --bidirectional-bars-in-logfc-cells no", CL)
      } else {
        CL = sprintf("%s --simple-logfc-format no --bidirectional-bars-in-logfc-cells no", CL)
      }

      if(Pars$white.background.in.Excel.reports){  CL = sprintf("%s --whitespace-mode yes", CL)
      } else                                       CL = sprintf("%s --whitespace-mode no",  CL)

      if(Pars$freeze.header.in.Excel.reports)  CL = sprintf("%s --freeze-rows 1", CL)
      if(Pars$freeze.gene.names.in.Excel.reports)  CL = sprintf("%s --freeze-cols 2", CL)
      
      if(Pars$Include.small.heatmaps.in.Excel.reports){
        CL = sprintf("%s --include-cpm-heatmaps yes --include-logfc-heatmaps yes", CL)
      }

      code = system(CL)
      if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
      if (!file.exists(out.Excel.file)){
        message(sprintf('Excel worksheet generation (inDetails DE results) for %s FAILED', Analysis.name))
      }
      file.copy(Verify.path(out.Excel.file), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, out.Excel.file)), overwrite = TRUE)
    }

    if(Pars$Create.rel.LogCPM.profiles.in.Excel.results){
      out.Excel.file = Verify.path(sprintf('inDetails - DE results, Log.rel.CPM profiles.xlsx'))
      CL = sprintf('%s "%s/Any.2.Excel.py" --in %s --out-excel "%s" --one-book yes --sheet-names %s --cpm-cells-formatting--logfc-max 2.5 --cpm-cells-formatting--cpm-max 500 --cpm-cells-formatting--gradient-color-schema 3 --cpm-cells-formatting--const-add 0.5 --predictor-rows-count %d %s --entry-type "Gene ID" --forced-cols "cpm_array:%s" "cpm:%s" %s --disable-coltype-autoassign yes --first-col-width 19 --first-col-italic yes --cpm-sparklines-groups %s --cpm-sparklines-mode %s --cpm-sparklines-start-col %d --logfc-sparklines-groups %s --logfc-sparklines-start-col %d %s --max-strings %s',
                     Pars$python.bin, Pars$suppl.data.dir,  paste(sprintf('"%s"', all.output.file.names), collapse = ' '), out.Excel.file, paste(sprintf('"%s"', all.group.Name_s), collapse = ' '),
                     dim(Analysis.Data$predictor.values.nr)[2],
                     forced_col_types_by_names..CL.insertion,
                     CPMs.cols.numbers.text, CPMs.cols.numbers.text,
                     paired.test...CL.insertion,
                     sapply(names(CPMs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(CPMs.col.start - 1 + CPMs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     'logfc', 7,
                     sapply(names(LogFCs.groups), function(x) sprintf('"%s:%s"', Analysis.Data$Main.Component.values..to..group.names[[toString(x)]] %>% gsub(":", "_", .), paste(LogFCs.col.start - 1 + LogFCs.groups[[x]], collapse = ','))) %>% paste(., collapse = ' '),
                     7,
                     ifelse('TBA LogFC' %in% colnames(ResTable), '--forced-col-widths-by-names 1:"&trimmed LogFC&,&Trimmed LogFC&"', ''),
                     as.character(Pars$Max.genes.in.Excel.results))

      if(Pars$bidirectional.bars.in.LogFC.cells.in.Excel.reports){
        CL = sprintf("%s --bidirectional-bars-in-logfc-cells yes --simple-logfc-format no", CL)
      } else if (Pars$colored.LogFC.cells.in.Excel.reports | Analysis.Data$perform.FC.associations.paired.test){
        CL = sprintf("%s --simple-logfc-format yes --bidirectional-bars-in-logfc-cells no", CL)
      } else {
        CL = sprintf("%s --simple-logfc-format no --bidirectional-bars-in-logfc-cells no", CL)
      }

      if(Pars$white.background.in.Excel.reports){  CL = sprintf("%s --whitespace-mode yes", CL)
      } else                                       CL = sprintf("%s --whitespace-mode no",  CL)

      if(Pars$freeze.header.in.Excel.reports)  CL = sprintf("%s --freeze-rows 1", CL)
      if(Pars$freeze.gene.names.in.Excel.reports)  CL = sprintf("%s --freeze-cols 2", CL)
      
      if(Pars$Include.small.heatmaps.in.Excel.reports){
        CL = sprintf("%s --include-cpm-heatmaps yes --include-logfc-heatmaps yes", CL)
      }

      code = system(CL)
      if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
      if (!file.exists(out.Excel.file)){
        message(sprintf('Excel worksheet generation (inDetails DE results) for %s FAILED', Analysis.name))
      }
      file.copy(Verify.path(out.Excel.file), Verify.path(sprintf('%s/%s', Analysis.Data$current.GLM.results.dir, out.Excel.file)), overwrite = TRUE)
    }
  }

  # CL = sprintf('%s "%s/DE.results.to.Excel.py"',Pars$python.bin,Pars$suppl.data.dir)
  # CL = gsub('/','\\',CL,fixed = T)
  # out.file.name = sprintf('inDetails - %s, DE results.xlsx',Analysis.name)
  # CL = paste(c(CL,sprintf('"%s"', Verify.path(out.file.name)),sprintf('"%s"', all.output.file.names)),collapse = ' ')
  # system(CL)
  # if (!file.exists(Verify.path(out.file.name))){
    #file.remove(output.file.names)
    # stop(sprintf('Excel worksheet generation FAILED [non-joint], %s',out.file.name))
  # }
  
  if (is.null(out.dir)) setwd(Analysis.Data$current.GLM.results.dir)
}


Collect.topGO.enrichment.result.files <- function(Startup.Data, GLM.models = NULL, forced.parameters = NULL,
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
                                  include.GLM.model.name.in.Result.names = FALSE){
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
        if (str_count(src.file.name.skeleton, "%s") == 2){
          src.file.name = sprintf(src.file.name.skeleton, Analysis.name, GO.type)
        } else if(str_count(src.file.name.skeleton, "%s") == 1) {
          src.file.name = sprintf(src.file.name.skeleton, GO.type)
        } else {
          stop('Incorrect skeleton')
        }
        
        src.file.full.name = sprintf('%s/%s/%s', current.GLM.results.dir, src.dir.name, src.file.name)

        if(file.exists(src.file.full.name)){
          src.files = c(src.files, src.file.full.name)
        } else if (file.exists(Verify.path(src.file.full.name))){
          src.files = c(src.files, Verify.path(src.file.full.name))
        } else {
          not.found.src.files = c(not.found.src.files, src.file.full.name)
          next
        }
        
        dest.file = Verify.path(sprintf('%s/[%s]    %s.xlsx', out.dir, printed.name, Analysis.name))
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
      if (str_count(src.file.name.skeleton, "%s") == 1){
        src.file.name = sprintf(src.file.name.skeleton, Analysis.name)
      } else if (str_count(src.file.name.skeleton, "%s") == 0){
        src.file.name = sprintf(src.file.name.skeleton)
      } else {
        stop('Incorrect skeleton')
      }
      
      src.file.full.name = sprintf('%s/%s/%s', current.GLM.results.dir, src.dir.name, src.file.name)
      
      if(file.exists(src.file.full.name)){
        src.files = c(src.files, src.file.full.name)
      } else if (file.exists(Verify.path(src.file.full.name))){
        src.files = c(src.files, Verify.path(src.file.full.name))
      } else {
        not.found.src.files = c(not.found.src.files, src.file.full.name)
        next
      }
      
      dest.file = Verify.path(sprintf('%s/[%s]    %s.xlsx', out.dir, printed.name, Analysis.name))
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

Summarize.results.by.Analyses.groups = function(Startup.Data,
                                                forced.GLM.Models = NULL, forced.group.Name = NULL,
                                                forced.parameters = NULL, out.dir = NULL, joint.GLM.results = TRUE,
                                                KEGG.ep = TRUE, KEGG.de = TRUE, Reactome.ep = TRUE, GO.ep = TRUE, topGO.ep = TRUE, bias = TRUE,
                                                collect.inDetails = TRUE,
                                                include.GLM.model.name.in.Result.names = FALSE, include.with.spaces = FALSE){
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
  
  if(is.null(forced.GLM.Models) & !is.null(forced.group.Name)){
    stop('forced.Models should be defined if forced.group.Name is defined too')
  }
  
  if(is.null(forced.GLM.Models)){
    group.names = names(Pars$Analyses.groups)
    Analyses.groups = Pars$Analyses.groups
  } else {
    if(is.null(forced.group.Name)){
      group.names = c('(group 1)')
    } else {
      group.names = c(forced.group.Name)
    }
    Analyses.groups = vector(mode ='list')
    Analyses.groups[[group.names[1]]] = forced.GLM.Models
  }
  
  group.name = group.names[1]
  for(group.name in group.names){
    GLM.models = Analyses.groups[[group.name]]
    
    current.group.out.dir = sprintf('%s/Essential results - %s', out.dir, group.name)
    dir.create(current.group.out.dir, showWarnings = FALSE)
    
    if(include.GLM.model.name.in.Result.names){
      Collect.GLM.results(Startup.Data, GLM.models, forced.parameters = Pars, out.dir = current.group.out.dir,
                          include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
      if(collect.inDetails){
        inDetails.dir = sprintf('%s/%s - inDetails DE summary', current.group.out.dir, group.name)
        dir.create(inDetails.dir, showWarnings = FALSE)
        Collect.inDetails.GLM.results(Startup.Data, GLM.models = GLM.models, out.dir = inDetails.dir, add.biotypes = Pars$inDetails.additional.biotypes)
      }
      if(!Startup.Data$miR.mode){
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
      }
    } else {
      Collect.GLM.results(Startup.Data, GLM.models, forced.parameters = Pars, out.dir = current.group.out.dir,
                          include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
      if(collect.inDetails){
        inDetails.dir = sprintf('%s/%s - inDetails DE summary', current.group.out.dir, group.name)
        dir.create(inDetails.dir, showWarnings = FALSE)
        Collect.inDetails.GLM.results(Startup.Data, GLM.models = GLM.models, out.dir = inDetails.dir, add.biotypes = Pars$inDetails.additional.biotypes)
      }
      if(!Startup.Data$miR.mode){
        Collect.Enrich.results(Startup.Data, GLM.models, forced.parameters = Pars, out.dir = current.group.out.dir,
                               database = 'topGO', printed.name = 'GO terms - topGO - classic enrichment', GO.types = c('BP','MF','CC'),
                               file.name.skeleton = 'GO Terms (%s) - topGO - Enrichment Analysis/GO Terms (%s) topGO enrichment.xlsx')
        Collect.Enrich.results(Startup.Data, GLM.models, forced.parameters = Pars, out.dir = current.group.out.dir,
                               database = 'GO', printed.name = 'GO terms classic enrichment', GO.types = c('BP','MF','CC'),
                               file.name.skeleton = 'GO Terms (%s) - classic Enrichment Analysis/GO Terms (%s) classic enrichment.xlsx')
        Collect.Enrich.results(Startup.Data, GLM.models, forced.parameters = Pars, out.dir = current.group.out.dir,
                               database = 'KEGG', printed.name = 'KEGG Pathways classic enrichment',
                               file.name.skeleton = 'KEGG Pathways - classic Enrichment Analysis/KEGG Pathways classic enrichment.xlsx')
        Collect.Enrich.results(Startup.Data, GLM.models, forced.parameters = Pars, out.dir = current.group.out.dir,
                               database = 'Reactome', printed.name = 'Reactome Pathways classic enrichment',
                               file.name.skeleton = 'Reactome Pathways - classic Enrichment Analysis/Reactome Pathways classic enrichment.xlsx')
      }
    }

    if(joint.GLM.results){
      cat('Creating DE joint summary reports....\n')
      Summarize.GLM.results(Startup.Data, GLM.models = GLM.models, forced.parameters = Pars,
                            out.file.name = sprintf('%s/%s - Summary expression info.xlsx', current.group.out.dir, group.name))
      inDetails.dir = sprintf('%s/%s - inDetails DE summary', current.group.out.dir, group.name)
      dir.create(inDetails.dir, showWarnings = FALSE)
      cat('Creating DE joint summary reports (inDetails)....\n')
      Summarize.inDetails.GLM.results(Startup.Data, GLM.models = GLM.models, out.dir = inDetails.dir)
    }

    if(Pars$Create.Venn_Euler.diagrams)
      cat('Creating Euler-Venn diagrams....\n')
      Create.Euler.Venn.diagrams(Startup.Data, GLM.models = GLM.models, forced.parameters = Pars,
        output.dir = current.group.out.dir, reported.analysis.group.name = group.name, bypass..too.large.comparisons = TRUE)
    
    if(!Startup.Data$miR.mode){
      if(KEGG.ep)
        Summarize.KEGG.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.models, out.dir = current.group.out.dir,
                                                out.file.prefix = sprintf('@ %s - ', group.name),
                                                include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names, include.with.spaces = include.with.spaces)
      if(KEGG.de)
        Summarize.KEGG.DE.info(Startup.Data, GLM.models = GLM.models, out.dir = current.group.out.dir,
                             out.file.prefix = sprintf('@ %s - ', group.name))
      
      if(Reactome.ep)
        Summarize.Reactome.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.models, out.dir = current.group.out.dir,
                                                      out.file.prefix = sprintf('@ %s - ', group.name),
                                                      include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names, include.with.spaces = include.with.spaces)
      if(GO.ep)
        Summarize.GO.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.models, out.dir = current.group.out.dir,
                                                out.file.prefix = sprintf('@ %s - ', group.name),
                                                include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names, include.with.spaces = include.with.spaces)

      if(topGO.ep)
        Summarize.topGO.Expression.Profiles.Custom(Startup.Data, GLM.models = GLM.models, out.dir = current.group.out.dir,
                                                   out.file.prefix = sprintf('@ %s - ', group.name),
                                                   include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names, include.with.spaces = include.with.spaces)
      
      if(bias){
        Summarize.Bias.Analysis.Results(Startup.Data, GLM.models = GLM.models, out.dir = current.group.out.dir,
                                      out.file.prefix = sprintf('@ %s - ', group.name), DB.suffix = '')
        Summarize.Bias.Analysis.Results(Startup.Data, GLM.models = GLM.models, out.dir = current.group.out.dir,
                                      out.file.prefix = sprintf('@ %s - final ', group.name), DB.suffix = '.adjusted')
      }
    }
  }
}


Bias.Factor.LogFC.Plots = function(approx_line__const_coeff, approx_line__slope, spearman_r, BF.values, LogFC.values, P.values, 
                             bf.name, Analysis.name, LogFC.col.name, out.dir, LogCPM.limit = NULL, logFC.range = 3, p.threshold = 0.05, BF.values.axis.limits = c(8, 14),
                             plots.suffix = '',
                             trimmed.approx_line__const_coeff = NULL, trimmed.approx_line__slope = NULL,
                             include.GLM.model.name.in.Result.names = FALSE){
  
  hist2d.cols <- colorRampPalette(c("#ffffff","#4bbc00","#e59400", "#f77700", "#f70000", "#f70000"))(512)
  hist2d.cols.v2 <- colorRampPalette(c("#ffffff",  '#a0c9e7',  '#789ae7', "#e59400", "#f77700", "#f70000", "#f70000"))(512)

  if(is.null(trimmed.approx_line__const_coeff) != is.null(trimmed.approx_line__slope)){
    stop('Both trimmed.approx_line__const_coeff and trimmed.approx_line__slope should be either defined or not defined')
  }

  if(is.null(trimmed.approx_line__const_coeff)){
    main.approx_line__const_coeff = approx_line__const_coeff
    main.approx_line__slope = approx_line__slope
    secondary.approx_line__const_coeff = NULL
    secondary.approx_line__slope = NULL
  } else {
    main.approx_line__const_coeff = trimmed.approx_line__const_coeff
    main.approx_line__slope = trimmed.approx_line__slope
    secondary.approx_line__const_coeff = approx_line__const_coeff
    secondary.approx_line__slope = approx_line__slope
  }
  
  equation.label = sprintf(
    'model: %s = %.3f * %s %s %.3f. Spearman r = %.2f', LogFC.col.name, main.approx_line__slope, bf.name,
    ifelse(main.approx_line__const_coeff > 0, "+", "-"), abs(main.approx_line__const_coeff),
    spearman_r)
  
  # logFC.range = 3
  # p.threshold = 0.05
  if(include.GLM.model.name.in.Result.names){
    filename = sprintf('%s/%s, 2D-hist, assoc. with %s%s%s.png', out.dir, Analysis.name, bf.name,
      ifelse(is.null(LogCPM.limit), "", sprintf(" (LogCPM gt. %g)", LogCPM.limit)), plots.suffix)
  } else {
    filename = sprintf('%s/Bias 2D-hist, assoc. with %s%s%s.png', out.dir, bf.name,
      ifelse(is.null(LogCPM.limit), "", sprintf(" (LogCPM gt. %g)", LogCPM.limit)), plots.suffix)
  }

  .pardefault <- par(no.readonly = TRUE)
  png.mod(filename = filename,
          width = 12, height = 10, pointsize = 12, res = 300)
  par(mar=c(6.5, 5, 4.5, 5))
  hist2d(x = BF.values, y = LogFC.values,
         nbins = 400, col = hist2d.cols, xlab = bf.name,
         ylab = LogFC.col.name, show = TRUE,
         xlim = BF.values.axis.limits, ylim = c(-logFC.range, logFC.range))
  title(main = sprintf('[%s]  %s vs %s  (slope = %.3f)%s', Analysis.name, LogFC.col.name, bf.name,
        main.approx_line__slope, ifelse(is.null(LogCPM.limit), "", sprintf(", LogCPM > %g", LogCPM.limit))), sub = equation.label)
  
  for(x in 1:(logFC.range-1)){
    abline(a = 0, b = 0, col = "gray60", lty = 3, h = x)
    abline(a = 0, b = 0, col = "gray60", lty = 3, h = -x)
  }
  abline(a = 0, b = 0, col = "black", lty = 2)
  if(!is.null(secondary.approx_line__const_coeff)){
    abline(a = secondary.approx_line__const_coeff, b = secondary.approx_line__slope, col = "gray80", lty = 2, lwd = 1.2)
  }
  abline(a = main.approx_line__const_coeff, b = main.approx_line__slope, col = "gray30", lty = 2, lwd = 1.5)
  dev.off()
  par(.pardefault)
  
  
  if(include.GLM.model.name.in.Result.names){
    filename = sprintf('%s/%s, bias plot, assoc. with %s%s%s.png', out.dir, Analysis.name, bf.name,
                             ifelse(is.null(LogCPM.limit), "", sprintf(" (LogCPM gt. %g)", LogCPM.limit)), plots.suffix)
  } else {
    filename = sprintf('%s/Bias plot, assoc. with %s%s%s.png', out.dir, bf.name,
                             ifelse(is.null(LogCPM.limit), "", sprintf(" (LogCPM gt. %g)", LogCPM.limit)), plots.suffix)
  }

  png.mod(filename = filename, width = 12, height = 10, pointsize = 12, res = 300)
  par(mar=c(6.5, 5, 4.5, 5))
  
  if(is.null(P.values)){
    plot(x = BF.values, y = LogFC.values,
           xlab = bf.name,
           ylab = LogFC.col.name,
           xlim = BF.values.axis.limits, ylim = c(-logFC.range, logFC.range), pch=1, cex = 0.4)
  } else {
    keep = P.values > p.threshold
    plot(x = BF.values[keep], y = LogFC.values[keep],
           xlab = bf.name,
           ylab = LogFC.col.name,
           xlim = BF.values.axis.limits, ylim = c(-logFC.range, logFC.range), pch=1, cex = 0.4)
    
    keep = P.values <= p.threshold & LogFC.values > 0
    points(x = BF.values[keep], y = LogFC.values[keep],
           pch=1, cex = 0.4, col = 'red')
    keep = P.values <= p.threshold & LogFC.values < 0
    points(x = BF.values[keep], y = LogFC.values[keep],
           pch=1, cex = 0.4, col = 'green')
  }  
  title(main = sprintf('[%s]  %s vs %s  (slope = %.3f)%s', Analysis.name, LogFC.col.name, bf.name,
                       main.approx_line__slope, ifelse(is.null(LogCPM.limit), "", sprintf(", LogCPM > %g", LogCPM.limit))), sub = equation.label)
  
  for(x in 1:(logFC.range-1)){
    abline(a = 0, b = 0, col = "gray60", lty = 3, h = x)
    abline(a = 0, b = 0, col = "gray60", lty = 3, h = -x)
  }
  abline(a = 0, b = 0, col = "black", lty = 2)
  if(!is.null(secondary.approx_line__const_coeff)){
    abline(a = secondary.approx_line__const_coeff, b = secondary.approx_line__slope, col = "gray80", lty = 2, lwd = 1.2)
  }
  abline(a = main.approx_line__const_coeff, b = main.approx_line__slope, col = "gray30", lty = 2, lwd = 1.5)
  dev.off()
  par(.pardefault)
}

Associations.LogFC.with.Bias.Factor__after.GLM = function(Startup.Data, Analysis.Data = NULL, GLM.model = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
                                       forced.DE.table = NULL, forced.Analysis.name = NULL, forced.parameters = NULL, gene.bias.factors.file = '{auto}',
                                       bias.factors = c('Avg. Transcript Length', 'Expression level'), out.dir = NULL, LogCPM.limits = c(0, 3, 5, 7), main.LogCPM.limit = 3,
                                       adjustable.bias.factor = 'Avg. Transcript Length', additional.plots = TRUE, script.dir = NULL, use.LogFC.col.name = '{auto}',
                                       db.suffix = '', plots.suffix = '', synchronize.pvalues = TRUE, include.GLM.model.name.in.Result.names = FALSE,
                                       LogFC.trimming = 0.05){
  

  if(LogFC.trimming >= 0.5){
    stop('LogFC.trimming cannot exceed 0.5')
  }

  if(!(main.LogCPM.limit %in% LogCPM.limits)){
    stop('main.LogCPM.limit should be present in LogCPM.limits')
  }

  allowed.bias.factors = c('Avg. Transcript Length', 'Avg. Transcript CG content', 'Expression level')
  absent.bias.factors = bias.factors[!(bias.factors %in% allowed.bias.factors)]
  if(length(absent.bias.factors) > 0){
    stop(sprintf('Unknown bias factors %s. Allowed only: %s', toString(absent.bias.factors), toString(allowed.bias.factors)))
  }

  if(!is.null(adjustable.bias.factor)){
    if(!adjustable.bias.factor %in% bias.factors){
      stop('adjustable.bias.factor should be present in bias.factors')
    }
  } else adjustable.bias.factor = 'none'

  if(gene.bias.factors.file == '{auto}'){
    if(is.null(script.dir)){
      tryCatch(expr = {
        script.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
      }, error = function (err){
        script.dir = '.'
        warning('rstudioapi failed to locate script directory')
      })
    }

    gene.bias.factors.file = sprintf('%s/Gene.Bias.Factors.tsv', script.dir)
    if(!file.exists(gene.bias.factors.file)){
      stop('Cannot find file with gene bias factors information (e.g. transcript length or CG content). Bias correction cannot not be applied. To enable bias correction place Gene.Bias.Factors.tsv (generated with PPLine) in the script directory')
    }
  }
  
  gene.bias.factors.table = read.table(gene.bias.factors.file, header = TRUE, check.names = FALSE, sep = '\t', quote = "")
  rownames(gene.bias.factors.table) = gene.bias.factors.table[,1]
  gene.bias.factors.table = gene.bias.factors.table[,-1]
  
  LogCPM.bf.present = 'Expression level' %in% bias.factors
  bias.factors = bias.factors[!(bias.factors %in% 'Expression level')]
    
  absent.bias.factors = bias.factors[!(bias.factors %in% colnames(gene.bias.factors.table))]
  if(length(absent.bias.factors) > 0){
    print(sprintf('%d bias factors are absent in %s file: %s', length(absent.bias.factors), gene.bias.factors.file, toString(absent.bias.factors)))
    bias.factors = bias.factors[(bias.factors %in% colnames(gene.bias.factors.table))]
  }
  
  if(LogCPM.bf.present) bias.factors = c(bias.factors, 'Expression level')


  if(is.null(forced.DE.table)){
    if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
      stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
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
      if(is.null(Analysis.Data.RDS.file))  Analysis.Data.RDS.file = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
      if(!file.exists(Analysis.Data.RDS.file)) stop(sprintf('RDS file %s is not found', Analysis.Data.RDS.file))
      Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
      Analysis.name = Analysis.Data$Analysis.name
      Analysis.Data$current.GLM.results.dir = GLM.working.dir
    } else {
      Analysis.name = Analysis.Data$Analysis.name
    }
  
  
    
    if(use.LogFC.col.name %in% c('auto', '{auto}')){
      LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part), preferred_columns = c('trimmed LogFC', 'LogFC', 'logFC'))
    } else {
      LogFC.col.name = use.LogFC.col.name
    }

    PValue.col.name = Get.PValue.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part))
    LogCPM.col.name = Get.LogCPM.col.name(colnames(Analysis.Data$ResTable.gene.part))
    common.genes = intersect(rownames(gene.bias.factors.table), rownames(Analysis.Data$ResTable.gene.part))
    DE.table = subset(Analysis.Data$ResTable.gene.part, select = c(LogFC.col.name, PValue.col.name, LogCPM.col.name))
    if(is.null(out.dir)){
      out.dir = Analysis.Data$current.GLM.results.dir
    }

  } else {
    DE.table = forced.DE.table
    if(dim(DE.table)[2] != 3){
      stop(sprintf('Incorrect dimensions of forced.DE.table (%d x %d). the second dimension should be 3', dim(forced.DE.table)[1], dim(forced.DE.table)[2]))
    }

    if(use.LogFC.col.name %in% c('auto', '{auto}')){
      LogFC.col.name = colnames(DE.table)[1]
    } else {
      LogFC.col.name = use.LogFC.col.name
    }

    
    PValue.col.name = colnames(DE.table)[2]
    LogCPM.col.name = colnames(DE.table)[3]
    common.genes = intersect(rownames(gene.bias.factors.table), rownames(DE.table))
    if(is.null(out.dir)) out.dir = getwd()
    Analysis.name = '(custom)'
  }


  if(!is.null(forced.Analysis.name)) Analysis.name = forced.Analysis.name
  cat(sprintf('[Associations.LogFC.with.Bias.Factor__after.GLM] Analyzing ~ %s...\n', Analysis.name))

  adjusted.DE.table = DE.table
  if(synchronize.pvalues){
    adjusted.DE.table = subset(adjusted.DE.table, select = c(colnames(adjusted.DE.table), LogFC.col.name, PValue.col.name, PValue.col.name))
    colnames(adjusted.DE.table)[4] = 'TBA LogFC'
    colnames(adjusted.DE.table)[5] = 'TBA p-value'
    colnames(adjusted.DE.table)[6] = 'TBA FDR'
  } else {
    adjusted.DE.table = subset(adjusted.DE.table, select = c(colnames(adjusted.DE.table), LogFC.col.name))
    colnames(adjusted.DE.table)[4] = 'TBA LogFC'
  }

  # DE.table = DE.table[common.genes, ]
  # gene.bias.factors.table = gene.bias.factors.table[common.genes, ]
  
  Bias.analysis.results = vector(mode = 'list')
  bf = bias.factors[1]
  for(bf in bias.factors){
    Bias.analysis.results[[bf]] = vector(mode = 'list')
    if(bf == 'Expression level'){
      c__LogCPM.limits = c(0)
    } else {
      c__LogCPM.limits = LogCPM.limits
    }
    current.LogCPM.limit = c__LogCPM.limits[1]
    for(current.LogCPM.limit in c__LogCPM.limits){
      if(bf == 'Expression level'){
        current.LogCPM.limit__char = '-'
        current.DE.table = DE.table[!(is.na(DE.table[,LogFC.col.name]) | is.na(DE.table[,PValue.col.name]) | is.na(DE.table[,LogCPM.col.name])), ]
        if(dim(current.DE.table)[1] < 4){
          cat(sprintf('[Associations.LogFC.with.Bias.Factor__after.GLM]  There is to few genes\n'))
          next
        }
        BF.values = current.DE.table[,LogCPM.col.name]
      } else {
        current.LogCPM.limit__char = as.character(current.LogCPM.limit)
        current.gene.bias.factors.table = gene.bias.factors.table[!is.na(gene.bias.factors.table[,bf]), ]

        current.DE.table = DE.table[!(is.na(DE.table[,LogFC.col.name]) | is.na(DE.table[,PValue.col.name]) | is.na(DE.table[,LogCPM.col.name])), ]
        current.DE.table = current.DE.table[current.DE.table[,LogCPM.col.name] > current.LogCPM.limit, ]

        if(dim(current.DE.table)[1] < 4){
          cat(sprintf('[Associations.LogFC.with.Bias.Factor__after.GLM]  There is to few genes with LogCPM > %g\n', current.LogCPM.limit))
          next
        }

        current.common.genes = intersect(rownames(current.gene.bias.factors.table), rownames(current.DE.table))

        if(length(current.common.genes) < 4){
          cat(sprintf('[Associations.LogFC.with.Bias.Factor__after.GLM]  There is to few common genes with LogCPM > %g\n', current.LogCPM.limit))
          next
        }

        current.gene.bias.factors.table = current.gene.bias.factors.table[current.common.genes, ]
        current.DE.table = current.DE.table[current.common.genes, ]

        BF.values = current.gene.bias.factors.table[,bf]
      }

      LogFC.values = current.DE.table[,LogFC.col.name]
      P.values = current.DE.table[,PValue.col.name]
      ct.spearman = suppressWarnings(cor.test(BF.values, LogFC.values, method = 'spearman'))
      ct.pearson = suppressWarnings(cor.test(BF.values, LogFC.values, method = 'pearson'))
      
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]] = vector(mode = 'list')
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['spearman.r']] = ct.spearman$estimate
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['spearman.p']] = ct.spearman$p.value
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['pearson.r']] = ct.pearson$estimate
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['pearson.p']] = ct.pearson$p.value
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['BF.values']] = BF.values
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['LogFC.values']] = LogFC.values
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['P.values']] = P.values
      
      sorted.LogFC.values = LogFC.values[order(LogFC.values)]
      min.LogFC.for.trimming = sorted.LogFC.values[min(length(LogFC.values), 1 + round(LogFC.trimming * length(LogFC.values)))]
      max.LogFC.for.trimming = sorted.LogFC.values[max(0, length(LogFC.values) - round(LogFC.trimming * length(LogFC.values)))]
      keep = LogFC.values >= min.LogFC.for.trimming & LogFC.values <= max.LogFC.for.trimming
      trimmed.LogFC.values = LogFC.values[keep]
      trimmed.BF.values = BF.values[keep]


      if (bf == 'Avg. Transcript Length'){
        fit = lm(LogFC.values ~ log2(BF.values))
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.interceipt']] = fit$coefficients[1]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.slope']] = fit$coefficients[2]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.pvalue']] = lmp(fit)

        trimmed.fit = lm(trimmed.LogFC.values ~ log2(trimmed.BF.values))
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.interceipt']] = trimmed.fit$coefficients[1]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.slope']] = trimmed.fit$coefficients[2]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.pvalue']] = lmp(trimmed.fit)

        #adjusted.LogFC.values = LogFC.values
        if(bf == adjustable.bias.factor & current.LogCPM.limit == main.LogCPM.limit){
          sub.gene.bias.factors.table = gene.bias.factors.table[!is.na(gene.bias.factors.table[,bf]) & (rownames(gene.bias.factors.table) %in% rownames(adjusted.DE.table)), ]
          
          adjusted.DE.table[rownames(sub.gene.bias.factors.table), 'TBA LogFC'] =
            adjusted.DE.table[rownames(sub.gene.bias.factors.table), LogFC.col.name] - trimmed.fit$coefficients[2] * log2(sub.gene.bias.factors.table[, bf]) - trimmed.fit$coefficients[1]

          if(synchronize.pvalues){
            all.abs.LogFCs  = abs(adjusted.DE.table[, LogFC.col.name])
            all.PValues = adjusted.DE.table[, PValue.col.name]

            abs.LogFC.bins = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.6, 2.9, 3.4, 4.0, 5.0, 1000)
            avg.PValue.per.bin = rep(1, length(abs.LogFC.bins) - 1)
            last.avg.PValue = 1
            for(bin_n in 1:(length(abs.LogFC.bins) - 1)){
              keep = all.abs.LogFCs >= abs.LogFC.bins[bin_n] & all.abs.LogFCs < abs.LogFC.bins[bin_n + 1]
              if(sum.mod(keep) == 0){
                avg.PValue.per.bin[bin_n] = last.avg.PValue
              } else {
                avg.PValue.per.bin[bin_n] = gm_mean(mean(all.PValues[keep]), gm_mean(all.PValues[keep]))
                last.avg.PValue = avg.PValue.per.bin[bin_n]
              }
            }
            for(x in 1:dim(adjusted.DE.table)[1]){
              current.LogFC = abs(adjusted.DE.table[x, 'TBA LogFC'])
              bin_n_found = NA
              for(bin_n in 1:(length(abs.LogFC.bins) - 1)){
                if(current.LogFC >= abs.LogFC.bins[bin_n] & current.LogFC < abs.LogFC.bins[bin_n + 1]){
                  bin_n_found = bin_n
                  break
                }
              }
              if(!is.na(bin_n_found)){  adjusted.DE.table[x, 'TBA p-value'] = max(adjusted.DE.table[x, 'TBA p-value'], avg.PValue.per.bin[bin_n_found])
              } else   adjusted.DE.table[x, 'TBA p-value'] = 1
            }
            adjusted.DE.table[x, 'TBA FDR'] = p.adjust(adjusted.DE.table[x, 'TBA p-value'], method = 'BH')
          }
        }
        
        if(additional.plots){
          if(current.LogCPM.limit == main.LogCPM.limit){
            Bias.Factor.LogFC.Plots(approx_line__const_coeff = fit$coefficients[1], approx_line__slope = fit$coefficients[2], spearman_r = ct.spearman$estimate,
                              BF.values = log2(BF.values), LogFC.values = LogFC.values, P.values = P.values, 
                              bf.name = 'Log[Transcript Length]', Analysis.name = Analysis.name, LogFC.col.name = LogFC.col.name,
                              out.dir = out.dir, LogCPM.limit = current.LogCPM.limit, logFC.range = 3, p.threshold = 0.05,
                              plots.suffix = plots.suffix,
                              include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                              trimmed.approx_line__const_coeff = trimmed.fit$coefficients[1], trimmed.approx_line__slope = trimmed.fit$coefficients[2])
          }
          current.out.dir = sprintf('%s/Bias Factor plots', out.dir)
          dir.create(current.out.dir, showWarnings = FALSE)
          Bias.Factor.LogFC.Plots(approx_line__const_coeff = fit$coefficients[1], approx_line__slope = fit$coefficients[2], spearman_r = ct.spearman$estimate,
                            BF.values = log2(BF.values), LogFC.values = LogFC.values, P.values = P.values, 
                            bf.name = 'Log[Transcript Length]', Analysis.name = Analysis.name, LogFC.col.name = LogFC.col.name,
                            out.dir = current.out.dir, LogCPM.limit = current.LogCPM.limit, logFC.range = 3, p.threshold = 0.05,
                            plots.suffix = plots.suffix,
                            include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                            trimmed.approx_line__const_coeff = trimmed.fit$coefficients[1], trimmed.approx_line__slope = trimmed.fit$coefficients[2])
        }

      } else if (bf == 'Expression level'){
        fit = lm(LogFC.values ~ BF.values)
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.interceipt']] = fit$coefficients[1]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.slope']] = fit$coefficients[2]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.pvalue']] = lmp(fit)

        trimmed.fit = lm(trimmed.LogFC.values ~ trimmed.BF.values)
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.interceipt']] = trimmed.fit$coefficients[1]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.slope']] = trimmed.fit$coefficients[2]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.pvalue']] = lmp(trimmed.fit)

        if(additional.plots)
          Bias.Factor.LogFC.Plots(approx_line__const_coeff = fit$coefficients[1], approx_line__slope = fit$coefficients[2], spearman_r = ct.spearman$estimate,
                          BF.values = BF.values, LogFC.values = LogFC.values, P.values = P.values, 
                          bf.name = 'Expression level', Analysis.name = Analysis.name, LogFC.col.name = LogFC.col.name,
                          out.dir = out.dir, LogCPM.limit = current.LogCPM.limit, logFC.range = 3, p.threshold = 0.05, BF.values.axis.limits = c(0, 15),
                          plots.suffix = plots.suffix,
                          include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                          trimmed.approx_line__const_coeff = trimmed.fit$coefficients[1], trimmed.approx_line__slope = trimmed.fit$coefficients[2])
        
      } else if (bf == 'Avg. Transcript CG content'){
        stop(sprintf('Avg. Transcript CG content is not supported now'))
      } else {
        stop(sprintf('Unknown bias bactor %s', bf))
      }
    }
    #current.gene.bias.factors.table[order(BF.values)[1:10],]
  }
  if(include.GLM.model.name.in.Result.names){
    saveRDS(Bias.analysis.results, file = sprintf("%s/%s, bias.analysis.results%s.rds", out.dir, Analysis.name, db.suffix))
  } else {
    saveRDS(Bias.analysis.results, file = sprintf("%s/bias.analysis.results%s.rds", out.dir, db.suffix))
  }
  return(invisible(adjusted.DE.table))
  
}


Summarize.Bias.Analysis.Results = function(Startup.Data, GLM.models = NULL,
                                           forced.parameters = NULL,
                                           out.dir = NULL, out.file.prefix = '', DB.suffix = ''){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters
  if (length(Startup.Data$topGO.Expression.Profiles___Custom.GO.terms) == 0) return()
  
  if (!Pars$Create.Excel.results){
    warning('Create.Excel.results parameter is turned OFF. "Summarize.Bias.Analysis.Results" will not run')
    return()
  }
  
  if(is.null(GLM.models)){
    GLM.models = Pars$Models.to.Test
  }
  
  if(is.null(out.dir)) out.dir = Pars$results.dir
  
  all.Bias.analysis.results = vector(mode = 'list')
  all.bias.factors = c()
  absent.file.names = c()
  for (GLM.model in GLM.models){
    GLM.model = Delete.prefix.in.model(GLM.model)
    Analysis.name  = Cleanup.Model.Name(GLM.model)
    file.name = sprintf('%s/~ %s, results/bias.analysis.results%s.rds', Pars$results.dir, Analysis.name, DB.suffix)
    if(file.exists(file.name)){
      all.Bias.analysis.results[[GLM.model]] = readRDS(file = file.name)
    } else {
      absent.file.names = c(absent.file.names, file.name)
    }
    
    all.bias.factors = c(all.bias.factors, names(all.Bias.analysis.results[[GLM.model]]))
  }
  
  if(length(absent.file.names) > 0){
    cat(sprintf('[Summarize.Bias.Analysis.Results] Total %d files with Bias analysis results do not exist: %s\n',
                length(absent.file.names),  toString(absent.file.names)))
    if(length(all.Bias.analysis.results) == 0){
      cat('No files found\n')
      return(invisible())
    }
  } else {
    cat(sprintf('[Summarize.Bias.Analysis.Results] All %d files with Bias analysis results exist\n', length(GLM.models)))
  }
  
  all.bias.factors %<>% .[!duplicated(.)]
  
  cat(sprintf('%d bias factor(s) found: %s\n', length(all.bias.factors), toString(all.bias.factors)))

  all.LogCPM.limits__by.bias.factor = vector(mode = 'list')
  for(BFA.an in names(all.Bias.analysis.results)){
    for(bias.factor in names(all.Bias.analysis.results[[BFA.an]])){
      if(!(bias.factor %in% names(all.LogCPM.limits__by.bias.factor))){
        all.LogCPM.limits__by.bias.factor[[bias.factor]] = c()
      }
      #all.Bias.analysis.results[[BFA.an]]
      for(LogCPM.limit in names(all.Bias.analysis.results[[BFA.an]][[bias.factor]])){
        if(!(LogCPM.limit %in% all.LogCPM.limits__by.bias.factor[[bias.factor]])){
          all.LogCPM.limits__by.bias.factor[[bias.factor]] = c(all.LogCPM.limits__by.bias.factor[[bias.factor]], LogCPM.limit)
        }
      }
    }
  }
  
  bias.factor = all.bias.factors[1]
  for(bias.factor in all.bias.factors){
    all.out.tsv.files = c()
    all.Excel.sheet.names = c()
    for(LogCPM.limit in all.LogCPM.limits__by.bias.factor[[bias.factor]]){
      BFA.table = data.frame(array(dim = c(length(all.Bias.analysis.results), 8)))
      colnames(BFA.table) = c('Spearman r', 'Spearman p-value', 'LM slope', 'LM interceipt', 'LM p-value', 'trimmed LM slope', 'trimmed LM interceipt', 'trimmed LM p-value')
      rownames(BFA.table) = names(all.Bias.analysis.results) #GLM.models[GLM.models %in% names(all.Bias.analysis.results)]
      BFA.an = names(all.Bias.analysis.results)[1]
      for(BFA.an in names(all.Bias.analysis.results)){
        current.res = all.Bias.analysis.results[[BFA.an]][[bias.factor]][[LogCPM.limit]]
        BFA.table[BFA.an, ] = c(
          current.res$spearman.r,
          current.res$spearman.p,
          current.res$lm.slope,
          current.res$lm.interceipt,
          current.res$lm.pvalue,
          current.res$trimmed.lm.slope,
          current.res$trimmed.lm.interceipt,
          current.res$trimmed.lm.pvalue
          )
      }
      
      out.tsv.file = sprintf('%s/%s(LogCPM gt. %s) %s.tsv', out.dir, out.file.prefix, LogCPM.limit, bias.factor)
      write.table.mod(BFA.table, out.tsv.file, quote = FALSE, sep = '\t')
      all.out.tsv.files = c(all.out.tsv.files, Verify.path(out.tsv.file))
      if(LogCPM.limit == '-' & length(all.LogCPM.limits__by.bias.factor[[bias.factor]]) == 1){
        all.Excel.sheet.names = c(all.Excel.sheet.names, bias.factor)
      } else {
        all.Excel.sheet.names = c(all.Excel.sheet.names, sprintf('(LogCPM>%s) %s', LogCPM.limit, bias.factor))
      }
    }
    out.Excel.file = Verify.path(sprintf('%s/%sassoc. with %s.xlsx', out.dir, out.file.prefix, bias.factor))
    CL = sprintf('%s "%s/Any.2.Excel.py" --in %s --out-excel "%s" --one-book yes --sheet-names %s --logfc-abs-max 0.7 --simple-logfc-format no --predictor-rows-count 0 --forced-col-types-by-names "LogFC:&LM slope&,&trimmed LM slope&" "p:&Spearman p-value&,&LM p-value&,&trimmed LM p-value&" --entry-type "Model" --first-col-width 19 --first-col-italic no',
                 Pars$python.bin, Pars$suppl.data.dir, paste(sprintf('"%s"', all.out.tsv.files), collapse = ' '),
                 out.Excel.file, paste(sprintf('"%s"', all.Excel.sheet.names), collapse = ' '))
    code = system(CL)
    if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
    
    if (!file.exists(out.Excel.file)){
      message(sprintf('[Summarize.Bias.Analysis.Results]  Excel worksheet generation (Bias factor analysis) for %s FAILED', bias.factor))
    } else {
      file.remove(all.out.tsv.files)
    }
  }
  

  #setwd(Pars$results.dir)

}

Bias.Factor.LogCPM.Plots = function(sample_name, equation.label, BF.values, y.Values, 
                             bf.name = 'Avg. Transcript Length' , out.dir, y.axis.label = 'LogCPM', CPM.limit = NULL,
                             approx.curve.points = NULL,
                             BF.values.axis.limits = c(8, 14), y.axis.limts = c(0, 12),
                             BF.values.axis.limits_ke2d = NULL, y.axis.limts_ke2d = NULL,
                             create.dotplot = TRUE, create_hist2D = TRUE, create_density2D = TRUE,
                             plots.suffix = '', nbins = 80, render.vertical.lines.at.x.positions = NULL){
  
  if(is.null(BF.values.axis.limits_ke2d))  BF.values.axis.limits_ke2d = BF.values.axis.limits
  if(is.null(y.axis.limts_ke2d))  y.axis.limts_ke2d = y.axis.limts

  hist2d.cols <- colorRampPalette(c("#ffffff","#4bbc00","#e59400", "#f77700", "#f70000", "#f70000"))(512)
  hist2d.cols.v2 <- colorRampPalette(c("#ffffff",  '#a0c9e7',  '#789ae7', "#e59400", "#f77700", "#f70000", "#f70000"))(512)
  
  if(create_hist2D){
    filename = sprintf('%s/[2D-hist] %s, assoc. %s - %s%s%s.png', out.dir, sample_name, y.axis.label, bf.name,
      ifelse(is.null(CPM.limit), "", sprintf(" (CPM gt. %g)", CPM.limit)), plots.suffix)

    .pardefault <- par(no.readonly = TRUE)
    png.mod(filename = filename,
            width = 12, height = 10, pointsize = 12, res = 300)
    par(mar=c(6.5, 5, 4.5, 5))
    hist2d(x = BF.values, y = y.Values,
           nbins = nbins, col = hist2d.cols.v2, xlab = bf.name,
           ylab = y.axis.label, show = TRUE,
           xlim = BF.values.axis.limits, ylim = y.axis.limts)
    if(!is.null(approx.curve.points))    lines(approx.curve.points, col='red', lwd=2)
    title(main = sprintf('%s,  %s vs %s%s', sample_name, y.axis.label, bf.name,
          ifelse(is.null(CPM.limit), "", sprintf(", CPM > %g", CPM.limit))), sub = equation.label)

    if(!is.null(render.vertical.lines.at.x.positions)){
      for(pos in render.vertical.lines.at.x.positions)   abline(v =pos, col = "gray60", lty = 2, lwd = 0.75)
    }
    
    # for(x in 1:(logFC.range-1)){
    #   abline(a = 0, b = 0, col = "gray60", lty = 3, h = x)
    #   abline(a = 0, b = 0, col = "gray60", lty = 3, h = -x)
    # }
    # abline(a = 0, b = 0, col = "black", lty = 2)
    # abline(a = approx_line__const_coeff, b = approx_line__slope, col = "gray30", lty = 2, lwd = 1.5)
    dev.off()
    par(.pardefault)
  }

  
  if(create.dotplot){
    filename = sprintf('%s/[dot plot] %s, assoc. with %s%s%s.png', out.dir, sample_name, bf.name,
                             ifelse(is.null(CPM.limit), "", sprintf(" (CPM gt. %g)", CPM.limit)), plots.suffix)

    .pardefault <- par(no.readonly = TRUE)
    png.mod(filename = filename, width = 12, height = 10, pointsize = 12, res = 300)
    par(mar=c(6.5, 5, 4.5, 5))
    
    plot(x = BF.values, y = y.Values,
           xlab = bf.name,
           ylab = y.axis.label,
           xlim = BF.values.axis.limits, ylim = y.axis.limts, pch=1, cex = 0.4)

    if(!is.null(approx.curve.points))    lines(approx.curve.points, col='red', lwd=2)
    
    title(main = sprintf('%s,  %s vs %s%s', sample_name, y.axis.label, bf.name,
                         ifelse(is.null(CPM.limit), "", sprintf(", CPM > %g", CPM.limit))), sub = equation.label)
    
    if(!is.null(render.vertical.lines.at.x.positions)){
      for(pos in render.vertical.lines.at.x.positions)   abline(v =pos, col = "gray60", lty = 2, lwd = 1.5)
    }

    # for(x in 1:(logFC.range-1)){
    #   abline(a = 0, b = 0, col = "gray60", lty = 3, h = x)
    #   abline(a = 0, b = 0, col = "gray60", lty = 3, h = -x)
    # }
    # abline(a = 0, b = 0, col = "black", lty = 2)
    # abline(a = approx_line__const_coeff, b = approx_line__slope, col = "gray30", lty = 2, lwd = 1.5)
    dev.off()
    par(.pardefault)
  }

  if(create_density2D){
    g = suppressWarnings(ggplot(as.data.frame(cbind(BF.values, y.Values)), aes(x=BF.values, y=y.Values) ) + stat_density_2d(aes(fill = ..level.., colour = "red"), geom = "polygon", bins = 15) + 
          xlim(BF.values.axis.limits_ke2d) + ylim(y.axis.limts_ke2d) + scale_fill_gradientn(colours = rev( brewer.pal( 7, "Spectral" ) )) +
          xlab(bf.name) + ylab(y.axis.label) + ggtitle(label = sprintf('%s,  %s vs %s%s', sample_name, y.axis.label, bf.name,
                                                                       ifelse(is.null(CPM.limit), "", sprintf(", CPM > %g", CPM.limit))),
                                                       subtitle = equation.label))
    if(!is.null(render.vertical.lines.at.x.positions)){
      g = g + geom_vline(xintercept = render.vertical.lines.at.x.positions, linetype=2, color = "gray60", size = 0.25)
    }

    filename = sprintf('%s/[2D-dens] %s, assoc. %s - %s%s%s.png', out.dir, sample_name, y.axis.label, bf.name,
      ifelse(is.null(CPM.limit), "", sprintf(" (CPM gt. %g)", CPM.limit)), plots.suffix)
    ggsave(filename = filename,  plot = g, width = 12, height = 10, units = 'in', dpi = 300, limitsize = FALSE)    
  }
}


Get.bins = function(bf){
  if(bf == 'Avg. Transcript Length'){
    bins = c(0, 8, 9.6, 10.0, 10.3, 10.5, 10.8, 11.0, 11.2, 11.5, 11.8, 15)
    bins = c(0, 8, 9.6, 10.0, 10.3, 10.5, 10.8, 11.1, 11.4, 11.8, 20)
    bins = c(0, 8, 9.6, 10.0, 10.35, 10.7, 11.05, 11.4, 11.8, 20)
  } else if(bf == 'Avg. Transcript CG Content'){
    bins = c(0, 30, 40, 50, 60, 70, 100)
  } else if(bf == 'Avg. Transcript CG Content (3\'-tail 1kb)'){
    bins = c(0, 30, 40, 50, 60, 70, 100)
  } else {
    stop(sprintf('bias factor %s is not supported now'))
  }
}      


analyze.LogCPM.each.Sample.vs.Bias.Factor.associations <- function(counts.table, Startup.Data, out.dir = NULL, forced.parameters = NULL, gene.bias.factors.file = '{auto}',
                                       bias.factors = c('Avg. Transcript Length'), adjustable.bias.factor = 'Avg. Transcript Length',
                                       min.CPM.for.association.analysis = 1.0,
                                       render.not.adjusted.accos.plots = TRUE, render.adjustment.CPM.density.plots.all.genes = TRUE,
                                       render.adjustment.CPM.density.plots = TRUE,
                                       script.dir = NULL, results.prefix = '', forced.Analysis.name = NULL,
                                       include.GLM.model.name.in.Result.names = FALSE,
                                       min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = 1.6,
                                       min.cpms.stdev_abs_value.in.bin.to.exclude = 0.007,
                                       max.bins.to.exclude.num = 2,
                                       max.bins.to.exclude.percent = 25){


  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters

  if('Expression level' %in% bias.factors)  bias.factors = bias.factors[!(bias.factors %in% 'Expression level')]
  allowed.bias.factors = c('Avg. Transcript Length', 'Avg. Transcript CG content')
  absent.bias.factors = bias.factors[!(bias.factors %in% allowed.bias.factors)]
  if(length(absent.bias.factors) > 0){
    stop(sprintf('Unknown bias factors %s. Allowed only: %s', toString(absent.bias.factors), toString(allowed.bias.factors)))
  }

  if(!is.null(adjustable.bias.factor)){
    if(!adjustable.bias.factor %in% bias.factors){
      stop('adjustable.bias.factor should be present in bias.factors')
    }
  } else {
    adjustable.bias.factor = ''
  }

  if(gene.bias.factors.file == '{auto}'){
    if(is.null(script.dir)){
      tryCatch(expr = {
        script.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
      }, error = function (err){
        script.dir = '.'
        warning('rstudioapi failed to locate script directory')
      })
    }

    gene.bias.factors.file = sprintf('%s/Gene.Bias.Factors.tsv', script.dir)
    if(!file.exists(gene.bias.factors.file)){
      stop('Cannot find file with gene bias factors information (e.g. transcript length or CG content). Bias correction cannot not be applied. To enable bias correction place Gene.Bias.Factors.tsv (generated with PPLine) in the script directory')
    }
  }
  
  gene.bias.factors.table = read.table(gene.bias.factors.file, header = TRUE, check.names = FALSE, sep = '\t', quote = "")
  rownames(gene.bias.factors.table) = gene.bias.factors.table[,1]
  gene.bias.factors.table = gene.bias.factors.table[,-1]
   
  absent.bias.factors = bias.factors[!(bias.factors %in% colnames(gene.bias.factors.table))]
  if(length(absent.bias.factors) > 0){
    print(sprintf('%d bias factors are absent in %s file: %s', length(absent.bias.factors), gene.bias.factors.file, toString(absent.bias.factors)))
    bias.factors = bias.factors[(bias.factors %in% colnames(gene.bias.factors.table))]
  }
  
  common.genes = intersect(rownames(gene.bias.factors.table), rownames(counts.table))
  if(is.null(out.dir)) out.dir = getwd()

  if(!is.null(forced.Analysis.name)){
    Analysis.name = forced.Analysis.name
  } else {
    Analysis.name = '(custom)'
  }
  cat(sprintf('[analyze.LogCPM.each.Sample.vs.Bias.Factor.associations] Analyzing ~ %s...\n', Analysis.name))
  adjusted_pseudo.counts.table = NULL
  CPM.table = cpm(counts.table)

  all.tsv.file.names = c()
  Bias.analysis.ES.results = vector(mode = 'list')
  bf = bias.factors[1]
  for(bf in bias.factors){
    Bias.analysis.ES.results[[bf]] = vector(mode = 'list')
    assoc.parameters = c('Spearman r', 'Spearman p-value', 'G-asymmetry', 'A-asymmetry')
    assoc.table = data.frame(array(dim = c(dim(counts.table)[2], length(assoc.parameters)) ))
    colnames(assoc.table) = assoc.parameters
    rownames(assoc.table) = colnames(counts.table)
    snn = 1
    for(snn in 1:ncol(counts.table)){
      sn = colnames(counts.table)[snn]
      cat(sprintf('\r%s - Analyzing bias and rendering plots. Sample %d of %d...      ', bf, snn, ncol(counts.table)))
      current.gene.bias.factors.table = gene.bias.factors.table[!is.na(gene.bias.factors.table[,bf]), ]
      current.CPMs = subset(CPM.table, select = sn, subset = CPM.table[,sn] > min.CPM.for.association.analysis)
      common.genes = intersect(rownames(current.CPMs), rownames(current.gene.bias.factors.table))
      BF.values = current.gene.bias.factors.table[common.genes, bf]
      current.CPMs = current.CPMs[common.genes,]

      if(length(common.genes) < 4){
        cat(sprintf('[analyze.LogCPM.each.Sample.vs.Bias.Factor.associations]  There is to few genes with LogCPM > %g for sample\n', min.CPM.for.association.analysis, sn))
        next
      }

      ct.spearman = suppressWarnings(cor.test(BF.values, current.CPMs, method = 'spearman'))
      # ct.pearson = suppressWarnings(cor.test(BF.values, current.CPMs, method = 'pearson'))
      
      Bias.analysis.ES.results[[bf]][[sn]] = vector(mode = 'list')
      Bias.analysis.ES.results[[bf]][[sn]][['spearman.r']] = ct.spearman$estimate
      Bias.analysis.ES.results[[bf]][[sn]][['spearman.p']] = ct.spearman$p.value
      assoc.table[sn, 'Spearman r'] = ct.spearman$estimate
      assoc.table[sn, 'Spearman p-value'] = ct.spearman$p.value
      assoc.table[sn, 'G-asymmetry'] = mean(log2(current.CPMs[BF.values < median(BF.values)])) - mean(log2(current.CPMs[BF.values > median(BF.values)]))
      assoc.table[sn, 'A-asymmetry'] = log2(mean(current.CPMs[BF.values < median(BF.values)]) / mean(current.CPMs[BF.values > median(BF.values)]))
      
      # approx.curve = loess(log2(current.CPMs) ~ log2(BF.values), span = 0.8, degree = 4)
      # approx.curve.points = predict(approx.curve)
      # ggplot(as.data.frame(cbind(x,y)), aes(x=x, y=y) ) + geom_density_2d(aes(color = ..level..), bins = 25) + xlim(8, 14) + ylim(0, 12) + 
      #   scale_fill_gradientn(colours = rev( brewer.pal( 7, "Spectral" ) ))
      # ggplot(as.data.frame(cbind(x,y)), aes(x=x, y=y) ) + stat_density_2d(aes(fill = ..level.., colour = "red"), geom = "polygon", bins = 15) + 
      #   xlim(8, 14) + ylim(-0.5, 12) + scale_fill_gradientn(colours = rev( brewer.pal( 7, "Spectral" ) ))
      

      if(render.not.adjusted.accos.plots){
        current.out.dir = sprintf('%s/%sBias Factor plots - per sample (%s)', out.dir, results.prefix, bf)
        dir.create(current.out.dir, showWarnings = FALSE)
        equation.label = sprintf('Spearman r = %.2f (p = %g), G-asymmetry = %.2f', ct.spearman$estimate, ct.spearman$p.value, assoc.table[sn, 'G-asymmetry'])

        if(bf == adjustable.bias.factor){  bins = Get.bins(bf)
        } else bins = NULL

        if (bf == 'Avg. Transcript Length'){

          Bias.Factor.LogCPM.Plots(sample_name = sn, equation.label = equation.label, BF.values = log2(BF.values), y.Values = log2(current.CPMs),
                             bf.name = 'Log [Avg. Transcript Length]' , out.dir = current.out.dir, y.axis.label = 'LogCPM', CPM.limit = NULL,
                             approx.curve.points = NULL,
                             BF.values.axis.limits = c(8, 14), y.axis.limts = c(0, 12),
                             BF.values.axis.limits_ke2d = c(8.5, 13), y.axis.limts_ke2d = c(-0.5, 10.5),
                             plots.suffix = '', render.vertical.lines.at.x.positions= bins)

        } else if (bf == 'Avg. Transcript CG content'){
          Bias.Factor.LogCPM.Plots(sample_name = sn, equation.label = equation.label, BF.values = BF.values, y.Values = log2(current.CPMs),
                             bf.name = bf, out.dir = current.out.dir, y.axis.label = 'LogCPM', CPM.limit = NULL,
                             approx.curve.points = NULL,
                             BF.values.axis.limits = c(0.2, 0.8), y.axis.limts = c(0, 12), plots.suffix = '', render.vertical.lines.at.x.positions= bins)
        } else if (bf == 'Avg. Transcript CG content (3\'-tail 1kb)'){
          Bias.Factor.LogCPM.Plots(sample_name = sn, equation.label = equation.label, BF.values = BF.values, y.Values = log2(current.CPMs),
                             bf.name = bf, out.dir = current.out.dir, y.axis.label = 'LogCPM', CPM.limit = NULL,
                             approx.curve.points = NULL,
                             BF.values.axis.limits = c(0.2, 0.8), y.axis.limts = c(0, 12), plots.suffix = '', render.vertical.lines.at.x.positions= bins)
        } else {
          stop(sprintf('Unknown bias bactor %s', bf))
        }
      }
    }
    cat(sprintf('\rCompleted                                                                       \n'))

    tsv.file.name = Verify.path(sprintf("%s/%sAssociations with %s, per sample.tsv", out.dir, results.prefix, bf))
    write.table.mod(x = assoc.table, file = tsv.file.name, quote = FALSE, sep = '\t')
    all.tsv.file.names = c(all.tsv.file.names, tsv.file.name)

    if(bf == adjustable.bias.factor){
      bins = Get.bins(bf)
      if(results.prefix == ''){
        current.out.dir = sprintf('%s/Adjusting for %s', out.dir, bf)
      } else  current.out.dir = sprintf('%s/%s, adjusting for %s', out.dir, results.prefix, bf)
      dir.create(current.out.dir, showWarnings = FALSE)

      d.all.genes = DGEList(counts = counts.table)
      d.all.genes = calcNormFactors(d.all.genes, method=Pars$RNA.Seq.norm.method)
      not.adjusted.cpm.all.genes = cpm(d.all.genes, normalized.lib.sizes = TRUE)

      if(bf == 'Avg. Transcript Length'){  printed.BF.name = 'Tr.Length'
      } else printed.BF.name = bf

      if(render.adjustment.CPM.density.plots.all.genes)
        Render.Density.Plots__color.by.sample(not.adjusted.cpm.all.genes, out.dir = current.out.dir,
          plot.title = sprintf('[before %s adj.] All genes, CPM density, normalized (%s)', printed.BF.name, Pars$RNA.Seq.norm.method), plot.subtitle = '',
          png.file.name = sprintf('[before %s adj.] All genes, Log2(CPM) density, normalized (%s).png', printed.BF.name, Pars$RNA.Seq.norm.method))

      read.counts.per.bin.per.sample = data.frame(array( dim = c(length(bins) - 1, ncol(counts.table))))
      rownames(read.counts.per.bin.per.sample) = 1:(length(bins) - 1)
      colnames(read.counts.per.bin.per.sample) = colnames(counts.table)

      norm.factors.per.bin.per.sample = data.frame(array( dim = c(length(bins) - 1, ncol(counts.table))))
      rownames(norm.factors.per.bin.per.sample) = 1:(length(bins) - 1)
      colnames(norm.factors.per.bin.per.sample) = colnames(counts.table)

      stdev.per.bin = c()

      adjusted_pseudo.counts.table.per.bin = vector(mode = 'list')

      sum_counts.table = sum(counts.table)


      
      bin_n = 1
      for(bin_n in 1:(length(bins) - 1)){
        current.gene.bias.factors.table = gene.bias.factors.table[!is.na(gene.bias.factors.table[,bf]), ]
        if(bf == 'Avg. Transcript Length'){ real.BF.values = log2(gene.bias.factors.table[,bf])
        } else  real.BF.values = gene.bias.factors.table[,bf]
        
        current.gene.bias.factors.table = gene.bias.factors.table[real.BF.values >= bins[bin_n] & real.BF.values < bins[bin_n + 1], ]
        common.genes = intersect(rownames(current.gene.bias.factors.table), rownames(counts.table))
        current.counts.table = counts.table[common.genes, ]
        sum_current.counts.table = sum(current.counts.table)
        
        if(length(common.genes) < 20){
          cat(sprintf('bin %d: to few genes in bin... Bypassing...\n', bin_n))
          adjusted_pseudo.counts.table.per.bin[[bin_n]] = counts.table[c(), ]
          stdev.per.bin = c(stdev.per.bin, 0)
          read.counts.per.bin.per.sample[bin_n, ] = 0
          norm.factors.per.bin.per.sample[bin_n, ] = 1
          next
        }
        

        d.bin = DGEList(counts = current.counts.table)
        d.bin = calcNormFactors(d.bin, method=Pars$RNA.Seq.norm.method)
        read.counts.per.bin.per.sample[bin_n, ] = d.bin$samples$lib.size
        norm.factors.per.bin.per.sample[bin_n, ] = d.bin$samples$norm.factors
        
        current.adjusted_pseudo.counts.table = t(t(current.counts.table) / d.bin$samples$lib.size * mean(d.bin$samples$lib.size) / d.bin$samples$norm.factors )
        adjusted_pseudo.counts.table.per.bin[[bin_n]] = current.adjusted_pseudo.counts.table
        # adjusted_pseudo.counts.table = rbind(adjusted_pseudo.counts.table, current.adjusted_pseudo.counts.table)

        adjusted.cpm.current.bin = cpm(d.bin, normalized.lib.sizes = TRUE) * sum_current.counts.table / sum_counts.table
        current.stdev = Get.Density.table.avg.StDev(adjusted.cpm.current.bin)
        stdev.per.bin = c(stdev.per.bin, current.stdev)

        if(bf == 'Avg. Transcript Length'){
          insertion = sprintf('transcript length from %.0f to %.0f bp', 2**bins[bin_n], 2**bins[bin_n + 1])
        } else {
          insertion = sprintf('%s from %.0f to %.0f', bf, bins[bin_n], bins[bin_n + 1])
        }
        cat(sprintf('bin %d, %s: %d genes, LogCPM StDev = %g\n', bin_n, insertion, length(common.genes), current.stdev))

        if(render.adjustment.CPM.density.plots){
          Render.Density.Plots__color.by.sample(not.adjusted.cpm.all.genes[common.genes, ], out.dir = current.out.dir,
            plot.title = sprintf('[before %s adj.] Bin %d (%s)\nLog2(CPM) density, normalized (%s)', printed.BF.name, bin_n, insertion, Pars$RNA.Seq.norm.method), plot.subtitle = '',
            png.file.name = sprintf('[before %s adj.] Bin %2.f - CPM density, normalized (%s).png', printed.BF.name, bin_n, Pars$RNA.Seq.norm.method))
  
          Render.Density.Plots__color.by.sample(adjusted.cpm.current.bin, out.dir = current.out.dir,
            plot.title = sprintf('[after %s adj.] Bin %d (%s)\nLog2(CPM) density, normalized (%s)', printed.BF.name, bin_n, insertion, Pars$RNA.Seq.norm.method), plot.subtitle = '',
            png.file.name = sprintf('[after %s adj.] Bin %2.f - CPM density, normalized (%s).png', printed.BF.name, bin_n, Pars$RNA.Seq.norm.method))
        }
      }


      bins_count = (length(bins) - 1)
      names(stdev.per.bin) = 1:bins_count
      stdev.per.bin_reordered = stdev.per.bin[rev(order(stdev.per.bin))]
      
      min.cpms.stdev..complex = max(min.cpms.stdev_abs_value.in.bin.to.exclude, mean(stdev.per.bin_reordered[(round(0.25 * bins_count) + 1) : bins_count]) * min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude)
      cat(sprintf('Setting LogCPM StDev threshold (per bin) as %g \n', min.cpms.stdev..complex))

      keep.bins = (stdev.per.bin_reordered <= min.cpms.stdev..complex)
      names(keep.bins) = names(stdev.per.bin_reordered)

      max.bins.to.exclude.real = min(max.bins.to.exclude.num, max.bins.to.exclude.percent / 100 * (length(bins) - 1))
      if(sum.mod(!keep.bins) > max.bins.to.exclude.real){
        for(x in (max.bins.to.exclude.real + 1) : length(keep.bins))  keep.bins[x] = TRUE
      }

      cat(sprintf('%d bins will be excluded:\n', sum.mod(!keep.bins)))

      adjusted_pseudo.counts.table = data.frame(array(dim = c(0, ncol(counts.table))))
      colnames(adjusted_pseudo.counts.table) = colnames(counts.table)
      for(bin_n in 1:(length(bins) - 1)){
        if(keep.bins[as.character(bin_n)]){
          adjusted_pseudo.counts.table = rbind(adjusted_pseudo.counts.table, adjusted_pseudo.counts.table.per.bin[[bin_n]])
        } else {
          if(bf == 'Avg. Transcript Length'){
            insertion = sprintf('transcript length from %.0f to %.0f bp', 2**bins[bin_n], 2**bins[bin_n + 1])
          } else {
            insertion = sprintf('%s from %.0f to %.0f', bf, bins[bin_n], bins[bin_n + 1])
          }
          genes.in.bin = nrow(adjusted_pseudo.counts.table.per.bin[[bin_n]])
          cat(sprintf('    (excluded) bin %d, %s: %d genes\n', bin_n, insertion, genes.in.bin))
        }
      }
      
      adjusted.cpm.all.genes = adjusted_pseudo.counts.table * 1e+6 / mean(colSums(adjusted_pseudo.counts.table))
      if(render.adjustment.CPM.density.plots.all.genes)
        Render.Density.Plots__color.by.sample(adjusted.cpm.all.genes, out.dir = current.out.dir,
          plot.title = sprintf('[after %s adj.] All genes, CPM density, normalized (%s)', printed.BF.name, Pars$RNA.Seq.norm.method), plot.subtitle = '',
          png.file.name = sprintf('[after %s adj.] All genes, Log2(CPM) density, normalized (%s).png', printed.BF.name, Pars$RNA.Seq.norm.method))

      write.table.mod(adjusted_pseudo.counts.table, sprintf("%s/%s-adjusted pseudo read counts, normalized.tsv", current.out.dir, printed.BF.name), quote = FALSE, sep = '\t')
      adjusted_pseudo.counts.table = round(adjusted_pseudo.counts.table)
      cs = colSums(adjusted_pseudo.counts.table)
      dummy.row = t(as.data.frame(max(cs) - cs))
      rownames(dummy.row) = 'dummy'
      adjusted_pseudo.counts.table = rbind(adjusted_pseudo.counts.table, dummy.row)

      tsv.file.name.rc = sprintf("%s/read counts per bin.tsv", current.out.dir)
      write.table.mod(t(read.counts.per.bin.per.sample), tsv.file.name.rc, quote = FALSE, sep = '\t')
      tsv.file.name.nf = sprintf("%s/norm factors per bin.tsv", current.out.dir)
      write.table.mod(t(norm.factors.per.bin.per.sample), tsv.file.name.nf, quote = FALSE, sep = '\t')
      tsv.file.name.cpm = sprintf("%s/cpm per bin.tsv", current.out.dir)
      write.table.mod(t(norm.factors.per.bin.per.sample) / colSums(norm.factors.per.bin.per.sample) * 1e+6, tsv.file.name.cpm, quote = FALSE, sep = '\t')
      tsv.file.name.stdev = sprintf("%s/stdev per bin.tsv", current.out.dir)
      write.table.mod(as.data.frame(stdev.per.bin), tsv.file.name.stdev, quote = FALSE, sep = '\t')
      
      Excel.file.name = Verify.path(sprintf('%s/bins stat.xlsx', current.out.dir))
      CL = sprintf('%s "%s/Any.2.Excel.py" --in "%s" "%s" "%s" "%s" --out-excel "%s" --one-book yes --sheet-names "raw read counts" "raw CPMs" "norm factors" "LogCPM density StDev" --predictor-rows-count 0 --first-col-italic no --entry-type "sample / bin" --forced-col-types-by-names "index:%s" --first-col-width 19',
                   Pars$python.bin, Pars$suppl.data.dir, tsv.file.name.rc, tsv.file.name.cpm, tsv.file.name.nf, tsv.file.name.stdev, Excel.file.name,
                   paste(sprintf('&%d&',1:length(bins)), collapse = ','))
      code = system(CL)
      if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
      if(file.exists(Excel.file.name)){
        tmp = file.remove(tsv.file.name.rc)
        tmp = file.remove(tsv.file.name.nf)
        tmp = file.remove(tsv.file.name.cpm)
        tmp = file.remove(tsv.file.name.stdev)
      } else {
        stop(sprintf('Creating Excel file %s FAILED', Excel.file.name))
      }
      
      
    }

  }
  
  Excel.file.name = Verify.path(sprintf('%s/%sAssociations with bias factors, per sample.xlsx', out.dir, results.prefix))
  CL = sprintf('%s "%s/Any.2.Excel.py" --in %s --out-excel "%s" --one-book yes --sheet-names %s --simple-logfc-format no --predictor-rows-count 0 --forced-col-types-by-names "LogFC:&A-asymmetry&,&G-asymmetry&" "p:&Spearman p-value&" --logfc-abs-max 2 --entry-type "Sample name" --first-col-width 21 --first-col-italic no',
               Pars$python.bin, Pars$suppl.data.dir, paste(sprintf('"%s"', all.tsv.file.names), collapse = ' '),
               Excel.file.name, paste(sprintf('"%s"', bias.factors), collapse = ' '))
  code = system(CL)
  if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
  if(file.exists(Excel.file.name)){
    tmp = file.remove(all.tsv.file.names)
  } else {
    stop(sprintf('Creating Excel file %s FAILED', Excel.file.name))
  }

  # if(include.GLM.model.name.in.Result.names){
  #   saveRDS(Bias.analysis.results, file = sprintf("%s/%s, single-sample.bias.analysis.results%s.rds", out.dir, Analysis.name, db.suffix))
  # } else {
  #   saveRDS(Bias.analysis.results, file = sprintf("%s/single-sample.bias.analysis.results%s.rds", out.dir, db.suffix))
  # }
  return(invisible(adjusted_pseudo.counts.table))
}

Associations.LogFC.with.Bias.Factor__before.GLM = function(Startup.Data, counts.table, Main.Component, Analysis.name, gene.bias.factors.file = '{auto}',
                                       bias.factors = c('Avg. Transcript Length', 'Expression level'), out.dir = NULL, LogCPM.limits = c(0, 3, 5, 7), main.LogCPM.limit = 3,
                                       adjustable.bias.factor = 'Avg. Transcript Length', additional.plots = TRUE, script.dir = NULL,
                                       db.suffix = '', plots.suffix = '', include.GLM.model.name.in.Result.names = FALSE,
                                       LogFC.trimming = 0.02, const.add.read.counts = 0.5, trim = 0){
  

  if(length(Main.Component[!duplicated(Main.Component)]) != 2){
    stop('Main.Component should have only two levels in order to use Associations.LogFC.with.Bias.Factor__before.GLM')
  }

  if(length(Main.Component) != dim(counts.table)[2]){
    stop('[Associations.LogFC.with.Bias.Factor__before.GLM]  Incompartable counts.table and Main.Component')
  }

  if(LogFC.trimming >= 0.5){
    stop('LogFC.trimming cannot exceed 0.5')
  }

  if(!(main.LogCPM.limit %in% LogCPM.limits)){
    stop('main.LogCPM.limit should be present in LogCPM.limits')
  }

  if(is.null(out.dir)) out.dir = getwd()
  
  allowed.bias.factors = c('Avg. Transcript Length', 'Avg. Transcript CG content', 'Expression level')
  absent.bias.factors = bias.factors[!(bias.factors %in% allowed.bias.factors)]
  if(length(absent.bias.factors) > 0){
    stop(sprintf('Unknown bias factors %s. Allowed only: %s', toString(absent.bias.factors), toString(allowed.bias.factors)))
  }

  if(!is.null(adjustable.bias.factor)){
    if(!adjustable.bias.factor %in% bias.factors){
      stop('adjustable.bias.factor should be present in bias.factors')
    }
  } else adjustable.bias.factor = 'none'

  if(gene.bias.factors.file == '{auto}'){
    if(is.null(script.dir)){
      tryCatch(expr = {
        script.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
      }, error = function (err){
        script.dir = '.'
        warning('rstudioapi failed to locate script directory')
      })
    }

    gene.bias.factors.file = sprintf('%s/Gene.Bias.Factors.tsv', script.dir)
    if(!file.exists(gene.bias.factors.file)){
      stop('Cannot find file with gene bias factors information (e.g. transcript length or CG content). Bias correction cannot not be applied. To enable bias correction place Gene.Bias.Factors.tsv (generated with PPLine) in the script directory')
    }
  }
  
  gene.bias.factors.table = read.table(gene.bias.factors.file, header = TRUE, check.names = FALSE, sep = '\t', quote = "")
  rownames(gene.bias.factors.table) = gene.bias.factors.table[,1]
  gene.bias.factors.table = gene.bias.factors.table[,-1]
  
  LogCPM.bf.present = 'Expression level' %in% bias.factors
  bias.factors = bias.factors[!(bias.factors %in% 'Expression level')]
    
  absent.bias.factors = bias.factors[!(bias.factors %in% colnames(gene.bias.factors.table))]
  if(length(absent.bias.factors) > 0){
    print(sprintf('%d bias factors are absent in %s file: %s', length(absent.bias.factors), gene.bias.factors.file, toString(absent.bias.factors)))
    bias.factors = bias.factors[(bias.factors %in% colnames(gene.bias.factors.table))]
  }
  
  if(LogCPM.bf.present) bias.factors = c(bias.factors, 'Expression level')
  
  cat(sprintf('[Associations.LogFC.with.Bias.Factor__before.GLM] Analyzing ~ %s...\n', Analysis.name))
  
  counts.table = counts.table[!(rownames(counts.table) %in% c('dummy')),]
  counts.table = counts.table[apply(counts.table, 1, sum) > 0, ]
  common.genes = intersect(rownames(gene.bias.factors.table), rownames(counts.table))

  if(length(common.genes) < 4){
    cat(sprintf('[Associations.LogFC.with.Bias.Factor__before.GLM]  There is to few genes\n'))
    return(invisible(counts.table))
  }

  counts.table = counts.table[common.genes,]
  gene.bias.factors.table = gene.bias.factors.table[common.genes,]

  adj.counts.table = counts.table
 
  levels = Main.Component[!duplicated(Main.Component)]
  level0 = levels[order(levels)][1]
  level1 = levels[order(levels)][2]
  group0.mask = Main.Component %in% level0
  group1.mask = Main.Component %in% level1

  LogFC.values = apply(counts.table + const.add.read.counts, 1, function(x) log2(mean(x[group1.mask], trim = trim) / mean(x[group0.mask], trim = trim)))

  LogFC.values = pmax(pmin(LogFC.values, 10), -10)
  names(LogFC.values) = rownames(counts.table)
  LogCPM.values = log2(apply(counts.table, 1, sum) / sum(counts.table) * 1e+6 )

  Bias.analysis.results = vector(mode = 'list')
  bf = bias.factors[1]
  for(bf in bias.factors){
    Bias.analysis.results[[bf]] = vector(mode = 'list')
    if(bf == 'Expression level'){
      c__LogCPM.limits = c(0)
    } else {
      c__LogCPM.limits = LogCPM.limits
    }
    current.LogCPM.limit = c__LogCPM.limits[1]
    for(current.LogCPM.limit in c__LogCPM.limits){
      if(bf == 'Expression level'){
        current.LogCPM.limit__char = '-'
        BF.values = LogCPM.values
        current.LogFC.values = LogFC.values
      } else {
        current.LogCPM.limit__char = as.character(current.LogCPM.limit)
        current.gene.bias.factors.table = gene.bias.factors.table[!is.na(gene.bias.factors.table[,bf]), ]

        current.common.genes = intersect(rownames(current.gene.bias.factors.table), names(LogCPM.values[LogCPM.values > current.LogCPM.limit]))
        if(length(current.common.genes) < 4){
          cat(sprintf('[Associations.LogFC.with.Bias.Factor__before.GLM]  There is to few common genes with LogCPM > %g\n', current.LogCPM.limit))
          next
        }

        current.LogFC.values = LogFC.values[current.common.genes]
        current.gene.bias.factors.table = current.gene.bias.factors.table[current.common.genes, ]
        BF.values = current.gene.bias.factors.table[,bf]
        names(BF.values) = rownames(current.gene.bias.factors.table)
      }

      ct.spearman = suppressWarnings(cor.test(BF.values, current.LogFC.values, method = 'spearman'))
      ct.pearson = suppressWarnings(cor.test(BF.values, current.LogFC.values, method = 'pearson'))
      
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]] = vector(mode = 'list')
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['spearman.r']] = ct.spearman$estimate
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['spearman.p']] = ct.spearman$p.value
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['pearson.r']] = ct.pearson$estimate
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['pearson.p']] = ct.pearson$p.value
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['BF.values']] = BF.values
      Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['LogFC.values']] = LogFC.values
      # Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['P.values']] = P.values
      
      sorted.LogFC.values = current.LogFC.values[order(current.LogFC.values)]
      min.LogFC.for.trimming = sorted.LogFC.values[min(length(current.LogFC.values), 1 + round(LogFC.trimming * length(current.LogFC.values)))]
      max.LogFC.for.trimming = sorted.LogFC.values[max(0, length(current.LogFC.values) - round(LogFC.trimming * length(current.LogFC.values)))]
      keep = current.LogFC.values >= min.LogFC.for.trimming & current.LogFC.values <= max.LogFC.for.trimming
      trimmed.LogFC.values = current.LogFC.values[keep]
      trimmed.BF.values = BF.values[keep]


      if (bf == 'Avg. Transcript Length'){
        fit = lm(current.LogFC.values ~ log2(BF.values))
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.interceipt']] = fit$coefficients[1]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.slope']] = fit$coefficients[2]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.pvalue']] = lmp(fit)

        trimmed.fit = lm(trimmed.LogFC.values ~ log2(trimmed.BF.values))
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.interceipt']] = trimmed.fit$coefficients[1]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.slope']] = trimmed.fit$coefficients[2]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.pvalue']] = lmp(trimmed.fit)

        #adjusted.LogFC.values = LogFC.values
        if(bf == adjustable.bias.factor & current.LogCPM.limit == main.LogCPM.limit){
          log2.BF.value.at.Zero.LogFC = - (trimmed.fit$coefficients[1] / trimmed.fit$coefficients[2])
          if(log2.BF.value.at.Zero.LogFC < quantile(log2(BF.values), 0.1) | log2.BF.value.at.Zero.LogFC > quantile(log2(BF.values), 0.8)){
            cat(sprintf('Approximating model "LogFC = %g + %g * Log[Avg. Transcript Length]" is not suitable for LogFC adjusting. This line crosses x axis at position = %g, whereas quartile range of Log2[Avg.Transcript Length] is %g-%g. LogFC adjustment will not be applied\n',
                        trimmed.fit$coefficients[1], trimmed.fit$coefficients[2], log2.BF.value.at.Zero.LogFC,
                        quantile(log2(BF.values), 0.25), quantile(log2(BF.values), 0.75)))
          } else {
            adjustable.genes.with.bias.factor.info = rownames(gene.bias.factors.table)[!is.na(gene.bias.factors.table[,bf])]
            common.genes = intersect(rownames(counts.table), adjustable.genes.with.bias.factor.info)
            adjustable.counts.table = counts.table[common.genes, ]
            adjustable.LogFC.values = LogFC.values[common.genes]
            adjustable.BF.values = gene.bias.factors.table[common.genes, bf]
  
            desired.LogFC.delta =  - log2(adjustable.BF.values) * trimmed.fit$coefficients[2] - trimmed.fit$coefficients[1]
            multipliers = 2**(desired.LogFC.delta / 2)
            multiplier.array.skeleton = rep(1, ncol(adjustable.counts.table))
            x = 1
            chunk = 0
            total.genes = nrow(counts.table)
            for(x in 1:total.genes){
              chunk = chunk + 1
              if(chunk > 50){
                chunk = 0
                cat(sprintf('\radjusting gene %d of %d...   ', x, total.genes))
              }
              current.multiplier.array = multiplier.array.skeleton
              current.multiplier.array[which(group1.mask)] = multipliers[x]
              current.multiplier.array[which(group0.mask)] = 1 / multipliers[x]
              adj.counts.table[x,] = counts.table[x,] * current.multiplier.array
            }
            adj.counts.table = round(adj.counts.table)
            cat(sprintf('\rcompleted                          \n'))
          }  
          #sub.gene.bias.factors.table = gene.bias.factors.table[!is.na(gene.bias.factors.table[,bf]) & (rownames(gene.bias.factors.table) %in% rownames(adjusted.DE.table)), ]
          
          #adjusted.DE.table[rownames(sub.gene.bias.factors.table), 'TBA LogFC'] =
            #adjusted.DE.table[rownames(sub.gene.bias.factors.table), LogFC.col.name] - fit$coefficients[2] * log2(sub.gene.bias.factors.table[, bf]) - fit$coefficients[1]
        }
        
        if(additional.plots){
          if(current.LogCPM.limit == main.LogCPM.limit){
            Bias.Factor.LogFC.Plots(approx_line__const_coeff = fit$coefficients[1], approx_line__slope = fit$coefficients[2], spearman_r = ct.spearman$estimate,
                              BF.values = log2(BF.values), LogFC.values = current.LogFC.values, P.values = NULL, 
                              bf.name = 'Log[Transcript Length]', Analysis.name = Analysis.name, LogFC.col.name = 'trimmed LogFC',
                              out.dir = out.dir, LogCPM.limit = current.LogCPM.limit, logFC.range = 3, p.threshold = 0.05,
                              plots.suffix = plots.suffix,
                              include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                              trimmed.approx_line__const_coeff = trimmed.fit$coefficients[1], trimmed.approx_line__slope = trimmed.fit$coefficients[2])
          }
          current.out.dir = sprintf('%s/Bias Factor plots', out.dir)
          dir.create(current.out.dir, showWarnings = FALSE)
          Bias.Factor.LogFC.Plots(approx_line__const_coeff = fit$coefficients[1], approx_line__slope = fit$coefficients[2], spearman_r = ct.spearman$estimate,
                            BF.values = log2(BF.values), LogFC.values = current.LogFC.values, P.values = NULL, 
                            bf.name = 'Log[Transcript Length]', Analysis.name = Analysis.name, LogFC.col.name = 'trimmed LogFC',
                            out.dir = current.out.dir, LogCPM.limit = current.LogCPM.limit, logFC.range = 3, p.threshold = 0.05,
                            plots.suffix = plots.suffix,
                            include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                            trimmed.approx_line__const_coeff = trimmed.fit$coefficients[1], trimmed.approx_line__slope = trimmed.fit$coefficients[2])
        }

      } else if (bf == 'Expression level'){
        fit = lm(current.LogFC.values ~ BF.values)
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.interceipt']] = fit$coefficients[1]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.slope']] = fit$coefficients[2]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['lm.pvalue']] = lmp(fit)

        trimmed.fit = lm(trimmed.LogFC.values ~ trimmed.BF.values)
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.interceipt']] = trimmed.fit$coefficients[1]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.slope']] = trimmed.fit$coefficients[2]
        Bias.analysis.results[[bf]][[current.LogCPM.limit__char]][['trimmed.lm.pvalue']] = lmp(trimmed.fit)

        if(additional.plots)
          Bias.Factor.LogFC.Plots(approx_line__const_coeff = fit$coefficients[1], approx_line__slope = fit$coefficients[2], spearman_r = ct.spearman$estimate,
                          BF.values = BF.values, LogFC.values = current.LogFC.values, P.values = NULL, 
                          bf.name = 'Expression level', Analysis.name = Analysis.name, LogFC.col.name = 'trimmed LogFC',
                          out.dir = out.dir, LogCPM.limit = current.LogCPM.limit, logFC.range = 3, p.threshold = 0.05, BF.values.axis.limits = c(0, 15),
                          plots.suffix = plots.suffix,
                          include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
                          trimmed.approx_line__const_coeff = trimmed.fit$coefficients[1], trimmed.approx_line__slope = trimmed.fit$coefficients[2])
        
      } else if (bf == 'Avg. Transcript CG content'){
        stop(sprintf('Avg. Transcript CG content is not supported now'))
      } else {
        stop(sprintf('Unknown bias bactor %s', bf))
      }
    }
    #current.gene.bias.factors.table[order(BF.values)[1:10],]
  }
  if(include.GLM.model.name.in.Result.names){
    saveRDS(Bias.analysis.results, file = sprintf("%s/%s, bias.analysis.results%s.rds", out.dir, Analysis.name, db.suffix))
  } else {
    saveRDS(Bias.analysis.results, file = sprintf("%s/bias.analysis.results%s.rds", out.dir, db.suffix))
  }
  return(invisible(adj.counts.table))
  
}


Get.CPM.bins = function(all.CPM.values, bins.count = 8, high.expression.bins.count = 3, std.to.HE.bins.size.ratio = 2.5){
  if(high.expression.bins.count >= bins.count){
    high.expression.bins.count = bins.count - 1
  }
  standart.bins.count = bins.count - high.expression.bins.count

  chunk = 1 / (std.to.HE.bins.size.ratio * standart.bins.count + high.expression.bins.count)

  bins = c(0)
  last.pos = 0

  for(bin_n in 1:standart.bins.count){
    bins = c(bins, last.pos + chunk * std.to.HE.bins.size.ratio)
    last.pos = last.pos + chunk * std.to.HE.bins.size.ratio
  }

  for(bin_n in 1:high.expression.bins.count){
    bins = c(bins, last.pos + chunk)
    last.pos = last.pos + chunk
  }

  cpms = quantile(all.CPM.values, bins)
  cpms[length(cpms)] = cpms[length(cpms)] + 1
  return(cpms)
}

Normalize.read.counts.by.CPM.bins <- function(counts.table, Startup.Data, out.dir = NULL, forced.parameters = NULL,
                                       bins.count = 7, high.expression.bins.count = 2, std.to.HE.bins.size.ratio = 2,
                                       render.adjustment.CPM.density.plots.all.genes = TRUE,
                                       render.adjustment.CPM.density.plots = TRUE,
                                       script.dir = NULL, results.prefix = '', forced.Analysis.name = NULL,
                                       include.GLM.model.name.in.Result.names = FALSE,
                                       min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = 1.6,
                                       min.cpms.stdev_abs_value.in.bin.to.exclude = 0.007,
                                       max.bins.to.exclude.num = 1,
                                       max.bins.to.exclude.percent = 25){


  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters

  if(!is.null(forced.Analysis.name)){
    Analysis.name = forced.Analysis.name
  } else {
    Analysis.name = '(custom)'
  }

  cat(sprintf('[Normalize.read.counts.by.CPM.bins] Analyzing ~ %s...\n', Analysis.name))
  adjusted_pseudo.counts.table = NULL
  CPM.table = cpm(counts.table)

  CPMs = apply(CPM.table, 1, mean)

  bins = Get.CPM.bins(CPMs, bins.count, high.expression.bins.count, std.to.HE.bins.size.ratio)

  if(results.prefix == ''){
    current.out.dir = sprintf('%s/Adjusting for expression level bias', out.dir)
  } else  current.out.dir = sprintf('%s/%s, adjusting for expression level bias', out.dir, results.prefix)
  dir.create(current.out.dir, showWarnings = FALSE)

  d.all.genes = DGEList(counts = counts.table)
  d.all.genes = calcNormFactors(d.all.genes, method=Pars$RNA.Seq.norm.method)
  not.adjusted.cpm.all.genes = cpm(d.all.genes, normalized.lib.sizes = TRUE)
  if(render.adjustment.CPM.density.plots.all.genes)
    Render.Density.Plots__color.by.sample(not.adjusted.cpm.all.genes, out.dir = current.out.dir,
      plot.title = sprintf('[before LogCPM-adj.] All genes, CPM density, normalized (%s)', Pars$RNA.Seq.norm.method), plot.subtitle = '',
      png.file.name = sprintf('[before LogCPM-adj.] All genes, Log2(CPM) density, normalized (%s).png', Pars$RNA.Seq.norm.method))

  read.counts.per.bin.per.sample = data.frame(array( dim = c(length(bins) - 1, ncol(counts.table))))
  rownames(read.counts.per.bin.per.sample) = 1:(length(bins) - 1)
  colnames(read.counts.per.bin.per.sample) = colnames(counts.table)

  norm.factors.per.bin.per.sample = data.frame(array( dim = c(length(bins) - 1, ncol(counts.table))))
  rownames(norm.factors.per.bin.per.sample) = 1:(length(bins) - 1)
  colnames(norm.factors.per.bin.per.sample) = colnames(counts.table)

  stdev.per.bin = c()
  adjusted_pseudo.counts.table.per.bin = vector(mode = 'list')
  sum_counts.table = sum(counts.table)

  
  bin_n = 1
  for(bin_n in 1:(length(bins) - 1)){
    current.counts.table = counts.table[CPMs >= bins[bin_n] & CPMs < bins[bin_n + 1], ]
    genes.in.bin = rownames(current.counts.table)
    sum_current.counts.table = sum(current.counts.table)
    
    if(nrow(current.counts.table) < 15){
      cat(sprintf('bin %d: to few genes in bin... Bypassing...\n', bin_n))
      adjusted_pseudo.counts.table.per.bin[[bin_n]] = counts.table[c(), ]
      stdev.per.bin = c(stdev.per.bin, 0)
      read.counts.per.bin.per.sample[bin_n, ] = 0
      norm.factors.per.bin.per.sample[bin_n, ] = 1
      next
    }
    
    d.bin = DGEList(counts = current.counts.table)
    d.bin = calcNormFactors(d.bin, method=Pars$RNA.Seq.norm.method)
    read.counts.per.bin.per.sample[bin_n, ] = d.bin$samples$lib.size
    norm.factors.per.bin.per.sample[bin_n, ] = d.bin$samples$norm.factors
    
    current.adjusted_pseudo.counts.table = t(t(current.counts.table) / d.bin$samples$lib.size * mean(d.bin$samples$lib.size) / d.bin$samples$norm.factors )
    adjusted_pseudo.counts.table.per.bin[[bin_n]] = current.adjusted_pseudo.counts.table
    # adjusted_pseudo.counts.table = rbind(adjusted_pseudo.counts.table, current.adjusted_pseudo.counts.table)

    adjusted.cpm.current.bin = cpm(d.bin, normalized.lib.sizes = TRUE) * sum_current.counts.table / sum_counts.table
    current.stdev = Get.Density.table.avg.StDev(adjusted.cpm.current.bin)
    stdev.per.bin = c(stdev.per.bin, current.stdev)

    insertion = sprintf('CPM from %.0f to %.0f', bins[bin_n], bins[bin_n + 1])

    cat(sprintf('bin %d, %s: %d genes, LogCPM StDev = %g\n', bin_n, insertion, length(genes.in.bin), current.stdev))

    if(render.adjustment.CPM.density.plots){
      Render.Density.Plots__color.by.sample(not.adjusted.cpm.all.genes[genes.in.bin, ], out.dir = current.out.dir,
        plot.title = sprintf('[before LogCPM-adj.] Bin %d (%s)\nLog2(CPM) density, normalized (%s)', bin_n, insertion, Pars$RNA.Seq.norm.method), plot.subtitle = '',
        png.file.name = sprintf('[before LogCPM-adj.] Bin %2.f - CPM density, normalized (%s).png', bin_n, Pars$RNA.Seq.norm.method))

      Render.Density.Plots__color.by.sample(adjusted.cpm.current.bin, out.dir = current.out.dir,
        plot.title = sprintf('[after LogCPM-adj.] Bin %d (%s)\nLog2(CPM) density, normalized (%s)', bin_n, insertion, Pars$RNA.Seq.norm.method), plot.subtitle = '',
        png.file.name = sprintf('[after LogCPM-adj.] Bin %2.f - CPM density, normalized (%s).png', bin_n, Pars$RNA.Seq.norm.method))
    }
  }


  bins_count = (length(bins) - 1)
  names(stdev.per.bin) = 1:bins_count
  stdev.per.bin_reordered = stdev.per.bin[rev(order(stdev.per.bin))]
  
  min.cpms.stdev..complex = max(min.cpms.stdev_abs_value.in.bin.to.exclude, mean(stdev.per.bin_reordered[(round(0.25 * bins_count) + 1) : bins_count]) * min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude)
  cat(sprintf('Setting LogCPM StDev threshold (per bin) as %g \n', min.cpms.stdev..complex))

  keep.bins = (stdev.per.bin_reordered <= min.cpms.stdev..complex)
  names(keep.bins) = names(stdev.per.bin_reordered)

  max.bins.to.exclude.real = min(max.bins.to.exclude.num, max.bins.to.exclude.percent / 100 * (length(bins) - 1))
  if(sum.mod(!keep.bins) > max.bins.to.exclude.real){
    for(x in (max.bins.to.exclude.real + 1) : length(keep.bins))  keep.bins[x] = TRUE
  }

  cat(sprintf('%d bins will be excluded:\n', sum.mod(!keep.bins)))

  adjusted_pseudo.counts.table = data.frame(array(dim = c(0, ncol(counts.table))))
  colnames(adjusted_pseudo.counts.table) = colnames(counts.table)
  for(bin_n in 1:(length(bins) - 1)){
    if(keep.bins[as.character(bin_n)]){
      adjusted_pseudo.counts.table = rbind(adjusted_pseudo.counts.table, adjusted_pseudo.counts.table.per.bin[[bin_n]])
    } else {
      insertion = sprintf('CPM from %.0f to %.0f', bins[bin_n], bins[bin_n + 1])
      genes.in.bin = nrow(adjusted_pseudo.counts.table.per.bin[[bin_n]])
      cat(sprintf('    (excluded) bin %d, %s: %d genes\n', bin_n, insertion, genes.in.bin))
    }
  }
  
  adjusted.cpm.all.genes = adjusted_pseudo.counts.table * 1e+6 / mean(colSums(adjusted_pseudo.counts.table))
  if(render.adjustment.CPM.density.plots.all.genes)
    Render.Density.Plots__color.by.sample(adjusted.cpm.all.genes, out.dir = current.out.dir,
      plot.title = sprintf('[after LogCPM-adj.] All genes, CPM density, normalized (%s)', Pars$RNA.Seq.norm.method), plot.subtitle = '',
      png.file.name = sprintf('[after LogCPM-adj.] All genes, Log2(CPM) density, normalized (%s).png', Pars$RNA.Seq.norm.method))

  write.table.mod(adjusted_pseudo.counts.table, sprintf("%s/LogCPM-adjusted pseudo read counts, normalized.tsv", current.out.dir), quote = FALSE, sep = '\t')
  adjusted_pseudo.counts.table = round(adjusted_pseudo.counts.table)
  cs = colSums(adjusted_pseudo.counts.table)
  dummy.row = t(as.data.frame(max(cs) - cs))
  rownames(dummy.row) = 'dummy'
  adjusted_pseudo.counts.table = rbind(adjusted_pseudo.counts.table, dummy.row)

  tsv.file.name.rc = sprintf("%s/read counts per bin.tsv", current.out.dir)
  write.table.mod(t(read.counts.per.bin.per.sample), tsv.file.name.rc, quote = FALSE, sep = '\t')
  tsv.file.name.nf = sprintf("%s/norm factors per bin.tsv", current.out.dir)
  write.table.mod(t(norm.factors.per.bin.per.sample), tsv.file.name.nf, quote = FALSE, sep = '\t')
  tsv.file.name.cpm = sprintf("%s/cpm per bin.tsv", current.out.dir)
  write.table.mod(t(norm.factors.per.bin.per.sample) / colSums(norm.factors.per.bin.per.sample) * 1e+6, tsv.file.name.cpm, quote = FALSE, sep = '\t')
  tsv.file.name.stdev = sprintf("%s/stdev per bin.tsv", current.out.dir)
  write.table.mod(as.data.frame(stdev.per.bin), tsv.file.name.stdev, quote = FALSE, sep = '\t')
  
  Excel.file.name = Verify.path(sprintf('%s/bins stat.xlsx', current.out.dir))
  CL = sprintf('%s "%s/Any.2.Excel.py" --in "%s" "%s" "%s" "%s" --out-excel "%s" --one-book yes --sheet-names "raw read counts" "raw CPMs" "norm factors" "LogCPM density StDev" --predictor-rows-count 0 --first-col-italic no --entry-type "sample / bin" --forced-col-types-by-names "index:%s" --first-col-width 19',
               Pars$python.bin, Pars$suppl.data.dir, tsv.file.name.rc, tsv.file.name.cpm, tsv.file.name.nf, tsv.file.name.stdev, Excel.file.name,
               paste(sprintf('&%d&',1:length(bins)), collapse = ','))
  code = system(CL)
  if(code > 0)   cat(sprintf('\nCommand FAILED with error code %d:\n%s\n', code, CL))
  if(file.exists(Excel.file.name)){
    tmp = file.remove(tsv.file.name.rc)
    tmp = file.remove(tsv.file.name.nf)
    tmp = file.remove(tsv.file.name.cpm)
    tmp = file.remove(tsv.file.name.stdev)
  } else {
    stop(sprintf('Creating Excel file %s FAILED', Excel.file.name))
  }
  

  # if(include.GLM.model.name.in.Result.names){
  #   saveRDS(Bias.analysis.results, file = sprintf("%s/%s, single-sample.bias.analysis.results%s.rds", out.dir, Analysis.name, db.suffix))
  # } else {
  #   saveRDS(Bias.analysis.results, file = sprintf("%s/single-sample.bias.analysis.results%s.rds", out.dir, db.suffix))
  # }
  return(invisible(adjusted_pseudo.counts.table))
}

        
Filter.counts.by.transcript.length = function(counts.table, gene.bias.factors.file = '{auto}', min.length = NULL, max.length = NULL, script.dir = NULL){
  if(gene.bias.factors.file == '{auto}'){
    if(is.null(script.dir)){
      tryCatch(expr = {
        script.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
      }, error = function (err){
        script.dir = '.'
        warning('rstudioapi failed to locate script directory')
      })
    }

    gene.bias.factors.file = sprintf('%s/Gene.Bias.Factors.tsv', script.dir)
    if(!file.exists(gene.bias.factors.file)){
      stop('Cannot find file with gene bias factors information (e.g. transcript length or CG content). Filter.counts.by.transcript.length cannot run. Place Gene.Bias.Factors.tsv (generated with PPLine) in the script directory')
    }
  }

  gene.bias.factors.table = read.table(gene.bias.factors.file, header = TRUE, check.names = FALSE, sep = '\t', quote = "")
  rownames(gene.bias.factors.table) = gene.bias.factors.table[,1]
  gene.bias.factors.table = gene.bias.factors.table[,-1]

  if(!('Avg. Transcript Length' %in% colnames(gene.bias.factors.table)))  stop("'Avg. Transcript Length' information is absent in gene.bias.factors.table")

  gene.bias.factors.table = gene.bias.factors.table[!is.na(gene.bias.factors.table[, 'Avg. Transcript Length']), ]

  keep = rep(TRUE, dim(gene.bias.factors.table)[1])
  if(!is.null(min.length)){
    current.keep = gene.bias.factors.table[, 'Avg. Transcript Length'] >= min.length
    keep = keep & current.keep
  }

  if(!is.null(max.length)){
    current.keep = gene.bias.factors.table[, 'Avg. Transcript Length'] <= max.length
    keep = keep & current.keep
  }

  return(counts.table[rownames(counts.table) %in% rownames(gene.bias.factors.table)[keep], ])
}


get.ranks = function(numeric.array){
  numeric.array.nna = numeric.array[!is.na(numeric.array)]
  return(sapply(numeric.array, function(x) {  sum.mod(x >= numeric.array.nna) / length(numeric.array.nna) }))
}

Add.Quantile.statistics__LogCPM = function(Stats.table = NULL, Quantile.Statistics.table = NULL){
  
  if(is.null(Quantile.Statistics.table)){
    if(is.null(Stats.table))   stop('Either Stats.table or Quantile.Statistics.table should be supplied')
    Quantile.Statistics.table = data.frame(array(dim = c(dim(Stats.table)[1],0)))
    rownames(Quantile.Statistics.table) = rownames(Stats.table)
  } else if(!is.null(Stats.table)){
    if(sum.mod(rownames(Stats.table) != rownames(Quantile.Statistics.table)))  stop('row names of Quantile.Statistics.table and Stats.table do not match')
  }

  LogCPM.col.name = Get.LogCPM.col.name(colnames(Stats.table))
  rank.column = get.ranks(Stats.table[, LogCPM.col.name]) * 100
  Quantile.Statistics.table = cbind(Quantile.Statistics.table, as.data.frame(rank.column))
  colnames(Quantile.Statistics.table)[dim(Quantile.Statistics.table)[2]] = 'LogCPM rank (%)'
  return(Quantile.Statistics.table)
}

Add.Quantile.statistics__BF = function(Stats.table = NULL, Quantile.Statistics.table = NULL, BF = 'Avg. Transcript Length',
                                       QStats.colname = 'Tr.Length rank (%)', gene.bias.factors.file = '{auto}', script.dir = NULL){
  
  if(is.null(Quantile.Statistics.table)){
    if(is.null(Stats.table))   stop('Either Stats.table or Quantile.Statistics.table should be supplied')
    Quantile.Statistics.table = data.frame(array(dim = c(dim(Stats.table)[1],0)))
    rownames(Quantile.Statistics.table) = rownames(Stats.table)
  } else if(!is.null(Stats.table)){
    if(sum.mod(rownames(Stats.table) != rownames(Quantile.Statistics.table)))  stop('row names of Quantile.Statistics.table and Stats.table do not match')
  }

  if(gene.bias.factors.file == '{auto}'){
    if(is.null(script.dir)){
      tryCatch(expr = {
        script.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
      }, error = function (err){
        script.dir = '.'
        warning('rstudioapi failed to locate script directory')
      })
    }

    gene.bias.factors.file = sprintf('%s/Gene.Bias.Factors.tsv', script.dir)
    if(!file.exists(gene.bias.factors.file)){
      stop('Cannot find file with gene bias factors information (e.g. transcript length or CG content). Filter.counts.by.transcript.length cannot run. Place Gene.Bias.Factors.tsv (generated with PPLine) in the script directory')
    }
  }

  gene.bias.factors.table = read.table(gene.bias.factors.file, header = TRUE, check.names = FALSE, sep = '\t', quote = "")
  rownames(gene.bias.factors.table) = gene.bias.factors.table[,1]
  gene.bias.factors.table = gene.bias.factors.table[,-1]      

  if(!(BF %in% colnames(gene.bias.factors.table)))  stop(sprintf('[Add.Quantile.statistics__BF]   bias factor %s is not found in file %s', BF, gene.bias.factors.file))

  BF.values = gene.bias.factors.table[rownames(Quantile.Statistics.table), BF]
  rank.column = get.ranks(BF.values) * 100
  Quantile.Statistics.table = cbind(Quantile.Statistics.table, as.data.frame(rank.column))
  colnames(Quantile.Statistics.table)[dim(Quantile.Statistics.table)[2]] = QStats.colname
  return(Quantile.Statistics.table)
}



# Summarize.GLM.results

# Create.Joint.DE.heatmaps = function(Startup.Data, GLM.models = NULL, GLM.working.dir = NULL, Analysis.Data.RDS.file = NULL,
#                                             GO.type = 'BP', profile.type = 'custom',
#                                             bypass.if.completed = FALSE, forced.parameters = NULL, forced.terms = NULL,
#                                             out.dir = NULL, forced.Analysis.name = NULL,GO.stats.file.name = NULL, sparklines.axis.limits = 2, min.term.size = NULL,
#                                             include.GLM.model.name.in.Result.names = FALSE){
#   if(sum(c(!is.null(Analysis.Data), !is.null(GLM.model), !is.null(Analysis.Data.RDS.file))) != 1) {
#     stop('Please specify either Analysis.Data, GLM.model or Analysis.Data.RDS.file arguments')
#   }
  
#   if(typeof(Analysis.Data) == 'character'){
#     GLM.model = Analysis.Data
#     Analysis.Data = NULL
#   }
  
#   if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
#   } else Pars = forced.parameters
  
#   if(is.null(Analysis.Data)){
#     Analysis.name  = Cleanup.Model.Name(GLM.model)
#     if(is.null(GLM.working.dir)) GLM.working.dir = sprintf("%s/~ %s, results", Pars$results.dir, Analysis.name)
#     if(is.null(Analysis.Data.RDS.file)){
#       cf1 = sprintf('%s/analysis data.rds', GLM.working.dir)
#       cf2 = sprintf('%s/%s, analysis data.rds', GLM.working.dir, Analysis.name)
#       if(file.exists(cf1)){
#         Analysis.Data.RDS.file = cf1
#       } else if(file.exists(cf2)){
#         Analysis.Data.RDS.file = cf2
#       } else {
#         stop(sprintf('Cannot find any *.rds RDS files: "%s" or "%s"', cf1, cf2))
#       }
#     }
#     Analysis.Data = readRDS(file = Analysis.Data.RDS.file)
#     Analysis.name = Analysis.Data$Analysis.name
#     Analysis.Data$current.GLM.results.dir = GLM.working.dir
#   } else {
#     Analysis.name = Analysis.Data$Analysis.name
#   }
  
#   if(!Analysis.Data$GLM.analysis.performed) stop('GO.Expression.Profiles(...) must be run after GLM/DE analysis. Run Analyze.GLM(...) first.')
  
#   # if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
#   # } else Pars = forced.parameters
  
#   if (is.null(forced.Analysis.name)) { Analysis.name = Analysis.Data$Analysis.name
#   } else Analysis.name = forced.Analysis.name
  
#   if(is.null(out.dir) != is.null(forced.terms))  stop('GO.Expression.Profiles: both forced.terms and out.dir parameters should be both NULL or both have a value')
  
#   if (bypass.if.completed & 
#       Read.Completed.steps.status(sprintf("%s.%s.%s.GO.expression.visualization",Analysis.Data$GLM.model,Startup.Data$Startup.hashmd5,Analysis.Data$step.hashmd5),Startup.Data$Completed.steps.file)){
#     message(sprintf('\nGO-centric expression profiles visualization for GLM model "%s" was previously completed. Bypassing...',Analysis.Data$GLM.model))
#     return()
#   }
  
#   LogFC.col.name = Get.LogFC.col.name(colnames(Analysis.Data$ResTable.gene.part))
#   PValue.col.name = Get.PValue.col.name(Analysis.Data$Main.test, colnames(Analysis.Data$ResTable.gene.part))
# }

Perform.WGCNA = function(Startup.Data, WGCNA.group, forced.parameters = NULL,
                       script.dir = NULL, Create.Excel.results = NULL,
                       include.GLM.model.name.in.Result.names = FALSE,
                       verbose = TRUE){
  if (is.null(forced.parameters)) { Pars = Startup.Data$Pars
  } else Pars = forced.parameters

  cat(sprintf('Processing WGCNA group "%s"...\n', WGCNA.group))

  if(is.null(script.dir))  script.dir = Get.script.path()

  setwd(Pars$results.dir)

  WGCNA.Data = new('RTrans.WGCNA.AnalysisData')

  WGCNA.Analysis.name = Cleanup.Model.Name(WGCNA.group)  
  WGCNA.Data$WGCNA.Analysis.name = WGCNA.Analysis.name

  current.WGCNA.results.dir = sprintf("%s/WGCNA - %s, results", Pars$results.dir, WGCNA.Analysis.name)
  WGCNA.Data$current.WGCNA.results.dir = current.WGCNA.results.dir
  dir.create(current.WGCNA.results.dir, showWarnings = FALSE)
  setwd(current.WGCNA.results.dir)



  if(Pars$Evaluate.bias.factors.Associations){
    if(Pars$gene.bias.factors.file == ''){
      gene.bias.factors.file = sprintf('%s/Gene.Bias.Factors.tsv', script.dir)
      if(file.exists(gene.bias.factors.file)){
        cat('Gene bias factors file is found.\n')
      } else {
        cat('Cannot find file with gene bias factors information (e.g. transcript length or CG content). Bias correction will not be applied. To enable bias correction place Gene.Bias.Factors.tsv (generated with PPLine) in the script directory\n')
        gene.bias.factors.file = NULL
      }
    } else if(!file.exists(gene.bias.factors.file)){
      stop(sprintf('Gene bias factors file %s is NOT found. Check parameter gene.bias.factors.file', gene.bias.factors.file))
    }
  } else {
    gene.bias.factors.file = NULL
  }

  WGCNA.schema = Startup.Data$WGCNA.schema
  if(!(WGCNA.group %in% colnames(WGCNA.schema)))
    stop(sprintf('WGCNA group "%s" is not found. The following WGCNA groups are available: %s', WGCNA.group, toString(colnames(WGCNA.schema)[-1])))

  if(Pars$DE.package != 'edger')  stop('WGCNA works only with edgeR package')
  
  #included.samples = rownames(WGCNA.schema)[!is.na(WGCNA.schema[,WGCNA.group])]

  WGCNA.schema = WGCNA.schema[!is.na(WGCNA.schema[,WGCNA.group]),, drop = FALSE]

  cat(sprintf('WGCNA group "%s" (%d samples):  loading and processing RNA-Seq data...\n', WGCNA.group, nrow(WGCNA.schema)))

  if(Pars$counts.dir != ''){
    if ("File names" %in% colnames(WGCNA.schema)){
      WGCNA.file.names = sprintf("%s/%s", Pars$counts.dir, WGCNA.schema$'File names')
    } else {
      WGCNA.file.names = sprintf("%s/%s%s", Pars$counts.dir, WGCNA.schema$'Sample names', Pars$counts.suffix)
    }

    if (sum.mod(!file.exists(WGCNA.file.names)) > 0){
      #warning(sprintf('Warning; The following files do not exist: %s',toString(file.names[!file.exists(file.names)])))
      cat(sprintf('\nTotal %d of %d *.counts files do not exist in the directory "%s": %s',
                  sum.mod(!file.exists(WGCNA.file.names)), length(WGCNA.file.names), Pars$counts.dir, toString(
                    sapply(strsplit(WGCNA.file.names[!file.exists(WGCNA.file.names)], split='/'), function(x) { tail(x,1)})
                    )))
      if (sum.mod(file.exists(WGCNA.file.names)) == 0) stop(sprintf('No *%s files found. Please fix "WGCNA groups" worksheet and check the directory "%s" ', Pars$counts.suffix, Pars$counts.dir))

      WGCNA.schema = WGCNA.schema[file.exists(WGCNA.file.names),]
      rownames(WGCNA.schema) = WGCNA.schema$'Sample names'
      WGCNA.file.names = WGCNA.file.names[file.exists(WGCNA.file.names)]
    }

  } else if (Pars$read.counts.table != ''){
    sep = Choose.Separator(Pars$read.counts.table)
    rt = read.table(Pars$read.counts.table, stringsAsFactors = FALSE, sep = sep, header = TRUE, check.names = FALSE)
    needed.samples = colnames(rt)[colnames(rt) %in% WGCNA.schema$'Sample names']
    rt = rt[,needed.samples]
    # print(rownames(rt))
    # stop('')
    absent.samples = WGCNA.schema$'Sample names'[!(WGCNA.schema$'Sample names' %in% colnames(rt))]
    if(length(absent.samples) > 0){
      cat(sprintf('The following samples are absent in file %s (total %d of %d): %s\n', Pars$read.counts.table, length(absent.samples),
        length(WGCNA.schema$'Sample names'), toString(absent.samples)))

      if(length(absent.samples) == length(WGCNA.schema$'Sample names'))  stop(sprintf('All samples are ABSENT in the table %s', Pars$read.counts.table))
    }

    sample.is.present = WGCNA.schema$'Sample names' %in% colnames(rt)
    WGCNA.schema = WGCNA.schema[sample.is.present, ]
    rownames(WGCNA.schema) = WGCNA.schema$'Sample names'
    rt = rt[, WGCNA.schema$'Sample names']

  } else stop('Both Pars$read.counts.table and Pars$counts.dir are not specified')

  Main.Component = WGCNA.schema[, WGCNA.group, drop=FALSE]
  
  if(!(Pars$DE.package %in% c('edger')) & !is.null(gene.bias.factors.file) & (Pars$Adjust.Transcript.length.bias__in.Read.counts | Pars$Adjust.Expression.level.bias__in.Read.counts)){
    stop('Bias factor adjustemnt is available only for edgeR')
  }

  # if(Pars$DE.package == 'edger'){
  if(Pars$counts.dir != ''){
    samples_counts=suppressMessages(readDGE(WGCNA.file.names))
    counts = samples_counts$counts
    # naming 'counts'
    colnames(counts) = WGCNA.schema$'Sample names'

  } else if (Pars$read.counts.table != ''){
    sep = Choose.Separator(Pars$read.counts.table)
    rt = read.table(Pars$read.counts.table, stringsAsFactors = FALSE, sep = sep, header = TRUE, check.names = FALSE, comment.char = '#')
    # if(all(!is.na(as.numeric(as.character(rownames(rt))))))
    #   stop(sprintf('Incorrect rownames in file "%s". Check the length of header (must be equal to the number of samples - 1)', Pars$read.counts.table))
    rt = rt[ , WGCNA.schema$'Sample names']
    samples_counts = DGEList(rt)
    counts = samples_counts$counts
    # naming 'counts'
    colnames(counts) = WGCNA.schema$'Sample names'

  } else stop('Both Pars$read.counts.table and Pars$counts.dir are not specified')

  
  # calculating stats and excluding meta-tags
  meta.tags = rownames(counts)[which(startsWith(rownames(counts), '__'))]
  # ("__no_feature","__ambiguous",'__too_low_aQual','__not_aligned','__alignment_not_unique')
  # meta.tags = meta.tags[meta.tags %in% rownames(counts)]
  if (length(meta.tags) > 0){
    mapping.stats.table = array(dim = c(dim(counts)[2],length(meta.tags) + 2))
    rownames(mapping.stats.table) = WGCNA.schema$'Sample names'
    colnames(mapping.stats.table) = c(meta.tags,'total hits','useful reads')
    for (mt in meta.tags){
      mapping.stats.table[WGCNA.schema$'Sample names',mt] = counts[mt,WGCNA.schema$'Sample names']
    }
    mapping.stats.table[WGCNA.schema$'Sample names','total hits'] = apply(counts,2,sum)
    
    noint = rownames(counts) %in% meta.tags
    counts = counts[!noint,]
    mapping.stats.table[WGCNA.schema$'Sample names','useful reads'] = apply(counts,2,sum)
    write.table.mod(mapping.stats.table,'mapping stats.tsv',sep='\t')
  }
  
  # excluding low expression genes
  cpms = cpm(counts)
  if (Pars$min.samples.with.sufficient.CPM.for.WGCNA.samples == '{auto}') {
    min.samples.with.sufficient.CPM.for.WGCNA.samples = Calc.min.high.CPM.samples.from.total.samples.count(dim(WGCNA.schema)[1])
  } else   min.samples.with.sufficient.CPM.for.WGCNA.samples = Pars$min.samples.with.sufficient.CPM.for.WGCNA.samples
  
  
  if(Pars$sufficient.CPM.to.analyze.for.WGCNA.samples %in% c('{auto}', 'auto')){
    min.CPM.single.sample = max(1, 3e+7 / mean(colSums(counts)) / 1.25 ) * Pars$sufficient.CPM..autoadjust..multiplier.for.WGCNA.samples
  } else   min.CPM.single.sample = Pars$sufficient.CPM.to.analyze.for.WGCNA.samples
  
  cat(sprintf('\nSetting CPM limit %.1f for at least %d samples', min.CPM.single.sample, min.samples.with.sufficient.CPM.for.WGCNA.samples))
  
  keep = rowSums(cpms > min.CPM.single.sample) >= min.samples.with.sufficient.CPM.for.WGCNA.samples
  cat(sprintf('\n%d of %d genes passed CPM threshold\n',sum.mod(keep),length(keep)))
  counts = counts[keep,]


    
  .pardefault <- par(no.readonly = TRUE)
  
  plot.subtitle = ''
  # if (length(levels(factor(Main.Component))) == 2) plot.subtitle = sprintf('color indicates %s (quantity; blue-orange)',Main.Predictor.name)
  # if (length(levels(factor(Main.Component))) == 3) plot.subtitle = sprintf('color indicates %s (quantity; blue-green-orange)',Main.Predictor.name)
  # if (length(levels(factor(Main.Component))) >= 4) plot.subtitle = sprintf('color indicates %s (quantity; blue-green-orange gradient)',Main.Predictor.name)
  
  # col.variety <- colorRampPalette(c("#0b00c4","#02babc","#25ad00","#cfcd00","#ff8400"))(512)
  # Main.Component.mod = Main.Component - min(Main.Component)
  # cl = col.variety[1 + as.integer(Main.Component.mod/max(Main.Component.mod)*511)]

  if(Pars$minimal.Transcript.length > 0){
    if(is.null(gene.bias.factors.file))  stop('parameter minimal.Transcript.length is set > 0, but gene.bias.factors.file is NULL')
    counts = Filter.counts.by.transcript.length(counts.table = counts, gene.bias.factors.file = gene.bias.factors.file,
                  min.length = Pars$minimal.Transcript.length, max.length = NULL)
  }
  
  adj.counts.table = NULL
  # if(!exclude.bins__in.Read.counts.adjustemnt___BF)  max.bins.to.exclude.num___BF = 0
  results.prefix = ''
  if(!is.null(gene.bias.factors.file) & length(Pars$Bias.factors.to.analyze) > 0){
    if(!Pars$Adjust.Transcript.length.bias__in.Read.counts){
      tmp = analyze.LogCPM.each.Sample.vs.Bias.Factor.associations(counts.table = counts, out.dir = current.WGCNA.results.dir,
        forced.parameters = Pars, gene.bias.factors.file = gene.bias.factors.file, bias.factors = Pars$Bias.factors.to.analyze,
        adjustable.bias.factor = NULL, min.CPM.for.association.analysis = 1.0,
        render.not.adjusted.accos.plots = Pars$Draw.bias.associations.plots..for.WGCNA.samples,
        render.adjustment.CPM.density.plots = Pars$Draw.CPM..to..bias.factor.density.plots..for.WGCNA.samples,
        render.adjustment.CPM.density.plots.all.genes = Pars$Draw.CPM..to..bias.factor.density.plots..for.WGCNA.samples,
        script.dir = script.dir, results.prefix = '',
        forced.Analysis.name = WGCNA.Analysis.name, include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
    
    } else {
      adj.counts.table = counts

      if(Pars$Adjust.Transcript.length.bias__in.Read.counts){
        aBF = 'Avg. Transcript Length'
        if(Pars$Adjust.Transcript.length.bias__in.Read.counts__exclude.bins){
          max.bins.to.exclude.num =                              Pars$Adjust.Transcript.length.bias__in.Read.counts__max.bins.to.exclude__count
          max.bins.to.exclude.percent =                          Pars$Adjust.Transcript.length.bias__in.Read.counts__max.bins.to.exclude__percentage
          min.cpms.stdev_abs_value.in.bin.to.exclude =           Pars$Adjust.Transcript.length.bias__in.Read.counts__StDev.abs.value__in.bin.to.exclude
          min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = Pars$Adjust.Transcript.length.bias__in.Read.counts__StDev.ratio__in.bin.to.exclude
        } else {
          max.bins.to.exclude.num = 0
          max.bins.to.exclude.percent = 0
          min.cpms.stdev_abs_value.in.bin.to.exclude  = 100
          min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = 100
        }

        adj.counts.table = adj.counts.table[!(rownames(adj.counts.table) %in% 'dummy'),]
        adj.counts.table = 
          analyze.LogCPM.each.Sample.vs.Bias.Factor.associations(counts.table = adj.counts.table, out.dir = current.WGCNA.results.dir,
            forced.parameters = Pars, gene.bias.factors.file = gene.bias.factors.file, bias.factors = Pars$Bias.factors.to.analyze,
            adjustable.bias.factor = aBF, min.CPM.for.association.analysis = 1.0,
            render.not.adjusted.accos.plots = Pars$Draw.bias.associations.plots..for.WGCNA.samples,
            render.adjustment.CPM.density.plots.all.genes = Pars$Draw.CPM..to..bias.factor.density.plots..for.WGCNA.samples,
            render.adjustment.CPM.density.plots = Pars$Draw.CPM..to..bias.factor.density.plots..for.WGCNA.samples,
            script.dir = script.dir, results.prefix = results.prefix,
            forced.Analysis.name = WGCNA.Analysis.name, include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names,
            min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude,
            min.cpms.stdev_abs_value.in.bin.to.exclude = min.cpms.stdev_abs_value.in.bin.to.exclude,
            max.bins.to.exclude.num = max.bins.to.exclude.num,
            max.bins.to.exclude.percent = max.bins.to.exclude.percent)
        results.prefix = sprintf('Adjusted to %s - %s', aBF, results.prefix)
      }

      tmp = analyze.LogCPM.each.Sample.vs.Bias.Factor.associations(counts.table = adj.counts.table, out.dir = current.WGCNA.results.dir,
        forced.parameters = Pars, gene.bias.factors.file = gene.bias.factors.file, bias.factors = Pars$Bias.factors.to.analyze,
        adjustable.bias.factor = NULL, min.CPM.for.association.analysis = 1.0,
        render.not.adjusted.accos.plots = Pars$Draw.bias.associations.plots..for.WGCNA.samples,
        render.adjustment.CPM.density.plots.all.genes = Pars$Draw.CPM..to..bias.factor.density.plots..for.WGCNA.samples,
        render.adjustment.CPM.density.plots = Pars$Draw.CPM..to..bias.factor.density.plots..for.WGCNA.samples,
        script.dir = script.dir, results.prefix = results.prefix,
        forced.Analysis.name = WGCNA.Analysis.name, include.GLM.model.name.in.Result.names = include.GLM.model.name.in.Result.names)
    }
  }

  if(Pars$Adjust.Expression.level.bias__in.Read.counts){
    if(!is.null(adj.counts.table)){
      src.pseudo.counts.table = adj.counts.table[!(rownames(adj.counts.table) %in% 'dummy'),]
    } else   src.pseudo.counts.table = counts
    
    if(Pars$Adjust.Expression.level.bias__in.Read.counts__exclude.bins){
      max.bins.to.exclude.num =                              Pars$Adjust.Expression.level.bias__in.Read.counts__max.bins.to.exclude__count
      max.bins.to.exclude.percent =                          Pars$Adjust.Expression.level.bias__in.Read.counts__max.bins.to.exclude__percentage
      min.cpms.stdev_abs_value.in.bin.to.exclude =           Pars$Adjust.Expression.level.bias__in.Read.counts__StDev.abs.value__in.bin.to.exclude
      min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = Pars$Adjust.Expression.level.bias__in.Read.counts__StDev.ratio__in.bin.to.exclude
    } else {
      max.bins.to.exclude.num = 0
      max.bins.to.exclude.percent = 0
      min.cpms.stdev_abs_value.in.bin.to.exclude  = 100
      min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = 100
    }

    adj.counts.table = Normalize.read.counts.by.CPM.bins(src.pseudo.counts.table, Startup.Data,
                               out.dir = current.WGCNA.results.dir, forced.parameters = Pars,
                               bins.count = 8, high.expression.bins.count = 3, std.to.HE.bins.size.ratio = 2.5,
                               render.adjustment.CPM.density.plots = Pars$Draw.CPM..to..bias.factor.density.plots..for.WGCNA.samples,
                               render.adjustment.CPM.density.plots.all.genes = Pars$Pars$Draw.CPM..to..bias.factor.density.plots..for.WGCNA.samples,
                               script.dir = script.dir, results.prefix = results.prefix, forced.Analysis.name = WGCNA.Analysis.name,
                               include.GLM.model.name.in.Result.names = FALSE,
                               min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude = min.cpms.stdev_to_mean.stdev.ratio.in.bin.to.exclude,
                               min.cpms.stdev_abs_value.in.bin.to.exclude = min.cpms.stdev_abs_value.in.bin.to.exclude,
                               max.bins.to.exclude.num = max.bins.to.exclude.num,
                               max.bins.to.exclude.percent = max.bins.to.exclude.percent)
  }

  
  if(Pars$Create.CPM.density.plots.for.WGCNA.samples){
    cpms = cpm(counts)
    Write.Density.table(cpms, tsv.file.name = 'CPM density, non-normalized.tsv')
    #Render.Density.Plots__color.by.condition(cpms, cl, out.dir = NULL, plot.title = 'Log2(CPM) density, non-normalized', plot.subtitle = plot.subtitle,
    #  png.file.name = 'CPM density, non-normalized, color by condition.png')
    Render.Density.Plots__color.by.sample(cpms, out.dir = NULL, plot.title = 'Log2(CPM) density, non-normalized', plot.subtitle = '',
      png.file.name = 'CPM density, non-normalized, color by sample.png')
  }    

  # creating DGEList object
  # Calculating Norm Factors
  d = DGEList(counts=counts)
  
  #RNA.Seq.norm.method = 'RLE'  #,"TMM","upperquartile","none","RLE"
  if(Pars$RNA.Seq.norm.method == 'RG'){
    d = Normalize..Using.Reference.Genes(d, Startup.Data$Reference.genes)
  } else {
    d = calcNormFactors(d, method=Pars$RNA.Seq.norm.method) #,"TMM","upperquartile","none","RLE"
  }
  WGCNA.Data$edgeR_d = d
  
  # plotting histograms - after normalization
  if(Pars$Create.CPM.density.plots.for.WGCNA.samples){
    cpms = cpm(d, normalized.lib.sizes = TRUE)
    Write.Density.table(cpms, tsv.file.name = sprintf('CPM density, normalized (%s).tsv', Pars$RNA.Seq.norm.method))
    Render.Density.Plots__color.by.sample(cpms, out.dir = NULL, plot.title = sprintf('Log2(CPM) density, normalized (%s)', Pars$RNA.Seq.norm.method), plot.subtitle = '',
      png.file.name = sprintf('CPM density, normalized (%s), color by sample.png', Pars$RNA.Seq.norm.method))
    #Render.Density.Plots__color.by.condition(cpms, cl, out.dir = NULL, plot.title = sprintf('Log2(CPM) density, normalized (%s)', Pars$RNA.Seq.norm.method), plot.subtitle = plot.subtitle,
    #  png.file.name = sprintf('CPM density, normalized (%s), color by condition.png', Pars$RNA.Seq.norm.method))
  }
  
  
  # plotting histograms - after normalization and ajustement
  if(!is.null(adj.counts.table) & Pars$Create.CPM.density.plots.for.WGCNA.samples){
    adj.counts.table.cl = adj.counts.table[!(rownames(adj.counts.table) %in% 'dummy'),]
    adj.cpms = adj.counts.table.cl * 1e+6 / mean(colSums(adj.counts.table.cl))
    Write.Density.table(adj.cpms, tsv.file.name = sprintf('CPM density, normalized (%s).tsv', Pars$RNA.Seq.norm.method))
    Render.Density.Plots__color.by.sample(adj.cpms, out.dir = NULL, plot.title = sprintf('Log2(CPM) density,\nadjusted and normalized (%s)', Pars$RNA.Seq.norm.method), plot.subtitle = '',
      png.file.name = sprintf('CPM density, adjusted and normalized (%s), color by sample.png', Pars$RNA.Seq.norm.method))
    #Render.Density.Plots__color.by.condition(adj.cpms, cl, out.dir = NULL, plot.title = sprintf('Log2(CPM) density,\nadjusted and normalized (%s)', Pars$RNA.Seq.norm.method), plot.subtitle = plot.subtitle,
    #  png.file.name = sprintf('CPM density, adjusted and normalized (%s), color by condition.png', Pars$RNA.Seq.norm.method))
    # cat(sprintf('\n\nUsing counts table that is adjusted to the following bias factors: %s%s\n\n', toString(bias.factors__to.adjust__read.counts),
    #   ifelse(adjust.to.Expression.Level, ", 'absolute' expression level (LogCPM)", "")))
    rm(adj.counts.table.cl)
  }
  par(.pardefault)
  # }

  # if(Pars$DE.package == 'edger'){
  C.Create.Distance.matrices = Pars$Create.Distance.matrices.for.WGCNA.samples
  if (Pars$Create.Distance.matrices.for.WGCNA.samples == TRUE)
    C.Create.Distance.matrices = (dim(cpms)[2] <= 70)
  
  if (C.Create.Distance.matrices)      Create.Distance.matrices.and.Dendrograms(Pars, cpms, WGCNA.schema, NULL, WGCNA.group, WGCNA.group, prefix = '', verbose = verbose)

  if(Pars$Create.MDS.plots.for.WGCNA.samples)
    Create.MDS.plots(d, Main.Component, WGCNA.Analysis.name, out.dir = NULL, prefix = '')


  if(!is.null(adj.counts.table)){
    d = DGEList(counts=adj.counts.table)
    if(Pars$RNA.Seq.norm.method == 'RG'){
      d = Normalize..Using.Reference.Genes(d, Startup.Data$Reference.genes)
    } else {
      d = calcNormFactors(d, method = 'none') #,"TMM","upperquartile","none","RLE"
    }
    WGCNA.Data$edgeR_d = d
    
    ##  for adjusted read counts
    if(Pars$Create.MDS.plots.for.WGCNA.samples)
      Create.MDS.plots(d, Main.Component, WGCNA.Analysis.name, out.dir = NULL, prefix = 'Bias-adjusted ')

    if (C.Create.Distance.matrices)      Create.Distance.matrices.and.Dendrograms(Pars, adj.cpms, WGCNA.schema, NULL, WGCNA.group, WGCNA.group, prefix = 'Bias adjusted ', verbose = verbose)
  }
  # }

  norm.CPMs = cpm(d, normalized.lib.sizes = TRUE)
  if(Pars$report.norm.read.counts.instead.of.CPM){
    lib.sizes = d$samples[,'lib.size']
    norm.factors = d$samples[,'norm.factors']
    norm.counts = t(t(d$counts) / lib.sizes / norm.factors * mean(lib.sizes))
  }
  n_genes = dim(norm.CPMs)[1]

  # print(dim(norm.CPMs))
  
  ##################################################
  ##### main WGCNA processing module
  ##################################################
  ##################################################

  options(stringsAsFactors = FALSE)
  datExpr = t(as.data.frame(cpm(d, normalized.lib.sizes = TRUE)))
  # removing genes containing too many NAs
  gsg = goodSamplesGenes(datExpr, verbose = 3);
  if (!gsg$allOK) {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  } else {
    printFlush('All genes are OK for WCGNA')
  }
  
  sampleTree = hclust(dist(datExpr), method = "average");
  sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  .pardefault <- par(no.readonly = TRUE)
  png(filename = sprintf('Clustering samples (average method), %d genes.png', n_genes),units="in",res=600,
      width=8.23,height=6.3,pointsize=12)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
       cex.axis = 1.5, cex.main = 2)
  dev.off()
  par(.pardefault)
  
  ### remove outliers
  # abline(h = 15, col = "red");
  # clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
  # table(clust)
  # keepSamples = (clust==1)
  # datExpr = datExpr0[keepSamples, ]
  # nGenes = ncol(datExpr)
  # nSamples = nrow(datExpr)
  
  if(Pars$WGCNA.threads %in% c('auto', '{auto}')){
    enableWGCNAThreads(nThreads = round(detectCores() / max(1, (dim(WGCNA.schema)[2] - 1))))
  } else enableWGCNAThreads(nThreads = Pars$WGCNA.threads)
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, blockSize = 20000)
  # Plot the results:
  png(filename = sprintf('Fit index vs power, %d genes.png', n_genes),units="in",res=600,
      width=13,height=7.5,pointsize=12)
  #sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  par(.pardefault)
  
  
  # png(filename = sprintf('Scale-free topology fit index, %d genes.png', n_genes),units="in",res=600,
  #     width=8.23,height=6.3,pointsize=12)
  # plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  #      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
  #      main = paste("Scale independence"));
  # text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  #      labels=powers,cex=cex1,col="red");
  # dev.off()
  # 
  # png(filename = sprintf('WGCNA - Mean connectivity, %d genes.png', n_genes),units="in",res=600,
  #     width=8.23,height=6.3,pointsize=12)
  # plot(sft$fitIndices[,1], sft$fitIndices[,5],
  #      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
  #      main = paste("Mean connectivity"))
  # text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  # dev.off()
  
  
  ## detecting co-expression networks
  
  if(Pars$WGCNA.power %in% c('auto', '{auto}')){
    x.values = sft$fitIndices[,1]
    y.values = -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
    
    set.seed(10)
    
    fit <- lm(y.values ~ x.values + I(x.values^2) + I(x.values^3))
    mc1 = model$coefficients[1]
    mc2 = model$coefficients[2]
    mc3 = model$coefficients[3]
    mc4 = model$coefficients[4]
    
    #Pars$WGCNA.fit.desired.slope = 1 / 20
    #mc1 + mc2*x + mc3*x*x + mc4*x*x*x
    #(mc2 - Pars$WGCNA.fit.desired.slope) + 2*mc3*x + 3*mc4*x*x = 0
    
    D = 4*mc3*mc3 - 4 * (mc2 - Pars$WGCNA.fit.desired.slope) * 3 * mc4
    if(D < 0){
      current.WGCNA.power = 12
    } else {
      current.WGCNA.power = (-2*mc3 - D^0.5) / 2 / 3 / mc4
      intercept = mc1 + mc2*current.WGCNA.power + mc3*(current.WGCNA.power^2) + mc4*(current.WGCNA.power^3) - Pars$WGCNA.fit.desired.slope * current.WGCNA.power
  
      src.data = data.frame(x.values, y.values)
      
      prd <- data.frame(hp = seq(from = min(x.values), to = max(x.values), length.out = 100))
      prd <- data.frame(hp = x.values)
      err <- predict(fit, newdata = prd, se.fit = TRUE)
      
      prd$lci <- err$fit - 1.96 * err$se.fit
      prd$fit <- err$fit
      prd$uci <- err$fit + 1.96 * err$se.fit
      
      g = ggplot(prd, aes(x = hp, y = fit)) +
        theme_bw() +
        geom_line() +
        geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity") + 
        geom_point(data = src.data, aes(x = x.values, y = y.values)) + 
        geom_vline(xintercept = current.WGCNA.power, linetype="dotted", color = "blue", size=1.5) +
        geom_abline(slope = Pars$WGCNA.fit.desired.slope, intercept = intercept, color = "orange") + 
        ylab('model fit') + xlab('Soft Threshold (power)') + ggtitle('Scale Free Topology Model Fit, signed R^2\nApproxomation with polynomial model')
        
  
      ggsave(filename = 'Scale-free topology fit index.png', plot = g, width = 12, height = 10, units = 'in', dpi = 300, limitsize = FALSE)    
    
    }

  } else   current.WGCNA.power = Pars$WGCNA.power
  
  net = blockwiseModules(datExpr, power = round(current.WGCNA.power),
                         TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = sprintf("%s_WGCNA", WGCNA.Analysis.name),
                         maxBlockSize = 30000,
                         verbose = 3)
  

  #table(net$colors) ## The label 0 is reserved for genes outside of all modules.

  # open a graphics window
  #sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath



  # plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
  #                     "Module colors",
  #                     dendroLabels = FALSE, hang = 0.03,
  #                     addGuide = TRUE, guideHang = 0.05)
  
  png(filename = sprintf('WGCNA - net cluster dendrogram, %d genes.png', n_genes), units="in",res=900,
      width=8.23, height=6.3, pointsize=8)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Net module color-ID",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  

  moduleLabels = net$colors  ### an array of the correspondence between gene and module
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs
  geneTree = net$dendrograms[[1]]
  
  ord = order(moduleLabels)
  datExpr.ord = t(datExpr)[ord, ]
  moduleLabels.ord = moduleLabels[ord]
  moduleColors.ord = moduleColors[ord]

  WGCNA.Data$WGCNA.performed = TRUE
  #WGCNA.Data$included.samples = included.samples
  WGCNA.Data$WGCNA.group = WGCNA.group
  WGCNA.Data$MEs = MEs
  WGCNA.Data$datExpr.ord = datExpr.ord
  WGCNA.Data$moduleLabels.ord = moduleLabels.ord
  WGCNA.Data$moduleColors.ord = moduleColors.ord
  WGCNA.Data$WGCNA.Analysis.name = WGCNA.Analysis.name

  saveRDS(WGCNA.Data, file = sprintf("WGCNA.data.rds"))

}