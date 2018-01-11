suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(GOSemSim))

get_KEGG_Env <- function() {
  if (! exists(".KEGG_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".KEGG_clusterProfiler_Env", new.env(), envir=envir)
  }
  get(".KEGG_clusterProfiler_Env", envir = .GlobalEnv)
}

download_KEGG <- function(species, keggType="KEGG", keyType="kegg") {
  KEGG_Env <- get_KEGG_Env()
  
  use_cached <- FALSE
  
  if (exists("organism", envir = KEGG_Env, inherits = FALSE) &&
      exists("_type_", envir = KEGG_Env, inherits = FALSE) ) {
    
    org <- get("organism", envir=KEGG_Env)
    type <- get("_type_", envir=KEGG_Env)
    
    if (org == species && type == keggType &&
        exists("KEGGPATHID2NAME", envir=KEGG_Env, inherits = FALSE) &&
        exists("KEGGPATHID2EXTID", envir=KEGG_Env, inherits = FALSE)) {
      
      use_cached <- TRUE
    }
  }
  
  if (use_cached) {
    KEGGPATHID2EXTID <- get("KEGGPATHID2EXTID", envir=KEGG_Env)
    KEGGPATHID2NAME <- get("KEGGPATHID2NAME", envir=KEGG_Env)
  } else {
    if (keggType == "KEGG") {
      kres <- download.KEGG.Path(species)
    } else {
      kres <- download.KEGG.Module(species)
    }
    
    KEGGPATHID2EXTID <- kres$KEGGPATHID2EXTID
    KEGGPATHID2NAME <- kres$KEGGPATHID2NAME
    
    assign("organism", species, envir=KEGG_Env)
    assign("_type_", keggType, envir=KEGG_Env)
    assign("KEGGPATHID2NAME", KEGGPATHID2NAME, envir=KEGG_Env)
    assign("KEGGPATHID2EXTID", KEGGPATHID2EXTID, envir=KEGG_Env)
  }
  
  if (keyType != "kegg") {
    need_idconv <- FALSE
    idconv <- NULL
    if (use_cached &&
        exists("key", envir=KEGG_Env, inherits = FALSE) &&
        exists("idconv", envir=KEGG_Env, inherits = FALSE)) {
      
      key <- get("key", envir=KEGG_Env)
      if (key == keyType) {
        idconv <- get("idconv", envir=KEGG_Env)
      } else {
        need_idconv <- TRUE
      }
    } else {
      neec_idconv <- TRUE
    }
    
    if (need_idconv || is.null(idconv)) {
      idconv <- KEGG_convert("kegg", keyType, species)
      assign("key", keyType, envir=KEGG_Env)
      assign("idconv", idconv, envir=KEGG_Env)
    }
    colnames(KEGGPATHID2EXTID) <- c("from", "kegg")
    KEGGPATHID2EXTID <- merge(KEGGPATHID2EXTID, idconv, by.x='kegg', by.y='from')
    KEGGPATHID2EXTID <- unique(KEGGPATHID2EXTID[, -1])
  }
  
  return(list(KEGGPATHID2EXTID = KEGGPATHID2EXTID,
              KEGGPATHID2NAME  = KEGGPATHID2NAME))
}

prepare_KEGG <- function(species, KEGG_Type="KEGG", keyType="kegg") {
  kegg <- download_KEGG(species, KEGG_Type, keyType)
  build_Anno(kegg$KEGGPATHID2EXTID,
             kegg$KEGGPATHID2NAME)
}

download.KEGG.Path <- function(species) {
  keggpathid2extid.df <- kegg_link(species, "pathway")
  if (is.null(keggpathid2extid.df))
    stop("'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'...")
  keggpathid2extid.df[,1] %<>% gsub("[^:]+:", "", .)
  keggpathid2extid.df[,2] %<>% gsub("[^:]+:", "", .)
  
  keggpathid2name.df <- kegg_list("pathway")
  keggpathid2name.df[,1] %<>% gsub("path:map", species, .)
  
  ## if 'species="ko"', ko and map path are duplicated, only keep ko path.
  ##
  ## http://www.kegg.jp/dbget-bin/www_bget?ko+ko00010
  ## http://www.kegg.jp/dbget-bin/www_bget?ko+map0001
  ##
  keggpathid2extid.df <- keggpathid2extid.df[keggpathid2extid.df[,1] %in% keggpathid2name.df[,1],]
  
  return(list(KEGGPATHID2EXTID=keggpathid2extid.df,
              KEGGPATHID2NAME=keggpathid2name.df))
}

download.KEGG.Module <- function(species) {
  keggmodule2extid.df <- kegg_link(species, "module")
  if (is.null(keggmodule2extid.df)) {
    stop("'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'...")
  }
  
  keggmodule2extid.df[,1] %<>% gsub("[^:]+:", "", .) %>% gsub(species, "", .) %>% gsub("^_", "", .)
  keggmodule2extid.df[,2] %<>% gsub("[^:]+:", "", .)
  
  keggmodule2name.df <- kegg_list("module")
  keggmodule2name.df[,1] %<>% gsub("md:", "", .)
  return(list(KEGGPATHID2EXTID=keggmodule2extid.df,
              KEGGPATHID2NAME =keggmodule2name.df))
}



browseKEGG <- function(x, pathID) {
  url <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?", pathID, '/', x[pathID, "geneID"])
  browseURL(url)
  invisible(url)
}


search_kegg_organism <- function(str, by="scientific_name", ignore.case=FALSE) {
  by <- match.arg(by, c("kegg_code", "scientific_name", "common_name"))
  kegg_species <- kegg_species_data()
  idx <- grep(str, kegg_species[, by], ignore.case = ignore.case)
  kegg_species[idx,]
}


kegg_species_data <- function() {
  utils::data(list="kegg_species", package="clusterProfiler")
  get("kegg_species", envir = .GlobalEnv)
}

get_kegg_species <- function() {
  pkg <- "XML"
  requireNamespace(pkg)
  readHTMLTable <- eval(parse(text="XML::readHTMLTable"))
  x <- readHTMLTable("http://www.genome.jp/kegg/catalog/org_list.html")
  
  y <- get_species_name(x[[2]], "Eukaryotes")
  y2 <- get_species_name(x[[3]], 'Prokaryotes')
  
  sci_name <- gsub(" \\(.*$", '', y[,2])
  com_name <- gsub("[^\\(]+ \\(([^\\)]+)\\)$", '\\1', y[,2])
  eu <- data.frame(kegg_code=unlist(y[,1]),
                   scientific_name = sci_name,
                   common_name = com_name,
                   stringsAsFactors = FALSE)
  pr <- data.frame(kegg_code=unlist(y2[,1]),
                   scientific_name = unlist(y2[,2]),
                   common_name = NA,
                   stringsAsFactors = FALSE)
  kegg_species <- rbind(eu, pr)
  save(kegg_species, file="kegg_species.rda")
  invisible(kegg_species)
}

get_species_name <- function(y, table) {
  idx <- get_species_name_idx(y, table)
  t(sapply(1:nrow(idx), function(i) {
    y[] = lapply(y, as.character)
    y[i, idx[i,]]
  }))
}


get_species_name_idx <- function(y, table='Eukaryotes') {
  table <- match.arg(table, c("Eukaryotes", "Prokaryotes"))
  t(apply(y, 1, function(x) {
    ii <- which(!is.na(x))
    n <- length(ii)
    if (table == "Eukaryotes") {
      return(ii[(n-2):(n-1)])
    } else {
      return(ii[(n-3):(n-2)])
    }
  }))
}


kegg_rest <- function(rest_url) {
  content <- tryCatch(suppressWarnings(readLines(rest_url)), error=function(e) NULL)
  if (is.null(content))
    return(content)
  
  content %<>% strsplit(., "\t") %>% do.call('rbind', .)
  res <- data.frame(from=content[,1],
                    to=content[,2])
  return(res)
}

## http://www.genome.jp/kegg/rest/keggapi.html
## kegg_link('hsa', 'pathway')
kegg_link <- function(target_db, source_db) {
  url <- paste0("http://rest.kegg.jp/link/", target_db, "/", source_db, collapse="")
  kegg_rest(url)
}


kegg_list <- function(db) {
  url <- paste0("http://rest.kegg.jp/list/", db, collapse="")
  kegg_rest(url)
}


ko2name <- function(ko) {
  p <- kegg_list('pathway')
  ko2 <- gsub("^ko", "path:map", ko)
  ko.df <- data.frame(ko=ko, from=ko2)
  res <- merge(ko.df, p, by = 'from', all.x=TRUE)
  res <- res[, c("ko", "to")]
  colnames(res) <- c("ko", "name")
  return(res)
}



idType <- function(OrgDb = "org.Hs.eg.db") {
  db <- load_OrgDb(OrgDb)
  keytypes(db)
}

bitr <- function(geneID, fromType, toType, OrgDb, drop=TRUE) {
  idTypes <- idType(OrgDb)
  msg <-  paste0("should be one of ", paste(idTypes, collapse=", "), ".")
  if (! fromType %in% idTypes) {
    stop("'fromType' ", msg)
  }
  if (! all(toType %in% idTypes)) {
    stop("'toType' ", msg)
  }
  
  geneID %<>% as.character %>% unique
  db <- load_OrgDb(OrgDb)
  res <- suppressWarnings(select(db,
                                 keys = geneID,
                                 keytype = fromType,
                                 columns=c(fromType, toType)))
  
  ii <- which(is.na(res[,2]))
  if (length(ii)) {
    n <- res[ii, 1] %>% unique %>% length
    if (n) {
      warning(paste0(round(n/length(geneID)*100, 2), "%"), " of input gene IDs are fail to map...")
    }
    if (drop) {
      res <- res[-ii, ]
    }
  }
  return(res)
}


bitr_kegg <- function(geneID, fromType, toType, organism, drop=TRUE) {
  id_types <- c("Path", "Module", "ncbi-proteinid", "ncbi-geneid", "uniprot", "kegg")
  fromType <- match.arg(fromType, id_types)
  toType <- match.arg(toType, id_types)
  
  if (fromType == toType)
    stop("fromType and toType should not be identical...")
  
  if (fromType == "Path" || fromType == "Module") {
    idconv <- KEGG_path2extid(geneID, organism, fromType, toType)
  } else if (toType == "Path" || toType == "Module") {
    idconv <- KEGG_extid2path(geneID, organism, toType, fromType)
  } else {
    idconv <- KEGG_convert(fromType, toType, organism)
  }
  
  res <- idconv[idconv[,1] %in% geneID, ]
  n <- sum(!geneID %in% res[,1])
  if (n > 0) {
    warning(paste0(round(n/length(geneID)*100, 2), "%"), " of input gene IDs are fail to map...")
  }
  
  if (! drop && n > 0) {
    misHit <- data.frame(from = geneID[!geneID %in% res[,1]],
                         to = NA)
    res <- rbind(res, misHit)
  }
  colnames(res) <- c(fromType, toType)
  rownames(res) <- NULL
  return(res)
}

KEGG_convert <- function(fromType, toType, species) {
  if (fromType == "kegg" || toType != "kegg") {
    turl <- paste("http://rest.kegg.jp/conv", toType, species, sep='/')
    tidconv <- kegg_rest(turl)
    if (is.null(tidconv))
      stop(toType, " is not supported for ", species, " ...")
    idconv <- tidconv
  }
  
  if (toType == "kegg" || fromType != "kegg") {
    furl <- paste("http://rest.kegg.jp/conv", fromType, species, sep='/')
    fidconv <- kegg_rest(furl)
    if (is.null(fidconv))
      stop(fromType, " is not supported for ", species, " ...")
    idconv <- fidconv
  }
  
  if (fromType != "kegg" && toType != "kegg") {
    idconv <- merge(fidconv, tidconv, by.x='from', by.y='from')
    idconv <- idconv[, -1]
  } else if (fromType != "kegg") {
    idconv <- idconv[, c(2,1)]
  }
  colnames(idconv) <- c("from", "to")
  
  idconv[,1] %<>% gsub("[^:]+:", "", .)
  idconv[,2] %<>% gsub("[^:]+:", "", .)
  return(idconv)
}


KEGG_path2extid <- function(keggID, species=sub("\\d+$", "", keggID),
                            keggType = "Path", keyType = "kegg") {
  path2extid <- KEGGPATHID2EXTID(species, keggType, keyType)
  path2extid[path2extid$from %in% keggID, ]
}

KEGG_extid2path <- function(geneID, species, keggType = "Path", keyType = "kegg") {
  path2extid <- KEGGPATHID2EXTID(species, keggType, keyType)
  res <- path2extid[path2extid$to %in% geneID, ]
  res <- res[, c(2,1)]
  colnames(res) <- colnames(path2extid)
  return(res)
}


KEGGPATHID2EXTID <- function(species, keggType = "Path", keyType = "kegg") {
  keggType <- match.arg(keggType, c("Path", "Module"))
  if (keggType == "Path") {
    keggType <- "KEGG"
  } else {
    keggType <- "MKEGG"
  }
  kegg <- download_KEGG(species, keggType, keyType)
  return(kegg$KEGGPATHID2EXTID)
}

build_Anno <- DOSE:::build_Anno
get_organism <- DOSE:::get_organism

# KEGG_DATA <- prepare_KEGG(species = 'mmu', keyType = 'kegg')
