download.kegg <-
  function (pathway.id = "00010", species = "hsa", kegg.dir = ".", file.type=c("xml", "png"), cached.dir = NULL)
  {
    if(!is.null(cached.dir)) dir.create(cached.dir, showWarnings = FALSE)
	
	npath=length(pathway.id)
    if(species!="ko") species=kegg.species.code(species, na.rm=T)
    nspec=length(species)

    if(npath!=1 | nspec!=1) {
      species=rep(species, npath)
      pathway.id=rep(pathway.id, each=nspec)
    }
    pathway.id <- paste(species, pathway.id, sep = "")
    uidx=!duplicated(pathway.id)
    pathway.id=pathway.id[uidx]
    species=species[uidx]
    npath=length(pathway.id)
    
    xml.fnames=paste(pathway.id, ".xml", sep="")
	  if (!is.null(cached.dir)) cached.xml.fnames = sprintf('%s/%s',cached.dir,xml.fnames)
    png.fnames=paste(pathway.id, ".png", sep="")
	  if (!is.null(cached.dir)) cached.png.fnames = sprintf('%s/%s',cached.dir,png.fnames)
    ##      xml.fmt="http://www.genome.jp/kegg-bin/download?entry=%s&format=kgml"
    ##      png.fmt="http://www.genome.jp/kegg/pathway/%s/%s"
    xml.fmt="http://rest.kegg.jp/get/%s/kgml"
    png.fmt="http://rest.kegg.jp/get/%s/image"
    all.status=rep("succeed", npath)
    names(all.status)=pathway.id
    warn.fmt.xml="Download of %s xml file failed!\nThis pathway may not exist!"
    warn.fmt.png="Download of %s png file failed!\nThis pathway may not exist!"
    

    if("xml" %in% file.type){
      for (i in 1:npath) {
        xml.url=sprintf(xml.fmt,  pathway.id[i])
        xml.target=sprintf("%s/%s", kegg.dir, xml.fnames[i])
    		do.download = T
    		if (!is.null(cached.dir)) {
    		  if (file.exists(cached.xml.fnames[i])){
    		    file.copy(from=cached.xml.fnames[i], to=xml.target)
    		    do.download = F
                msg=sprintf("xml file for %s, %d/%d pathways is found in cache..", pathway.id[i], i, length(pathway.id))
                message("Info: ", msg)
    		  }
    		}
    		if (do.download){
          msg=sprintf("Downloading xml files for %s, %d/%d pathways..", pathway.id[i], i, length(pathway.id))
          message("Info: ", msg)
          xml.status=try(download.file(xml.url, xml.target, quiet=T), silent=T)

          if(xml.status!=0) { all.status[i]="failed"
    		  } else if (!is.null(cached.dir)) file.copy(to=cached.xml.fnames[i], from=xml.target)
          if(class(xml.status)=="try-error"){
            warn.msg=sprintf(warn.fmt.xml, pathway.id[i])
            message("Warning: ", warn.msg)
            unlink(xml.target)
          }
    		}
      }
    }
    
    if("png" %in% file.type){
      for (i in 1:npath) {
        png.url=sprintf(png.fmt,  pathway.id[i])#species[i], png.fnames[i])
        png.target=sprintf("%s/%s", kegg.dir, png.fnames[i])
    		do.download = T
    		if (!is.null(cached.dir)) {
    		  if (file.exists(cached.png.fnames[i])){
    		    file.copy(from=cached.png.fnames[i], to=png.target)
    		    do.download = F
            msg=sprintf("png file for %s, %d/%d pathways is found in cache..", pathway.id[i], i, length(pathway.id))
            message("Info: ", msg)
    		  }
    		}
    		if (do.download){
          msg=sprintf("Downloading png files for %s, %d/%d pathways..", pathway.id[i], i, length(pathway.id))
          message("Info: ", msg)

          png.status=suppressWarnings(try(download.file(png.url, png.target, quiet=T, mode="wb"), silent=T))

          if(png.status!=0) { all.status[i]="failed"
    		  } else if (!is.null(cached.dir)) file.copy(to=cached.png.fnames[i], from=png.target)
          if(class(png.status)=="try-error"){
            warn.msg=sprintf(warn.fmt.png, pathway.id[i])
            message("Warning: ", warn.msg)
            unlink(png.target)
    		  }
        }
      }
    }

    return(all.status)
  }

