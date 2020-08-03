node.map <-
function(mol.data=NULL, node.data, node.types=c("gene", "ortholog", "compound")[1], node.sum =c("sum","mean", "median", "max", "max.abs", "random")[1], entrez.gnodes=TRUE,cpm.data=NULL){
	#### GK insertion 7
	if (!is.null(cpm.data)) node.sum='cpm.aware'

	## checking if cpm.data is data.matrix and converting it to named vector
	if (!is.vector(cpm.data) & !is.null(cpm.data)){
		cpm.data.mod = as.vector(cpm.data)
		names(cpm.data.mod) = as.character(rownames(cpm.data))
		cpm.data = cpm.data.mod
	}
	#### END

	type.sel=node.data$type %in% node.types
	if(sum(type.sel)<1){
	  message("Note: ", "No specified node types in the pathway!")
	  plot.data=NULL
	  return(plot.data)
	}
	node.data=lapply(node.data, "[", type.sel)
	n.nodes=length(node.data$kegg.names)
	spacials=as.matrix(as.data.frame(node.data[c("type", "x", "y", "width", "height")]))
	if(node.types[1]=="gene"){
	kng=node.data$kegg.names[node.data$type=="gene"]
	kng.char=gsub("[0-9]", "", unlist(kng))
	if(any(kng.char>"")) entrez.gnodes=FALSE
	}

	na.plot.data=function(){
		sapply(1:n.nodes, function(i){
		  kns=node.data$kegg.names[[i]]
		if(node.types[1]=="gene" & entrez.gnodes) items=as.numeric(kns)
		else items=kns
		ord=order(items)
		items=items[ord]
		 kns=kns[ord]
		  return(c(kns[1],"", spacials[i,], NA))
		})
	  }

	if(is.null(mol.data)){
	  plot.data=na.plot.data()
	} else{
	  
	#map gene data  
		if(is.character(mol.data)){
		  gd.names=mol.data
		  mol.data=rep(1, length(mol.data))
		  names(mol.data)=gd.names
		}
		mol.data=cbind(mol.data)

		if(is.null(colnames(mol.data))) colnames(mol.data)=paste("ge", 1:ncol(mol.data),sep="")
		mapped.mols <- intersect(unlist(node.data$kegg.names), row.names(mol.data))
		if(length(mapped.mols)==0){
		  message("Note: ", paste("None of the genes or compounds mapped to the pathway!",
						"Argument gene.idtype or cpd.idtype may be wrong.", sep="\n"))
		  plot.data=na.plot.data()
		} else{
		  if(node.types[1]=="gene" & entrez.gnodes) mapped.mols =as.numeric(mapped.mols)


		  plot.data=sapply(1:n.nodes, function(i){
			kns=node.data$kegg.names[[i]]
			if(node.types[1]=="gene" & entrez.gnodes) items=as.numeric(kns)
			else items=kns
			ord=order(items)
			items=items[ord]
			kns=kns[ord]
			hit=items %in% mapped.mols 
			if(sum(hit)==0) {
			  return(c(kns[1], "", spacials[i,], rep(NA, ncol(mol.data))))
			} else if(sum(hit)==1) {
			  edata=mol.data[as.character(items[hit]),]
			  return(c(kns[hit], kns[hit], spacials[i,], edata))

	##### GK insertion 6
			} else if (node.sum=='cpm.aware') {
			  edata=apply(cbind(mol.data[as.character(items[hit]),]), 2, function(x){

				x = x[!is.na(x)]
				if(length(x)<1) { return(NA) 
				} else {

					## calculating weighted sum LogFC
					weighted.sum = 0       ## sum of LogFC1*CPM1 + LogFC2*CPM2 + ...
					sum.of.weights = 0     ## sum of CPM1 + CPM2 + CPM3 + ...
					
					for (src.id in names(x)){
					  weighted.sum = weighted.sum + cpm.data[src.id]*x[src.id]
					  sum.of.weights = sum.of.weights + cpm.data[src.id]
					}
					return(weighted.sum/sum.of.weights)    ## weighted avg
				}
			  })
			  return(c(kns[hit][1], paste(kns[hit],collapse=","), spacials[i,], edata))
			} else {
	#####
			  node.sum=eval(as.name(node.sum))
											#      edata=apply(cbind(mol.data[as.character(items[hit]),]), 2, node.sum, na.rm=T)
			  edata=apply(cbind(mol.data[as.character(items[hit]),]), 2, function(x){
				x=x[!is.na(x)]
				if(length(x)<1) { return(NA) 
				} else  {
					return(node.sum(x, na.rm=F))
				}
			  })
			  return(c(kns[hit][1], paste(kns[hit],collapse=","), spacials[i,], edata))
			}    
		  })
		}
	  }

	colnames(plot.data)=names(node.data$kegg.names)
	plot.data=as.data.frame(t(plot.data), stringsAsFactors = F)
	  plot.data$labels=node.data$labels
	  ncs=ncol(plot.data)
	  plot.data=plot.data[,c(1,ncs,2:(ncs-1))]
	if(is.null(mol.data)) cns="mol.data" else cns=colnames(mol.data)
	colnames(plot.data)[c(1,3,9:ncs)]=c("kegg.names","all.mapped",cns)#c(1,8:ncs)
	for(ic in (1:ncol(plot.data))[-c(1:4)]) plot.data[,ic]=as.numeric(plot.data[,ic])#-c(1:3)

	return(plot.data)
	}

