mol.sum <-
function(mol.data, id.map, gene.annotpkg="org.Hs.eg.db", cpm.data = NULL,debug.info=F,sum.method=c("sum","mean", "median", "max", "max.abs", "random","cpm.aware")[1]){
  #### GK insertion 1
  if (is.null(cpm.data) & sum.method=='cpm.aware')	stop("CPM-aware mapping needs the expression data")
  #### END
  if(is.character(mol.data)){
    gd.names=mol.data
    mol.data=rep(1, length(mol.data))
    names(mol.data)=gd.names
    ng=length(mol.data)
  } else if(!is.null(mol.data)){
    if(length(dim(mol.data))==2){
    gd.names=rownames(mol.data)
    ng=nrow(mol.data)
    } else if(is.numeric(mol.data) & is.null(dim(mol.data))){
    gd.names=names(mol.data)
    ng=length(mol.data)
    } else stop("wrong mol.data format!")
  } else stop("NULL mol.data!")

  if(is.character(id.map) & length(id.map)==1){
    id.map=id2eg(gd.names, category=id.map, pkg.name=gene.annotpkg)
  }
  
  sel.idx=id.map[,2]>"" & !is.na(id.map[,2])
  id.map=id.map[sel.idx,]
  eff.idx=gd.names %in% id.map[,1]
  mapped.ids=id.map[match(gd.names[eff.idx], id.map[,1]),2]
  
  #### GK insertion 2
  if(sum.method == "cpm.aware"){
	## checking if cpm.data is data.matrix and converting it to named vector
	if (!is.vector(cpm.data)){
		cpm.data.mod = as.vector(cpm.data)
		names(cpm.data.mod) = as.character(rownames(cpm.data))
		cpm.data = cpm.data.mod
	}
	
    mapped.data=apply(cbind(cbind(mol.data)[eff.idx,]),2,function(x){
		names(mapped.ids) = names(x)

		## removing NAs
		x = x[!is.na(mapped.ids)]
		mapped.ids = mapped.ids[!is.na(mapped.ids)]

		## calculating weighted sum LogFC
		mapped.ids.list = levels(factor(mapped.ids))
		mapped.ids.count = length(mapped.ids.list)
		weighted.sums = rep(0,mapped.ids.count)      ## sum of LogFC1*CPM1 + LogFC2*CPM2 + ...
		names(weighted.sums) = mapped.ids.list
		sums.of.weights = rep(0,mapped.ids.count)    ## sum of CPM1 + CPM2 + CPM3 + ...
		names(sums.of.weights) = mapped.ids.list
		for (src.id in names(x)){
		  cur.mapped.id = mapped.ids[src.id]
		  weighted.sums[cur.mapped.id] = weighted.sums[cur.mapped.id] + cpm.data[src.id]*x[src.id]
		  sums.of.weights[cur.mapped.id] = sums.of.weights[cur.mapped.id] + cpm.data[src.id]
		}

		return(weighted.sums/sums.of.weights)   ## weighted avg
    })
  } else if(sum.method %in% c("sum","mean")){
  #### END
	sum.method=eval(as.name(sum.method))
    mapped.data=apply(cbind(cbind(mol.data)[eff.idx,]),2,function(x){
      sum.res=tapply(x, mapped.ids, sum.method, na.rm=T)
      return(sum.res)
    })
  } else{
    sum.method=eval(as.name(sum.method))
    mol.data=cbind(cbind(mol.data)[eff.idx,])
    if(all(mol.data>=0) | all(mol.data<=0)){
      vars=apply(cbind(mol.data), 1, IQR)
    } else vars=apply(cbind(mol.data), 1, sum, na.rm=T)
    
    sel.rn=tapply(1:sum(eff.idx), mapped.ids, function(x){
      if(length(x)==1) return(x)
      else return(x[which.min(abs(vars[x]-sum.method(vars[x], na.rm=T)))])
    })
    mapped.data=cbind(mol.data[sel.rn,])
    rownames(mapped.data)=names(sel.rn)
  }
  return(mapped.data)
}

