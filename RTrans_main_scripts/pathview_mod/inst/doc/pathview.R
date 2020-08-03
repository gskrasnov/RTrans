### R code from vignette source 'pathviewmod.Rnw'

###################################################
### code chunk number 1: synopsis1 (eval = FALSE)
###################################################
## library(pathviewmod)
## data(gse16873.d)
## pv.out <- pathviewmod(gene.data = gse16873.d[, 1], pathway.id = "04110",
##                    species = "hsa", out.suffix = "gse16873")


###################################################
### code chunk number 2: install (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("pathviewmod")


###################################################
### code chunk number 3: <install (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("Rgraphviz", "png", "KEGGgraph", "org.Hs.eg.db"))


###################################################
### code chunk number 4: install (eval = FALSE)
###################################################
## install.packages("pathviewmod",repos="http://R-Forge.R-project.org")


###################################################
### code chunk number 5: install (eval = FALSE)
###################################################
## install.packages("/your/local/directory/pathviewmod_1.0.0.tar.gz", 
##     repos = NULL, type = "source")


###################################################
### code chunk number 6: <install (eval = FALSE)
###################################################
## install.packages("/your/local/directory/XML_3.95-0.2.zip", repos = NULL)


###################################################
### code chunk number 7: start
###################################################
options(width=80)


###################################################
### code chunk number 8: start
###################################################
library(pathviewmod)


###################################################
### code chunk number 9: start (eval = FALSE)
###################################################
## library(help=pathviewmod)


###################################################
### code chunk number 10: start (eval = FALSE)
###################################################
## help(pathviewmod)
## ?pathviewmod


###################################################
### code chunk number 11: dataPrep.gse16873
###################################################
data(gse16873.d)


###################################################
### code chunk number 12: readin (eval = FALSE)
###################################################
## filename=system.file("extdata/gse16873.demo", package = "pathviewmod")
## gse16873=read.delim(filename, row.names=1)
## gse16873.d=gse16873[,2*(1:6)]-gse16873[,2*(1:6)-1]


###################################################
### code chunk number 13: dataPrep.demo.paths
###################################################
data(demo.paths)


###################################################
### code chunk number 14: dataPreppaths.hsa
###################################################
data(paths.hsa)
head(paths.hsa,3)


###################################################
### code chunk number 15: kegg.native
###################################################
i <- 1
pv.out <- pathviewmod(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873", kegg.native = T)
list.files(pattern="hsa04110", full.names=T)
str(pv.out)
head(pv.out$plot.data.gene)


###################################################
### code chunk number 16: kegg.native_2layer
###################################################
pv.out <- pathviewmod(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.2layer", kegg.native = T,
                   same.layer = F)


###################################################
### code chunk number 17: graphviz
###################################################
pv.out <- pathviewmod(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873", kegg.native = F,
                   sign.pos = demo.paths$spos[i])
#pv.out remains the same
dim(pv.out$plot.data.gene)
head(pv.out$plot.data.gene)


###################################################
### code chunk number 18: graphviz.2layer
###################################################
pv.out <- pathviewmod(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i], 
    species = "hsa", out.suffix = "gse16873.2layer", kegg.native = F, 
    sign.pos = demo.paths$spos[i], same.layer = F)


###################################################
### code chunk number 19: split
###################################################
pv.out <- pathviewmod(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i], 
    species = "hsa", out.suffix = "gse16873.split", kegg.native = F, 
    sign.pos = demo.paths$spos[i], split.group = T)
dim(pv.out$plot.data.gene)
head(pv.out$plot.data.gene)


###################################################
### code chunk number 20: expanded
###################################################
pv.out <- pathviewmod(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i], 
    species = "hsa", out.suffix = "gse16873.split.expanded", kegg.native = F,
    sign.pos = demo.paths$spos[i], split.group = T, expand.node = T)
dim(pv.out$plot.data.gene)
head(pv.out$plot.data.gene)


###################################################
### code chunk number 21: dataPrep.sim.cpd
###################################################
sim.cpd.data=sim.mol.data(mol.type="cpd", nmol=3000)
data(cpd.simtypes)


###################################################
### code chunk number 22: gene_cpd.data
###################################################
i <- 3
print(demo.paths$sel.paths[i])
pv.out <- pathviewmod(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data, 
    pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873.cpd", 
    keys.align = "y", kegg.native = T, key.pos = demo.paths$kpos1[i])
str(pv.out)
head(pv.out$plot.data.cpd)


###################################################
### code chunk number 23: graphviz.gene_cpd.data
###################################################
pv.out <- pathviewmod(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data, 
    pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873.cpd", 
    keys.align = "y", kegg.native = F, key.pos = demo.paths$kpos2[i], 
    sign.pos = demo.paths$spos[i], cpd.lab.offset = demo.paths$offs[i])


###################################################
### code chunk number 24: sim.cpd.data2
###################################################
set.seed(10)
sim.cpd.data2 = matrix(sample(sim.cpd.data, 18000, 
    replace = T), ncol = 6)
rownames(sim.cpd.data2) = names(sim.cpd.data)
colnames(sim.cpd.data2) = paste("exp", 1:6, sep = "")
head(sim.cpd.data2, 3)


###################################################
### code chunk number 25: multisample.gene_cpd.data
###################################################
#KEGG view
pv.out <- pathviewmod(gene.data = gse16873.d[, 1:3], 
    cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i], 
    species = "hsa", out.suffix = "gse16873.cpd.3-2s", keys.align = "y", 
    kegg.native = T, match.data = F, multi.state = T, same.layer = T)
head(pv.out$plot.data.cpd)
#KEGG view with data match
pv.out <- pathviewmod(gene.data = gse16873.d[, 1:3], 
    cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i], 
    species = "hsa", out.suffix = "gse16873.cpd.3-2s.match", 
    keys.align = "y", kegg.native = T, match.data = T, multi.state = T, 
    same.layer = T)
#graphviz view
pv.out <- pathviewmod(gene.data = gse16873.d[, 1:3], 
    cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i], 
    species = "hsa", out.suffix = "gse16873.cpd.3-2s", keys.align = "y", 
    kegg.native = F, match.data = F, multi.state = T, same.layer = T,
    key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i])


###################################################
### code chunk number 26: multisample.gene_cpd.data.seperate (eval = FALSE)
###################################################
## #plot samples/states separately
## pv.out <- pathviewmod(gene.data = gse16873.d[, 1:3], 
##     cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i], 
##     species = "hsa", out.suffix = "gse16873.cpd.3-2s", keys.align = "y", 
##     kegg.native = T, match.data = F, multi.state = F, same.layer = T)


###################################################
### code chunk number 27: multisample.gene_cpd.data.2layer (eval = FALSE)
###################################################
## pv.out <- pathviewmod(gene.data = gse16873.d[, 1:3], 
##     cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i], 
##     species = "hsa", out.suffix = "gse16873.cpd.3-2s.2layer", 
##     keys.align = "y", kegg.native = T, match.data = F, multi.state = T, 
##     same.layer = F)


###################################################
### code chunk number 28: discrete.gene_cpd.data
###################################################
require(org.Hs.eg.db)
gse16873.t <- apply(gse16873.d, 1, function(x) t.test(x, 
    alternative = "two.sided")$p.value)
sel.genes <- names(gse16873.t)[gse16873.t < 0.1]
sel.cpds <- names(sim.cpd.data)[abs(sim.cpd.data) > 0.5]
pv.out <- pathviewmod(gene.data = sel.genes, cpd.data = sel.cpds, 
    pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "sel.genes.sel.cpd", 
    keys.align = "y", kegg.native = T, key.pos = demo.paths$kpos1[i], 
    limit = list(gene = 5, cpd = 2), bins = list(gene = 5, cpd = 2), 
    na.col = "gray", discrete = list(gene = T, cpd = T))
pv.out <- pathviewmod(gene.data = sel.genes, cpd.data = sim.cpd.data, 
    pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "sel.genes.cpd",
    keys.align = "y", kegg.native = T, key.pos = demo.paths$kpos1[i], 
    limit = list(gene = 5, cpd = 1), bins = list(gene = 5, cpd = 10), 
    na.col = "gray", discrete = list(gene = T, cpd = F))


###################################################
### code chunk number 29: gene.ensprot_cpd.cas
###################################################
cpd.cas <- sim.mol.data(mol.type = "cpd", id.type = cpd.simtypes[2], 
    nmol = 10000)
gene.ensprot <- sim.mol.data(mol.type = "gene", id.type = gene.idtype.list[4], 
    nmol = 50000)
pv.out <- pathviewmod(gene.data = gene.ensprot, cpd.data = cpd.cas, 
    gene.idtype = gene.idtype.list[4], cpd.idtype = cpd.simtypes[2], 
    pathway.id = demo.paths$sel.paths[i], species = "hsa", same.layer = T, 
    out.suffix = "gene.ensprot.cpd.cas", keys.align = "y", kegg.native = T, 
    key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
    limit = list(gene = 3, cpd = 3), bins = list(gene = 6, cpd = 6))


###################################################
### code chunk number 30: gene.ensprot_cpd.cas.manual.map
###################################################
id.map.cas <- cpdidmap(in.ids = names(cpd.cas), in.type = cpd.simtypes[2], 
    out.type = "KEGG COMPOUND accession")
cpd.kc <- mol.sum(mol.data = cpd.cas, id.map = id.map.cas)
id.map.ensprot <- id2eg(ids = names(gene.ensprot), 
    category = gene.idtype.list[4], org = "Hs")
gene.entrez <- mol.sum(mol.data = gene.ensprot, id.map = id.map.ensprot)
pv.out <- pathviewmod(gene.data = gene.entrez, cpd.data = cpd.kc, 
    pathway.id = demo.paths$sel.paths[i], species = "hsa", same.layer = T, 
    out.suffix = "gene.entrez.cpd.kc", keys.align = "y", kegg.native = T, 
    key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
    limit = list(gene = 3, cpd = 3), bins = list(gene = 6, cpd = 6))


###################################################
### code chunk number 31: korg
###################################################
data(korg)
head(korg)
#number of species which use Entrez Gene as the default ID
sum(korg[,"entrez.gnodes"]=="1",na.rm=T)
#number of species which use other ID types or none as the default ID
sum(korg[,"entrez.gnodes"]=="0",na.rm=T)
#species which do not have Entrez Gene annotation at all
na.idx=is.na(korg[,"ncbi.geneid"])
korg[na.idx,]


###################################################
### code chunk number 32: bods_gene.idtype.list
###################################################
data(bods)
bods
data(gene.idtype.list)
gene.idtype.list


###################################################
### code chunk number 33: eco.dat.kegg
###################################################
eco.dat.kegg <- sim.mol.data(mol.type="gene",id.type="kegg",species="eco",nmol=3000)
head(eco.dat.kegg)
pv.out <- pathviewmod(gene.data = eco.dat.kegg, gene.idtype="kegg",
    pathway.id = "00640", species = "eco", out.suffix = "eco.kegg",
    kegg.native = T, same.layer=T)


###################################################
### code chunk number 34: eco.dat.kegg
###################################################
eco.dat.entrez <- sim.mol.data(mol.type="gene",id.type="entrez",species="eco",nmol=3000)
head(eco.dat.entrez)
pv.out <- pathviewmod(gene.data = eco.dat.entrez, gene.idtype="entrez",
    pathway.id = "00640", species = "eco", out.suffix = "eco.entrez",
    kegg.native = T, same.layer=T)


###################################################
### code chunk number 35: eco.dat.symbol (eval = FALSE)
###################################################
## egid.eco=eg2id(names(eco.dat.entrez), category="symbol", pkg="org.EcK12.eg.db")
## eco.dat.symbol <- eco.dat.entrez
## names(eco.dat.symbol) <- egid.eco[,2]
## head(eco.dat.kegg)
## pv.out <- pathviewmod(gene.data = eco.dat.symbol, gene.idtype="symbol",
##     pathway.id = "00640", species = "eco", out.suffix = "eco.symbol.2layer",
##     kegg.native = T, same.layer=F)


###################################################
### code chunk number 36: gene.ensprot_cpd.cas.manual.map
###################################################
ko.data=sim.mol.data(mol.type="gene.ko", nmol=5000)
pv.out <- pathviewmod(gene.data = ko.data, pathway.id = "04112",
                   species = "ko", out.suffix = "ko.data", kegg.native = T)


###################################################
### code chunk number 37: GAGE.Pathview.pipeline (eval = FALSE)
###################################################
## library(gage)
## data(gse16873)
## cn <- colnames(gse16873)
## hn <- grep('HN',cn, ignore.case =TRUE)
## dcis <- grep('DCIS',cn, ignore.case =TRUE)
## data(kegg.gs)
## #pathway analysis using gage
## gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs,
##     ref = hn, samp = dcis)
## #prepare the differential expression data
## gse16873.d <- gagePrep(gse16873, ref = hn, samp = dcis)
## #equivalently, you can do simple subtraction for paired samples
## gse16873.d <- gse16873[,dcis]-gse16873[,hn]
## #select significant pathways and extract their IDs
## sel <- gse16873.kegg.p$greater[, "q.val"] < 0.1 & !is.na(gse16873.kegg.p$greater[, 
##     "q.val"])
## path.ids <- rownames(gse16873.kegg.p$greater)[sel]
## path.ids2 <- substr(path.ids[c(1, 2, 7)], 1, 8)
## #pathviewmod visualization
## pv.out.list <- sapply(path.ids2, function(pid) pathviewmod(gene.data = gse16873.d[, 
##     1:2], pathway.id = pid, species = "hsa")) 


