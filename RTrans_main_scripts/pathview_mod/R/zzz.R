.onLoad <- function(libname, pkgname) {
  pnames=rownames(installed.packages())
  if("pathviewmod" %in% pnames){
    data(gene.idtype.list, package ="pathviewmod")
    data(gene.idtype.bods, package ="pathviewmod")
    data(cpd.simtypes, package ="pathviewmod")
  }
disclaimer="##############################################################################\nPathview is an open source software package distributed under GNU General Public License version 3 (GPLv3). Details of GPLv3 is available at http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to formally cite the original Pathview paper (not just mention it) in publications or products. For details, do citation(\"pathviewmod\") within R. \n\nThe pathviewmod downloads and uses KEGG data. Non-academic uses may require a KEGG license agreement (details at http://www.kegg.jp/kegg/legal.html).\n##############################################################################\n\n"
packageStartupMessage(wordwrap(disclaimer, 80))
}
