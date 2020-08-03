parsePathwayInfo <- function(root) {
  attrs <- xmlAttrs(root)
  ## required: name, org, number
  name <- attrs[["name"]]
  org <- attrs[["org"]]
  number <- attrs[["number"]]
  ## implied: title, image, link
  title <- getNamedElement(attrs, "title")
  image <- getNamedElement(attrs, "image")
  link <- getNamedElement(attrs, "link")

  return(new("KEGGPathwayInfo",
             name=name,
             org=org,
             number=number,
             title=title,
             image=image,
             link=link))
}

parseGraphics <- function(graphics) {
  if(is.null(graphics))
    return(new("KEGGGraphics"))
  attrs <- xmlAttrs(graphics)
  g <- new("KEGGGraphics",
           name=getNamedElement(attrs,"name"),
           x=as.integer(getNamedElement(attrs,"x")),
           y=as.integer(getNamedElement(attrs,"y")),
           type=getNamedElement(attrs,"type"),
           width=as.integer(getNamedElement(attrs, "width")),
           height=as.integer(getNamedElement(attrs,"height")),
           fgcolor=getNamedElement(attrs, "fgcolor"),
           bgcolor=getNamedElement(attrs, "bgcolor")
           )
  return(g)
  
}

parseEntry <- function(entry) {
  attrs <- xmlAttrs(entry)

  ## required: id, name,type
  entryID <- attrs[["id"]]
  name <- unname(unlist(strsplit(attrs["name"]," ")))
  type <- attrs[["type"]]

  ## implied: link, reaction, map
  link <- getNamedElement(attrs,"link")
  reaction <- getNamedElement(attrs, "reaction")
  map <- getNamedElement(attrs, "map")

  ## graphics
  graphics <- xmlChildren(entry)$graphics
  g <- parseGraphics(graphics)
  
  ## types: ortholog, enzyme, gene, group, compound and map
  if(type != "group") {
    newNode <- new("KEGGNode",
                   entryID=entryID,
                   name=name,
                   type=type,
                   link=link,
                   reaction=reaction,
                   map=map,
                   graphics=g)
  } else if(type=="group") {
    children <- xmlChildren(entry)
    children <- children[names(children) == "component"]
    if(length(children)==0) {
      component <- as.character(NA)
    } else {
      component <- sapply(children, function(x) {
        if(xmlName(x) == "component") {
          return(xmlAttrs(x)["id"])
        } else {
          return(as.character(NA))
        }     
      })
    }
    component <- unname(unlist(component))
    newNode <- new("KEGGGroup",
                   component=component,
                   entryID=entryID,
                   name=name,
                   type=type,
                   link=link,
                   reaction=reaction,
                   map=map,
                   graphics=g
                   )
  }
  return(newNode)
}

parseSubType <- function(subtype) {
  attrs <- xmlAttrs(subtype)
  name <- attrs[["name"]]
  value <- attrs[["value"]]
  return(new("KEGGEdgeSubType",name=name, value=value))
}
parseRelation <- function(relation) {
  attrs <- xmlAttrs(relation)

  ## required: entry1, entry2, type
  entry1 <- attrs[["entry1"]]
  entry2 <- attrs[["entry2"]]
  type <- attrs[["type"]]

  subtypeNodes <- xmlChildren(relation)
  subtypes <- sapply(subtypeNodes, parseSubType)
  newEdge <- new("KEGGEdge",
                 entry1ID=entry1,
                 entry2ID=entry2,
                 type=type,
                 subtype=subtypes
                 )                     
  return(newEdge)
}

.xmlChildrenWarningFree <- function(xmlNode) {
    if(is.null(xmlNode$children))
        return(NULL)
    return(XML::xmlChildren(xmlNode))
}

parseReaction <- function(reaction) {
  attrs <- xmlAttrs(reaction)

  ## required: name,type
  name <- attrs[["name"]]
  type <- attrs[["type"]]

  children <- xmlChildren(reaction)

  ## more than one substrate/product possible
  childrenNames <- names(children)
  substrateIndices <- grep("^substrate$", childrenNames)
  productIndices <- grep("^product$", childrenNames)
  substrateName <- substrateAltName <- vector("character", length(substrateIndices))
  productName <- productAltName <- vector("character", length(productIndices))  
  
  for (i in seq(along=substrateIndices)) {
    ind <- substrateIndices[i]
    substrate <- children[[ind]]
    substrateName[i] <- xmlAttrs(substrate)[["name"]]
    substrateAltName[i] <- as.character(NA)
    
    substrateChildren <- .xmlChildrenWarningFree(substrate)
    if (!is.null(substrateChildren)) {
        substrateAlt <- substrateChildren$alt
        substrateAltName[i] <- xmlAttrs(substrateAlt)[["name"]]
    }

  }

  for(i in seq(along=productIndices)) {
    ind <- productIndices[i]
    product <- children[[ind]]
    productName[i] <- xmlAttrs(product)[["name"]]
    productChildren <- .xmlChildrenWarningFree(product)
    productAltName[i] <- as.character(NA)
    if(!is.null(productChildren)) {
      productAlt <- productChildren$alt
      productAltName[i] <- xmlAttrs(productAlt)[["name"]]
    }
  }

  new("KEGGReaction",
      name = name,
      type = type,
      substrateName = substrateName,
      substrateAltName = substrateAltName,
      productName = productName,
      productAltName = productAltName)
}
