# gtfAttributeParser
# author : Frederic JEHL
# 
# Improved version of the attributeParser function
#
# INPUT: a dataframe corresponding to a GTF file, with at least the ninth column named "attribute"
# OUTPUT: the attribute column in the form of a dataframe, with as many rows as the GTF, and as many columns as the total number of uniques fields found
# in the attribute colum. In lines where a field was missing, its value is set to NA


gtfAttributeParser <- function(gtfFile){
  attribute <- gtfFile$attribute
  cat("Treating the GTF file composed of", nrow(gtfFile), "lines ... ")
  attribute <- gsub(" \\(.*\\)", "_bracketsRemoved", attribute)
  list.splitAttribute <- strsplit(attribute, "; ")
  list.splitAttribute <- lapply(list.splitAttribute, function(x){gsub(";", "", x)})
  namesOfParsedAttributes <- unlist(list.splitAttribute)
  namesOfParsedAttributes <- gsub("\\s.*", "", namesOfParsedAttributes)
  namesOfParsedAttributes <- unique(namesOfParsedAttributes)
  
  splitAttribute <- matrix(nrow = length(list.splitAttribute), ncol = length(namesOfParsedAttributes))
  colnames(splitAttribute) <- namesOfParsedAttributes  
  list.splitAttribute <- lapply(list.splitAttribute, function(x){unlist(strsplit(x, " "))})
  totalSize <- length(list.splitAttribute)
  by_percent <- 10
  slice_size <- totalSize / by_percent
  
  for (i in 1:length(list.splitAttribute)){
    
    length_list.splitAttribute_i <- length(list.splitAttribute[[i]])
    
    names(list.splitAttribute[[i]]) <- list.splitAttribute[[i]][rep(seq(from = 1, to = length_list.splitAttribute_i, by = 2), each = 2)]
    list.splitAttribute[[i]] <- list.splitAttribute[[i]][- rep(seq(from = 1, to = length_list.splitAttribute_i, by = 2), each = 2)]
    
    splitAttribute[i, names(list.splitAttribute[[i]])] <- unname(list.splitAttribute[[i]])
  }
  
  cat("DONE\n")
  
  return(splitAttribute)
}
