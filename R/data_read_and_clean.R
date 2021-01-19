#' @title Read, clean and annotate MS TMT data.
#' @description Read, clean and annotate MS TMT data with Maxquant process data and meta file.
#' @details Input MaxQuang processed MS TMT data and meta file, then return a data set with annotation and TMT quantification results.
#' @param filename A string of processed file path.
#' @param metafile A string of metafile path.
#' @param total A boolean of considering zero row or not.
#' @return A DataFrame.
#' @export


read_maxquant_prot <- function(filename,metafile,total=FALSE){
  prot_raw <- read.table(filename,header = TRUE,sep = '\t')
  prot_dat_clean <- data_clean(prot_raw)
  prot_dat_extract <- data_extract(prot_dat_clean,metafile,total)
  prot_annotation <- data_annotation(prot_dat_extract)
  rownames(prot_annotation) <- prot_annotation$id
  return(prot_annotation)
}

data_clean <- function(prot_raw){
  temp <- (prot_raw$Potential.contaminant!='+') & (prot_raw$Reverse!='+') & (prot_raw$Only.identified.by.site!='+')
  prot_dat_clean <- prot_raw[temp,]
  return(prot_dat_clean)
}

data_extract <- function(prot_dat,metafile,total=FALSE){
  meta <- read.table(metafile,header = TRUE, sep = '\t')

  if (total == FALSE){
    temp <- TRUE
    for (t in meta$channel[!is.na(meta$reference)]) {
      temp <- temp & (prot_dat[,t] != 0)
    }
    prot_dat <- prot_dat[temp,]
  }

  prot_dat_tmt <- prot_dat[,meta$channel]
  colnames(prot_dat_tmt) <- meta$sample
  target <- meta$sample[is.na(meta$reference)]
  prot_dat_tmt <- prot_dat_tmt[,target]

  prot_dat_extract <- cbind(prot_dat_tmt,prot_dat[,c('id','Fasta.headers')])

  return(prot_dat_extract)
}

data_annotation <- function(prot_dat_extract){
  prot_annotation <- prot_dat_extract[,c('id','Fasta.headers')]
  prot_annotation$ACCID <- apply(prot_annotation,1,function(x){unlist(strsplit(x[2],'\\|'))[2]})
  prot_annotation$PROTEIN <- apply(prot_annotation,1,function(x){unlist(strsplit(x[2],'\\|'))[3]})
  prot_annotation$PROTEIN <- apply(prot_annotation,1,function(x){unlist(strsplit(x[4],'_HUMAN'))[1]})
  prot_annotation$FULLNAME <- apply(prot_annotation,1,function(x){unlist(strsplit(x[2],'_HUMAN'))[2]})
  prot_annotation$FULLNAME <- apply(prot_annotation,1,function(x){unlist(strsplit(x[5],'OS=Homo'))[1]})

  prot_annotation <- cbind(prot_dat_extract,prot_annotation[,3:5])
  prot_annotation <- prot_annotation[,!(colnames(prot_annotation)=='Fasta.headers')]

  return(prot_annotation)
}

