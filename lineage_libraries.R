require(ShortRead)
require(fastmatch)
require(Biostrings)
require(dplyr)
require(plyr)


#####################
# SCRIPT functions  #
#####################

# Rcpp progress bar style
progress_for <- function(n, tot) {
  message("0%   10   20   30   40   50   60   70   80   90   100%")
  message("[----|----|----|----|----|----|----|----|----|----|")
  # n:tot = nstars:50 -> nstars = (n*50)/tot
  nstars = floor((n*50)/tot)
  if(nstars>0)
    for (i in 1:nstars) {
      message("*", appendLF = FALSE)
      utils::flush.console()
    }
  message("|")
}

stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
tsmessage <- function(..., domain = NULL, appendLF = TRUE, verbose = TRUE,time_stamp = TRUE) {
  if (verbose) {
    msg <- ""
    if (time_stamp) {
      msg <- paste0(stime(), " ")
    }
    message(msg, ..., domain = domain, appendLF = appendLF)
    utils::flush.console()
  }
}

encode_seqs = function(seqs){
  dictionary = 0:4
  names(dictionary) = c("N","A","C","G","T")
  if(!is.character(seqs)) seqs = t(sapply(seqs, function(x) strsplit(as.character(x), "")[[1]]))
  else seqs = t(sapply(seqs, function(x) strsplit(x, "")[[1]]))
  enc_seq = matrix(data = dictionary[seqs], nrow = nrow(seqs),ncol = ncol(seqs),byrow = F)
  return(enc_seq)
}

correct_bc = function(seqs, bctype,bc14,bc30,threads=6){
  seqs2corr = encode_seqs(seqs[,3])
  if(bctype=='bc14') {
    seqsRef = encode_seqs(bc14$BC14_SEQUENCE)
    th = 5
    names = bc14$BC14_SEQUENCE
  } else {
    seqsRef = encode_seqs(bc30$BC30_SEQUENCE)
    th = 9
    names = bc30$BC30_SEQUENCE
  }
  ## Compute distance matrix
  seqs$id = paste0(seqs$group,":",seqs$start)
  rownames(seqs) = seqs$id
  hd = hammingDistAll(m = seqs2corr, m2 = seqsRef, ncores = threads, verbose = F, th=th)
  colnames(hd) = names
  rownames(hd) = seqs$id
  
  ## Search the lowest and unique value
  idx = apply(hd, 1, function(x) {ifelse(sum(x==sort(x,method="quick")[1])==1, which.min(x), 0)})
  idx = idx[!idx==0]
  corrseqs =  colnames(hd)[idx]
  if (bctype=="bc14"){
      corrseqs = sapply(corrseqs, function(x,y=bc14) y$BC14_ID[fastmatch::fmatch(x,y$BC14_SEQUENCE)], USE.NAMES = F); rm(hd)
    } else {
      corrseqs = sapply(corrseqs, function(x,y=bc30) y$BC30_ID[fastmatch::fmatch(x,y$BC30_SEQUENCE)], USE.NAMES = F); rm(hd)
    }
  
  corrseqs = data.frame(group = seqs[names(idx), ]$group, start = seqs[names(idx), ]$start, id = corrseqs,row.names = names(idx),stringsAsFactors = F)
  
  return(corrseqs)
}


readBC_QSQ = function(fqin2,sample,max.mismatch,condition,outDir,fqchunks = (10^6)/2)
{
  strm2 <- FastqStreamer(fqin2,n = fqchunks)
  
  tsmessage("Estimating the number of reads..")
  total = countLines(fqin2)/4
  tsmessage(paste0(total," reads found"))
  tsmessage("Estimating Cell lineage barcode..")
  
  
  
  count = 0
  progress_for(count,total)
  res = NULL
  repeat {
    # read FASTQ chunks and exit if needed
    fq2 <- yield(strm2)
    count = count + length(fq2)
    if (length(fq2) == 0) {
      progress_for(100,100)
      break
    }
    # find reads in which the anchor TGGT occur 
    r = as.data.frame(Biostrings::vmatchPattern(pattern = "TGGT",
                                                subject = fq2@sread,
                                                max.mismatch = max.mismatch,
                                                min.mismatch = 0,
                                                with.indels = F))
    
    # take only the reads in which the anchor occur and can contain both the barcodes
    r = subset(r,start>=15 & end <= 71)
    rownames(r) = paste0(r$group,":",r$start)
    
    ## Get bc14 and bc30 sequences
    r$bc14 = Biostrings::substr(x = fq2@sread[r$group],start = r$start-14,stop = r$start-1)
    r$bc30 = Biostrings::substr(x = fq2@sread[r$group],start = r$end+1,stop = r$end+30)
    
    ## Search in whitelist
    r$bc14idx = bc14$BC14_ID[fastmatch::fmatch(r$bc14,bc14$BC14_SEQUENCE)]
    r$bc30idx = bc30$BC30_ID[fastmatch::fmatch(r$bc30,bc30$BC30_SEQUENCE)]
    
    ## Correct barcodes 
    bc14_corr = subset(r,is.na(bc14idx)& !is.na(bc30idx))[, c("group", "start", "bc14")] ## correct bc14
    if(nrow(bc14_corr)>0) {
      bc14_corr = correct_bc(bc14_corr, bctype='bc14',bc14 = bc14,bc30 = bc30,threads = 12)
      r[rownames(bc14_corr),"bc14idx"] = bc14_corr$id
    }
    
    bc30_corr = subset(r,!is.na(bc14idx)& is.na(bc30idx))[, c("group", "start", "bc30")] ## correct bc30
    if(nrow(bc30_corr)>0) {
      bc30_corr = correct_bc(bc30_corr, bctype='bc30',bc14 = bc14,bc30 = bc30,threads = 12)
      r[rownames(bc30_corr),"bc30idx"] = bc30_corr$id
    }
    
    ## Retain reads in which there are both bc14 and bc30
    r = subset(r,!is.na(bc14idx) & !is.na(bc30idx))
    
    ## Add read 1 information (i.e. BC cell and UMI)
    res = rbind(res,r[,c("bc14idx","bc30idx")])
    progress_for(n = count,total)
  }
  
  tsmessage("Estimating Cell lineage barcode.. FINISHED!")
  
  # close connections
  close(strm2)
  
  res$linageBC = paste0(res$bc14,":",res$bc30)

  df = ddply(res,c("linageBC"),summarise,n=length(linageBC))
  df = df[order(df$n,decreasing = T),]
  df$freq = df$n/sum(df$n)
  df$sample = sample
  df$condition = condition
  saveRDS(df,file = outDir)
  return(df)
}


readBC_scRNA = function(fqin,fqin2,sample,max.mismatch,condition,outDir,fqchunks = (10^6)/2)
{
  strm1 <- FastqStreamer(fqin,n = fqchunks)
  strm2 <- FastqStreamer(fqin2,n = fqchunks)
  
  tsmessage("Estimating the number of reads..")
  total = countLines(fqin)/4
  tsmessage(paste0(total," reads found"))
  tsmessage("Estimating Cell lineage barcode..")
  
  
  count = 0
  progress_for(count,total)
  res = NULL
  repeat {
    # read FASTQ chunks and exit if needed
    fq <- yield(strm1)
    fq2 <- yield(strm2)
    count = count + length(fq)
    if (length(fq) == 0) {
      progress_for(100,100)
      break
    }
    # find reads in which the anchor TGGT occur 
    r = as.data.frame(Biostrings::vmatchPattern(pattern = "TGGT",
                                                subject = fq2@sread,
                                                max.mismatch = max.mismatch,
                                                min.mismatch = 0,
                                                with.indels = F))
    
    # take only the reads in which the anchor occur and can contain both the barcodes
    r = subset(r,start>=15 & end <= 54)
    rownames(r) = paste0(r$group,":",r$start)
    
    ## Get bc14 and bc30 sequences
    r$bc14 = Biostrings::substr(x = fq2@sread[r$group],start = r$start-14,stop = r$start-1)
    r$bc30 = Biostrings::substr(x = fq2@sread[r$group],start = r$end+1,stop = r$end+30)
    
    ## Search in whitelist
    r$bc14idx = bc14$BC14_ID[fastmatch::fmatch(r$bc14,bc14$BC14_SEQUENCE)]
    r$bc30idx = bc30$BC30_ID[fastmatch::fmatch(r$bc30,bc30$BC30_SEQUENCE)]
    
    ## Correct barcodes 
    bc14_corr = subset(r,is.na(bc14idx)& !is.na(bc30idx))[, c("group", "start", "bc14")] ## correct bc14
    if(nrow(bc14_corr)>0) {
      bc14_corr = correct_bc(bc14_corr, bctype='bc14',bc14 = bc14,bc30 = bc30,threads = 12)
      r[rownames(bc14_corr),"bc14idx"] = bc14_corr$id
    }
    
    bc30_corr = subset(r,!is.na(bc14idx)& is.na(bc30idx))[, c("group", "start", "bc30")] ## correct bc30
    if(nrow(bc30_corr)>0) {
      bc30_corr = correct_bc(bc30_corr, bctype='bc30',bc14 = bc14,bc30 = bc30,threads = 12)
      r[rownames(bc30_corr),"bc30idx"] = bc30_corr$id
    }
    
    ## Retain reads in which there are both bc14 and bc30
    r = subset(r,!is.na(bc14idx) & !is.na(bc30idx))
    
    ## Add read 1 information (i.e. BC cell and UMI)
    r$cell = Biostrings::substr(x = fq@sread[r$group],start = 1,stop = 16)
    r$umi = Biostrings::substr(x = fq@sread[r$group],start = 17, stop = 28)
    res = rbind(res,r[,c("cell","umi","bc14idx","bc30idx")])
    progress_for(n = count,total)
  }
  
  tsmessage("Estimating Cell lineage barcode.. FINISHED!")
  
  # close connections
  close(strm1)
  close(strm2)
  
  res$linageBC = paste0(res$bc14,":",res$bc30)
  res = unique(res) # use cell-UMI trick to remove duplicated reads
  
  
  df = ddply(res,c("cell","linageBC"),summarise,n=length(linageBC))
  df = df[order(df$n,decreasing = T),]
  df$sample = sample
  df$condition = condition

  saveRDS(df,file = outDir)
  return(df)
}

