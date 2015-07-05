library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(snpStats)
library(doMC)
registerDoMC(4)

## function to convert snpid to target region
targetRegion <- function(target.snp, up = 200000, down = 200000){
  target.region <- target.snp
  start(target.region) <- start(target.snp) - abs(up)
  end(target.region) <- end(target.snp) + abs(down)
  target.region
}

## function to download the nessessary files for analysis
vcfxtract <- function(target.region, population){
  chr <- gsub("chr","",as.character(seqnames(target.region)))
  start <- start(target.region)
  end <- end(target.region)
  region <- paste(paste(chr,start,sep=":"),end, sep = "-")
  panel.file <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
  vcf.file <- sprintf("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", chr) 
  sprintf("./vcf_to_ped_convert.pl -vcf %s -sample_panel %s -region %s -population %s -output_ped %s -output_info %s", 
          vcf.file, panel.file, region, population, 
          sprintf("%s_%s.ped", region, population),
          sprintf("%s_%s.info", region, population))
}


## function to generate decent gene models from txdb
txdbCollapse <- function(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, tr){
  require(biovizBase)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(Homo.sapiens)
  tx <- crunch(TxDb.Hsapiens.UCSC.hg38.knownGene,which = tr)
  print(tx)
  tx <- tx[!nchar(as.character(tx$gene_id)) == 0]
  print(tx)
  if(length(tx) >0){
    tx <- unlist(lapply(split(tx,tx$gene_id), function(x) unlist(reduce(split(x,x$type)))))
    tx[nchar(names(tx)) == 0] <- NULL
    tx.names <- select(Homo.sapiens,keys = names(tx), columns = "SYMBOL", keytype = "GENEID")$SYMBOL
    names(tx) <- tx.names
    lapply(seq_along(tx),function(x){
      mcols(tx[[x]]) <<- data.frame("type" = names(tx[[x]]))
    })
    return(ggplot() + geom_alignment(GRangesList(tx)))
  }else{
    return(NULL)
  }
}


ld.analysis <- function(data, snp){
  depth <- ncol(data$genotypes@.Data)
#   ld.stats <- with(data, ld(genotypes[,snp], genotypes,
#                             stats=c("D.prime", "R.squared"), depth=depth))
  ld.stats <- with(data, ld(genotypes[,which.min(start(snp)-map$V2)[1]], genotypes,
                            stats=c("D.prime", "R.squared"), depth=depth))
  ld.stats <- with(ld.stats, data.frame(t(D.prime), 
                                        t(R.squared),
                                        data$map$V2, 
                                        data$map$snp.names, 
                                        stringsAsFactors=F,
                                        row.names = NULL))
  colnames(ld.stats) <- c("D.prime", "R.squared","position","snp.names")
  #ld.stats <- ld.stats[!ld.stats[,"snp.names"] %in% snp,]
  ld.stats
}


p63fun <- function(target.region, p63){
  target.p63 <- subsetByOverlaps(p63, target.region)
  if(length(target.p63) == 0){
    return(NULL)
  }else{
    return(target.p63)
  }
}


analysis <- function(snp,populations){
  require(ggplot2)
  require(ggbio)
  require(plyr)
  heights <- c(1,1,1,1)
  target.snp <- snpid2grange(SNPlocs.Hsapiens.dbSNP141.GRCh38, snp)
  target.snp <- GRanges(Rle(gsub("ch","chr",as.character(seqnames(target.snp)))),
                        IRanges(start = start(target.snp),end = end(target.snp)))
  target.region <- targetRegion(target.snp)
  pop <- llply(populations, function(p){
    response <- system(vcfxtract(target.region,p),intern = T)
    data.files <- regmatches(response,gregexpr("\\d+:\\d+-\\d+_\\w{3}\\.{1}(ped|info)",response,perl = T))[[1]]
    dat <- read.pedfile(data.files[2], snps = data.files[1])
    pop <- ld.analysis(dat,target.snp)
    mycols <- pop$R.squared >= 0.75
    mycols[which(mycols == T)] <- "red"
    mycols[which(mycols == F)] <- "black"
    pop$mycols <- mycols
    ggplot(pop) + ylim(c(0,1)) + geom_point(aes(position,round(D.prime,digits = 2),size = R.squared,fill = mycols,colour = mycols), alpha = 0.2) + scale_color_manual(values = c("black","black","red")) + geom_vline(xintercept = start(target.snp)) + geom_text(x=start(target.snp),y=0.5, label=snp)  },.parallel = T)
  names(pop) <- populations
  txs <- txdbCollapse(TxDb.Hsapiens.UCSC.hg38.knownGene,
                      tr = target.region)
  pop[["p63"]] <- p63fun(target.region,reduce(p63_hg38))
  if(!is.null(pop[["p63"]])){
    heights <- c(heights,0.2)
  }
  pop[["txs"]] <- txs
  if(!is.null(pop[["txs"]])){
    heights <- c(heights,0.75)
  }
  return(list(pop,heights))
}
