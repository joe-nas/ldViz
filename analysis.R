library(snpStats)
library(ggplot2)
library(ggbio)
library(rtracklayer)
library(plyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

p63_1 <- with(read.table("GSM439776_p63_1_peaks.bed.gz",header = F),GRanges(Rle(V1), IRanges(start=V2,end=V3),score = V4))
p63_2 <- with(read.table("GSM439777_p63_2_peaks.bed.gz",header = F),GRanges(Rle(V1), IRanges(start=V2,end=V3),score = V4))
p63_3 <- with(read.table("GSM538930_p63_3_peaks.bed.gz",header = F),GRanges(Rle(V1), IRanges(start=V2,end=V3),score = V4))
p63<- reduce(c(p63_1,p63_2,p63_3))

ch <- import.chain(ChainFile("hg18ToHg38.over.chain"))
p63_hg38 <- unlist(liftOver(x = p63,chain = ch))
nchar(getSeq(BSgenome.Hsapiens.UCSC.hg38, p63_hg38)[223])

names(p63_hg38) <- paste("seq",names(p63_hg38),sep="")
myseqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, p63_hg38)
file.remove("tmp.fa")
writeXStringSet(myseqs, "tmp.fa")
myDat <- ldply(system("p63scan.py -i tmp.fa", intern=T),
                 function(x) unlist(strsplit(x,split = "\t")))
myDat$V1<-myDat$V1
myDat$V4<-as.integer(myDat$V4)+1
myDat$V5<-as.integer(myDat$V5)
myDat$V9<-gsub(" ","",myDat$V9)


myres <- laply(seq_along(myseqs[myDat$V1]),function(i)grep(myDat$V9[i],myseqs[myDat$V1][i]),.parallel = T)

#subseq access
subseqfun <- function(gr_){
  seqs_ <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_)
  file.remove("tmp.fa")
  writeXStringSet(seqs_, "tmp.fa")
  scanDat <- ldply(system("p63scan.py -i tmp.fa", intern=T),
                   function(x) unlist(strsplit(x,split = "\t")))
  
  ids <- as.integer(scanDat$V1)
  starts <- as.integer(scanDat$V4) + 1
  ends <- as.integer(scanDat$V5)
  for(i in 1:length(ids)){
    print(ids[i])
    print(scanDat[i,])
    print(starts[i])
    print(subseq(seqs_[ids[i]],starts[i],ends[i]))
  }
  # list(seqs_,scanDat)
}
test <- subseqfun(p63_hg38)


conservedFun <- function(gr_,scan_,i){
  tmp_gr <- gr_[as.integer(scan_[,1])]
  new_starts <- start(tmp_gr) + as.numeric(scan_[,4]) + i
  tmp_ir <- IRanges(start = new_starts, width = 1)
  ranges(tmp_gr) <- tmp_ir
  tmp_gr
}

# conserved 3,6,13,16
disruptedBySnp <- ""


snps <- c("rs642961","rs987525","rs17085106")

snps<-c("rs2166975","rs3771494","rs3821261","rs426081","rs3771475","rs748044","rs6755375","rs277547","rs210910","rs729287","rs7940667","rs1546124","rs4783099","rs16974880","rs2326398","rs8061351","rs1466094","rs1055333","rs686606","rs1962113")

test <- analysis(snps[2],c("CEU","YRI","JPT","CHB"))
tracks(test[[1]],heights = test[[2]])





