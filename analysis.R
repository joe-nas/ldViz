library(snpStats)
library(ggplot2)
library(ggbio)
library(rtracklayer)

p63_1 <- with(read.table("GSM439776_p63_1_peaks.bed.gz",header = F),GRanges(Rle(V1), IRanges(start=V2,end=V3),score = V4))
p63_2 <- with(read.table("GSM439777_p63_2_peaks.bed.gz",header = F),GRanges(Rle(V1), IRanges(start=V2,end=V3),score = V4))
p63_3 <- with(read.table("GSM538930_p63_3_peaks.bed.gz",header = F),GRanges(Rle(V1), IRanges(start=V2,end=V3),score = V4))
p63<- c(p63_1,p63_2,p63_3)

ch <- import.chain(ChainFile("hg18ToHg38.over.chain"))
p63_hg38 <- unlist(liftOver(x = p63,chain = ch))


snps <- c("rs642961","rs987525","rs17085106")

snps<-c("rs2166975","rs3771494","rs3821261","rs426081","rs3771475","rs748044","rs6755375","rs277547","rs210910","rs729287","rs7940667","rs1546124","rs4783099","rs16974880","rs2326398","rs8061351","rs1466094","rs1055333","rs686606","rs1962113")

test <- analysis(snps[2],c("CEU","YRI","JPT","CHB"))
tracks(test[[1]],heights = test[[2]])





