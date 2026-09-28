library(dplyr)

args <- commandArgs(trailingOnly = T)
k <- args[1]
type1 <- args[2]
type2 <- args[3]
flist <- "chromosome_list.txt"

df<-read.delim(flist,head=F)

j<-0
for(i in df$V1){
	chr <- i
	file <- paste0("chromosomes/",chr,"_",type1,".txt")
	df0 <- read.delim(file, head=F)
	df0$R <- df0$V2/sum(df0$V2)

	if (j == 0){
		df_c <- df0[,c(1,2)]
		df_r <- df0[,c(1,3)]
		colnames(df_c) <- c("Kmer",chr)
		colnames(df_r) <- c("Kmer",chr)
	}else{
		df1 <- df0[,c(1,2)]
		df2 <- df0[,c(1,3)]
		colnames(df1) <- c("Kmer",chr)
		colnames(df2) <- c("Kmer",chr)
		df_c <- df_c %>% full_join(df1 ,by="Kmer")
		df_r <- df_r %>% full_join(df2 ,by="Kmer")
	}
	j <- j+1
}
df_c[is.na(df_c)] <- 0
df_r[is.na(df_r)] <- 0

write.table(df_c,file=paste0(type2,"_freq.txt"),sep="\t",quote=F,row.names=F)
write.table(df_r,file=paste0(type2,"_relfreq.txt"),sep="\t",quote=F,row.names=F)
