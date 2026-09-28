#library(maptools)
require(ca)
library(tibble)

args <- commandArgs(trailingOnly = T)
file <- args[1]
flist <- "chromosome_list.txt"

outf1 <- sub(".txt", ".coord.txt", file)
outf2 <- sub(".txt", ".contr.txt", file)

df <- read.delim(file, row.names=1)
ca1 = ca(t(df))

df1 <-data.frame(ca1$colcoord[,1:10])
df2 <- df1[order(df1$Dim1, decreasing=T),] %>% rownames_to_column("Oligo")

df3 <- data.frame(
  Dim = seq_along(ca1$sv),
  Eigenvalue = ca1$sv^2,
  Prop = (ca1$sv^2)/sum(ca1$sv^2)*100,
  Cum = cumsum((ca1$sv^2)/sum(ca1$sv^2)*100),
  check.names = FALSE
)



write.table(df2,file=outf1,sep="\t",quote=F,row.names=F)
write.table(df3,file=outf2,sep="\t",quote=F,row.names=F)
