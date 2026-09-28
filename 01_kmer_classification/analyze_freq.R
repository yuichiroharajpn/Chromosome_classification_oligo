#library(maptools)
library(ape)
library(phangorn)
library(dplyr)
library(tibble)
require(ca)

args <- commandArgs(trailingOnly = T)
file <- args[1]
flist <- "chromosome_list.txt"

outf1 <- sub(".txt", ".dendro.pdf", file)
outf2 <- sub(".txt", ".ca.pdf", file)
outf3 <- sub(".txt", ".dim.txt", file)
outf4 <- sub(".txt", ".var.txt", file)
outf5 <- sub(".txt", ".dendro.grp.txt", file)
outf51 <- sub(".txt", ".dendro.grp.h.txt", file)

df <- read.delim(file, row.names=1)
df_l <- read.delim(flist, head=F, row.names=1)
df_g <- data.frame(len=df_l[,2], row.names = rownames(df_l))
df_g$Chr <- rownames(df_g)
df_h <- data.frame(len=df_l[,2], row.names = rownames(df_l))
df_h$Chr <- rownames(df_h)

xd <- dist(t(df), method="euclidean")
clus1 <- hclust(xd, method="ward.D2")
pdf(outf1,width=10,height=8)
plot(clus1)
dev.off()

tree <- as.phylo(clus1)
tree <- makeNodeLabel(tree, method = "number", prefix = "N")

df_g0 <- data.frame(matrix(rep(0, nrow(df_g)*5), ncol=5))
colnames(df_g0) <- c("L1", "L2", "L3", "L4", "L5")
rownames(df_g0) <- tree$tip.label
C <- c(0, 0, 0, 0, 0)

root <- ncol(df)+1
a1 <- Descendants(tree, root, type="children")
for( i1 in a1){
	C[1] <- C[1]+1
	a10 <- Descendants(tree, i1, type="tips")
	#df_g0[df_g0$L1 %in% a10[[1]],"L1"] <- C[1]
	df_g0[a10[[1]],"L1"] <- C[1]

	a2 <- Descendants(tree, i1, type="children")
	for( i2 in a2){
		C[2] <- C[2]+1
		a20 <- Descendants(tree, i2, type="tips")
		df_g0[a20[[1]],"L2"] <- C[2]

		a3 <- Descendants(tree, i2, type="children")
		for( i3 in a3){
			C[3] <- C[3]+1
			a30 <- Descendants(tree, i3, type="tips")
			df_g0[a30[[1]],"L3"] <- C[3]

			a4 <- Descendants(tree, i3, type="children")
			for( i4 in a4){
				C[4] <- C[4]+1
				a40 <- Descendants(tree, i4, type="tips")
				df_g0[a40[[1]],"L4"] <- C[4]
				#cat(C[4] ,"|", a40[[1]], "\n")

				a5 <- Descendants(tree, i4, type="children")
				for( i5 in a5){
					C[5] <- C[5]+1
					a50 <- Descendants(tree, i5, type="tips")
					df_g0[a50[[1]],"L5"] <- C[5]
				}
			}
		}
	}
}

for(i in which(df_g0$L5==0)){
	C[5] <- C[5]+1
	df_g0$L5[i] <- C[5]
}

for(i in which(df_g0$L4==0)){
	C[4] <- C[4]+1
	df_g0$L4[i] <- C[4]
}



df_g <- df_g %>% left_join(df_g0 %>% rownames_to_column("Chr"), by="Chr")

df_g$Length <- df_g$len
df_g <- df_g[,-1]

for( i in 2:7 ){
	df_h0 <-data.frame(cutree(clus1,k=i))
	colnames(df_h0) <- c(paste0("K",i))
	df_h <- cbind(df_h,df_h0)
}
df_h$Length <- df_h$len
df_h <- df_h[,-1]


ca1 = ca(t(df))

pdf(outf2,width=8,height=8)
plot.ca(ca1, what=c("all","none"), col=c("#ff8c00","black"), col.lab=c("black","black"), bty="n")
dev.off()

#write.table(df_g,"",sep="\t",quote=F,row.names=F)
#write.table(df_h,"",sep="\t",quote=F,row.names=F)
write.table(df_g,file=outf51,sep="\t",quote=F,row.names=F)
write.table(df_h,file=outf5,sep="\t",quote=F,row.names=F)
write.table(ca1$rowcoord,file=outf3,sep="\t",quote=F,row.names=F)
write.table(summary(ca1)$scree,file=outf4,sep="\t",quote=F,row.names=F)
