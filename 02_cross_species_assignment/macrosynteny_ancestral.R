# Ancestral linkage groups on sp1; target chromosomes on sp2.
library(macrosyntR)
library(dplyr)
args <- commandArgs(trailingOnly = T)

one2one <- args[1]
sp1 <- args[2]
sp2 <- args[3]
anc <- args[4]

df_ortho <- read.delim(one2one, head=F)
df_anc <- read.delim(anc)
df_anc$sp1.ID <- df_anc$Chicken.gene
df_anc <- df_anc %>% select(anc_jv, sp1.ID)
df_sp1 <- read.delim(sp1, head=F)
df_sp2 <- read.delim(sp2, head=F)
colnames(df_ortho) <- c("sp1", "sp2")
colnames(df_sp1) <- c("chr", "st", "ed", "pid")
colnames(df_sp2) <- c("chr", "st", "ed", "pid")

df_sp1 <- df_sp1 %>% inner_join(df_anc, by = c("pid" = "sp1.ID"))

df2 <- df_sp2 %>% group_by(chr) %>% summarize(count = n())
df1 <- df_sp1 %>% group_by(anc_jv) %>% summarize(count = n())
df2 <- df2 %>% left_join(df_sp2 %>% filter(pid %in% df_ortho$sp2) %>% group_by(chr) %>% summarize(ocount = n()), by="chr") %>% replace(is.na(.), 0)
df1 <- df1 %>% left_join(df_sp1 %>% filter(pid %in% df_ortho$sp1) %>% group_by(anc_jv) %>% summarize(ocount = n()), by="anc_jv") %>% replace(is.na(.), 0)

my_orthologs <- load_orthologs(orthologs_table = one2one, sp1_bed = sp1, sp2_bed = sp2)

my_orthologs2 <- my_orthologs %>% inner_join(df_anc, by="sp1.ID")
my_orthologs2$sp1.Chr <- as.factor(my_orthologs2$anc_jv)
my_orthologs2 <- my_orthologs2 %>% select(1:10)

my_macrosynteny <- compute_macrosynteny(my_orthologs2)
res <- my_macrosynteny[my_macrosynteny$significant=="yes",]

ocount1<-as.numeric(df1$ocount)
names(ocount1) <- df1$anc_jv
ocount2<-as.numeric(df2$ocount)
names(ocount2) <- df2$chr
res$sp1.Chr.ratio <- format(res$orthologs/ocount1[as.character(res$sp1.Chr)], nsmall = 2)
res$sp2.Chr.ratio <- format(res$orthologs/ocount2[as.character(res$sp2.Chr)], nsmall = 2)

write.table(res,file="",sep="\t",quote=F,row.names=F)
