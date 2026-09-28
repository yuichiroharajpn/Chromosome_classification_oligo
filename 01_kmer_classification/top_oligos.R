library(dplyr)
library(readr)
library(stringr)

# 1. 引数の取得
args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
out1 <- sub(".txt", ".top100.Dim1.txt", input_file)
out2 <- sub(".txt", ".top100.Dim1.plus_revcomp.txt", input_file)


# 2. データの読み込み
df <- read_table(input_file, show_col_types = FALSE)

# 3. オリジナルのTop 100を抽出して保存
top_df <- df %>%
  arrange(desc(Dim1)) %>%
  slice_head(n = 100)

top_oligos <- top_df$Oligo
write_lines(top_oligos, out1)

# 4. 相補・逆相補配列の生成
# chartr を使うのが R で最も確実で高速な塩基置換方法です
get_rev_comp <- function(seqs) {
  # 相補鎖(Complement)の作成: A->T, T->A, C->G, G->C
  comp <- chartr("ATCGatcg", "TAGCtagc", seqs)

  # 逆転(Reverse)させて「逆相補鎖」にする
  rev_comp <- sapply(strsplit(comp, ""), function(x) paste(rev(x), collapse = ""))
  return(rev_comp)
}

rev_comp_oligos <- get_rev_comp(top_oligos)

# 5. マージして重複削除
# オリジナル(CGTTTT等)と逆相補(AAAACG等)を合わせる
final_list <- unique(c(top_oligos, rev_comp_oligos))

write_lines(final_list, out2)

# 検証用出力
if("CGTTTT" %in% top_oligos) {
  message("Check: CGTTTT found in top100")
  message(paste("Check: Reverse Complement is", get_rev_comp("CGTTTT")))
}
