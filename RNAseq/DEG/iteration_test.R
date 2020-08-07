src_dir <- c("/Users/hmkim/data/quant_data/LJK")
setwd("/Users/hmkim/data/quant_data/LJK")
src_files <- list.files(src_dir)
src_files

for (i in src_files){
  print(i)
}
