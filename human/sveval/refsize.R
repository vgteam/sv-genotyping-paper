library(sveval)
library(dplyr)
args = commandArgs(TRUE)

df = lapply(args[-1], function(vcf.file){
  message(vcf.file)
  sv = readSVvcf(vcf.file)
  data.frame(vcf=vcf.file, size=sv$size, type=sv$type, stringsAsFactors=TRUE)
})
df = do.call(rbind, df)

df = df %>% mutate(
              vcf=basename(vcf),
              size=cut(size, breaks=c(50,100,200,300,400,600,800,
                                      1000,2500,5000,Inf))) %>%
  group_by(vcf, size, type) %>% summarize(n=n())

write.table(df, file=args[1], quote=FALSE, sep='\t', row.names=FALSE)
