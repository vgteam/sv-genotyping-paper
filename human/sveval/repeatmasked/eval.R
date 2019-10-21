library(sveval)
library(VariantAnnotation)
library(ggplot2)

truth.vcf = readSVvcf('hgsvc-truth-baseline.norm.rmsk.vcf.gz', vcf.object=TRUE)
calls.vcf = readSVvcf('hgsvc-vg-HG00514.norm.rmsk.vcf.gz', vcf.object=TRUE)

min.cov = .8
truth.reps = ifelse(info(truth.vcf)$RMSKCOV>min.cov, info(truth.vcf)$RMSKCLASS, NA)
truth.ids = paste(as.character(seqnames(rowRanges(truth.vcf))),
                  start(rowRanges(truth.vcf)),
                  info(truth.vcf)$SVTYPE, info(truth.vcf)$SIZE)

calls.reps = ifelse(info(calls.vcf)$RMSKCOV>min.cov, info(calls.vcf)$RMSKCLASS, NA)
calls.ids = paste(as.character(seqnames(rowRanges(calls.vcf))),
                  start(rowRanges(calls.vcf)),
                  info(calls.vcf)$SVTYPE, info(calls.vcf)$SIZE)

load('../rdata/sveval-hgsvc-vg-HG00514-all-geno.RData')
svs.geno = eval.o$svs
load('../rdata/sveval-hgsvc-vg-HG00514-all-call.RData')
svs.call = eval.o$svs
svs = list(geno=svs.geno, call=svs.call)

reps = table(calls.reps)
reps = names(reps)[reps>100]
eval.df = lapply(names(svs), function(eval){
  svs = svs[[eval]]
  eval.df = lapply(reps, function(repc){
    ## For each SV type
    eval.df = lapply(names(svs), function(svtype){
      svs = svs[[svtype]]
      ## For each class of variant
      df = lapply(c('TP', 'TP.baseline', 'FP', 'FN'), function(metric){
        svs = svs[[metric]]
        sv.ids = paste(as.character(seqnames(svs)), start(svs), svs$type, svs$size)
        if(metric %in% c('TP.baseline', 'FN')){
          svs = svs[which(sv.ids %in% truth.ids[which(truth.reps==repc)])]
        } else {
          svs = svs[which(sv.ids %in% calls.ids[which(calls.reps==repc)])]
        }
        data.frame(rep=repc, type=svtype, metric=metric, n=length(svs), eval=eval,
                   stringsAsFactors=FALSE)
      })
      do.call(rbind, df)
    })
    eval.df = do.call(rbind, eval.df)
  })
  eval.df = do.call(rbind, eval.df)
})
eval.df = do.call(rbind, eval.df)

## Reformat into one row per size class/type with columns TP, FP, etc
eval.df = tidyr::spread(eval.df, 'metric', 'n', fill=0)

## Precision, recall and F1
eval.df = prf(eval.df)

write.table(eval.df, file='eval-rmsk-hgsvc-vg-HG00514-call-geno.tsv', sep='\t', quote=FALSE, row.names=FALSE)
