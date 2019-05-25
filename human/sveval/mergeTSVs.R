## Merge persize files
tsvs = list.files('tsv', 'persize.tsv')
tsvs = grep('merged', tsvs, value=TRUE, invert=TRUE)
df = lapply(tsvs, function(ff){
  df = read.table(paste0('tsv/', ff), as.is=TRUE, header=TRUE)
  df$exp = gsub('(.*)-.*-.*-.*-.*-persize.tsv', '\\1', ff)
  df$method = gsub('.*-(.*)-.*-.*-.*-persize.tsv', '\\1', ff)
  df$sample = gsub('.*-.*-(.*)-.*-.*-persize.tsv', '\\1', ff)
  df$region = gsub('.*-.*-.*-(.*)-.*-persize.tsv', '\\1', ff)
  df$eval = gsub('.*-.*-.*-.*-(.*)-persize.tsv', '\\1', ff)
  df
})
df = do.call(rbind, df)
write.table(df, file='tsv/human-merged-persize.tsv', quote=FALSE, sep='\t', row.names=FALSE)

## Merge PR curve files
tsvs = list.files('tsv', 'prcurve.tsv')
tsvs = grep('merged', tsvs, value=TRUE, invert=TRUE)
df = lapply(tsvs, function(ff){
  df = read.table(paste0('tsv/', ff), as.is=TRUE, header=TRUE)
  df$exp = gsub('(.*)-.*-.*-.*-.*-prcurve.tsv', '\\1', ff)
  df$method = gsub('.*-(.*)-.*-.*-.*-prcurve.tsv', '\\1', ff)
  df$sample = gsub('.*-.*-(.*)-.*-.*-prcurve.tsv', '\\1', ff)
  df$region = gsub('.*-.*-.*-(.*)-.*-prcurve.tsv', '\\1', ff)
  df$eval = gsub('.*-.*-.*-.*-(.*)-prcurve.tsv', '\\1', ff)
  df
})
df = do.call(rbind, df)
write.table(df, file='tsv/human-merged-prcurve.tsv', quote=FALSE, sep='\t', row.names=FALSE)
