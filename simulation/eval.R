## Reads all the VCFs and evaluate the genotyping accuracy
## using the sveval package (https://github.com/jmonlong/sveval).

library(sveval)
library(magrittr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

winsor <- function(x, u){
  if(any(x>u)) x[x>u] = u
  x
}

## Should cached results be used (is present)
CACHE=TRUE

## Read VCFs
if(!file.exists('vcfs.RData') | !CACHE){
  ## Read truth set for each sample
  samps = c('s0', 's1', 's2')
  truth = lapply(samps, function(samp){
    gr = readSVvcf('truth.symb.vcf', sample.name=samp, check.inv=TRUE)
    gr$sample = samp
    gr
  })
  truth = do.call(c, truth)

  message('Importing VCFs')
  ## Read each call set and create a list of GRanges
  covs = list.files('.', 'cov')
  graphs = c('calls', 'truth')
  methods = c('toilvg', 'svtyper', 'delly', 'bayestyper', 'paragraph')
  calls = lapply(covs, function(cov){
    lapply(samps, function(samp){
      lapply(graphs, function(graph){
        lapply(methods, function(meth){
          message(paste(cov, samp, graph, meth))
          vcf.file = paste0(cov, '/', graph, '.', samp, '.', meth, '.vcf')
          if(!file.exists(vcf.file)){
            vcf.file = paste0(vcf.file, '.gz')
          }
          res = readSVvcf(vcf.file, check.inv=TRUE)
          if(meth=='delly'){
            vcf.file = paste0(cov, '/', graph, '.', samp, '.', meth, '.inv.vcf')
            if(!file.exists(vcf.file)){
              vcf.file = paste0(vcf.file, '.gz')
            }
            res.inv = readSVvcf(vcf.file, check.inv=TRUE)
            res = c(res, res.inv)
          }
          if(length(res)==0) return(NULL)
          res$method = meth
          res$graph = graph
          res$sample = samp
          res$depth = as.numeric(gsub('cov','', cov))
          res
        })
      })
    })
  })
  calls = unlist(calls)

  calls.symb = lapply(covs, function(cov){
    lapply(samps, function(samp){
      lapply(graphs, function(graph){
        meth = 'toilvg'
        ## message(paste(cov, samp, graph, meth))
        vcf.file = paste0(cov, '/', graph, '.symb.', samp, '.', meth, '.vcf')
        if(!file.exists(vcf.file)){
          vcf.file = paste0(vcf.file, '.gz')
        }
        res = readSVvcf(vcf.file, check.inv=TRUE)
        if(length(res)==0) return(NULL)
        res$method = paste0(meth, '-symb')
        res$graph = graph
        res$sample = samp
        res$depth = as.numeric(gsub('cov','', cov))
        res
      })
    })
  })

  calls = c(calls, unlist(calls.symb))

  save(truth, calls, file='vcfs.RData')
} else {
  load('vcfs.RData')
}



## SV evaluation
message('Evaluation')

#### "Calling" evaluation (presence of SV)
eval.l = lapply(calls, function(call.gr){
  samp = call.gr$sample[1]
  message(paste(samp, call.gr$method[1], call.gr$graph[1], call.gr$depth[1]))
  truth = truth[which(truth$sample==samp)]
  res = suppressMessages(svevalOl(call.gr, truth, max.ins.dist=20, min.cov=.5, min.del.rol=.1, min.size=50, check.inv=TRUE))
  eval.df = res$eval
  eval.pr = res$curve
  eval.pr$method = call.gr$method[1]
  eval.pr$graph = call.gr$graph[1]
  eval.pr$sample = call.gr$sample[1]
  eval.pr$depth = call.gr$depth[1]
  return(list(eval=eval.df, curve=eval.pr))
})
nonna = unlist(lapply(eval.l, function(e) any(colnames(e$eval)=='precision')))
eval.pr = do.call(rbind, lapply(eval.l[which(nonna)], function(e) e$curve))
write.table(eval.pr, file='simerror-prcurve.tsv', row.names=FALSE, quote=FALSE, sep='\t')

#### Genotyping evaluation (actual genotype)
eval.geno.l = lapply(calls, function(call.gr){
  samp = call.gr$sample[1]
  message(paste(samp, call.gr$method[1], call.gr$graph[1], call.gr$depth[1]))
  truth = truth[which(truth$sample==samp)]
  res = suppressMessages(svevalOl(call.gr, truth, max.ins.dist=20, min.cov=.5, min.del.rol=.1, min.size=50, check.inv=TRUE,
                                  geno.eval=TRUE, stitch.hets=TRUE, merge.hets=TRUE))
  eval.df = res$eval
  eval.pr = res$curve
  eval.pr$method = call.gr$method[1]
  eval.pr$graph = call.gr$graph[1]
  eval.pr$sample = call.gr$sample[1]
  eval.pr$depth = call.gr$depth[1]
  return(list(eval=eval.df, curve=eval.pr))
})
nonna = unlist(lapply(eval.geno.l, function(e) any(colnames(e$eval)=='precision')))
eval.geno.pr = do.call(rbind, lapply(eval.geno.l[which(nonna)], function(e) e$curve))
write.table(eval.geno.pr, file='simerror-geno-prcurve.tsv', row.names=FALSE, quote=FALSE, sep='\t')
