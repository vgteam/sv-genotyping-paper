## Arguments
## 1: call vcf
## 2: truth vcf
## 3: sample name
## 4: look for inversions
## 5: regions of interest (NA means whole genome)
## 6: genotype evaluation?
## 7: output Rdata file
## 8: minimum coverage to match variants (optional, default:0.5)
args = commandArgs(TRUE)

if(length(args)==7){
  args = c(args, '0.5')
}

library(sveval)
bed = NULL
if(args[5] != 'NA'){
  bed = args[5]
}
eval.o = svevalOl(args[1], args[2], sample.name=args[3], check.inv=as.logical(args[4]),
                  bed.regions=bed, min.cov=as.numeric(args[8]), 
                  geno.eval=args[6], stitch.hets=as.logical(args[6]),
                  merge.hets=args[6], min.size=50)
save(eval.o, file=args[7])
