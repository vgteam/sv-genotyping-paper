library(sveval)
args = commandArgs(TRUE)
load(args[1])
df = plot_persize(eval, plot=FALSE,
                  size.breaks=c(50,100,200,300,400,600,800,
                                1000,2500,5000,Inf))
write.table(df, file=args[2], sep='\t', quote=FALSE, row.names=FALSE)
