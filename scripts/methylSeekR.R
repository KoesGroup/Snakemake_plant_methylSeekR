library("BSgenomeGenome", lib.loc="localR")
library("mhsmm", lib.loc="localR")
library("MethylSeekR", lib.loc="localR")
library(optparse)

# arguments to provide
option_list = list(
  make_option(c("-a", "--CGin"), type="character", default="None", help="methyl calling output for CG", metavar="character"),
  make_option(c("-b", "--CCGin"), type="character", default="None", help="methyl calling output for CCG", metavar="character"),
  make_option(c("-c", "--CWGin"), type="character", default="None", help="methyl calling output for CWG", metavar="character"),
  make_option(c("-d", "--CHHin"), type="character", default="None", help="methyl calling output for CHH", metavar="character"),
  make_option(c("-e", "--CHGin"), type="character", default="None", help="methyl calling output for CHG", metavar="character"),
  make_option(c("-w", "--CGout"), type="character", default="None", help="methyl calling output for CG", metavar="character"),
  make_option(c("-x", "--CCGout"), type="character", default="None", help="methyl calling output for CCG", metavar="character"),
  make_option(c("-y", "--CWGout"), type="character", default="None", help="methyl calling output for CWG", metavar="character"),
  make_option(c("-z", "--CHHout"), type="character", default="None", help="methyl calling output for CHH", metavar="character"),
  make_option(c("-u", "--CHGout"), type="character", default="None", help="methyl calling output for CHG", metavar="character"),
  make_option(c("-l", "--upperLMR"), type="double", default=50.0, help="upperlimit to be LMR", metavar="double")
  make_option(c("-n", "--CpGnumber"), type="integer", default=4, help="upperlimit to be LMR", metavar="double")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


m.sel = opt$upperLMR/100     # threshold for % methylation
n.sel = opt$CpGnumber       # threshold for # of CpGs
core.num = 5    # # of cores to use
set.seed(123)

genome <- BSgenomeGenome
sLengths=seqlengths(genome)


meth.gr <- readMethylome(opt$CGin, seqLengths=sLengths)
UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, myGenomeSeq=genome, seqLengths=sLengths, pdfFilename=paste0("Bvul_CG_MTR_m", m.sel, "_n", n.sel, ".pdf"))
saveUMRLMRSegments(segs=UMRLMRsegments.gr, TableFilename=opt$CGout)

meth.gr <- readMethylome(opt$CCGin, seqLengths=sLengths)
UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, myGenomeSeq=genome, seqLengths=sLengths, pdfFilename=paste0("Bvul_CG_MTR_m", m.sel, "_n", n.sel, ".pdf"))
saveUMRLMRSegments(segs=UMRLMRsegments.gr, TableFilename=opt$CCGout)

meth.gr <- readMethylome(opt$CWGin, seqLengths=sLengths)
UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, myGenomeSeq=genome, seqLengths=sLengths, pdfFilename=paste0("Bvul_CG_MTR_m", m.sel, "_n", n.sel, ".pdf"))
saveUMRLMRSegments(segs=UMRLMRsegments.gr, TableFilename=opt$CWGout)

meth.gr <- readMethylome(opt$CHHin, seqLengths=sLengths)
UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, myGenomeSeq=genome, seqLengths=sLengths, pdfFilename=paste0("Bvul_CG_MTR_m", m.sel, "_n", n.sel, ".pdf"))
saveUMRLMRSegments(segs=UMRLMRsegments.gr, TableFilename=opt$CHHout)

meth.gr <- readMethylome(opt$CHGin, seqLengths=sLengths)
UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, myGenomeSeq=genome, seqLengths=sLengths, pdfFilename=paste0("Bvul_CG_MTR_m", m.sel, "_n", n.sel, ".pdf"))
saveUMRLMRSegments(segs=UMRLMRsegments.gr, TableFilename=opt$CHGout)
