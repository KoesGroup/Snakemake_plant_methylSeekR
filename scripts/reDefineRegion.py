from optparse import OptionParser 

# get the command line arguments
parser = OptionParser(description="script written for a pipeline to forge a BSgenome data package")
parser.add_option('-a', '--CG', 
                    type=str,
                    default="None",
                    metavar="",
                    help = "Name of MethylseekR output file for CG methylation regions")
parser.add_option('-b', '--CCG', 
                    type=str,
                    default="None",
                    metavar="",
                    help = "Name of MethylseekR output file for CCG methylation regions")
parser.add_option('-c', '--CWG', 
                    type=str,
                    default="None",
                    metavar="",
                    help = "Name of MethylseekR output file for CWG methylation regions")
parser.add_option('-d', '--CHH', 
                    type=str,
                    default="None",
                    metavar="",
                    help = "Name of MethylseekR output file for CHH methylation regions")
parser.add_option('-w', '--CGp', 
                    type=str,
                    default="None",
                    metavar="",
                    help = "Name of output file containing regions adjusted to plants for CG methylation regions")
parser.add_option('-x', '--CCGp', 
                    type=str,
                    default="None",
                    metavar="",
                    help = "Name of output file containing regions adjusted to plants for CCG methylation regions")
parser.add_option('-y', '--CWGp', 
                    type=str,
                    default="None",
                    metavar="",
                    help = "Name of output file containing regions adjusted to plants for CWG methylation regions")
parser.add_option('-z', '--CHHp', 
                    type=str,
                    default="None",
                    metavar="",
                    help = "Name of output file containing regions adjusted to plants for CHH methylation regions")
parser.add_option('-m', '--UMRmax', 
                    type=float,
                    default=10,
                    metavar="",
                    help = "maximum percentage of methylation accepted to be defined unmethylated")    

(options, args) = parser.parse_args()


bed = open(options.CG)
uit = open(options.CGp, "w")

for l in bed:
    if l.startswith("chr\tstart"):
        uit.write(l)
    else:
        l = l.split("\t")
        if float(l[-1].rstrip()) > options.UMRmax:
            l[3] = "LMR"
        else:
            l[3] = "UMR"
        l = "\t".join(l)
        uit.write(l)



uit.close()
bed.close()