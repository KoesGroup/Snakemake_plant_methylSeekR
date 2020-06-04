library(BSgenome)
require(Biostrings)
library(optparse)

# arguments to provide
option_list = list(
  make_option(c("-s", "--seed_file"), type="character", default="BSgenome_seed", help="name of the seed file to be used for forging the data package", metavar="character")
) 

# parse the command-line arguments and pass them to a list called 'opt'
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## run forge BSgenome on the seed file
forgeBSgenomeDataPkg(opt$seed_file)
