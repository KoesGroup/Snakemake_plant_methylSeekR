import os 
import glob
from datetime import date
from optparse import OptionParser 

# get the command line arguments
parser = OptionParser(description="script written for a pipeline to forge a BSgenome data package")
parser.add_option('-g', '--genome', 
                    type=str,
                    default="genome.fasta",
                    metavar="",
                    help = "name of of the genome fasta file, default = genome.fasta")
parser.add_option('-o', '--organism', 
                    type=str,
                    default="Unknown",
                    metavar="",
                    help = "orgasm the genome originates from, default =Unknown")
parser.add_option('-c', '--common_name', 
                    type=str,
                    default="unknow species",
                    metavar="",
                    help = "common name of the organism , default = unknown species")
parser.add_option('-p', '--package_name', 
                    type=str,
                    default="BSgenomeGenome",
                    metavar="",
                    help = "Name of the data package to be created, default = BSgenomeGenome")
parser.add_option('-d', '--genomeDir', 
                    type=str,
                    default="genome",
                    metavar="",
                    help = "name of the directory contianing the genome fasta file, default = genome")
parser.add_option('-b', '--BSgenomeObjname', 
                    type=str,
                    default="Chromosomes",
                    metavar="",
                    help = "the BSgenome object name. Usualy the name of the species used, default = Chromosomes")


(options, args) = parser.parse_args()

## check if genome in zipped, if so, unzip.

genomeFasta = options.genome.replace("genome/","")
print(genomeFasta)

os.chdir("genome/")

print(os.getcwd())
print(glob.glob("*.gz"))
if genomeFasta in glob.glob("*.gz"):
    print("unzipping......")
    os.system("gunzip " + genomeFasta)


## create seperate files for each of the chromosomes.
genomeFasta = options.genome.replace(".gz", "")
genome = open(genomeFasta)
go = False
chroms = []
for l in genome:
    if l.startswith(">"):
        if go:
            chromosome.close()
        name = l[1:].rstrip() + ".fa" 
        chromosome = open(name, "w")
        chromosome.write(l)
        chroms.append(l[1:].rstrip())
        go = True
    elif go:
        chromosome.write(l)
if go:
    chromosome.close()
else:
    print("No fasta header found!.")

os.chdir("..")
## create seed file
seedFile = open("BSgenome_seed", "w")
d1 = date.today().strftime("%d/%m/%Y")   ## get date
chroms = ",".join(['"'+x+'"' for x in chroms])  ## list of the chromosome names

seedText = "Package: " + options.package_name + "\n\
Title: BS_gen\n\
Description: BSgenome data package needed to run the methylseekR script.\n\
Version: 1.0\n\
organism: " + options.organism + "\n\
common_name: " + options.common_name + "\n\
provider: NA\n\
provider_version: NA\n\
release_date: " + d1 + "\n\
release_name: Plant_methylSeekR_BSgenome\n\
organism_biocview: NA\n\
BSgenomeObjname: " + options.BSgenomeObjname + "\n\
seqs_srcdir: " + options.genomeDir + "\n\
seqfiles_suffix: .fa\n\
seqnames: c(" + chroms + ")\n"

seedFile.write(seedText)
seedFile.close()

## create and install package.
command = "rscript Forge.R"
print(command)
#os.system(command)
command = "R CMD build " + options.package_name + " --no-manual"
print(command)
#os.system(command)
command = "R CMD check " + options.package_name + " --no-manual"
print(command)
#os.system(command)
command = "R CMD INSTALL " + options.package_name + " --no-manual"
print(command)
#os.system(command)
