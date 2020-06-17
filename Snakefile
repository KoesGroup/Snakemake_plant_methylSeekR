######################################################
# Snakemake pipeline for bisulfite sequencing analysis
######################################################


###########
# Libraries
###########
import pandas as pd

###############
# Configuration
###############
configfile: "config.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]

# fetch URL to transcriptome multi fasta from configfile
genome_url = config["refs"]["genome"]

########################
# Samples and conditions
########################

# read the tabulated separated table containing the sample, condition and fastq file information∂DE
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# create lists containing the sample names and conditions
SAMPLES = units.index.get_level_values('sample').unique().tolist()
samples = pd.read_csv(config["units"], dtype=str,index_col=0,sep="\t")
#CONDITIONS = list(pd.read_table(config["units"])["condition"])
samplefile = config["units"]


###########################
# Input functions for rules
###########################

def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    if "fq2" not in samples.columns:
        return True
    else:
        return pd.isnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
    """ This function checks if the sample has paired end or single end reads
    and returns 1 or 2 names of the fastq files """
    if sample_is_single_end(wildcards.sample):
        return samples.loc[(wildcards.sample), ["fq1"]].dropna()
    else:
        return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_trimmed(wildcards):
    """ This function checks if sample is paired end or single end
    and returns 1 or 2 names of the trimmed fastq files """
    if sample_is_single_end(wildcards.sample):
        return WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz"
    else:
        return [WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz", WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"]

def get_filenames(wildcards):
    """ This function checks if sample is paired end or single end
    and returns 1 or 2 names of the trimmed fastq files """
    if sample_is_single_end(wildcards.sample):
        return "-i " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz"
    else:
        return "-1 " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz -2 " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"


#################
# Desired outputs
#################
rule all:
    input:
        #expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES),
        #expand(WORKING_DIR + "result/{sample}.ATCGmap.gz", sample = SAMPLES),
        #expand(WORKING_DIR + "result/{sample}_CG.msr", sample = SAMPLES),
        #WORKING_DIR + "BSgenome_seed",
        #WORKING_DIR + "BSgenomeGenome/DESCRIPTION",
        #WORKING_DIR + "BSgenomeGenome/NAMESPACE",
        expand(WORKING_DIR + "results/{sample}_CG_UMRLMR_plant.bed", sample = SAMPLES),
        expand(WORKING_DIR + "results/{sample}_CCG_UMRLMR_plant.bed", sample = SAMPLES),
        expand(WORKING_DIR + "results/{sample}_CWG_UMRLMR_plant.bed", sample = SAMPLES),
        expand(WORKING_DIR + "results/{sample}_CHH_UMRLMR_plant.bed", sample = SAMPLES),
        expand(WORKING_DIR + "results/{sample}_Active_regions.bed", sample = SAMPLES),

    message:
        "Job done!\n\n\t#=========================#\n\t|       tijs bliek        |\n\t| University of Amsterdam |\n\t#=========================#\n"

#######
# Rules
#######


#####################
# Download references
#####################

rule get_genome_fasta:
    output:
        WORKING_DIR + "genome/genome.fasta.gz"
    message:
        "downloading the required genomic fasta file"
#    conda:
#        "envs/wget.yaml"
    shell:
        "wget -O {output} {genome_url}"


##################################
# Fastp
##################################

rule fastp:
    input:
        get_fastq
    output:
        fq1  = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html"
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    run:
        if sample_is_single_end(params.sampleName):
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --in1 {input} --out1 {output} \
            2> {log}; \
            touch {output.fq2}")
        else:
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2}; \
            2> {log}")
###########################
# BSseeker2 read alignement
###########################

rule BSseeker_build_index:
    input:
        WORKING_DIR + "genome/genome.fasta.gz"
    output:
        WORKING_DIR + "genome/genome.check"
        ##[WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome"
    params:
        WORKING_DIR + "genome/genome"
    conda:
        "envs/python2.yaml"
    threads: 10
    shell:
        "python2 BSseeker2/bs_seeker2-build.py -f {input} && touch {output}"

rule BSseeker_mapping:
    input:
        get_trimmed,
        check      = WORKING_DIR + "genome/genome.check",
        genome     = WORKING_DIR + "genome/genome.fasta.gz"
    output:
        bams  = WORKING_DIR + "mapped/{sample}.bam",
        log   = WORKING_DIR + "mapped/{sample}.bam.bs_seeker2_log"
    params:
        indexName  = WORKING_DIR + "genome/genome",
        sampleName = "{sample}",
        shellInput = get_filenames
    conda:
        "envs/python28.yaml"
    message:
        "mapping reads to genome to bam files."
    threads: 10
    shell:
       "python2 BSseeker2/bs_seeker2-align.py --bt-p {threads} {params.shellInput} -g {input.genome} -o {output.bams}"
#        print(shellInput)
#        if sample_is_single_end(params.sampleName):
#            shell("python2 BSseeker2/bs_seeker2-align.py --bt-p {threads} -i {input[0]} -g {input.genome} -o {output.bams}")
#        else:
#            shell("python2 BSseeker2/bs_seeker2-align.py --bt-p {threads} -1 {input[0]} -2 {input[1]} -g {input.genome} -o {output.bams}")


#####################################################
# Calling of the methylation levels and split to type
#####################################################

rule methylation_calling:
    input:
        bams  = WORKING_DIR + "mapped/{sample}.bam",
        index = WORKING_DIR + "genome/genome.check",
    output:
        atcgmap = WORKING_DIR + "results/{sample}.ATCGmap.gz"
    message:
        "calling methylation levels"
    params:
        prefix = "results/{sample}",
        index  = WORKING_DIR + "BSseeker2/bs_utils/reference_genomes/genome.fasta.gz_bowtie/"
    conda:
        "envs/python28.yaml"
    threads: 10
    shell:
        "python2 BSseeker2/bs_seeker2-call_methylation.py -i {input.bams} -o {params.prefix} -d {params.index}"


rule split_methylation_types:
    input:
        atcgmap = WORKING_DIR + "results/{sample}.ATCGmap.gz"
    output:
        CG  = WORKING_DIR + "results/{sample}_CG.msr",
        CCG = WORKING_DIR + "results/{sample}_CCG.msr",
        CWG = WORKING_DIR + "results/{sample}_CWG.msr",
        CHH = WORKING_DIR + "results/{sample}_CHH.msr",
        CHG = WORKING_DIR + "results/{sample}_CHG.msr"
    shell:"""
zless {input.atcgmap} | awk '{{if ($4 == "CG") print}}' | awk '{{ if ($2=="C") print $1 "\\t" $3 "\\t" $7+$8 "\\t" $8; if ($2=="G") print $1 "\\t" $3 "\\t" $11+$14 "\\t" $14}}' > {output.CG}
zless {input.atcgmap} | awk '{{if ($4 == "CHG" && $5 == "CC") print}}' | awk '{{ if ($2=="C") print $1 "\\t" $3 "\\t" $7+$8 "\\t" $8; if ($2=="G") print $1 "\\t" $3 "\\t" $11+$14 "\\t" $14}}' > {output.CCG}
zless {input.atcgmap} | awk '{{if ($4 == "CHG" && $5 ~ /C[AT]/) print}}' | awk '{{ if ($2=="C") print $1 "\\t" $3 "\\t" $7+$8 "\\t" $8; if ($2=="G") print $1 "\\t" $3 "\\t" $11+$14 "\\t" $14}}' > {output.CWG}
zless {input.atcgmap} | awk '{{if ($4 == "CHH") print}}' | awk '{{ if ($2=="C") print $1 "\\t" $3 "\\t" $7+$8 "\\t" $8; if ($2=="G") print $1 "\\t" $3 "\\t" $11+$14 "\\t" $14}}' > {output.CHH}
zless {input.atcgmap} | awk '{{if ($4 == "CHG") print}}' | awk '{{ if ($2=="C") print $1 "\\t" $3 "\\t" $7+$8 "\\t" $8; if ($2=="G") print $1 "\\t" $3 "\\t" $11+$14 "\\t" $14}}' > {output.CHG}
"""


###########################################
# determin regions of low or no methylation 
###########################################

rule forge_genome_data_package:
    input:
        genome = WORKING_DIR + "genome/genome.fasta.gz"
    output:
        seed        = WORKING_DIR + "BSgenome_seed",
        discription = WORKING_DIR + "BSgenomeGenome/DESCRIPTION",
        namespace   = WORKING_DIR + "BSgenomeGenome/NAMESPACE"
    message:
        "forging BSgenome"
    params:
        package_name    = config["forge"]["package_name"],
        organism        = config["forge"]["organism"],
        common_name     = config["forge"]["common_name"],
        genome_dir      = config["forge"]["genome_dir"],
        BSgenomeObjname = config["forge"]["BSgenomeObjname"]
    conda:
        "envs/forge_genome.yaml"
    shell:
        "python scripts/ForgeBSgenome.py "
        "-g {input.genome} "
        "-p {params.package_name} "
        "-o {params.organism} "
        "-c {params.common_name} "
        "-d {params.genome_dir} "
        "-b {params.BSgenomeObjname}"

rule methylSeekR:
    input:
        CG          = WORKING_DIR + "results/{sample}_CG.msr",
        CCG         = WORKING_DIR + "results/{sample}_CCG.msr",
        CWG         = WORKING_DIR + "results/{sample}_CWG.msr",
        CHH         = WORKING_DIR + "results/{sample}_CHH.msr",
        CHG         = WORKING_DIR + "results/{sample}_CHG.msr",
        seed        = WORKING_DIR + "BSgenome_seed",
        discription = WORKING_DIR + "BSgenomeGenome/DESCRIPTION",
        namespace   = WORKING_DIR + "BSgenomeGenome/NAMESPACE"
    output:
        CG  = WORKING_DIR + "results/{sample}_CG_UMRLMR.bed",
        CCG = WORKING_DIR + "results/{sample}_CCG_UMRLMR.bed",
        CWG = WORKING_DIR + "results/{sample}_CWG_UMRLMR.bed",
        CHH = WORKING_DIR + "results/{sample}_CHH_UMRLMR.bed",
        CHG = WORKING_DIR + "results/{sample}_CHG_UMRLMR.bed",
    params:
        LMR = config["UMRLMR"]["LMR"],
    conda:
        "envs/methylseekr.yaml"
    message:
        "running R-sript methylSeekR.R"        
    shell:
        "Rscript scripts/methylSeekR.R "
        "-l {params.LMR} "
        "--CGin {input.CG} "
        "--CCGin {input.CCG} "
        "--CWGin {input.CWG} "
        "--CHHin {input.CHH} "
        "--CHGin {input.CHG} "
        "--CGout {output.CG} "
        "--CCGout {output.CCG} "
        "--CWGout {output.CWG} "
        "--CHHout {output.CHH} "
        "--CHGout {output.CHG}"

###################################
# Redefine regions to plants sample
###################################

rule reDefine_regions:
    input:
        CG  = WORKING_DIR + "results/{sample}_CG_UMRLMR.bed",
        CCG = WORKING_DIR + "results/{sample}_CCG_UMRLMR.bed",
        CWG = WORKING_DIR + "results/{sample}_CWG_UMRLMR.bed",
        CHH = WORKING_DIR + "results/{sample}_CHH_UMRLMR.bed",
        CHG = WORKING_DIR + "results/{sample}_CHG_UMRLMR.bed",
    output:
        CG  = WORKING_DIR + "results/{sample}_CG_UMRLMR_plant.bed",
        CCG = WORKING_DIR + "results/{sample}_CCG_UMRLMR_plant.bed",
        CWG = WORKING_DIR + "results/{sample}_CWG_UMRLMR_plant.bed",
        CHH = WORKING_DIR + "results/{sample}_CHH_UMRLMR_plant.bed",
        CHG = WORKING_DIR + "results/{sample}_CHG_UMRLMR_plant.bed",
    params:
        UMR = config["UMRLMR"]["UMR"],
    message:
        "running reDefineRegions.py"        
    shell:
        "python scripts/reDefineRegion.py "
        "-m {params.UMR} "
        "--CG {input.CG} "
        "--CCG {input.CCG} "
        "--CWG {input.CWG} "
        "--CHH {input.CHH} "
        "--CHG {input.CHG} "
        "--CGp {output.CG} "
        "--CCGp {output.CCG} "
        "--CWGp {output.CWG} "
        "--CHHp {output.CHH} "
        "--CHGp {output.CHG}"

#############################
# find active genomic regions
#############################

rule get_active_regions:
    input:
        CG  = WORKING_DIR + "results/{sample}_CG_UMRLMR_plant.bed",
        CHG = WORKING_DIR + "results/{sample}_CHG_UMRLMR_plant.bed",
    output:
        CG  = WORKING_DIR + "results/{sample}_CG_UMR_plant.bed",
        CHG = WORKING_DIR + "results/{sample}_CHG_UMR_plant.bed",
        AR  = WORKING_DIR + "results/{sample}_Active_regions.bed"
    message:
        "running reDefineRegions.py"
    shell: """
grep UMR {input.CG} > {output.CG}
grep UMR {input.CHG} > {output.CHG}
bedtools intersect -wao -a {output.CG} -b {output.CHG} > {output.AR}
"""
