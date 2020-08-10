configfile: "{}/ichorCNA/config_hg38.yaml".format(workflow.basedir)
configfile: "{}/ref.yaml".format(workflow.basedir)

include: "./pre_pro_af_merge_alt_bed.snakefile"

all_samples = config['tumors'] + config['normals']

rule all:
  input:
#   expand("results/ichorCNA/{tumor}/{tumor}.cna.seg", tumor = config["pairings"]),
    expand("results/ichorCNA/{tumor}/{tumor}.cna.seg", tumor = config['tumors']),
#    expand("results/readDepth/{samples}.bin{binSize}.wig", samples = all_samples, binSize = str(config["binSize"]))
#    "results/PoN/PoN_median.rds"

rule read_counter:
  input:
    lambda wildcards: "GATK_runs/{samples}/ApplyBQSR/{samples}_recal.bam"
  output:
    "results/readDepth/{samples}.bin{binSize}.wig"
  conda:
    "envs_dir/ichorCNA_env_pon.yaml"
	params:
		readCounter = config["readCounterScript"],
		binSize = config["binSize"],
		qual = "20",
		chrs = config["chrs"]
	resources:
		mem_mb = 4000
	log:
		"logs/readDepth/{samples}.bin{binSize}.log"
	shell:
		"{params.readCounter} {input} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"

rule install_ichorCNA:
  output: "signal.txt"
  conda:
    "envs_dir/ichorCNA_env_pon.yaml"
  resources:
    mem_mb = 4000
  log:
    "logs/install_ichorCNA.log"
  shell:
    """
      Rscript -e "library(devtools);install_github(\\"broadinstitute/ichorCNA\\")" &> {log}
      echo "installed ichorCNA successfully" > {output}
    """


rule make_PoN:
  input:
    wigs = expand("results/readDepth/{samples}.bin{binSize}.wig", samples = config['normals'], binSize = str(config["binSize"])),
    sig = rules.install_ichorCNA.output
  output:
    "results/PoN/PoN_median.rds"
  conda:
    "envs_dir/ichorCNA_env_pon.yaml"
  params:
    gcwig = config["ichorCNA_gcWig"],
    chrs = config["ichorCNA_chrs"],
    chrs_norm = config["ichorCNA_chrTrain"],
    mapwig = config["ichorCNA_mapWig"],
    centromere = config["ichorCNA_centromere"],
    exons = f"--exons.bed {config['alt_bed']}" if config["alt_bed"] is not None else "",
    rscript = config['bioinfo_workflows_path'] + '/ichorCNA/createPanelOfNormals.R',
    outfile_pref = "results/PoN/PoN"
  resources:
    mem_mb = 4000
  log:
    "logs/ichorCNA/PoN.log"
  shell:
    """
      if [[ !(-d results/PoN) ]]; then
        mkdir -p results/PoN/
      fi
      echo "{input.wigs}" | sed 's/ /\\n/g' > results/PoN/norm_files.txt
      {params.rscript} \
        --filelist results/PoN/norm_files.txt \
        --gcWig {params.gcwig} \
        --mapWig {params.mapwig} \
        {params.exons} \
        --centromere {params.centromere} \
        --outfile {params.outfile_pref} &> {log}
    """

rule ichorCNA:
  input:
    tum = "results/readDepth/{tumor}.bin" + str(config["binSize"]) + ".wig",
    pon = rules.make_PoN.output,
#   norm = lambda wildcards: "results/readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
    sig = rules.install_ichorCNA.output
  output:
    #corrDepth = "results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
    #param = "results/ichorCNA/{tumor}/{tumor}.params.txt",
    cna = "results/ichorCNA/{tumor}/{tumor}.cna.seg",
    #segTxt = "results/ichorCNA/{tumor}/{tumor}.seg.txt",
    #seg = "results/ichorCNA/{tumor}/{tumor}.seg",
    #rdata = "results/ichorCNA/{tumor}/{tumor}.RData",
    #outDir = "results/ichorCNA/{tumor}/",
  conda:
    "envs_dir/ichorCNA_env_pon.yaml"
  params:
    outDir = "results/ichorCNA/{tumor}/",
    rscript = config['bioinfo_workflows_path'] + '/new_ichorCNA/ichorCNA/scripts/runIchorCNA.R',
    id = "{tumor}",
    ploidy = config["ichorCNA_ploidy"],
    normal = config["ichorCNA_normal"],
    gcwig = config["ichorCNA_gcWig"],
    mapwig = config["ichorCNA_mapWig"],
    normalpanel = config["ichorCNA_normalPanel"],
    estimateNormal = config["ichorCNA_estimateNormal"],
    estimatePloidy = config["ichorCNA_estimatePloidy"],
    estimateClonality = config["ichorCNA_estimateClonality"],
    scStates = config["ichorCNA_scStates"],
    maxCN = config["ichorCNA_maxCN"],
    includeHOMD = config["ichorCNA_includeHOMD"],
    chrs = config["ichorCNA_chrs"],
    chrTrain = config["ichorCNA_chrTrain"],
    genomeStyle = config["ichorCNA_genomeStyle"],
    centromere = config["ichorCNA_centromere"],
    exons = f"--exons.bed {config['alt_bed']}" if config["alt_bed"] is not None else "",
    txnE = config["ichorCNA_txnE"],
    txnStrength = config["ichorCNA_txnStrength"],
    fracReadsChrYMale = "0.001",
    plotFileType = config["ichorCNA_plotFileType"],
    plotYlim = config["ichorCNA_plotYlim"],
    libdir = config["ichorCNA_libdir"]
  resources:
    mem_mb = 4000
  log:
    "logs/ichorCNA/{tumor}.log"
  shell:
    "Rscript {params.rscript} \
            --id {params.id}  \
            --libdir {params.libdir} \
            --WIG {input.tum} \
            --gcWig {params.gcwig} \
            --mapWig {params.mapwig} \
            --ploidy \"{params.ploidy}\" \
            {params.exons}  \
            --normalPanel {input.pon}  \
            --maxCN {params.maxCN} \
            --includeHOMD {params.includeHOMD} \
            --chrs \"{params.chrs}\" \
            --chrTrain \"{params.chrTrain}\" \
            --genomeStyle {params.genomeStyle} \
            --estimateNormal {params.estimateNormal} \
            --estimatePloidy {params.estimatePloidy} \
            --estimateScPrevalence {params.estimateClonality} \
            --scStates \"{params.scStates}\" \
            --centromere {params.centromere} \
            --txnE {params.txnE} \
            --txnStrength {params.txnStrength} \
            --fracReadsInChrYForMale {params.fracReadsChrYMale} \
            --plotFileType {params.plotFileType} \
            --plotYLim \"{params.plotYlim}\" \
            --outDir {params.outDir} > {log} 2> {log}"
#
##"Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --fracReadsInChrYForMale {params.fracReadsChrYMale} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"
#
