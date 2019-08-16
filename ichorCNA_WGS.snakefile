configfile: "{}/ichorCNA/config_hg38.yaml".format(workflow.basedir)

include: "pre_pro_af_merge.snakefile"

rule all:
  input:
    expand("results/ichorCNA/{tumor}/{tumor}.cna.seg", tumor=config["pairs"]),
#   expand("results/ichorCNA/{tumor}/{tumor}.cna.seg", tumor=config["samples"]),
    expand("results/readDepth/{samples}.bin{binSize}.wig", samples=config["merge_libs"], binSize=str(config["binSize"]))

rule read_counter:
  input:
    "GATK_runs/{samples}/ApplyBQSR/{samples}_recal.bam"
  output:
    "results/readDepth/{samples}.bin{binSize}.wig"		
  conda:    
    "envs_dir/ichorCNA_env.yaml"
	params:
		readCounter=config["readCounterScript"],
		binSize=config["binSize"],
		qual="20",
		chrs=config["chrs"]
	resources:
		mem_mb=4000
	log:
		"logs/readDepth/{samples}.bin{binSize}.log"
	shell:
		"{params.readCounter} {input} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"

rule ichorCNA:
  input:
    tum="results/readDepth/{tumor}.bin" + str(config["binSize"]) + ".wig",
    norm=lambda wildcards: "results/readDepth/" + config["pairs"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
  output:
    corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
    #param="results/ichorCNA/{tumor}/{tumor}.params.txt",
    cna="results/ichorCNA/{tumor}/{tumor}.cna.seg",
    #segTxt="results/ichorCNA/{tumor}/{tumor}.seg.txt",
    #seg="results/ichorCNA/{tumor}/{tumor}.seg",
    #rdata="results/ichorCNA/{tumor}/{tumor}.RData",
    #outDir="results/ichorCNA/{tumor}/",
  conda:    
    "envs_dir/ichorCNA_env.yaml"
  params:
    outDir="results/ichorCNA/{tumor}/",
    rscript=config["ichorCNA_rscript"],
    id="{tumor}",
    ploidy=config["ichorCNA_ploidy"],
    normal=config["ichorCNA_normal"],
    gcwig=config["ichorCNA_gcWig"],
    mapwig=config["ichorCNA_mapWig"],
    normalpanel=config["ichorCNA_normalPanel"],
    estimateNormal=config["ichorCNA_estimateNormal"],
    estimatePloidy=config["ichorCNA_estimatePloidy"],
    estimateClonality=config["ichorCNA_estimateClonality"],
    scStates=config["ichorCNA_scStates"],
    maxCN=config["ichorCNA_maxCN"],
    includeHOMD=config["ichorCNA_includeHOMD"],
    chrs=config["ichorCNA_chrs"],
    chrTrain=config["ichorCNA_chrTrain"],
    genomeStyle=config["ichorCNA_genomeStyle"],
    centromere=config["ichorCNA_centromere"],
    exons=config["ichorCNA_exons"],
    txnE=config["ichorCNA_txnE"],
    txnStrength=config["ichorCNA_txnStrength"],
    fracReadsChrYMale="0.001",
    plotFileType=config["ichorCNA_plotFileType"],
    plotYlim=config["ichorCNA_plotYlim"],
    libdir=config["ichorCNA_libdir"]
  resources:
    mem_mb=4000
  log:
    "logs/ichorCNA/{tumor}.log"	
  shell:
    "Rscript {params.rscript} \
            --id {params.id}  \
            --libdir {params.libdir} \
            --WIG {input.tum} \
            --gcWig {params.gcwig} \
            --NORMWIG {input.norm} \
            --mapWig {params.mapwig} \
            --ploidy \"{params.ploidy}\" \
            --normal \"{params.normal}\" \
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

#"Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --fracReadsInChrYForMale {params.fracReadsChrYMale} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"

