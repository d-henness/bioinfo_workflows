configfile: "{}/TitanCNA/scripts/snakemake/config/config.yaml".format(workflow.basedir)

include: "ichorCNA_alt.snakefile"
include: "getAlleleCounts.snakefile"
import os.path

CLUST = {1:[1], 2:[1,2], 3:[1,2,3], 4:[1,2,3,4], 5:[1,2,3,4,5], 6:[1,2,3,4,5,6], 7:[1,2,3,4,5,6,7], 8:[1,2,3,4,5,6,7,8], 9:[1,2,3,4,5,6,7,8,9], 10:[1,2,3,4,5,6,7,8,9,10]}
PLOIDY = {2:[2], 3:[2,3], 4:[2,3,4]}


rule all_Titan:
	input:
		expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairs"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		#expand("results/titan/hmm/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		"results/titan/hmm/optimalClusterSolution.txt"

#rule makeOutDir:
#  output:
#    "results/titan/hmm/titanCNA_ploidy{ploidy}/"
#  resources:
#    mem_mb = 1
#  shell:
#    "mkdir -p {output}"

rule runTitanCNA:
  input:
    alleleCounts="results/titan/tumCounts/{tumor}.tumCounts.txt",
    corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt"
  output:
#   outRoot="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}/",
    titan="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
    param="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
    segTxt="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt",
    seg="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.seg"
  params:
    outRoot="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}/",
    titanRscript=config["TitanCNA_rscript"],
    libdir=config["TitanCNA_libdir"],
    numCores=config["TitanCNA_numCores"],
    normal=config["TitanCNA_normalInit"],
    chrs=config["TitanCNA_chrs"],
    genomeStyle=config["genomeStyle"],
    genomeBuild=config["genomeBuild"],
    cytobandFile=config["cytobandFile"],
    estimatePloidy=config["TitanCNA_estimatePloidy"],
    estimateClonality=config["TitanCNA_estimateClonality"],
    estimateNormal=config["TitanCNA_estimateNormal"],
    centromere=config["centromere"],
    alphaK=config["TitanCNA_alphaK"],
    alphaKmax=config["TitanCNA_alphaK"],
    #alphaR=config["TitanCNA_alphaR"],
    #alleleModel=config["TitanCNA_alleleModel"],
    txnExpLen=config["TitanCNA_txnExpLen"],
    plotYlim=config["TitanCNA_plotYlim"],
  conda:
    "envs_dir/Titan.yaml"
  resources:
    mem_mb = 16384
  log:
    "logs/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.log"
  benchmark:
    "results/benchmarks/runTitanCNA_titanCNA_ploidy{ploidy}_{tumor}_cluster{clustNum}.log"
  shell:
    "Rscript {workflow.basedir}/TitanCNA/scripts/R_scripts/titanCNA.R --hetFile {input.alleleCounts} --cnFile {input.corrDepth} --outFile {output.titan} --outSeg {output.segTxt} --outParam {output.param} --outIGV {output.seg} --outPlotDir {params.outRoot} --libdir {params.libdir} --id {wildcards.tumor} --numClusters {wildcards.clustNum} --numCores {params.numCores} --normal_0 {params.normal} --ploidy_0 {wildcards.ploidy} --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --cytobandFile {params.cytobandFile} --chrs \"{params.chrs}\" --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateClonality {params.estimateClonality}  --centromere {params.centromere} --alphaK {params.alphaK} --alphaKHigh {params.alphaK} --txnExpLen {params.txnExpLen} --plotYlim \"{params.plotYlim}\" > {log} 2> {log}"
#   "Rscript $CONDA_PREFIX/bin/{params.titanRscript} --hetFile {input.alleleCounts} --cnFile {input.corrDepth} --outFile {output.titan} --outSeg {output.segTxt} --outParam {output.param} --outIGV {output.seg} --outPlotDir {params.outRoot} --libdir {params.libdir} --id {wildcards.tumor} --numClusters {wildcards.clustNum} --numCores {params.numCores} --normal_0 {params.normal} --ploidy_0 {wildcards.ploidy} --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --cytobandFile {params.cytobandFile} --chrs \"{params.chrs}\" --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateClonality {params.estimateClonality}  --centromere {params.centromere} --alphaK {params.alphaK} --txnExpLen {params.txnExpLen} --plotYlim \"{params.plotYlim}\" > {log} 2> {log}"

#--alleleModel {params.alleleModel} --alphaR {params.alphaR}


rule selectSolution:
  input:
#   ploidyDirs=expand("results/titan/hmm/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
    resultFiles=expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairs"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]])
  output:
    "results/titan/hmm/optimalClusterSolution.txt"
  params:
    ploidyDirs="results/titan/hmm/titanCNA_ploidy2",
    solutionRscript=config["TitanCNA_selectSolutionRscript"],
    threshold=config["TitanCNA_solutionThreshold"],
  conda:
    "envs_dir/Titan.yaml"
  resources:
    mem_mb = 4096
  log:
    "logs/titan/selectSolution.log"
  shell:
    """
    if [ -d results/titan/hmm/titanCNA_ploidy3/ ]; then
      ploidyRun3=results/titan/hmm/titanCNA_ploidy3/
    else
      ploidyRun3=NULL
    fi
    if [ -d results/titan/hmm/titanCNA_ploidy4/ ]; then
      ploidyRun4=results/titan/hmm/titanCNA_ploidy4/
    else
      ploidyRun4=NULL
    fi
    Rscript {workflow.basedir}/TitanCNA/scripts/R_scripts/selectSolution.R --ploidyRun2 {params.ploidyDirs} --ploidyRun3 $ploidyRun3 --ploidyRun4 $ploidyRun4 --threshold {params.threshold} --outFile {output} &> {log}
    """
#   Rscript $CONDA_PREFIX/bin/{params.solutionRscript} --ploidyRun2 {params.ploidyDirs} --ploidyRun3 $ploidyRun3 --ploidyRun4 $ploidyRun4 --threshold {params.threshold} --outFile {output} &> {log}
#	run:
#		if "results/titan/hmm/titanCNA_ploidy3" in input:
#			ploidyRun3 = input[1]
#		else:
#			ploidyRun3 = "NULL"
#		if "results/titan/hmm/titanCNA_ploidy4" in input:
#			ploidyRun4 = input[2]
#		else:
#			ploidyRun4 = "NULL"
#		os.system("Rscript params.solutionRscript --ploidyRun2 input[0] --ploidyRun3 ploidyRun3 --ploidyRun4 ploidyRun4 --threshold params.threshold --outFile output")


