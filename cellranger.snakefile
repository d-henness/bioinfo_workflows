rule cellranger_all:
    input: expand("{sample}_cellranger/outs/possorted_genome_bam.bam", sample = config['samples'])

rule cellranger:
    input:
        fq_dir = "{sample}"
    output:
        maf = "{sample}_cellranger/outs/possorted_genome_bam.bam",
    log: "log/{sample}_cellranger.log"
    benchmark: "benchmark/{sample}_cellranger.benchmark"
    params:
        transcriptome = "$HOME/cellranger/refdata-gex-GRCh38-2020-A",
        cellranger_path = "$HOME/cellranger/cellranger-10.0.0/bin/cellranger",
        outdir = "{sample}_cellranger",
        sample = lambda wildcards: config['samples'][wildcards.sample]
    threads: 8
    resources:
        mem_mb = 64 * 1024,
        mem_gb = 64 # keep up to date with mem_mb
    shell:
        """
            rm -rf {params.outdir}

            {params.cellranger_path} count --id={wildcards.sample} \
                --transcriptome={params.transcriptome} \
                --fastqs={input.fq_dir} \
                --sample={params.sample} \
                --create-bam=true \
                --localcores={threads} \
                --localmem={resources.mem_gb} \
                --output-dir={params.outdir} &> {log}
        """
