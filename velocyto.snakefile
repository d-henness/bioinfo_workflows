rule velocyto_all:
    input: expand("velocyto_output_dir/{sample}_cellranger/{sample}.loom", sample = config['samples'])

rule velocyto:
    input:
        bam = "{sample}_cellranger/outs/possorted_genome_bam.bam",
        barcodes = "{sample}_cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        loom = "velocyto_output_dir/{sample}_cellranger/{sample}.loom",
        cellsorted_bam = "{sample}_cellranger/outs/cellsorted_possorted_genome_bam.bam"
    log: "log/{sample}_velocyto.log"
    benchmark: "benchmark/{sample}_velocyto.benchmark"
    params:
        gtf = "/project/def-rgni0001/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
        rmsk = "/project/def-rgni0001/ref_data/hg38/hg38_rmsk.gtf",
        container = "/project/def-rgni0001/dhenness/containers/velocyto.sif",
        outdir = "velocyto_output_dir/{sample}_cellranger",
        sample_id = "{sample}"
    threads: 4
    resources:
        mem_mb = 128 * 1024 * 2,
        mem_gb = 128 * 2, # keep up to date with mem_mb
        runtime = 12 * 60 # minutes
    shell:
        """
            module load apptainer
            mkdir -p {params.outdir}
            apptainer exec \
                --bind /project/def-rgni0001/ \
                --bind /etc/pki/tls/certs/ca-bundle.crt:/etc/pki/tls/certs/ca-bundle.crt:ro \
                --bind /etc/ssl/certs:/etc/ssl/certs:ro \
                --bind /home/dhenness/.cache:/home/dhenness/.cache \
                --bind /etc/localtime:/etc/localtime:ro \
                --bind /usr/share/zoneinfo:/usr/share/zoneinfo:ro \
                --env TZ=America/Edmonton \
                --env LANG=C.UTF-8,LC_ALL=C.UTF-8 \
                --env BIOCFILECACHE_CACHE=/home/dhenness/.cache/R/BiocFileCache \
                --env SSL_CERT_FILE=/etc/pki/tls/certs/ca-bundle.crt \
                --env CURL_CA_BUNDLE=/etc/pki/tls/certs/ca-bundle.crt \
                {params.container} \
                    velocyto run \
                        -b {input.barcodes} \
                        -o {params.outdir} \
                        -m {params.rmsk} \
                        -e {params.sample_id} \
                        -@ {threads} \
                        {input.bam} \
                        {params.gtf} &> {log}
        """
