import yaml
import argparse
import os

parser = argparse.ArgumentParser(description = 'make sample_annotation.tsv for DROP')
parser.add_argument('config', help = 'yaml file used for muctect2 and rsem')
args = parser.parse_args()

with open(args.config, 'r') as data:
    yaml_data = yaml.safe_load(data)

with open('sample_annotation.tsv', 'w') as data:
    data.write('RNA_ID\tRNA_BAM_FILE\tDNA_VCF_FILE\tDNA_ID\tDROP_GROUP\tPAIRED_END\tCOUNT_MODE\tCOUNT_OVERLAPS\tSTRAND\tANNOTATION\tGENE_COUNTS_FILE\n')
    for sample in yaml_data['rna2dna']:
        vcf_loc = os.getcwd() + f'/GATK_runs/{yaml_data["rna2dna"][sample]}/FilterByOrientationBias/{yaml_data["rna2dna"][sample]}.vcf'
        rna_bam_loc = os.getcwd() + f'/rsem/{sample}/{sample}.transcript-sorted.bam'
        assert(os.path.exists(vcf_loc)), f'could not find {rna_bam_loc}'
        assert(os.path.exists(rna_bam_loc)), f'could not find {rna_bam_loc}'
        assert(yaml_data['stranded'][sample] is not None), f'specify standedness of {sample} in {args.config}'
        assert(yaml_data['stage'][sample] is not None), f'specify stage of {sample} in {args.config}'
        data.write('\t'.join([sample, rna_bam_loc, vcf_loc, yaml_data['rna2dna'][sample], yaml_data['stage'][sample], 'True', yaml_data['stranded'][sample]]) + '\n')
