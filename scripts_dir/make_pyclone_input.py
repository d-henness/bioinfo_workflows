import sys

vcf_snp = sys.argv[1]
vcf_indel = sys.argv[2]
cna = sys.argv[3]

vcf_data = {}
with open(vcf_snp, 'r') as data:
    for line in data:
        if line[:len('chr')] == 'chr':
            id_string = f"{line.split()[0]}_{line.split()[1]}"
            ref_depth, alt_depth = line.split()[9].split(':')[1].split(',')
            vcf_data[id_string] = [ref_depth, alt_depth]
with open(vcf_indel, 'r') as data:
    for line in data:
        if line[:len('chr')] == 'chr':
            id_string = f"{line.split()[0]}_{line.split()[1]}"
            ref_depth, alt_depth = line.split()[9].split(':')[1].split(',')
            vcf_data[id_string] = [ref_depth, alt_depth]

cna_data = {
        'chr1':{},
        'chr2':{},
        'chr3':{},
        'chr4':{},
        'chr5':{},
        'chr6':{},
        'chr7':{},
        'chr8':{},
        'chr9':{},
        'chr10':{},
        'chr11':{},
        'chr12':{},
        'chr13':{},
        'chr14':{},
        'chr15':{},
        'chr16':{},
        'chr17':{},
        'chr18':{},
        'chr19':{},
        'chr20':{},
        'chr21':{},
        'chr22':{},
        'chrX':{},
        'chrY':{},
        }
with open(cna, 'r') as data:
    for i, line in enumerate(data):
        if i == 0:
            fields = line.strip().split()
            minor_cn_field = fields.index('MinorCN')
            major_cn_field = fields.index('MajorCN')
        if i > 0:
            split_line = line.split()
            cna_data[split_line[1]][int(split_line[2])] = [int(split_line[3]), split_line[minor_cn_field], split_line[major_cn_field]]

print(f"mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn")
for mut in vcf_data:
    chrom, pos = mut.split('_')[0], int(mut.split('_')[1])
    if (chrom != 'chrX') and (chrom != 'chrY'):
        for region in sorted(cna_data[chrom]):
            if (pos >= region) and (pos <= cna_data[chrom][region][0]) and (cna_data[chrom][region][2] != '0'):
                print(f"{mut}\t{vcf_data[mut][0]}\t{vcf_data[mut][1]}\t2\t{cna_data[chrom][region][1]}\t{cna_data[chrom][region][2]}")

