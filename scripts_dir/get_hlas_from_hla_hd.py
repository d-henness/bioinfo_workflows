import sys
import re

type_1_hlas = ['A', 'B', 'C']
type_2_hlas = ['DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 'DRB3', 'DRB4', 'DRB5']

def clean_hla1_allele_helper(allele):
    allele_sub_type_list = allele.split(':')           
    allele_cleaned = ':'.join(allele_sub_type_list[:min(len(allele_sub_type_list), 2)])
    allele_cleaned = re.sub(r'\*', '', allele_cleaned)
    return allele_cleaned

def clean_hla1_allele(allele1, allele2):
    if allele1 == 'Not typed':
        return
    else:
        allele1_cleaned = clean_hla1_allele_helper(allele1)

    if allele2 != '-':
        allele2_cleaned = clean_hla1_allele_helper(allele2)
        return f'{allele1_cleaned},{allele2_cleaned}'

    return allele1_cleaned

def clean_hla2_allele_alpha_beta_helper(allele):
    if allele == '-' or allele == 'Not typed':
        return
    cleaned_allele = re.sub(r'HLA-', '', allele)
    cleaned_allele = re.sub(r'\*', '', cleaned_allele)
    cleaned_allele = ''.join(cleaned_allele.split(':')[:2])
    return cleaned_allele

def clean_hla2_allele_alpha_beta(alpha_line, beta_line):
    alpha_line_split = alpha_line.strip().split('\t')
    alpha_alleles = []
    for allele in alpha_line_split[1:]:
        alpha_alleles.append(clean_hla2_allele_alpha_beta_helper(allele))

    beta_line_split = beta_line.strip().split('\t')
    beta_alleles = []
    for allele in beta_line_split[1:]:
        beta_alleles.append(clean_hla2_allele_alpha_beta_helper(allele))

    allele_list = []
    for alpha in alpha_alleles:
        if alpha is not None:
            for beta in beta_alleles:
                if beta is not None:
                    allele_list.append(f'HLA-{alpha}-{beta}')

    return allele_list

def clean_hla2_alleles_DRB_helper(allele):
    if allele == '-' or allele == 'Not typed':
        return
    cleaned_allele = re.sub(r'HLA-', '', allele)
    cleaned_allele = re.sub(r'\*', '_', cleaned_allele)
    cleaned_allele = ''.join(cleaned_allele.split(':')[:2])
    return cleaned_allele

def clean_hla2_alleles_DRB(line):
    split_line = line.strip().split('\t')

    allele_list = []
    for allele in split_line[1:]:
        allele_list.append(clean_hla2_alleles_DRB_helper(allele))

    allele_list = [allele for allele in allele_list if allele is not None]

    return allele_list


def clean_hla2_alleles(hla2_dict):
    dqa_line = hla2_dict['DQA1']
    dqb_line = hla2_dict['DQB1']
    dq_list = clean_hla2_allele_alpha_beta(dqa_line, dqb_line)

    dpa_line = hla2_dict['DPA1']
    dpb_line = hla2_dict['DPB1']
    dp_list = clean_hla2_allele_alpha_beta(dpa_line, dpb_line)

    drb1 = clean_hla2_alleles_DRB(hla2_dict['DRB1'])
    drb3 = clean_hla2_alleles_DRB(hla2_dict['DRB3'])
    drb4 = clean_hla2_alleles_DRB(hla2_dict['DRB4'])
    drb5 = clean_hla2_alleles_DRB(hla2_dict['DRB5'])

    return dq_list + dp_list + drb1 + drb3 + drb4 + drb5
    

filenm = sys.argv[1]
hla1_file = sys.argv[2]
hla2_file = sys.argv[3]

hla1_dict = {}
hla2_dict = {}

with open(filenm, 'r') as data:
    for line in data:
        split_line = line.strip().split('\t')
        locus = split_line[0]
        allele1 = split_line[1]
        allele2 = split_line[2]

        if locus in type_1_hlas:
            allele_string = clean_hla1_allele(allele1, allele2)
            if allele_string != None:
                hla1_dict[locus] = allele_string

        elif locus in type_2_hlas:
            hla2_dict[locus] = line

hla1_string = ','.join(hla1_dict[locus] for locus in hla1_dict.keys() if locus in type_1_hlas)

with open(hla1_file, 'w') as data:
    data.write(hla1_string)

hla2_list = clean_hla2_alleles(hla2_dict)
hla2_string = ','.join(hla2_list)

with open(hla2_file, 'w') as data:
    data.write(hla2_string)
