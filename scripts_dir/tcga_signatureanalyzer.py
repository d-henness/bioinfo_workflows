import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('sig_tsv', help = 'path to tsv with signatures')
parser.add_argument('tcga_data', help = 'path to tcda data')
parser.add_argument('new_file', help = 'path to new file to make')
args = parser.parse_args()

fields = [
        'patient.bcr_patient_barcode',
        'patient.age_at_initial_pathologic_diagnosis',
        'patient.gender',
        'patient.melanoma_origin_skin_anatomic_site',
        'patient.race_list.race',
        'patient.samples.sample-2.sample_type',
        'patient.samples.sample.sample_type',
        'patient.sites_of_primary_melanomas.site.tumor_tissue_site',
        'patient.stage_event.pathologic_stage',
        'patient.stage_event.tnm_categories.pathologic_categories.pathologic_m',
        'patient.stage_event.tnm_categories.pathologic_categories.pathologic_n',
        'patient.stage_event.tnm_categories.pathologic_categories.pathologic_t',
        'patient.tissue_source_site',
        'patient.tumor_tissue_site',
        'patient.vital_status',
]



kept_data = {}
with open(args.tcga_data, 'r') as data:
    for line in data:
        split_line = line.strip().split('\t')
        if split_line[0] in fields:
            kept_data[split_line[0]] = split_line[1:]

new_lines = []
with open(args.sig_tsv, 'r') as data:
    for i, line in enumerate(data):
        split_line = line.strip().split('\t')
        signatures = '\t'.join(split_line[1:])
        if i == 0:
            wanted_info = '\t'.join(fields)
            new_line = f'{split_line[0]}\t{wanted_info}\t{signatures}\n'
            new_lines.append(new_line)
        else:
            parsed_sample = '-'.join(split_line[0].lower().split('_')[:3])
            patient_index = kept_data['patient.bcr_patient_barcode'].index(parsed_sample)
            wanted_info = '\t'.join(kept_data[field][patient_index] for field in fields)
            new_line = f'{split_line[0].lower().replace("_", "-")}\t{wanted_info}\t{signatures}\n'
            new_lines.append(new_line)

with open(args.new_file, 'w') as data:
    data.writelines(new_lines)
