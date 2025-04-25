import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('count_file')
    parser.add_argument('mapping_file')
    parser.add_argument('count_output_file')
    parser.add_argument('phenotype_output_file')
    args = parser.parse_args()

    mapping_input = pd.read_csv(args.mapping_file)

    conditions = mapping_input["condition"].unique()
    conditions_map = {condition: i for i, condition in enumerate(conditions)}

    mapping_input["condition_map"] = mapping_input["condition"].map(conditions_map)
    mapping_input = mapping_input.sort_values(by = "condition_map")
    new_sample_order = mapping_input["sample"]

    deseq2_input = pd.read_csv(args.count_file, sep = '\t')
    nas = ["na" for _ in range(len(deseq2_input))]
    deseq2_input["DESCRIPTION"] = nas

    deseq2_input = deseq2_input[["DESCRIPTION"] + list(new_sample_order)]
    deseq2_input.index.name = "NAME"
    deseq2_input = deseq2_input[pd.notna(deseq2_input.index)]
    deseq2_input.to_csv(args.count_output_file, sep = '\t')

    with open(args.phenotype_output_file, 'w') as data:
        data.write(f"{len(mapping_input)} {len(conditions)} 1\n")
        condition_string = " ".join(conditions)
        data.write(f"# {condition_string}\n")
        condition_string = " ".join(mapping_input["condition"])
        data.write(f"{condition_string}\n")







if __name__ == '__main__':
    main()



