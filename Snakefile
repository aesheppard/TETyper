import json
import sys
import os
import pandas as pd
import subprocess

def load_config(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Check that number of arguments is the same.
def check_args(sample_data):
    data = [len(value) if value is not None else 0 for key, value in sample_data.items()]
    if (data[0] == data[1] and data[0] == data[2]) or (data[0] == data[3]):
        return 0
    return 1

# Read csv file and convert to df.
def convert_dataframe(file_path):
    df = pd.read_csv(file_path, sep='\t')
    columns = ["outprefix", "fq1", "fq2", "bam"]
    output_data = {"--" + col: df[col].tolist() for col in columns}
    for col in columns[1:]:
        output_data["--" + col] = [value if pd.notna(value) else None for value in output_data["--" + col]]
    return output_data
    lengths = [len(sample_data[col]) for col in sample_data]
    return all(length == lengths[0] for length in lengths)

# Generate dict that contains all arguments for TETyper runs.
def generate_parsing_dict(sample_data):
    sample_name_list = sample_data["--outprefix"]
    parsing_dict = {}
    for name in sample_name_list:
        indices = [i for i, val in enumerate(sample_name_list) if val == name]
        return_string = ""
        for index in indices:
            for key, val in sample_data.items():
                if val is not None and val[index] is not None:
                    return_string += f"{key} {val[index]} "
        parsing_dict[name] = return_string.strip()
    return parsing_dict

file = "config.json"
file_path = "sample_data.txt"
data = load_config(file)

if os.path.isfile(file_path):
    sample_data = convert_dataframe(file_path)
    data[0] = sample_data
sample_data = data[0]
arg_data = data[1]

if check_args(sample_data):
    print("Error: number of sample names, and fq1/2 or bam files must be the same.")
    sys.exit(1)
parsing_dict = generate_parsing_dict(sample_data)

# Generates summary file after TETyper has finished.
rule all:
    input:
        all_sums = ["{sample_name}_summary.txt".format(sample_name = sample_name) for sample_name in sample_data["--outprefix"]],
        spades_output = "correction.txt"
    params:
        summary_args = sample_data["--outprefix"]
    output:
        total_sum = "all_summary.txt"
    shell:
        "./generate_summary.sh {params.summary_args}"

# Ensure spades has been fixed.
rule correct_spades:
    output:
        spades_output = "correction.txt"
    shell:
        """
        echo "Please run spades_corrector.sh using 'source spades_corrector.sh. If this issue persists, see README for fix.'"
        exit 1
        """

# Run TETyper for each cluster of arguments in generated dict.
rule run_TETyper:
    input:
        script = "TETyper.py"
    params:
        sample = lambda wildcards: parsing_dict[wildcards.sample_name],
        args = [f" {key} {value}" for key, value in arg_data.items() if value is not None]
    output:
        summary = "{sample_name}_summary.txt",
    shell:
        "python3 {input.script} {params.sample} {params.args}" 