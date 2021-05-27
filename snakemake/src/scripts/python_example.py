# I first set all variables
output_file = snakemake.output.python_output
input_file = snakemake.input[0]
samples = snakemake.params.sample_list

with open(output_file, "w") as out_file:
    with open(input_file, "r") as in_file:
        for line in in_file:
           out_file.write(line)
        for sample in samples:
            out_file.write(sample + "\n")