import os
import csv

def split_vcf_files(vcf_folder, csv_folder, output_folder):
    vcf_files = os.listdir(vcf_folder)

    for vcf_file in vcf_files:
        if vcf_file.endswith(".vcf"):
            sample_id = vcf_file.split(".")[0]

            # Find matching CSV files
            clonal_csv_file = os.path.join(csv_folder, sample_id + ".clonal.csv")
            subclonal_csv_file = os.path.join(csv_folder, sample_id + ".subclonal.csv")

            # Check if both CSV files exist
            if os.path.isfile(clonal_csv_file) and os.path.isfile(subclonal_csv_file):
                # Create output file paths
                clonal_vcf_file = os.path.join(output_folder, sample_id + ".clonal.vcf")
                subclonal_vcf_file = os.path.join(output_folder, sample_id + ".subclonal.vcf")

                # Read mutation positions from CSV files
                clonal_positions = set()
                subclonal_positions = set()

                with open(clonal_csv_file, "r") as clonal_csv:
                    clonal_reader = csv.reader(clonal_csv)
                    header = next(clonal_reader)
                    if len(header) == 2 and header[0] == "CHR" and header[1] == "POS":
                        for row in clonal_reader:
                            clonal_positions.add((row[0], row[1]))

                with open(subclonal_csv_file, "r") as subclonal_csv:
                    subclonal_reader = csv.reader(subclonal_csv)
                    header = next(subclonal_reader)
                    if len(header) == 2 and header[0] == "CHR" and header[1] == "POS":
                        for row in subclonal_reader:
                            subclonal_positions.add((row[0], row[1]))

                # Split VCF file based on mutation positions
                with open(os.path.join(vcf_folder, vcf_file), "r") as vcf:
                    with open(clonal_vcf_file, "w") as clonal_vcf:
                        with open(subclonal_vcf_file, "w") as subclonal_vcf:
                            for line in vcf:
                                if line.startswith("#"):
                                    clonal_vcf.write(line)
                                    subclonal_vcf.write(line)
                                else:
                                    fields = line.split("\t")
                                    position = (fields[0], fields[1])
                                    if position in clonal_positions:
                                        clonal_vcf.write(line)
                                    elif position in subclonal_positions:
                                        subclonal_vcf.write(line)

# Specify the paths
vcf_folder = "/home/alhendi/gamble-lung_precancer/Lukas_WES/Summaries/CONIPHER/combined_multisite/mutational_sig/VCF_multiSite/"
csv_folder = "/home/alhendi/gamble-lung_precancer/Lukas_WES/Summaries/CONIPHER/combined_multisite/mutational_sig/clonality_to_use_for_split/"
output_folder = "/home/alhendi/gamble-lung_precancer/Lukas_WES/Summaries/CONIPHER/combined_multisite/mutational_sig/VCF_multiSite_split_by_clonality/"

# Create the output directory if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Call the function
split_vcf_files(vcf_folder, csv_folder, output_folder)
