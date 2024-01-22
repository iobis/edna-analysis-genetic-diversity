#! /usr/bin/env python3
# Script written by Chandra Earl

from src.functions import *


if __name__ == "__main__":
    raw_data_folder = "raw_data"
    
    downloaddata("https://obis-edna-results.s3.amazonaws.com/output.zip", raw_data_folder)
    
    input_data = importdata(raw_data_folder + "/output")

    cleaned_input_data = cleaninputdata(input_data)

    #fasta_list = buildfastas(cleaned_input_data)
    #aligned_fasta_list = alignseqs(fasta_list)
    #tree_list = createphylogeny(aligned_fasta_list)
    tree_list = ["MiFish_MiMammal_aligned.tre", "mlCOIintF_aligned.tre", "teleo_F_L1848_aligned.tre", "Vert-16S-eDNA-F1_aligned.tre"]
    calculatePD(tree_list, cleaned_input_data)