#Script written by Chandra Earl

import os, sys
import urllib.request
from zipfile import ZipFile
import progressbar
import pandas as pd
import glob
import shutil
from Bio import Phylo

pbar = None

def show_progress(block_num, block_size, total_size):
    global pbar
    if pbar is None:
        pbar = progressbar.ProgressBar(maxval=total_size)
        pbar.start()

    downloaded = block_num * block_size
    if downloaded < total_size:
        pbar.update(downloaded)
    else:
        pbar.finish()
        pbar = None

def downloaddata(url, destination_folder="."):
    # Create destination folder if it doesn't exist
    os.makedirs(destination_folder, exist_ok=True)

    # Set the paths for the zip file and the destination folder
    zip_file_path = os.path.join(destination_folder, "output.zip")
    output_folder_path = os.path.join(destination_folder, "output")

    # Check if the unzipped files already exist
    if os.path.exists(output_folder_path):
        print("Files already exist. Skipping download and unzip.")
        return

    # Check if the zip file already exists
    if not os.path.exists(zip_file_path):
        # Download the file using urllib if it doesn't exist
        print("Downloading raw data.")
        urllib.request.urlretrieve(url, zip_file_path, show_progress)

    # Unzip the downloaded file
    with ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(destination_folder)

    # Remove the zip file after extraction
    os.remove(zip_file_path)

def importdata(folder_path):
    # Get a list of all files in the folder
    pattern = os.path.join(folder_path, '*_DNADerivedData.tsv')
    files = glob.glob(pattern)
    
    # Initialize an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Iterate over the files in the folder
    for file_name in files:
        site_name = os.path.basename(file_name).replace("_DNADerivedData.tsv", "")
        print("Importing " + site_name)
        
		#get Occurrence file
        file_path = os.path.join(folder_path, site_name + "_Occurrence.tsv")

        #read data
        DNADeriveddata = pd.read_csv(file_name, sep='\t', low_memory=False)
        Occurrencedata = pd.read_csv(file_path, sep='\t', low_memory=False)
        mergeddata = pd.merge(DNADeriveddata, Occurrencedata, on='occurrenceID', how='inner')
        

		# Combine the data based on the 'occurrenceID' field
        combined_data = pd.concat([combined_data, mergeddata], ignore_index=True)

    return combined_data

def cleaninputdata(input_data):
    clean_input_data = input_data[(~input_data['kingdom'].isin(['Bacteria', 'Fungi', 'undef_Bacteria_bacteria', 'Archaea']) & \
                            ~input_data['phylum'].isin(['Acidobacteria', 'Bacillariophyta', 'Bacteria incertae sedis', 'Cyanobacteria', 'Deferribacteres', 'Firmicutes', 'Proteobacteria', 'Verrucomicrobia', 'Tenericutes', 'Spirochaetes']) & \
                        ~input_data['scientificName'].isin(['Bacteria'])
    )]
    
    #maybe remove incertae sedis
    
    return clean_input_data


def remove_duplicate_seqs(df):
    df['asv_pattern'] = df['occurrenceID'].str.split("_", n=1).str[0]
    df_filtered = df.drop_duplicates(subset='asv_pattern', keep='first')
    return df_filtered

def write_fasta(primer_name, file_name, df):
    os.makedirs("outputs/unaligned_fasta", exist_ok=True)
    with open("outputs/unaligned_fasta/"+file_name, 'a') as f:
        for index, row in df.iterrows():
            fasta_identifier = f">{row['asv_pattern']}_{row['kingdom']}_{row['phylum']}_{row['class']}_{row['order']}_{row['family']}_{row['scientificName']}_{primer_name}\n"
            fasta_identifier = fasta_identifier.replace(" ", "_")
            sequence = row['DNA_sequence'] + '\n'
            f.write(fasta_identifier)
            f.write(sequence)
	
def buildfastas(cleaned_input_data):
    if os.path.exists("outputs"):
        shutil.rmtree("outputs")
        
    outputfiles = []
    primer_names = cleaned_input_data['pcr_primer_name_forward'].unique()
    
    for primer_name in primer_names:
        primer_df = cleaned_input_data[cleaned_input_data['pcr_primer_name_forward'] == primer_name]
        primer_df_asv = remove_duplicate_seqs(primer_df)
        print(primer_df_asv)
        if primer_name == "MiFish-UE-F" or primer_name == "MiMammal-UEB-F":
            file_name = 'MiFish_MiMammal.fa'
        else:
            file_name = f'{primer_name}.fa'
        
        write_fasta(primer_name, file_name, primer_df_asv)

        print(f"Fasta file '{file_name}' created.")
        outputfiles.append(file_name)
        
    return outputfiles
        
def alignseqs(fasta_list):
    outputfiles = []
    os.makedirs("outputs/aligned_fasta", exist_ok=True)
    mafft_command = "mafft --thread -1 --adjustdirection --auto {} > {}"
    
    for input_file in fasta_list:
        output_file = f'{os.path.splitext(os.path.basename(input_file))[0]}_aligned.fasta'
        os.system(mafft_command.format("outputs/unaligned_fasta/"+input_file, "outputs/aligned_fasta/"+output_file))
        outputfiles.append(output_file)
    return outputfiles

def createphylogeny(aligned_fasta_list):
    outputfiles = []
    os.makedirs("outputs/fasttree", exist_ok=True)
    mafft_command = "fasttree -nt {} > {}"
    
    for input_file in aligned_fasta_list:
        output_file = f'{os.path.splitext(os.path.basename(input_file))[0]}.tre'
        os.system(mafft_command.format("outputs/aligned_fasta/"+input_file, "outputs/fasttree/"+output_file))	
        outputfiles.append(output_file)
    return outputfiles    
    
def faith_pd_subset(tree, subset_leaves):
    """
    Calculate Faith's Phylogenetic Diversity (PD) for a subset of leaves in a given tree.

    Parameters:
    - tree (Bio.Phylo.BaseTree.Tree): Phylogenetic tree object.
    - subset_leaves (list): List of leaf labels for the subset.

    Returns:
    - float: Faith's PD value for the subset.
    """
    subset_leaves = set(subset_leaves)
    subset_tree = tree.copy()
    
    # Prune the tree to include only the subset of leaves
    subset_tree.prune(subset_leaves)
    
    # Calculate Faith's PD for the subset
    pd_value = sum(subset_tree.distance(node) for node in subset_tree.get_terminals())
    return pd_value

def get_tip_list(site, cleaned_input_data):
    site_data = cleaned_input_data[cleaned_input_data['pcr_primer_name_forward'] == site]
    print(site_data)
    site_data['asv_pattern'] = site_data['occurrenceID'].str.split("_", n=1).str[0]
    columns_to_concatenate = ['asv_pattern', 'kingdom', 'phylum', 'class', 'order', 'family', 'scientificName', 'pcr_primer_name_forward']
    tip_list = site_data.apply(lambda row: '_'.join(map(str, row[columns_to_concatenate])), axis=1).tolist()
    print(tip_list)

    
def calculatePD(input_tree_files, cleaned_input_data):
    site_names = cleaned_input_data['higherGeography'].unique()
    
    for tree_file in input_tree_files:
        for site in site_names:
            get_tip_list(site, cleaned_input_data)
        
    
    
    subset_of_leaves = ['leaf1', 'leaf2', 'leaf3']  # Replace with the labels of your subset
    
    # Read the phylogenetic tree from a Newick file
    phylo_tree = Phylo.read(tree_file, 'newick')
    
    # Calculate Faith's PD for the subset of leaves
    pd_value_subset = faith_pd_subset(phylo_tree, subset_of_leaves)