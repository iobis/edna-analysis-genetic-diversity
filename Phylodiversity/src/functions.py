#Script written by Chandra Earl

import os, sys
import urllib.request
from zipfile import ZipFile
import progressbar
import pandas as pd
import glob
import shutil
from Bio import Phylo
import csv

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

        combined_data = pd.concat([combined_data, mergeddata], ignore_index=True)

    return combined_data

def cleaninputdata(input_data):
    clean_input_data = input_data[(~input_data['kingdom'].isin(['Bacteria', 'Fungi', 'undef_Bacteria_bacteria', 'Archaea']) & \
                                   ~input_data['phylum'].isin(['Acidobacteria', 'Bacteria incertae sedis', 'Cyanobacteria', 'Deferribacteres', 'Firmicutes', 'Proteobacteria', 'Verrucomicrobia', 'Tenericutes', 'Spirochaetes']) & \
                                   ~input_data['scientificName'].isin(['Bacteria'])
    )]
    
    #maybe remove incertae sedis
    return clean_input_data


def remove_duplicate_seqs(df):

    df[['asv_number', 'EE_number', 'rest']] = df['occurrenceID'].str.split('_', n=2, expand=True)
    df_grouped = df.groupby(['asv_number', 'DNA_sequence'])['EE_number'].agg(lambda x: '_'.join(x)).reset_index()
    result_df = pd.merge(df, df_grouped, on=['asv_number', 'DNA_sequence'])
    return result_df

def write_fasta(file_name, df):
    current_asv = ''
    current_EE = ''
    os.makedirs("outputs/unaligned_fasta", exist_ok=True)
    with open("outputs/unaligned_fasta/"+file_name, 'a') as f:
        for index, row in df.iterrows():
            if current_asv == row['asv_number'] and current_EE == row['EE_number_y']:
                continue
            fasta_identifier = f">{row['asv_number']}_{row['EE_number_y']}_{row['rest']}_{row['kingdom']}_{row['phylum']}_{row['class']}_{row['order']}_{row['family']}_{row['scientificName']}\n"
            fasta_identifier = fasta_identifier.replace(" ", "_").replace("[", "").replace("]", "")
            sequence = row['DNA_sequence'] + '\n'
            f.write(fasta_identifier)
            f.write(sequence)
            current_asv = row['asv_number']
            current_EE = row['EE_number_y']
	
def buildfastas(cleaned_input_data):
    if os.path.exists("outputs"):
        shutil.rmtree("outputs")
        
    outputfiles = []
    primer_names = cleaned_input_data['pcr_primer_name_forward'].unique()
    
    for primer_name in primer_names:
        primer_df = cleaned_input_data[cleaned_input_data['pcr_primer_name_forward'] == primer_name]
        primer_df = primer_df.sort_values('occurrenceID')
        
        primer_df = remove_duplicate_seqs(primer_df)
        if primer_name == "MiFish-UE-F" or primer_name == "MiMammal-UEB-F":
            file_name = 'MiFish_MiMammal.fa'
        else:
            file_name = f'{primer_name}.fa'

        write_fasta(file_name, primer_df)
        
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
    mafft_command = "fasttree -nt -gtr -fastest {} > {}"
    
    for input_file in aligned_fasta_list:
        output_file = f'{os.path.splitext(os.path.basename(input_file))[0]}.tre'
        os.system(mafft_command.format("outputs/aligned_fasta/"+input_file, "outputs/fasttree/"+output_file))	
        outputfiles.append(output_file)
    return outputfiles    
    
def faith_pd_subset(tree, subset_leaves):
    pd_value = 0
    errors = 0
    for node in subset_leaves:
        if tree.find_clades({"name": node}):
            try:
                pd_value += tree.distance(node)
            except:
                try:
                    pd_value += tree.distance("_R_" + node) 
                except:
                    print("error: "+node)
                    errors += 1
    return pd_value, errors

def get_tip_list(tree_file, site, cleaned_input_data):
    #gets list of tips that are from a site for a specific primer
    site_data = cleaned_input_data[cleaned_input_data['higherGeography'] == site]
    if tree_file.replace("_aligned.tre", "") == "MiFish_MiMammal":
        
        site_primer_data = site_data[
            (site_data['pcr_primer_name_forward'] == 'MiFish-UE-F') | 
            (site_data['pcr_primer_name_forward'] == 'MiMammal-UEB-F')
        ]
    else:
        primer = tree_file.replace("_aligned.tre", "")
        site_primer_data = site_data[site_data['pcr_primer_name_forward'] == primer]
    columns_to_concatenate = ['occurrenceID', 'kingdom', 'phylum', 'class', 'order', 'family', 'scientificName']
    tip_list = site_primer_data.apply(lambda row: '_'.join(map(str, row[columns_to_concatenate])).replace(' ', '_').replace("[", "").replace("]", ""), axis=1).tolist()
    return tip_list

    
def calculatePD(input_tree_files, cleaned_input_data):
    site_names = cleaned_input_data['higherGeography'].unique()
    site_names = site_names[~pd.isna(site_names)]
    
    with open("PD_raw_out.csv", 'w') as csv_file:
        with open("PD_perspecies_out.csv", 'w') as csv_file2:
            csv_writer = csv.writer(csv_file, delimiter=',')
            csv_writer.writerow([''] + input_tree_files)
            
            csv_writer2 = csv.writer(csv_file2, delimiter=',')
            csv_writer2.writerow([''] + input_tree_files)
                
            for site in site_names:
                row_data = [site]
                row_data2 = [site]
                for tree_file in input_tree_files:
                    print(tree_file)
                    subset_of_leaves = get_tip_list(tree_file, site, cleaned_input_data)
                    phylo_tree = Phylo.read("outputs/fasttree/"+tree_file, 'newick')
                    pd_value_subset, errors = faith_pd_subset(phylo_tree, subset_of_leaves)
                    row_data.append(str(pd_value_subset))
                    row_data2.append(str(pd_value_subset/len(subset_of_leaves)))
                    print(f'site: {site}, tree/locus: {tree_file}, PD: {pd_value_subset}, errors: {errors}')
                csv_writer.writerow(row_data)
                csv_writer2.writerow(row_data2)
            
    
    
    
    