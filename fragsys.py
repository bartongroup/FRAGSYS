### IMPORTS ###

import os
import re
import Bio
import math
import scipy
import shutil
import Bio.SeqIO
import importlib
import prointvar
import statistics
import subprocess
import Bio.AlignIO
import numpy as np
import configparser
import pandas as pd
import seaborn as sns
from Bio.Seq import Seq
from Bio import pairwise2
from scipy import cluster
import varalign.alignments
import scipy.stats as stats
import varalign.align_variants
import matplotlib.pyplot as plt
from prointvar.pdbx import PDBXreader
from prointvar.pdbx import PDBXwriter
import matplotlib.patches as mpatches
from prointvar.sifts import SIFTSreader
from prointvar.config import config as cfg
from prointvar.dssp import DSSPrunner, DSSPreader
from scipy.spatial.distance import squareform, pdist
from prointvar.fetchers import download_sifts_from_ebi
from prointvar.fetchers import download_structure_from_pdbe

### DICTIONARIES AND LISTS ###

arpeggio_suffixes = [
        "atomtypes", "bs_contacts", "contacts", "specific.sift", 
        "sift","specific.siftmatch", "siftmatch", "specific.polarmatch",
        "polarmatch", "ri", "rings", "ari", "amri", "amam", "residue_sifts"
]

pdb_clean_suffixes = ["break_residues", "breaks"]

simple_ions = [
    "ZN", "MN", "CL", "MG", "CD", "NI", "NA", "IOD", "CA", "BR", "XE"
]

acidic_ions = [
    "PO4", "ACT", "SO4", "MLI", "CIT", "ACY", "VO4"
]

non_relevant_ligs_manual = [
    "DMS", "EDO", "HOH", "TRS", "GOL", "OGA", "FMN", "PG4", "PGR",
    "MPD", "TPP", "MES", "PLP", "HYP", "CSO", "UNX", "EPE", "PEG",
    "PGE", "DOD", "SUI"
]

non_relevant = non_relevant_ligs_manual + simple_ions + acidic_ions

pdb_resnames = [
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",       ### ONE OF THESE MIGHT BE ENOUGH ###
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"
]

aas = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',                   ### ONE OF THESE MIGHT BE ENOUGH ###
        'THR', 'TRP', 'TYR', 'VAL', 'GLX', 'GLI', 'NLE', 'CYC'
]

bbone = ["N", "CA", "C", "O"]

aa_code = {
        "ALA" : 'A', "CYS" : 'C', "ASP" : 'D', "GLU" : 'E',
        "PHE" : 'F', "GLY" : 'G', "HIS" : 'H', "ILE" : 'I',
        "LYS" : 'K', "LEU" : 'L', "MET" : 'M', "ASN" : 'N',
        "PRO" : 'P', "GLN" : 'Q', "ARG" : 'R', "SER" : 'S',
        "THR" : 'T', "VAL" : 'V', "TRP" : 'W', "TYR" : 'Y',
        "PYL" : 'O', "SEC" : 'U', "HYP" : 'P', "CSO" : 'C', # WEIRD ONES
        "SUI" : 'D',
}

### CONFIG FILE READING AND VARIABLE SAVING ###

config = configparser.ConfigParser()
config.read("fragsys_config.txt")

dssp_bin = config["binaries"].get("dssp_bin")
stamp_bin = config["binaries"].get("stamp_bin")
transform_bin = config["binaries"].get("transform_bin")
clean_pdb_python_bin = config["binaries"].get("clean_pdb_python_bin")
clean_pdb_bin = config["binaries"].get("clean_pdb_bin")
arpeggio_python_bin = config["binaries"].get("arpeggio_python_bin")
arpeggio_bin = config["binaries"].get("arpeggio_bin")
gnomad_vcf = config["dbs"].get("gnomad_vcf")
swissprot = config["dbs"].get("swissprot")
stampdir = config["other"].get("stampdir")
oc_dist = float(config["clustering"].get("oc_dist"))
oc_metric = config["clustering"].get("oc_metric")
oc_method = config["clustering"].get("oc_method")
mes_sig_t = float(config["other"].get("mes_sig_t"))

### PRE-PROCESSING DATA CODE ###

def setup_fragsys_start(main_dir, prot, df):
    """
    given a main directory, protein accession and structures dataframe,
    downloads all biological assemblies of such structures if they have
    not already been downloaded and classified in subdirectories
    """
    wd = os.path.join(main_dir, prot)
    unsupp_cifs_dir = os.path.join(wd, "unsupp_cifs")
    if not os.path.isdir(wd):
        os.mkdir(wd)
    else:
        pass
    if not os.path.isdir(unsupp_cifs_dir):
        os.mkdir(unsupp_cifs_dir)
    else:
        pass
    prot_df = df[df.entry_uniprot_accession == prot]
    prot_strucs = prot_df.pdb_id.unique().tolist()
    struc_files = []
    unsupp_cifs_subirs = os.listdir(unsupp_cifs_dir)
    if len(unsupp_cifs_subirs) == 0: #it's empty, need to download
        pass
    else: #not empty
        for element in unsupp_cifs_subirs:
            element_path = os.path.join(unsupp_cifs_dir, element)
            if os.path.isfile(element_path):
                pass
            elif os.path.isdir(element_path):
                element_files = os.listdir(element_path)
                for element_file in element_files:
                    prot_strucs.remove(element_file[:4])
                    unsupp_file_path = os.path.join(unsupp_cifs_dir, element_file)
                    if os.path.isfile(unsupp_file_path):
                        os.remove(unsupp_file_path)
                
    for struc in prot_strucs: #only those that are downloaded or have yet to be (not those classified already)
        input_struc = os.path.join(cfg.db_root, cfg.db_pdbx, struc + "_bio.cif")
        input_struc_moved = os.path.join(unsupp_cifs_dir, os.path.basename(input_struc))
        if os.path.isfile(input_struc) or os.path.isfile(input_struc_moved):
            pass
        else:
            download_structure_from_pdbe(struc, bio = True) # downloads only files missing
        struc_files.append(input_struc)
        
    for file in struc_files:
        if os.path.isfile(file):
            shutil.move(file, os.path.join(unsupp_cifs_dir, os.path.basename(file))) #moves downloaded cif files
            
    pdbx_files = os.listdir(os.path.join(cfg.db_root, cfg.db_pdbx))
    if len(pdbx_files) != 0:
        for pdbx in pdbx_files:
            os.remove(os.path.join(cfg.db_root, cfg.db_pdbx, pdbx)) # removes any leftover files?

def setup_dirs(dirs):
    """
    sets up the directories for the program to tun
    """
    for dirr in dirs:
        if os.path.isdir(dirr):
            continue
        else:
            os.mkdir(dirr)

def classify_pdbs(prot, main_dir, h = 0.2, show = False):
    """
    Classifies a group of cif files into different
    folders according to their amino acid sequence
    """
    z, un_seqs, seqs_df = cluster_protein_sequences(prot, main_dir, h, show)
    
    seqs_df = assign_seqs_clusters(z, un_seqs, seqs_df)
    
    create_cluster_dirs(main_dir, prot, seqs_df)

def cluster_protein_sequences(prot, main_dir, h = 0.2, show = False):
    """
    Extracts sequences from PDB structures of a certain UniProt accession,
    calculates distances between them, genereates matrix, dendrogram and
    clustermap.
    """
    wd = os.path.join(main_dir, prot)
    unsupp_cifs_dir = os.path.join(wd, "unsupp_cifs")
    cif_paths = []
    for file in os.listdir(unsupp_cifs_dir):
        if os.path.isfile(os.path.join(unsupp_cifs_dir, file)): #cifs without classifying
            cif_paths.append(os.path.join(unsupp_cifs_dir, file))
        elif os.path.isdir(os.path.join(unsupp_cifs_dir, file)): #checking for directories, might have classified PDBs already
            for file2 in os.listdir(os.path.join(unsupp_cifs_dir, file)):
                cif_paths.append(os.path.join(unsupp_cifs_dir, file, file2))
        else:
            pass
    seqs_df_path = os.path.join(wd, "{}_seqs_df.csv".format(prot))
    if os.path.isfile(seqs_df_path):
        seqs_df = pd.read_csv(seqs_df_path)
    else:
        seqs = {}
        for cif_path in cif_paths:
            seq = get_seq_from_pdb(cif_path, pdb_fmt = "mmcif")
            cif_id = cif_path.split("/")[-1][:4]
            seqs[cif_id] = seq
        seqs_df = pd.DataFrame(list(seqs.items()), index = None, columns = ["pdb_id", "seq"])
        seqs_df.to_csv(seqs_df_path, index = False)
        
    un_seqs = seqs_df.drop_duplicates("seq").seq.tolist()
    msa_df_out = os.path.join(wd, "{}_msa_dists.csv".format(prot))
    if os.path.isfile(msa_df_out):
        msa_df_norm = pd.read_csv(msa_df_out)
    else:
        msa_mat_norm = get_msa_matrix(un_seqs)
        msa_df_norm = pd.DataFrame(msa_mat_norm)
        msa_df_norm.to_csv(msa_df_out, index = False)
    if len(un_seqs) == 1:
        return np.array([]), un_seqs, seqs_df
    else:
        pass
        
    plot_clustermap(msa_df_norm, prot, wd, show)
    
    z = scipy.cluster.hierarchy.linkage(scipy.spatial.distance.squareform(msa_df_norm), method = "complete")
    
    plot_dendrogram(z, msa_df_norm, wd, prot, h, show)
    
    return z, un_seqs, seqs_df

def get_msa_matrix(seqs_list, norm = True):
    """
    Given a list of sequences, calculates all pairwise global sequence alignments
    and from the scores, generates a distance from 0-1.
    """
    scores = {i: {} for i in range(len(seqs_list))}
    for i in range(len(seqs_list)):
        scores[i][i] = pairwise2.align.globalxx(seqs_list[i], seqs_list[i])[0][2]
        for j in range(i+1, len(seqs_list)):
            scores[i][j] = pairwise2.align.globalxx(seqs_list[i], seqs_list[j])[0][2]
            scores[j][i] = scores[i][j]
    if norm == True:
        scores_norm = {i: {} for i in range(len(seqs_list))}
        for i in range(len(seqs_list)):
            scores_norm[i][i] = abs(scores[i][i] - scores[i][i])/scores[i][i]
            for j in range(i+1, len(seqs_list)):
                scores_norm[i][j] = abs(scores[i][j] - scores[i][i])/scores[i][i]
                scores_norm[j][i] = abs(scores[j][i] - scores[i][i])/scores[i][i]
        return scores_norm
    else:
        return scores

def plot_clustermap(msa_df_norm, prot, wd, show = False):
    """
    Plots a clustermap given a matrix of distances between sequences
    """
    plt.rcParams['figure.dpi'] = 300
    g = sns.clustermap(
        msa_df_norm, #square = True,
        cmap = "viridis_r", linewidths = 1, linecolor = "k",
        vmin = 0, vmax = 1
    )
    
    fig_out = os.path.join(wd, "{}_unique_seqs_clustermap.png".format(prot))
    
    if os.path.isfile(fig_out):
        pass
    
    else:
        plt.savefig(fig_out)
        
    if show == True:
        plt.show()
    plt.close()

def plot_dendrogram(z, msa_df_norm, wd, prot, h, show = False):
    """
    Plots a dendrogram of the complete linkage hierarchical clustering
    calculated from the matrix of distances between sequences
    """
    plt.figure(figsize = (7.5, 4), dpi = 300)

    d = scipy.cluster.hierarchy.dendrogram(
        z, labels = msa_df_norm.index, color_threshold = h
    )
    plt.title("Complete linkage clustering of unique sequences", fontsize = 15, pad = 15)
    plt.xlabel("Unique sequence ID", fontsize = 10, labelpad = 10)
    plt.ylabel("Alignment score distance", fontsize = 10, labelpad = 10)
    plt.axhline(y = 0.2, linestyle = "--", color = "k", linewidth = 1)
    fig_out = os.path.join(wd, "{}_unique_seqs_dendrogram.png".format(prot))
    if os.path.isfile(fig_out):
        pass
    else:
        plt.savefig(fig_out)
    if show == True:
        plt.show()
    plt.close()

def assign_seqs_clusters(z, unique_seqs, seqs_df, h = 0.2):
    """
    Given a dendrogram and group of sequences,
    assigns a cluster ID to those sequences
    """
    if z.size == 0:
        print("This condition is being met")  # only one sequence
        seqs_df["seq_id"] = [0 for i in range(len(seqs_df))]
        seqs_df["cluster_id"] = [0 for i in range(len(seqs_df))]
    else:
        cutree = scipy.cluster.hierarchy.cut_tree(z, height = h) # clusters coordinates
        cluster_ids = [int(cut) for cut in cutree]
        un_seqs_dict = {el: i for i, el in enumerate(unique_seqs)}
        cluster_id_dict = {i: el for i, el in enumerate(cluster_ids)}
        seqs_df["seq_id"] = seqs_df.seq.map(un_seqs_dict)
        seqs_df["cluster_id"] = seqs_df.seq_id.map(cluster_id_dict)
    return seqs_df

def create_cluster_dirs(main_dir, prot, seqs_df):
    """
    Once sequences have been assigned a cluster ID, the .cif files
    are divided into subdirectories
    """
    unsupp_cifs_dir = os.path.join(main_dir, prot, "unsupp_cifs")
    cluster_ids = sorted([str(cluster_id) for cluster_id in seqs_df.cluster_id.unique().tolist()])
    unsupp_files = sorted(os.listdir(unsupp_cifs_dir))
    if unsupp_files == cluster_ids:
        return
    for cluster_id in cluster_ids:
        cluster_unsupp_dir = os.path.join(unsupp_cifs_dir, str(cluster_id))
        if not os.path.isdir(cluster_unsupp_dir):
            os.mkdir(cluster_unsupp_dir)
        cluster_pdbs = seqs_df[seqs_df.cluster_id == int(cluster_id)].pdb_id.tolist()
        cluster_pdb_paths = ["{}_bio.cif".format(cluster_pdb) for cluster_pdb in cluster_pdbs]
        for cluster_pdb_path in cluster_pdb_paths:
            if os.path.isfile(os.path.join(unsupp_cifs_dir, cluster_pdb_path)):
                shutil.move(
                os.path.join(unsupp_cifs_dir, cluster_pdb_path),
                os.path.join(cluster_unsupp_dir, cluster_pdb_path)
                )

def get_swissprot(): 
    """
    Retrieves sequences and their data from Swiss-Prot

    :param db: absolute path to a fasta file containing sequences, Swiss-Prot database by default
    :type db: str
    :returns: dictionary containing the sequence id, description and sequence for all proteins in Swiss-Prot
    :rtpe: dict
    """
    swissprot_dict = Bio.SeqIO.parse(swissprot, "fasta")
    proteins = {}
    for protein in swissprot_dict:
        acc = protein.id.split("|")[1]
        proteins[acc] = {}
        proteins[acc]["id"] = protein.id
        proteins[acc]["desc"] = protein.description
        proteins[acc]["seq"] = protein.seq
    return proteins

### STAMP SUPERIMPOSITION CODE ###

def cif2pdb(cif_in, pdb_out):
    """
    converts cif to pdb format, needed to run STAMP
    """
    cif_df = PDBXreader(inputfile = cif_in).atoms(format_type = "mmcif", excluded=())
    w = PDBXwriter(outputfile = pdb_out)
    id_equiv_dict = w.run(cif_df, format_type = "pdb", category = 'auth')
    return id_equiv_dict

def get_chain_dict(cif_path):
    """
    Creates an equivalence dict between asymemtric unit chain IDs and biological assembly chain IDs
    """
    cif_df = PDBXreader(inputfile = cif_path).atoms(format_type = "mmcif", excluded=())
    orig_chains = cif_df[cif_df.group_PDB == "ATOM"][['orig_auth_asym_id', 'auth_asym_id']].drop_duplicates().orig_auth_asym_id.tolist()
    new_chains = cif_df[cif_df.group_PDB == "ATOM"][['orig_auth_asym_id', 'auth_asym_id']].drop_duplicates().auth_asym_id.tolist()
    chain_dict = {new_chains[i]: orig_chains[i] for i in range(len(orig_chains))}
    return chain_dict

def generate_domains(pdbs_dir, domains_out, roi = "ALL"):
    """
    Genereates domains file, needed to run STAMP
    """
    with open(domains_out, "w+") as fh:
        for pdb in os.listdir(pdbs_dir):
            fh.write("{} {} {{{}}}\n".format(os.path.join(pdbs_dir, pdb), pdb[:-4], roi))

def stamp(domains, prefix, out):
    """
    runs stamp using the domains file as input
    """
    if "STAMPDIR" not in os.environ:
        os.environ["STAMPDIR"] = stampdir
    args = [
        stamp_bin, "-l", domains, "-rough", "-n",
        str(2), "-prefix", prefix, ">", out
    ]
    exit_code = os.system(" ".join(args))
    if exit_code != 0:
        print(" ".join(args))
    return exit_code

def transform(matrix):
    """
    runs transform to obtain set of transformed coordinates
    """
    if "STAMPDIR" not in os.environ:
        os.environ["STAMPDIR"] = stampdir
    #print(stampdir)
    args = [transform_bin, "-f", matrix, "-het"]
    #print(" ".join(args))
    exit_code = os.system(" ".join(args))
    if exit_code != 0:
        print(" ".join(args))
    return exit_code

def move_supp_files(unsupp_pdbs_dir, supp_pdbs_dir):
    """
    moves set of supperimposed coordinate files to appropriate directory
    """
    struc_files = os.listdir(unsupp_pdbs_dir)
    cwd = os.getcwd()
    for file in struc_files:
        if os.path.isfile(os.path.join(cwd, file)):
            shutil.move(
                os.path.join(cwd, file),
                os.path.join(supp_pdbs_dir, file)
            )

def move_stamp_output(wd, prefix, stamp_out_dir):
    """
    moves STAMP output files to appropriate directory
    """
    cwd = os.getcwd()
    stamp_files = sorted([file for file in os.listdir(cwd) if prefix in file]) + ["stamp_rough.trans"]
    for file in stamp_files:
        filepath = os.path.join(cwd, file)
        if os.path.isfile(filepath):
            shutil.move(filepath, os.path.join(stamp_out_dir, file))
    out_from = os.path.join(wd, prefix + ".out")
    out_to = os.path.join(stamp_out_dir, prefix + ".out")
    doms_from = os.path.join(wd, prefix + ".domains")
    doms_to = os.path.join(stamp_out_dir, prefix + ".domains")
    if os.path.isfile(out_from):
        shutil.move(out_from, out_to)
    if os.path.isfile(doms_from):
        shutil.move(doms_from, doms_to)

### DSSP + SIFTS + ARPEGGIO ###

def get_lig_data(supp_pdbs_dir, ligs_df_path):
    """
    From a directory containing a set of structurally superimposed pdbs,
    writes a csv file indicating the name, chain and residue number of the
    ligand(s) of interest in every pdb
    """
    ligs_df = pd.DataFrame([])
    for struc in os.listdir(supp_pdbs_dir):
        struc_path = os.path.join(supp_pdbs_dir, struc)
        df = PDBXreader(inputfile = struc_path).atoms(format_type = "pdb", excluded=())
        hetatm_df = df[df.group_PDB == "HETATM"]
        ligs = hetatm_df.label_comp_id.unique().tolist()
        lois = [lig for lig in ligs if lig not in non_relevant]
        for loi in lois:
            loi_df = hetatm_df[hetatm_df.label_comp_id == loi]
            lois_df_un = loi_df.drop_duplicates(["label_comp_id", "label_asym_id"])[["label_comp_id", "label_asym_id", "auth_seq_id"]]
            lois_df_un["struc_name"] = struc
            ligs_df = ligs_df.append(lois_df_un)
    ligs_df = ligs_df[["struc_name","label_comp_id", "label_asym_id", "auth_seq_id"]]
    ligs_df.to_csv(ligs_df_path, index = False)
    return ligs_df

def run_dssp(struc, supp_pdbs_dir, dssp_dir):
    """
    runs DSSP, saves and return resulting output dataframe
    """
    dssp_csv = os.path.join(dssp_dir, "dssp_" + struc.replace("pdb", "csv")) # output csv filepath
    dssp_out = os.path.join(dssp_dir, struc.replace("pdb", "dssp"))
    struc_in = os.path.join(supp_pdbs_dir, struc)
    DSSPrunner(inputfile = struc_in, outputfile = dssp_out).write()            # runs DSSP
    dssp_data = DSSPreader(inputfile = dssp_out).read()            # reads DSSP output
    dssp_data = dssp_data.rename(index = str, columns = {"RES": "PDB_ResNum"})
    dssp_data.PDB_ResNum = dssp_data.PDB_ResNum.astype(str)
    dssp_cols = ["PDB_ResNum", "SS", "ACC", "KAPPA", "ALPHA", "PHI", "PSI", "RSA"]    # selects subset of columns
    dssp_data.to_csv(dssp_csv, index = False)
    return dssp_data[dssp_cols]

def process_sifts_data(input_sifts, sifts_dir, pdb_id):
    """
    processes SIFTS table, saves and returns the processed tables
    """
    sifts_data = SIFTSreader(inputfile = input_sifts).read() # read sifts data
    sifts_mapping = sifts_data[[
            "UniProt_dbResNum", "UniProt_dbResName", "PDB_dbResName",
            "PDB_dbResNum", "PDB_dbChainId"
        ]] # subsetting sifts table
        
    sifts_mapping = sifts_mapping.rename(index = str, columns = {
            "UniProt_dbResNum":"UniProt_ResNum", "UniProt_dbResName":"UniProt_ResName",
            "PDB_dbResName":"PDB_ResName", "PDB_dbResNum":"PDB_ResNum",
            "PDB_dbChainId":"PDB_ChainID"
        }) # renaming table columns
    try:
        sifts_mapping = sifts_mapping[sifts_mapping.UniProt_ResNum != "null"]
    except:
        pass
    try:
        sifts_mapping = sifts_mapping[~sifts_mapping.UniProt_ResNum.isnull()]
    except:
        pass
    
    sifts_mapping = sifts_mapping[sifts_mapping.PDB_ResNum != "null"]
    sifts_mapping.UniProt_ResNum =  sifts_mapping.UniProt_ResNum.astype(int)
    sifts_mapping.PDB_ResNum =  sifts_mapping.PDB_ResNum.astype(str)
    
    sifts_csv = os.path.join(sifts_dir, "sifts_" + pdb_id + ".csv") # sifts output csv file
    sifts_data.to_csv(sifts_csv, index = False) # write sifts table to csv
    
    sifts_csv2 = os.path.join(sifts_dir, "sifts_mapping_" + pdb_id + ".csv")
    sifts_mapping.to_csv(sifts_csv2, index = False)
    
    shutil.move(input_sifts, os.path.join(sifts_dir, os.path.basename(input_sifts)))
    return sifts_data, sifts_mapping

### ARPEGGIO ###

def run_clean_pdb(pdb_path):
    """
    runs pdb_clean.py to prepare files for arpeggio
    """
    args = [
        clean_pdb_python_bin, clean_pdb_bin, pdb_path
    ]
    #print(args)
    exit_code = os.system(" ".join(args))
    if exit_code != 0:
        print(" ".join(args))
    return exit_code

def run_arpeggio(pdb_path, lig_name, dist = 5):
    """
    runs Arpeggio
    """
    args = [
        arpeggio_python_bin, arpeggio_bin, pdb_path, "-s",
        "RESNAME:{}".format(lig_name), "-i", str(dist), "-v"
    ]
    exit_code = os.system(" ".join(args))
    if exit_code != 0:
        print(" ".join(args))
    return exit_code

def move_arpeggio_output(wd, subdir, strucs, supp_pdbs_subdir, clean_pdbs_subdir, struc2ligs):
    """
    moves pdb_clean.py and Arpeggio output files to appropriate directory
    """
    
    for struc in strucs:
        clean_struc = struc.replace(".pdb", ".clean.pdb")
        clean_struc_path_from = os.path.join(supp_pdbs_subdir, clean_struc)
        clean_struc_path_to = os.path.join(clean_pdbs_subdir, clean_struc)
        if os.path.isfile(clean_struc_path_from):
            shutil.move(clean_struc_path_from, clean_struc_path_to)
        for arpeggio_suff in arpeggio_suffixes:
            for the_lig in struc2ligs[struc]:
                arpeggio_file_from_supp = os.path.join(supp_pdbs_subdir, struc[:-3] + "clean_{}".format(the_lig) + "." + arpeggio_suff)
                arpeggio_file_to = os.path.join(wd, "results/arpeggio/", subdir, struc[:-3] + "clean_{}".format(the_lig) + "." + arpeggio_suff)
                arpeggio_file_from_clean = os.path.join(clean_pdbs_subdir, struc[:-3] + "clean_{}".format(the_lig) + "." + arpeggio_suff)
                if os.path.isfile(arpeggio_file_from_supp):
                    shutil.move(arpeggio_file_from_supp, arpeggio_file_to)
                elif os.path.isfile(arpeggio_file_from_clean):
                    shutil.move(arpeggio_file_from_clean, arpeggio_file_to)
        for pdb_clean_suff in pdb_clean_suffixes:
            pdb_clean_file_from = os.path.join(supp_pdbs_subdir, struc + "." + pdb_clean_suff)
            pdb_clean_file_to = os.path.join(wd, "results/pdb_clean/", subdir, struc + "." + pdb_clean_suff)
            if os.path.isfile(pdb_clean_file_from):
                shutil.move(pdb_clean_file_from, pdb_clean_file_to)
            elif os.path.isfile(os.path.join(wd, "results/pdb_clean", struc + "." + pdb_clean_suff)):
                shutil.move(os.path.join(wd, "results/pdb_clean", struc + "." + pdb_clean_suff), pdb_clean_file_to)

def process_arpeggio(struc, all_ligs, clean_pdbs_dir, arpeggio_dir, sifts_dir, bio2asym_chain_dict, cif2pdb_chain_dict):
    """
    processes arpeggio output to generate the two tables that will
    be used later in the analysis
    """
    lig_cons_splits = []
    arpeggio_lig_conss = []
        
    for lig in all_ligs:
        arpeggio_out_path_lig = os.path.join(arpeggio_dir, struc.replace("pdb", "clean_{}.bs_contacts".format(lig)))
        
        all_cons_lig = pd.read_table(arpeggio_out_path_lig, header = None)

        lig_cons_split_lig = reformat_arpeggio(all_cons_lig)

        lig_cons_split_lig = add_resnames_to_arpeggio_table(struc, clean_pdbs_dir, lig_cons_split_lig)

        lig_cons_split_lig = ligand_to_atom2(lig_cons_split_lig, lig)
        
        lig_cons_split_lig["contact_type"] = lig_cons_split_lig.apply(lambda row: contact_type(row), axis = 1)

        arpeggio_lig_cons_lig = lig_cons_split_lig.sort_values(by = ["Chain (Atom1)", "ResNum (Atom1)"])
        arpeggio_lig_cons_lig = arpeggio_lig_cons_lig[["ResNum (Atom1)","Chain (Atom1)", 'ResName (Atom1)']] 
        arpeggio_lig_cons_lig = arpeggio_lig_cons_lig.drop_duplicates(subset = ["ResNum (Atom1)", "Chain (Atom1)"])
        arpeggio_lig_cons_lig = arpeggio_lig_cons_lig.rename(index = str, columns = {"ResNum (Atom1)": "PDB_ResNum", "Chain (Atom1)": "PDB_ChainID"}) 
        arpeggio_lig_cons_lig = arpeggio_lig_cons_lig.astype({"PDB_ResNum": int})
        
        lig_cons_splits.append(lig_cons_split_lig)
        arpeggio_lig_conss.append(arpeggio_lig_cons_lig)
    
    lig_cons_split = pd.concat(lig_cons_splits)
    arpeggio_lig_cons = pd.concat(arpeggio_lig_conss)
    
    ########################### ADDED TO AVOID INCORRECT BS DEFINITION DUE TO PDB RESNUMS ###########################

    pdb_id = struc[:4]
    sifts_mapping = pd.read_csv(os.path.join(sifts_dir, "sifts_mapping_{}.csv".format(pdb_id)))

    sifts_dict = {}
    for chain, chain_df in sifts_mapping.groupby("PDB_ChainID"):
        for i, row in chain_df.iterrows():
            sifts_dict[(chain, str(row["PDB_ResNum"]))] = row["UniProt_ResNum"]
    
    sifts_dict = extend_sifts_mapping_dict(sifts_dict, bio2asym_chain_dict, cif2pdb_chain_dict)
    
    lig_cons_split['ResNum (Atom1)'] = lig_cons_split['ResNum (Atom1)'].astype(str)
    arpeggio_lig_cons['PDB_ResNum'] = arpeggio_lig_cons['PDB_ResNum'].astype(str)
    lig_cons_split["UniProt_Resnum"] = lig_cons_split.set_index(['Chain (Atom1)', 'ResNum (Atom1)']).index.map(sifts_dict.get)
    arpeggio_lig_cons["UniProt_Resnum"] = arpeggio_lig_cons.set_index(['PDB_ChainID', 'PDB_ResNum']).index.map(sifts_dict.get)

    ########################### ADDED TO AVOID INCORRECT BS DEFINITION DUE TO PDB RESNUMS ###########################

    old_len1 = len(arpeggio_lig_cons)
    old_len2 = len(lig_cons_split)
    lig_cons_split = lig_cons_split[~lig_cons_split.UniProt_Resnum.isnull()]
    arpeggio_lig_cons = arpeggio_lig_cons[~arpeggio_lig_cons.UniProt_Resnum.isnull()]
    new_len1 = len(arpeggio_lig_cons)
    new_len2 = len(lig_cons_split)
    if new_len1 != old_len1:
        print("WARNING!!! {} residues lacked mapping from PDB to UniProt at arpeggio_lig_cons for {}".format(new_len1-old_len1, pdb_id))
    if new_len2 != old_len2:
        print("WARNING!!! {} residues lacked mapping from PDB to UniProt at lig_cons_split for {}".format(new_len2-old_len2, pdb_id))

    ########################### ADDED TO AVOID INCORRECT BS DEFINITION DUE TO LACKING MAPPING TO PDB RESNUMS ###########################

    lig_cons_split.to_csv(os.path.join(arpeggio_dir,  "arpeggio_all_cons_split_" + struc[:4] + ".csv"), index = False)
    arpeggio_lig_cons.to_csv(os.path.join(arpeggio_dir,  "arpeggio_lig_cons_" + struc[:4] + ".csv"), index = False)
    
    return lig_cons_split, arpeggio_lig_cons

def reformat_arpeggio(arpeggio_df):
    """
    starts formatting arpeggio table
    """
    arpeggio_df.columns = [
        'Atom_1', 'Atom_2', 'Clash', 'Covalent', 'VdW Clash', 'Vdw', 'Proximal', 'Hydrogen Bond',
        'Weak Hydrogen Bond', 'Halogen bond',  'Ionic', 'Metal Complex', 'Aromatic', 'Hydrophobic',
        'Carbonyl', 'Polar', 'Weak Polar', 'Atom proximity', 'Vdw proximity', 'Interacting entities'
    ]
    lig_cons = arpeggio_df.loc[arpeggio_df["Interacting entities"] == "INTER"] # Selecting only the interactions between our specified atoms (i.e ligands) and other selections
    lig_cons = lig_cons.sort_values(by = ["Atom_1"])
    
    split_atom1 = lig_cons.Atom_1.str.split("/", expand = True) # splits atom1 column into three new columns: chain, resnum and atom
    split_atom1.columns = ["Chain (Atom1)", "ResNum (Atom1)", "Atom (Atom1)"]
    split_atom2 = lig_cons.Atom_2.str.split("/", expand = True)  # splits atom2 column into three new columns: chain, resnum and atom
    split_atom2.columns = ["Chain (Atom2)", "ResNum (Atom2)", "Atom (Atom2)"]
    
    lig_cons_split = pd.merge(split_atom2, lig_cons, left_index = True, right_index = True) # Making a table of the contacts, but with the atom identifier split into chain, resnum, and atom
    lig_cons_split = pd.merge(split_atom1, lig_cons_split, left_index = True, right_index = True) # Making a table of the contacts, but with the atom identifier split into chain, resnum, and atom
    lig_cons_split = lig_cons_split.drop(axis = 1, labels = ["Atom_1", "Atom_2"])
    lig_cons_split["ResNum (Atom1)"] = lig_cons_split["ResNum (Atom1)"].astype(int)
    lig_cons_split["ResNum (Atom2)"] = lig_cons_split["ResNum (Atom2)"].astype(int)
    return lig_cons_split

def add_resnames_to_arpeggio_table(structure, clean_pdbs_dir, arpeggio_cons_split):
    """
    adds residue names to arpeggio table, needed for later table mergings
    """
    structure_path = os.path.join(clean_pdbs_dir, structure.replace("pdb", "clean.pdb"))
    pdb_structure = PDBXreader(inputfile = structure_path).atoms(format_type = "pdb", excluded=())
    resnames_dict = {(row.label_asym_id, int(row.label_seq_id)): row.label_comp_id for index, row in pdb_structure.drop_duplicates(['label_asym_id', 'label_seq_id']).iterrows()}
    arpeggio_cons_split["ResName (Atom1)"] = arpeggio_cons_split.set_index(["Chain (Atom1)", "ResNum (Atom1)"]).index.map(resnames_dict.get)
    arpeggio_cons_split["ResName (Atom2)"] = arpeggio_cons_split.set_index(["Chain (Atom2)", "ResNum (Atom2)"]).index.map(resnames_dict.get)
    return arpeggio_cons_split

def ligand_to_atom2(lig_cons_split, lig):
    """
    formats arpeggio table so that the ligand atoms are always Atom2
    """
    ordered_cols = [
        'Chain (Atom1)', 'ResNum (Atom1)', 'ResName (Atom1)', 'Atom (Atom1)',
        'Chain (Atom2)', 'ResNum (Atom2)', 'ResName (Atom2)', 'Atom (Atom2)',
        'Clash', 'Covalent', 'VdW Clash', 'Vdw', 'Proximal', 'Hydrogen Bond',
        'Weak Hydrogen Bond', 'Halogen bond', 'Ionic', 'Metal Complex', 'Aromatic',
        'Hydrophobic', 'Carbonyl', 'Polar', 'Weak Polar','Atom proximity',
        'Vdw proximity', 'Interacting entities'
    ]
    lig_is_atom1 = lig_cons_split[lig_cons_split["ResName (Atom1)"] == lig].sort_values("ResNum (Atom2)")
    lig_is_atom2 = lig_cons_split[lig_cons_split["ResName (Atom2)"] == lig].sort_values("ResNum (Atom1)")

    lig_is_atom1.rename(columns = {
        'Chain (Atom1)': 'Chain (Atom2)', 'Chain (Atom2)': 'Chain (Atom1)',
        'ResNum (Atom1)': 'ResNum (Atom2)', 'ResNum (Atom2)': 'ResNum (Atom1)',
        'Atom (Atom1)': 'Atom (Atom2)', 'Atom (Atom2)': 'Atom (Atom1)',
        'ResName (Atom1)': 'ResName (Atom2)', 'ResName (Atom2)': 'ResName (Atom1)'
    }, inplace = True)

    lig_cons_split_rf = pd.concat([lig_is_atom1[ordered_cols], lig_is_atom2[ordered_cols]]) # new dataframe so ligand is always atom2
    lig_cons_split_rf = lig_cons_split_rf[lig_cons_split_rf["ResName (Atom1)"].isin(aas)]
    return lig_cons_split_rf

def extend_sifts_mapping_dict(sifts_dict, bio2asym_chain_dict, cif2pdb_chain_dict):
    """
    extends PDB-UniProt residue number mapping so ALL chains in preferred biological assembly
    are present and a UniProt residue number can be obtained
    """
    unique_sifts_chains = list(set([k[0] for k in list(sifts_dict.keys())]))
    for bio_chain, asym_chain in bio2asym_chain_dict.items():
        if bio_chain not in unique_sifts_chains:
            asym_chain_ks = [k for k in list(sifts_dict.keys()) if k[0]  == asym_chain]
            for asym_chain_k in asym_chain_ks:
                sifts_dict[(cif2pdb_chain_dict[bio_chain], asym_chain_k[1])] = sifts_dict[asym_chain_k]
    return sifts_dict

def contact_type(row):
    """
    determines whether a row is backbone or sidechain contact
    """
    if row["Atom (Atom1)"] in bbone:
        return "backbone"
    else:
        return "sidechain"

### BINDING SITE DEFINITION ###

def def_bs_oc(results_dir, pdb_files, prot, subdir, lig_names, bs_def_out, attr_out, chimera_script_out, arpeggio_dir, metric = oc_metric, dist = oc_dist, method = oc_method, alt_fmt = False):
    """
    given a set of pdb structures, and other arguments, clusters ligands in space,
    defines binding sites and writes chimera attribute files and chimera script to
    format the superimposed structures to facilitate visualisation
    
    alt_fmt is a boolean I added so it works with slightly different input. Still PDB
    files, but some which coordinates were transformed using PDBe-KB transformation
    matrices. They have different nomenclature and therefore indices to get pdb_files_dict
    must be different
    """
    lig_data_df, labs = generate_ligs_res_df(arpeggio_dir, alt_fmt = alt_fmt)

    if len(lig_data_df) == 1: #should only happen in the case of only one LOI
        lig_data_df["binding_site"] = 0
    else:
        dis_out = os.path.join(results_dir, "{}_{}_{}.dis".format(prot, subdir, metric))
        get_dis_file(lig_data_df, labs, dis_out, metric = metric)
        
        ocout, ec = oc(dis_out, method = method, cut_t = dist)
        oc_dict = oc2dict(ocout)
        cluster_id_dict = {}
        for k, v in oc_dict.items():
            for member in v["members"]:
                cluster_id_dict[member] = v["new_id"]
                
        sample_colors = list(itertools.islice(rgbs(), 200)) # new_colours
        
        lig_data_df["lab"] = labs 
        lig_data_df["binding_site"] = lig_data_df.lab.map(cluster_id_dict)

    if alt_fmt == True:
        pdb_files_dict = {f[-16:-10]: f.split("/")[-1] for f in pdb_files}
        lig_data_df["pdb_id2"] = lig_data_df.pdb_id + "_" + lig_data_df.lig_chain
        lig_data_df["pdb_path"] = lig_data_df.pdb_id2.map(pdb_files_dict)
    else:
        pdb_files_dict = {f[-18:-14]: f.split("/")[-1] for f in pdb_files}
        lig_data_df["pdb_path"] = lig_data_df.pdb_id.map(pdb_files_dict)
    
    write_bs_files(lig_data_df, bs_def_out, attr_out, chimera_script_out)
    
    return lig_data_df, labs

def generate_ligs_res_df(arpeggio_dir, alt_fmt = False):
    """
    Given a directory containing processed Arpeggio output files,
    returns a dataset containing information about the ligands of
    interest binding the protein and their labels.
    """
    lig_files = sorted([file for file in os.listdir(arpeggio_dir) if file.startswith("arpeggio_all_cons_split")])
    pdbs, ligs, resnums, chains, ligs_ress = [[], [], [], [], []]
    for file in lig_files:
        if alt_fmt == True:
            pdb_id = file[-10:-6]
        else:
            pdb_id = file[-8:-4]
        df = pd.read_csv(os.path.join(arpeggio_dir, file))
        for lig, lig_df in df.groupby(["ResName (Atom2)", "ResNum (Atom2)", "Chain (Atom2)"]):
            lig_ress = lig_df["UniProt_Resnum"].unique().tolist()
            pdbs.append(pdb_id)
            ligs.append(lig[0])
            resnums.append(lig[1])
            chains.append(lig[2])
            ligs_ress.append(lig_ress)
    lig_data_df = pd.DataFrame(list(zip(pdbs, ligs, resnums, chains, ligs_ress)), columns = ["pdb_id", "lig_name", "lig_resnum", "lig_chain", "binding_res"])
    labs = [pdbs[i] + "_" + str(ligs[i]) + "_" + str(resnums[i]) + "_" + str(chains[i]) for i in range(len(ligs))]
    return lig_data_df, labs

def get_dis_file(lig_data_df, labs, out, metric = oc_metric):
    """
    creates dis file to be fed to OC
    """
    lig_res = lig_data_df.binding_res.tolist() #this is a list of lists, each list contains residue numbers interacting with ligand
    if metric == "i_rel":
        intersect_dict = get_intersect_rel_matrix(lig_res)
    elif metric == "i_abs":
        intersect_dict = get_intersect_matrix(lig_res)
    n_ligs = len(lig_res)
    with open(out, "w+") as fh:
        fh.write(str(n_ligs) + "\n")
        for lab in labs: #labs contain unique identifier for a ligand
            fh.write(lab + "\n")
        for i in range(n_ligs):
            for j in range(i+1, n_ligs):
                fh.write(str(intersect_dict[i][j]) + "\n")

def get_intersect_rel_matrix(binding_ress):
    """
    Given a set of ligand binding residues, calcualtes a
    similarity matrix between all the different sets of ligand
    binding residues.
    """
    inters = {i: {} for i in range(len(binding_ress))}
    for i in range(len(binding_ress)):
        inters[i][i] = intersection_rel(binding_ress[i], binding_ress[i])
        for j in range(i+1, len(binding_ress)):
            inters[i][j] = intersection_rel(binding_ress[i], binding_ress[j])
            inters[j][i] = inters[i][j]
    return inters

def intersection_rel(l1, l2):
    """
    Calculates relative intersection.
    """
    len1 = len(l1)
    len2 = len(l2)
    I_max = min([len1, len2])
    I = len(list(set(l1).intersection(l2)))
    return I/I_max

def get_intersect_matrix(binding_ress):
    """
    Given a set of ligand binding residues, calcualtes a
    similarity matrix between all the different sets of ligand
    binding residues.
    """
    inters = {i: {} for i in range(len(binding_ress))}
    for i in range(len(binding_ress)):
        inters[i][i] = intersection(binding_ress[i], binding_ress[i])
        for j in range(i+1, len(binding_ress)):
            inters[i][j] = intersection(binding_ress[i], binding_ress[j])
            inters[j][i] = inters[i][j]
    return inters

def intersection(l1, l2):
    """
    Calculates intersection.
    """
    I = len(list(set(l1).intersection(l2)))
    return I

def oc(oc_in, type_mat = "sim", method = oc_method, cut_t = oc_dist):
    """
    runs OC and returns exit code, should be 0 if all is ok
    """
    oc_out = oc_in.replace(".dis", "_{}_{}_cut_{}.ocout".format(type_mat, method, cut_t))
    args = [
        "/homes/2394007/oc", type_mat, method, "id", "cut", str(cut_t),
        "ps", oc_out[:-6], "<", oc_in, ">", oc_out
    ]
    exit_code = os.system(" ".join(args))
    return oc_out, exit_code

def oc2dict(ocout):
    """
    parses OC output to generate dict
    """
    oc_dict = {}
    re_cluster = re.compile("""##\s+(\d+)\s(\d+\.*\d*)\s(\d+)\s*""")
    with open (ocout, "r") as fh:
        s = 0
        n_singleton = 0
        n_cluster = 0
        for line in fh:
            if s == 0:
                if line.startswith("##"):
                    m = re_cluster.match(line)
                    cluster_id, score, cluster_size = m.groups()
                    oc_dict[cluster_id] = {"score": float(score), "size": int(cluster_size), "new_id": n_cluster}
                    n_cluster += 1
                elif line.startswith(" "):
                    oc_dict[cluster_id]["members"] = line.strip().split(" ")
                else:
                    if line.strip() == "":
                        pass
                    elif line.strip() == "UNCLUSTERED ENTITIES":
                        s = 1
            else:
                oc_dict["singleton_{}".format(n_singleton)] = {"score": "", "size": 1, "members": [line.strip(),], "new_id": n_cluster}
                n_singleton += 1
                n_cluster += 1
    return oc_dict

def write_bs_files(frag_mean_coords, bs_def_out, attr_out, chimera_script_out):
    """
    doc
    """
    frag_mean_coords = frag_mean_coords.dropna()
    frag_mean_coords.binding_site = frag_mean_coords.binding_site.astype(int)
    frag_mean_coords.lig_resnum = frag_mean_coords.lig_resnum.astype(int)
    chimera_atom_spec = (':'+ frag_mean_coords.lig_resnum.astype(str) +
                     '.'+ frag_mean_coords.lig_chain +
                     '&#/name==' + frag_mean_coords.pdb_path)
    frag_mean_coords = frag_mean_coords.assign(chimera_atom_spec = chimera_atom_spec)  
    frag_mean_coords.to_csv(bs_def_out, index = False) # saves table to csv
    write_bs_attribute_file(frag_mean_coords, attr_out)
    bs_labs = frag_mean_coords.binding_site.unique().tolist()
    write_chimera_script(chimera_script_out, bs_labs)

def write_bs_attribute_file(clustered_fragments, attr_out):
    """
    writes Chimera attribute file to later colour ligands
    according to the binding site they bind to
    """
    with open(attr_out, "w") as out:
        out.write("attribute: binding_site\n")
        out.write("match mode: 1-to-1\n")
        out.write("recipient: residues\n")
        out.write("\n".join("\t" + clustered_fragments.chimera_atom_spec.values + "\t" + clustered_fragments.binding_site.astype(str)))

def write_chimera_script(chimera_script_out, bs_labels):
    """
    writes Chimera script that will format the superimposed structures
    as well as colour ligands according to their binding site
    """
    sample_colors = list(itertools.islice(rgbs(), 200)) # new_colours
    
    chimera_rgb_string = [','.join(map(str, rgb)) for rgb in sample_colors]
    cmds = [
        "~rib", "rib #0", "ksdssp", "set silhouette", "set silhouettewidth 3",
        "background solid white", "~dis", "sel ~@/color=white", "dis sel", "namesel lois",
        "~sel"
    ]
    with open(chimera_script_out, 'w') as out:
        out.write('# neutral colour for everything not assigned a cluster\n')
        out.write('colour white\n')
    
        out.write('# colour each binding site\n')
        for i in range(0, len(bs_labels)):
            out.write('colour {} :/binding_site=={}\n'.format(','.join(list(map(str, list(sample_colors[i])))), i))
        out.write("### SOME FORMATTING ###\n")
        out.write("\n".join(cmds))
    print("Chimera script successfully created!")

### CONSERVATION AND VARIATION ###

def create_alignment_from_struc(example_struc, fasta_path, pdb_fmt = "pdb", n_it = 3, seqdb = swissprot):
    """
    given an example structure, creates and reformats an MSA
    """
    create_fasta_from_seq(get_seq_from_pdb(example_struc, pdb_fmt = pdb_fmt), fasta_path) # CREATES FASTA FILE FROM PDB FILE
    hits_out = fasta_path.replace("fa", "out")
    hits_aln = fasta_path.replace("fa", "sto")
    hits_aln_rf = fasta_path.replace(".fa", "_rf.sto")
    jackhmmer(fasta_path, hits_out, hits_aln , n_it = n_it, seqdb = seqdb,) # RUNS JACKHAMMER USING AS INPUT THE SEQUENCE FROM THE PDB AND GENERATES ALIGNMENT
    add_acc2msa(hits_aln, hits_aln_rf)

def get_seq_from_pdb(pdb_path, pdb_fmt = "mmcif"): # MIGHT NEED TO BE MORE SELECTIVE WITH CHAIN, ETC
    """
    generates aa sequence string from a pdb coordinates file
    """
    struc = PDBXreader(pdb_path).atoms(format_type = pdb_fmt, excluded=())
    return "".join([aa_code[aa] for aa in struc[struc.group_PDB == "ATOM"].drop_duplicates(["auth_seq_id"]).auth_comp_id.tolist()])

def create_fasta_from_seq(seq, out):
    """
    saves input sequence to fasta file to use as input for jackhmmer
    """
    with open(out, "w+") as fh:
        fh.write(">query\n{}\n".format(seq))

def jackhmmer(seq, hits_out, hits_aln, n_it = 3, seqdb = swissprot):
    """
    runs jackhmmer on an input seq for a number of iterations and returns exit code, should be 0 if all is ok
    """
    args = ["jackhmmer", "-N", str(n_it), "-o", hits_out, "-A", hits_aln, seq, seqdb]
    exit_code = os.system(" ".join(args))
    return exit_code

def add_acc2msa(aln_in, aln_out, fmt_in = "stockholm"):
    """
    modifies AC field of jackhmmer alignment in stockholm format.
    
    :param aln_in: path of input alignment
    :type aln_in: str, required
    :param aln_out: path of output alignment
    :type aln_in: str, required
    :param fmt_in: input and output MSA format
    :type aln_in: str, defaults to stockholm
    """
    aln = Bio.SeqIO.parse(aln_in, fmt_in)
    recs = []
    for rec in aln:
        if rec.id == "query":
            continue
        else:
            rec.annotations["accession"] = rec.id.split("|")[1]
            recs.append(rec)
    Bio.SeqIO.write(recs, aln_out, fmt_in)

def get_human_subset_msa(aln_in, human_msa_out, fmt_in = "stockholm"):
    """
    creates a subset MSA containing only human sequences
    """
    msa = Bio.AlignIO.read(aln_in, fmt_in)
    h = 0 # human seqs
    n = 0 # total seqs
    human_recs = []
    for rec in msa:
        if "HUMAN" in rec.name:
            human_recs.append(rec)
            h += 1
        n += 1
    Bio.SeqIO.write(human_recs, human_msa_out, "stockholm")
    print(aln_in.split("/")[-1][:-4], h, n)

def get_target_prot_cols(msa_in, msa_fmt = "stockholm"): 
    """
    returns list of MSA col idx that are popualted on the protein target
    """
    seqs = [str(rec.seq) for rec in Bio.SeqIO.parse(msa_in, msa_fmt) if "query" in rec.id]
    occupied_cols = [i+1 for seq in seqs for i, el in enumerate(seq) if el != "-"]
    return sorted(list(set(occupied_cols)))

def calculate_shenkin(aln_in, aln_fmt, out = None):
    """
    given an MSA, calculates Shenkin ans occupancy, gap
    percentage for all columns
    """
    cols = in_columns(aln_in, aln_fmt)
    scores = []
    occ = []
    gaps = []
    occ_pct = []
    gaps_pct = []
    for k, v in cols.items():
        scores.append(get_shenkin(v))
        stats = (get_stats(v))
        occ.append(stats[0])
        gaps.append(stats[1])
        occ_pct.append(stats[2])
        gaps_pct.append(stats[3])
    df = pd.DataFrame(list(zip(list(range(1,len(scores)+1)),scores, occ,gaps, occ_pct, gaps_pct)), columns = ['col','shenkin','occ','gaps','occ_pct','gaps_pct'])
    if out != None:
        df.to_csv(out, index = False)
    return df

def in_columns(aln_in, infmt):
    """
    returns dictionary in which column idx are the key
    and a list containing all aas aligned to that column
    is the value
    """
    aln = Bio.AlignIO.read(aln_in, infmt)
    n_cols = len(aln[0])
    cols = {}
    for col in range(1,n_cols+1):
        cols[col] = []
    for row in aln:
        seq = str(row.seq)
        for i in range(0,len(seq)):
            cols[i+1].append(seq[i])
    return cols

def get_shenkin(col):
    """
    calculates Shenkin score for an MSA column
    """
    S = get_entropy(get_freqs(col))
    return (2**S)*6

def get_freqs(col):
    """
    calculates amino acid frequences for a given MSA column
    """
    abs_freqs = {
        'A':0, 'R':0, 'N':0, 'D':0,'C':0,'Q':0, 'E':0,'G':0,'H':0, 'I':0,
        'L':0, 'K':0,'M':0, 'F':0, 'P':0,'S':0,'T':0, 'W':0, 'Y':0, 'V':0, '-':0
    }
    for aa in col:
        aa = aa.upper()
        if col.count('-') == len(col):
            abs_freqs['-'] = 1
            return abs_freqs
        if aa not in ['-', 'X', 'B']:
            abs_freqs[aa] += 1
    rel_freqs = {k: v/(len(col)-(col.count('-')+col.count('X')+col.count('B'))) for k, v in abs_freqs.items()}
    return rel_freqs

def get_entropy(freqs):
    """
    calculates Shannon's entropy from a set of aa frequencies
    """
    S = 0
    for f in freqs.values():
        if f != 0:
            S += f*math.log2(f)
    return -S

def get_stats(col):
    """
    for a given MSA column, calculates some basic statistics
    such as column residue occupancy ang gaps
    """
    n_seqs = len(col)
    gaps = col.count('-')
    occ = n_seqs - gaps
    occ_pct = occ/n_seqs
    gaps_pct = 1 - occ_pct
    return occ, gaps, occ_pct, gaps_pct

def get_and_format_shenkin(shenkin, prot_cols, out = None):
    """
    DOC
    """
    shenkin_filt = shenkin[shenkin.col.isin(prot_cols)]#.copy()
    shenkin_filt.index = range(1, len(shenkin_filt) + 1) # CONTAINS SHENKIN SCORE, OCCUPANCY/GAP PROPORTION OF CONSENSUS COLUMNS
    min_shenkin = min(shenkin_filt.shenkin)
    max_shenkin = max(shenkin_filt.shenkin)
    shenkin_filt["rel_norm_shenkin"] = 100*(shenkin_filt.shenkin - min_shenkin)/(max_shenkin - min_shenkin) # ADDING NEW COLUMNS WITH DIFFERENT NORMALISED SCORES
    shenkin_filt["abs_norm_shenkin"] = 100*(shenkin_filt.shenkin - 6)/(120 - 6)
    if out != None:
        shenkin_filt.to_csv(out, index = False)
    return shenkin_filt

def format_variant_table(df, col_mask, vep_mask = ["missense_variant"], tab_format = "gnomad"):
    """
    formats variant table, by gettint rid of empty rows that are not human sequences,
    changning column names and only keeping those variants that are of interest and
    are present in columns of interest
    """
    df_filt = df.copy(deep = True)
    df_filt.reset_index(inplace=True)
    if tab_format == "gnomad":
        df_filt.columns = [' '.join(col).strip() for col in df_filt.columns.tolist()]
    df_filt.columns = [col.lower().replace(" ", "_") for col in df_filt.columns.tolist()]
    df_filt = df_filt[df_filt.source_id.str.contains("HUMAN")]
    df_filt = df_filt.dropna(subset=["vep_consequence"])
    df_filt = df_filt[df_filt.vep_consequence.isin(vep_mask)]
    df_filt = df_filt[df_filt.alignment_column.isin(col_mask)]
    return df_filt

def get_missense_df(aln_in, aln_fmt, variants_df, shenkin_aln, prot_cols, aln_out = None, get_or = True):
    """
    doc
    """
    variants_aln = generate_subset_aln(aln_in, aln_fmt, variants_df, aln_out)
    variants_aln_info = calculate_shenkin(variants_aln, aln_fmt)
    variants_aln_info = variants_aln_info[variants_aln_info.col.isin(prot_cols)]
    vars_df = pd.DataFrame(variants_df.alignment_column.value_counts().reindex(prot_cols, fill_value = 0).sort_index()).reset_index()
    vars_df.index = range(1, len(prot_cols)+1)
    vars_df.columns = ["col", "variants"]
    merged = pd.merge(variants_aln_info, vars_df, on = "col", how = 'left')
    merged.index = range(1, len(vars_df)+1)
    merged.variants = merged.variants + 1 #pseudo count
    merged.occ = merged.occ + 1 #pseudo count
    merged["shenkin"] = shenkin_aln["shenkin"]
    merged["rel_norm_shenkin"] = shenkin_aln["rel_norm_shenkin"] 
    merged["abs_norm_shenkin"] = shenkin_aln["abs_norm_shenkin"]
    
    if get_or == True:
        merged_or = get_OR(merged)
        return merged_or
    else:
        return merged

def get_OR(df, variant_col = "variants"):
    """
    calculates OR, ln(OR) and associated p-value and CI,
    given a missense dataframe with variants and occupancy
    """
    tot_occ = sum(df.occ)
    tot_vars = sum(df[variant_col])
    idx = df.index.tolist()
    for i in idx:
        i_occ = df.loc[i,"occ"]
        i_vars = df.loc[i,variant_col]
        rest_occ = tot_occ - i_occ
        rest_vars = tot_vars - i_vars
        oddsr, pval = stats.fisher_exact([[i_vars, rest_vars], [i_occ, rest_occ]])
        vals = [i_vars,rest_vars,i_occ, rest_occ]
        se_logor = 1.96*(math.sqrt(sum(list(map((lambda x: 1/x),vals)))))
        logor = math.log(oddsr)
        df.loc[i,'oddsratio'] = oddsr
        df.loc[i,'log_oddsratio'] = logor
        df.loc[i,'pvalue'] = pval
        df.loc[i,'ci_dist'] = se_logor
    return df

def generate_subset_aln(aln_in, aln_fmt, df, aln_out = None):
    """
    creates a subset MSA containing only human sequences that present
    missense variants and returns the path of such MSA
    """
    seqs_ids = df.source_id.unique().tolist()
    aln = Bio.SeqIO.parse(aln_in, aln_fmt)
    variant_seqs = [rec for rec in aln if rec.id in seqs_ids]
    if aln_out == None:
        pref, fmt = aln_in.split(".")
        aln_out =  pref + "_variant_seqs." + fmt
    Bio.SeqIO.write(variant_seqs, aln_out, aln_fmt)
    return aln_out

def add_miss_class(df, miss_df_out = None, cons_col = "shenkin", thresholds = [30, 70], colours = ["royalblue", "green", "grey", "firebrick", "orange"]):
    """
    adds two columns to missense dataframe. These columns will put columns
    into classes according to their divergence and missense enrichment score.
    It also adds a column that will help colour MSA columns according to their
    classifications
    """
    for i in df.index:
        if df.loc[i, cons_col] <= thresholds[0] and df.loc[i, "log_oddsratio"] < 0:
            df.loc[i,"miss_class"] = "CMD"
        elif df.loc[i, cons_col] <= thresholds[0] and df.loc[i, "log_oddsratio"] > 0:
            df.loc[i,"miss_class"] = "CME"
        elif df.loc[i, cons_col] >= thresholds[1] and df.loc[i, "log_oddsratio"] < 0:
            df.loc[i,"miss_class"] = "UMD"
        elif df.loc[i, cons_col] >= thresholds[1] and df.loc[i, "log_oddsratio"] > 0:
            df.loc[i,"miss_class"] = "UME"
        else:
            df.loc[i,"miss_class"] = "None"
           
    coloring = {
        "CMD": colours[0],
        "CME": colours[1],
        "UMD": colours[3],
        "UME": colours[4],
        "None": colours[2]
    }
    df["miss_color"] =  df.miss_class.map(coloring)
    
    if miss_df_out != None:
            df.to_csv(miss_df_out, index = False)
    return df

### MERGING TABLES ###

def merge_shenkin_df_and_mapping(shenkin_df, mapping_df, aln_ids):
    shenkin_df = shenkin_df.rename(index = str,columns = {"col": "alignment_column"}) # renaming columns to be consistent with other StruVarPi dataframes
    # Narrowing the alignment down to just our protein of interest...
    prot_pfam_alns = []
    for aln_id in aln_ids:
        prot_pfam_alns.append(mapping_df.loc[aln_id])
    prot_pfam_alignment = pd.concat(prot_pfam_alns)
    prot_pfam_alignment = prot_pfam_alignment.reset_index()
    prot_pfam_alignment = prot_pfam_alignment.rename(index = str, columns = {"Alignment": "Pfam_column", "Protein_position": "UniProt_ResNum"})
    prot_pfam_alignment.columns = prot_pfam_alignment.columns.droplevel(1)
        
    # Merging the VarAlign data to the Pfam alignment to gain conservation and variation data for the whole family...
    mapped_data = pd.merge(prot_pfam_alignment, shenkin_df, left_on = "Pfam_column", right_on = "alignment_column")
    return mapped_data

def get_bs_residues(bs_definition, sifts_dir, arpeggio_dir):
    """
    This function creates a dictionary that indicates which residue
    numbers form each ligand binding site for each structure.
    """
    binding_site_res = {}
    for index, row in bs_definition.iterrows():
        if row.pdb_id not in binding_site_res:
            binding_site_res[row.pdb_id] = {}
        if row.binding_site not in binding_site_res[row.pdb_id]:
            binding_site_res[row.pdb_id][row.binding_site] = []
        sifts_mapping_dict, prot_chains = get_sifts_mapping_dict(row.pdb_id, sifts_dir)
        fragment_cons = pd.read_csv(os.path.join(arpeggio_dir, "arpeggio_all_cons_split_{}.csv".format(row.pdb_id)))
        binding_res = fragment_cons[(fragment_cons["Chain (Atom2)"] == row.lig_chain) & (fragment_cons["ResNum (Atom2)"] == row.lig_resnum) & (fragment_cons["Chain (Atom1)"].isin(prot_chains)) & (fragment_cons["ResName (Atom1)"].isin(pdb_resnames))].drop_duplicates(subset = ["Chain (Atom1)", "ResNum (Atom1)", "Chain (Atom2)", "ResNum (Atom2)"]).sort_values(by = ["ResNum (Atom1)"])["ResNum (Atom1)"].tolist()
        binding_res_chain = fragment_cons[(fragment_cons["Chain (Atom2)"] == row.lig_chain) & (fragment_cons["ResNum (Atom2)"] == row.lig_resnum) & (fragment_cons["Chain (Atom1)"].isin(prot_chains)) & (fragment_cons["ResName (Atom1)"].isin(pdb_resnames))].drop_duplicates(subset = ["Chain (Atom1)", "ResNum (Atom1)", "Chain (Atom2)", "ResNum (Atom2)"]).sort_values(by = ["ResNum (Atom1)"])["Chain (Atom1)"].tolist()
        ress = [str(binding_res[i]) + binding_res_chain[i] for i in range(0, len(binding_res))]
        binding_site_res[row.pdb_id][row.binding_site].extend([sifts_mapping_dict[res] for res in ress])
    return binding_site_res

def get_sifts_mapping_dict(pdb_id, sifts_dir):
    """
    This function creates a sifts mapping dictionary for the
    ligand binding residues: PDB_ResNum + PDB_ChainID: UniProt_ResNum.
    It also returns the PDB chain IDs that correspond to the protein.
    """
    mapping_df = pd.read_csv(os.path.join(sifts_dir, "sifts_mapping_{}.csv".format(pdb_id)))
    mapping_df.PDB_ResNum = mapping_df.PDB_ResNum.astype(str)
    binding_site_res = {str(row.PDB_ResNum) + row.PDB_ChainID : row.UniProt_ResNum for index, row in mapping_df.iterrows()}
    prot_chains = mapping_df.PDB_ChainID.unique().tolist()
    return binding_site_res, prot_chains

def get_bs_sig_cols(strucs, bs_ids, binding_site_res, all_contact_variations):
    """
    This function creates a dictionary that indicates to which binding
    site the ligand interacting residues for each structure belong to.
    """
    bs_sig_cols = {}
    for struc in [struc[:4] for struc in strucs]:
        if struc not in binding_site_res:
            print("{} does not have any LOIs, can't get its ligand binding residues!".format(struc))
            continue
        bs_sig_cols[struc] = {}
        for i in bs_ids:
            bs_sig_cols[struc][i] = []
            if i not in binding_site_res[struc]:
                bs_sig_cols[struc][i].extend([0 for i in range(0, len(all_contact_variations[all_contact_variations.structure == struc]))])
            else:
                fragment_binding_res = all_contact_variations[all_contact_variations.structure == struc].UniProt_ResNum.tolist()
                for binding_res in fragment_binding_res:
                    if binding_res in binding_site_res[struc][i]:
                        bs_sig_cols[struc][i].append(1)
                    else:
                        bs_sig_cols[struc][i].append(0)
    return bs_sig_cols

def add_bs_info2df(bs_sig_cols, all_contact_variations, out = None):
    """
    This function adds the binding site fingerprint
    for every ligand binding residue in the dataset
    """
    df_sig_cols = {}
    for k1, v1 in bs_sig_cols.items(): #k1 is structure id
        for k2, v2 in v1.items(): #k2 is bs id
            if k2 not in df_sig_cols:
                df_sig_cols[k2] = []
            df_sig_cols[k2].extend(v2)
    l = ["BS"+str(i) for i in df_sig_cols.keys()]
    bs_fingerprint = pd.DataFrame(df_sig_cols)
    bs_fingerprint.columns = l
    contact_vars_bs_sig_df = pd.concat([all_contact_variations.reset_index(), bs_fingerprint], axis=1)
    if out != None:
        contact_vars_bs_sig_df.to_csv(out, index = False)
    return contact_vars_bs_sig_df

def get_totals(mapped_data, prot, sifts_dir):
    """
    returns total variants and human occupancy of all
    protein structure positions with residues aligned in MSA
    """
    coords = get_seq_range(sifts_dir, prot)
    mapped_data_filt = mapped_data[mapped_data.UniProt_ResNum.isin(range(coords[0], coords[1]+1))]
    mapped_data_filt = mapped_data_filt.drop_duplicates("alignment_column")
    total_vars = mapped_data_filt.variants.sum()
    total_occ = mapped_data_filt.human_occ.sum()
    return (total_vars, total_occ)

def get_seq_range(sifts_dir, prot):
    """
    returns the sequence range of an example structure
    """
    sifts_path = os.path.join(sifts_dir, [file for file in os.listdir(sifts_dir) if file.startswith("sifts_mapping")][0])
    example_mapping = pd.read_csv(sifts_path)
    mapped_res = sorted(example_mapping.UniProt_ResNum.unique().tolist()) #without sorted, resulted in wrong coordinates for 5 examples
    return [mapped_res[0], mapped_res[-1]]

def create_binding_site_df(struvarpi_results):
    """
    Creates binding site dataframe from StruVarPi table.
    Contains information about variants, MES, shenkin,
    and other statistics for each binding site.
    """
    l_sites, l_vars, l_occ, l_vars_occ, l_mes, l_pvals, l_shenkin, l_shenkin_ci, l_mes_ci, l_bs_res = ([] for i in range(10))
    total_bs = 0
    struvarpi_df = struvarpi_results[0]                                             #first element of list contains struvarpi result df for the protein
    bs_labels = [col for col in struvarpi_df.columns if col.startswith("BS")]  #list of dataframe columns that contain info about BSs
    total_vars, total_occ = struvarpi_results[1]                               #variants an occupancy of whole protein structure
    for bs_label in bs_labels: # loops through all defined binding sites in protein structure
        total_bs += 1
        bs_occ = struvarpi_df[(struvarpi_df[bs_label] == 1)].drop_duplicates(subset = ["alignment_column"]).human_occ.sum() # human alignment occupancy of columns forming bs
        bs_vars = struvarpi_df[(struvarpi_df[bs_label] == 1)].drop_duplicates(subset = ["alignment_column"]).variants.sum() # human missense variants in columns forming bs
        shenkins = struvarpi_df[(struvarpi_df[bs_label] == 1)].drop_duplicates(subset = ["UniProt_ResNum"]).rel_norm_shenkin.tolist() #shenkins of columns forming bs
        rest_occ = total_occ - bs_occ
        rest_vars = total_vars - bs_vars
        if math.isnan(bs_occ):
            print("CAREFUL!!! {}\t{}\t{} --> NaN fields".format(bs_label, bs_vars, bs_occ)) # patch for cases where BS is not covered at all by MSA
            continue
        l_bs_res.append(len(shenkins))
        l_sites.append(bs_label)
        l_vars.append(bs_vars)
        l_occ.append(bs_occ)
        if bs_occ == 0:
            print("CAREFUL!!! {}\t{}\t{}".format(bs_label, bs_vars, bs_occ))
            l_vars_occ.append("")
            l_mes.append("")
            l_pvals.append("")
            l_shenkin.append("")
            l_shenkin_ci.append("")
            l_mes_ci.append("")
        else:
            if len(shenkins) == 1:
                shenkin_ci = 0
                avg_shenkin = statistics.mean(shenkins)
            else:
                avg_shenkin = statistics.mean(shenkins)
                shenkin_ci = 1.96*(statistics.stdev(shenkins)/math.sqrt(len(shenkins)))
            vals = [bs_vars, rest_vars, bs_occ, rest_occ]
            se_logor = 1.96*(math.sqrt(sum(list(map((lambda x: 1/x),vals))))) #crashes if there is only one binding site
            oddsr, pval = stats.fisher_exact([[bs_vars, rest_vars], [bs_occ, rest_occ]])
            l_vars_occ.append(bs_vars/bs_occ)
            l_mes.append(math.log(oddsr))
            l_pvals.append(pval)
            l_shenkin.append(avg_shenkin)
            l_shenkin_ci.append(shenkin_ci)
            l_mes_ci.append(se_logor)
    
    mes_sgc_df = pd.DataFrame(list(zip(
    l_sites, l_vars, l_occ, l_vars_occ, l_mes, l_pvals, l_shenkin, l_shenkin_ci, l_mes_ci, l_bs_res))
             ,columns = ['bs_id', "vars", "occ", "vars_per_occ", "MES", "p", "norm_shenkin_rel", "shenkin_ci", "MES_ci", "number_bs_res"])
    return mes_sgc_df

def add_bs_msa_coverage(mes_sgc_df, struvarpi_df, binding_site_res):
    """
    Returns a dictionary with the MSA coverage for each BS,
    i.e., of all residues that bind ligands, how many are
    covered
    """
    bs_coverage_data = {}
    bs_coverage_dict = {}
    bs_labels = mes_sgc_df.bs_id.unique().tolist()
    for bs_label in bs_labels:
        bs_coverage_data[bs_label] = {}
        bs_res = struvarpi_df[(struvarpi_df[bs_label] == 1)].drop_duplicates(subset = ["UniProt_ResNum"]).UniProt_ResNum.tolist()
        bs_res_pdb = struvarpi_df[(struvarpi_df[bs_label] == 1)].drop_duplicates(subset = ["PDB_ResNum"]).PDB_ResNum.tolist()
        bs_coverage_data[bs_label]["table_bs_res"] = sorted(bs_res)
    bs_unique_res = get_unique_bs_res(binding_site_res)
    for k, v in bs_unique_res.items():
        if "BS"+ str(k) in bs_coverage_data:
            bs_coverage_data["BS"+ str(k)]["real_bs_res"] = v
            bs_coverage_data["BS"+ str(k)]["bs_cov"] = len(bs_coverage_data["BS"+ str(k)]["table_bs_res"])/len(bs_coverage_data["BS"+ str(k)]["real_bs_res"])
            bs_coverage_dict["BS"+ str(k)] = len(bs_coverage_data["BS"+ str(k)]["table_bs_res"])/len(bs_coverage_data["BS"+ str(k)]["real_bs_res"])
    
    mes_sgc_df["bs_cov"] = mes_sgc_df["bs_id"].map(bs_coverage_dict)
    return mes_sgc_df, bs_unique_res

def get_unique_bs_res(binding_site_res):
    """
    returns a dictionary with the unique residues of each binding site
    """
    bs_unique_res = {} # contains unique UniProt residues forming each BS
    for k1, v1 in binding_site_res.items():
        for k2, v2 in v1.items():
            if k2 not in bs_unique_res:
                bs_unique_res[k2] = []
            bs_unique_res[k2].extend(v2)

    for k, v in bs_unique_res.items():
            bs_unique_res[k] = sorted(list(set(bs_unique_res[k])))
    return bs_unique_res

### PLOTTING ###

def plot_binding_site(df, cons_col = "shenkin", color = "blue", thresholds = [30,70], out = None, label = True, error = True, ylims = None, pltitle = None, show = True):
    """
    doc
    """
    fig = plt.figure(figsize = (17, 12))
    unsig_df = df[df.pvalue > mes_sig_t]
    sig_df = df[df.pvalue <= mes_sig_t]
    
    plt.scatter(unsig_df[cons_col], unsig_df.log_oddsratio, c = unsig_df.pvalue, cmap = plt.cm.Greys_r, s = 250, edgecolor = 'black', linewidth = 1.5, marker = ".")
    
    plt.scatter(sig_df[cons_col], sig_df.log_oddsratio, c = sig_df.pvalue, cmap = plt.cm.Greys_r, s = 250, edgecolor = 'black', linewidth = 1.5, marker = "D")
    
    cbar = plt.colorbar()
    plt.clim(0,1)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('p-value', fontsize = 25)# rotation = 270)
    cbar.ax.tick_params(labelsize=20)
    cbar.set_ticks([n/10 for n in range(0, 11)])
    cbar.set_ticklabels([n/10 for n in range(0, 11)])
    if error == True:
        plt.errorbar(df[cons_col], df.log_oddsratio, yerr = df.ci_dist, linewidth = 1.5, linestyle = "None", capsize = 4.5, capthick = 1.9)
    if label == True:
        for x, y, z in zip(df[cons_col], df.log_oddsratio, df.col):
            label = "{:d}".format(z)
            plt.annotate(label, # this is the text
                            (x, y), # this is the point to label
                            textcoords = "offset points", # how to position the text
                            xytext = (-20, 0), # distance from text to points (x,y)
                            ha = 'right',
                            fontsize = 17.5) # horizontal alignment can be left, right or center
    plt.xlabel("Normalised Shenkin divergence score", fontsize = 22.5, labelpad = 15)
    plt.ylabel("MES", fontsize = 22.5, labelpad = 15)
    
    if pltitle != None:
        plt.title(pltitle, fontsize = 25, pad = 20)
        
    plt.tick_params(axis = 'both' , which = 'major', pad = 7.5, width = 1.9, length = 6.25, labelsize = 20)
    ylims_round = [round(ylim, 1) for ylim in ylims]
    lower_yticks = list(np.arange(ylims_round[0], 0, 0.25))
    higher_yticks = list(np.arange(lower_yticks[-1], ylims_round[1], 0.25))
    all_yticks = [round(n,2) for n in lower_yticks + [0, ] + higher_yticks[1:]]
    plt.yticks(all_yticks)
    plt.xticks(np.arange(0, 110, 10))
    
    if ylims != None:
        plt.ylim(ylims[0]-0.1, ylims[1]+0.1)
    plt.xlim(-5, 105)
    
    plt.axhline(color = "black", linewidth = 1, linestyle = '--')
    plt.axvline(x = thresholds[0], color = "black", linewidth = 1, linestyle = '--')
    plt.axvline(x = thresholds[1], color = "black", linewidth = 1, linestyle = '--')
    
    if out != None:
        if os.path.isfile(out):
            pass
        else:
            plt.savefig(out)
        
    if show == False:
        plt.close(fig)
    else:
        plt.show()

def plot_prot_bss(df, prot, bs_color_dict, out = None, show = False, override = False):
    """
    doc
    """
    bs_ids = df.bs_id.unique().tolist()
    fig = plt.figure(figsize = (15, 7.5))
    plt.scatter(df.norm_shenkin_rel, df.MES, s = 100, linewidth = 1, edgecolor = "black", color = df.bs_color)
    plt.errorbar(df.norm_shenkin_rel, df.MES, yerr = df.MES_ci, linestyle = "None", color = "black", linewidth= 0.75, capsize  = 5, capthick = 0.75)
    plt.axhline(y = 0, linewidth = 1, linestyle = "--", color = "black")
    plt.axvline(x = 25, linewidth = 1, linestyle = "--", color = "black")
    plt.axvline(x = 75, linewidth = 1, linestyle = "--", color = "black")
    plt.xlim(0, 100)
    plt.xticks(range(0, 105, 5))
    plt.xlabel("Normalised Shenkin divergence score", fontsize = 20)
    plt.ylabel("MES", fontsize = 20)
    plt.title(prot, fontsize = 25, pad = 15)
    plt.tick_params(axis= 'both' , which = 'major', pad = 6, width = 1.5, length = 5, labelsize = 14)
    fig_leg = []
    for bs_id in bs_ids:
        fig_leg.append(mpatches.Patch(color = bs_color_dict[bs_id], label = bs_id))
    legend = plt.legend(handles = fig_leg, loc='center left', bbox_to_anchor=(1.0125, 0.5),fontsize = 12)
    legend.get_frame().set_linewidth(1.5)
    legend.get_frame().set_edgecolor("black")
    for legobj in legend.legendHandles:
        legobj.set_linewidth(1)
        legobj.set_edgecolor("black")
    if out != None:
        if os.path.isfile(out) and override == False:
            pass
        else:
            plt.savefig(out)
    if show == False:
        plt.close(fig)
    else:
        plt.show()

def get_overall_stats(wd, prot, subdir, lig_data_path, bs_ids, varalign_dir, totals):
    """
    doc
    """
    out = os.path.join(wd, "{}_{}_fragsys_report.txt".format(prot, subdir))
    ligdat = pd.read_csv(lig_data_path)
    n_strucs = len(ligdat.struc_name.unique().tolist())
    n_ligs = len(ligdat.label_comp_id.tolist())
    n_unique_ligs = len(ligdat.label_comp_id.unique().tolist())
    full_msa_path = os.path.join(varalign_dir, "{}_{}_rf.sto".format(prot, subdir))
    human_msa_path = os.path.join(varalign_dir, "{}_{}_rf_human.sto".format(prot, subdir))
    variants_msa_path = os.path.join(varalign_dir, "{}_{}_rf_human_missense_variants_seqs.sto".format(prot, subdir))
    full_msa_len = len([rec.id for rec in Bio.SeqIO.parse(full_msa_path, "stockholm")])
    human_msa_len = len([rec.id for rec in Bio.SeqIO.parse(human_msa_path, "stockholm")])
    variants_msa_len = len([rec.id for rec in Bio.SeqIO.parse(variants_msa_path, "stockholm")])
    sum_df_out = os.path.join(wd, "{}_{}_results.csv").format(prot, subdir)
    if os.path.isfile(sum_df_out):
        pass
    else:
        sum_df = pd.DataFrame(
            [
                [
                    prot, subdir, n_strucs, n_ligs, n_unique_ligs, len(bs_ids),
                    full_msa_len, human_msa_len, variants_msa_len, totals[0], totals[1]
                ]
            ],
            columns = [
                "acc", "group", "n_strucs", "n_ligs", "n_un_ligs", "n_bs",
                "n_seqs", "n_human_seqs", "n_var_seqs", "n_vars", "n_human_res"
            ]
        )
        sum_df.to_csv(sum_df_out, index = False)
    if os.path.isfile(out):
        pass
    else:
        with open(out, "w+") as fh:
            fh.write("This is a report for the fragment screening analysis of {}\n".format(prot))
            fh.write("There were {} structures for this protein\n".format(n_strucs))
            fh.write("There were a total of {} ligands, {} of which were unique\n".format(n_ligs, n_unique_ligs))
            fh.write("{} binding sites were defined on this protein\n".format(len(bs_ids)))
            fh.write("The full alignment has {} sequences, {} of which are human and {} of which present missense variants\n".format(full_msa_len, human_msa_len, variants_msa_len))
            fh.write("{} human missense variants were found across {} human residues\n".format(totals[0], totals[1]))
        print("Overall stats were printed to report!")

### COLORS ###

# This code I did not write, I grabbed it from the following URL:

# https://stackoverflow.com/questions/470690/how-to-automatically-generate-n-distinct-colors

import colorsys
import itertools
from fractions import Fraction
from typing import Iterable, Tuple

def zenos_dichotomy() -> Iterable[Fraction]:
    """
    http://en.wikipedia.org/wiki/1/2_%2B_1/4_%2B_1/8_%2B_1/16_%2B_%C2%B7_%C2%B7_%C2%B7
    """
    for k in itertools.count():
        yield Fraction(1,2**k)

def fracs() -> Iterable[Fraction]:
    """
    [Fraction(0, 1), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4), Fraction(1, 8), Fraction(3, 8), Fraction(5, 8), Fraction(7, 8), Fraction(1, 16), Fraction(3, 16), ...]
    [0.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 0.0625, 0.1875, ...]
    """
    yield Fraction(0)
    for k in zenos_dichotomy():
        i = k.denominator # [1,2,4,8,16,...]
        for j in range(1,i,2):
            yield Fraction(j,i)

# can be used for the v in hsv to map linear values 0..1 to something that looks equidistant
# bias = lambda x: (math.sqrt(x/3)/Fraction(2,3)+Fraction(1,3))/Fraction(6,5)

HSVTuple = Tuple[Fraction, Fraction, Fraction]
RGBTuple = Tuple[float, float, float]

def hue_to_tones(h: Fraction) -> Iterable[HSVTuple]:
    for s in [Fraction(6,10)]: # optionally use range
        for v in [Fraction(8,10),Fraction(5,10)]: # could use range too
            yield (h, s, v) # use bias for v here if you use range

def hsv_to_rgb(x: HSVTuple) -> RGBTuple:
    return colorsys.hsv_to_rgb(*map(float, x))

flatten = itertools.chain.from_iterable

def hsvs() -> Iterable[HSVTuple]:
    return flatten(map(hue_to_tones, fracs()))

def rgbs() -> Iterable[RGBTuple]:
    return map(hsv_to_rgb, hsvs())

def rgb_to_css(x: RGBTuple) -> str:
    uint8tuple = map(lambda y: int(y*255), x)
    return "rgb({},{},{})".format(*uint8tuple)

def css_colors() -> Iterable[str]:
    return map(rgb_to_css, rgbs())
        
### UNUSED ###

def get_dist_mat(lig_data_df, labs, arpeggio_dir, pdb_files):
    """
    doc
    """
    pdb_ids = [pdb_file.split("/")[-1][:4] for pdb_file in pdb_files]
    lig_data_df_filt = lig_data_df[lig_data_df.pdb_id.isin(pdb_ids)]
    labs_filt = [lab for lab in labs if lab[:4] in pdb_ids]
    print("lig_data_df == lig_data_df_filt is {}".format(lig_data_df.equals(lig_data_df_filt)))
    print("labs == labs_filt is {}".format(labs == labs_filt))
    sims = get_sims_matrix(lig_data_df_filt.binding_res.tolist())
    sims_df = pd.DataFrame(sims)
    dists_df = 1 -(sims_df)
    return dists_df, labs_filt, lig_data_df_filt

def get_sims_matrix(binding_ress):
    """
    Given a set of ligand binding residues, calcualtes a
    similarity matrix between all the different sets of ligand
    binding residues.
    """
    sims = {i: {} for i in range(len(binding_ress))}
    for i in range(len(binding_ress)):
        sims[i][i] = 1
        for j in range(i+1, len(binding_ress)):
            sims[i][j] = jaccard_sim(binding_ress[i], binding_ress[j])
            sims[j][i] = sims[i][j]
    return sims

def jaccard_sim(l1, l2):
    """
    Calculates Jaccard Similiarity index.
    """
    I = len(list(set(l1).intersection(l2)))
    U = (len(set(l1)) + len(set(l2))) - I
    return float(I) / U

### THE END ###
