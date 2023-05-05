### this file contains the main fragsys function ### 

from fragsys import *

### FRAGSYS ###

def main(main_dir, prot, input_df):
    """
    This is the main fragsys function. It carries
    out a series of calculations and operations for
    a protein given an input pandas dataframe and saves
    all results to a subdirectory of the main directory
    """
    
    ### DIRECTORIES SETUP ###

    wd = os.path.join(main_dir, prot)
    unsupp_cifs_dir = os.path.join(wd, "unsupp_cifs")
    unsupp_pdbs_dir = os.path.join(wd, "unsupp_pdbs")
    supp_pdbs_dir = os.path.join(wd, "supp_pdbs")
    clean_pdbs_dir = os.path.join(wd, "clean_pdbs")
    stamp_out_dir = os.path.join(wd, "stamp_out")
    results_dir = os.path.join(wd, "results")
    sifts_dir = os.path.join(results_dir, "sifts")
    dssp_dir = os.path.join(results_dir, "dssp")
    pdb_clean_dir = os.path.join(results_dir, "pdb_clean")
    arpeggio_dir = os.path.join(results_dir, "arpeggio")
    varalign_dir = os.path.join(results_dir, "varalign")
    figs_dir = os.path.join(results_dir, "figs")

    dirs = [
            unsupp_cifs_dir, unsupp_pdbs_dir, supp_pdbs_dir,
            clean_pdbs_dir, stamp_out_dir, results_dir, pdb_clean_dir,
            arpeggio_dir, sifts_dir, dssp_dir, varalign_dir, figs_dir
        ]
    
    setup_fragsys_start(main_dir, prot, input_df)
    
    setup_dirs(dirs)
    
    classify_pdbs(prot, main_dir)
    
    ### CONVERSION TO PDB AND STAMPING PROCESS ###
    
    subdirs = list(map(str, sorted([int(el) for el in os.listdir(unsupp_cifs_dir)])))
    
    sp = get_swissprot()

    cfg.dssp_bin = dssp_bin
    aln_fmt = "stockholm"
    
    for direc in dirs[1:]:
        if direc == results_dir:
            continue
        for subdir in subdirs:
            try:
                os.mkdir(os.path.join(direc, subdir))
            except:
                pass #directory exists

    for subdir in subdirs:

        print("Starting to process group {} of {} in {}".format(subdir, subdirs, prot))
        
        unsupp_pdbs_subdir = os.path.join(unsupp_pdbs_dir, subdir)
        unsupp_cifs_subdir = os.path.join(unsupp_cifs_dir, subdir)
        supp_pdbs_subdir = os.path.join(supp_pdbs_dir, subdir)
        stamp_out_subdir = os.path.join(stamp_out_dir, subdir)
        clean_pdbs_subdir = os.path.join(clean_pdbs_dir, subdir)
        arpeggio_subdir = os.path.join(arpeggio_dir, subdir)
        dssp_subdir = os.path.join(dssp_dir, subdir)
        sifts_subdir = os.path.join(sifts_dir, subdir)
        varalign_subdir = os.path.join(varalign_dir, subdir)
        figs_subdir = os.path.join(figs_dir, subdir)
        
        print("Starting CIF to PDB conversion section!")
         
        subdir_cifs = os.listdir(unsupp_cifs_subdir)
        cif2pdb_chain_dict = {}
        bio2asym_chain_dict = {}
        for cif in subdir_cifs:
            cif_in = os.path.join(unsupp_cifs_subdir, cif)
            pdb_out = os.path.join(unsupp_pdbs_subdir, cif[:-3] + "pdb")
            if os.path.isfile(pdb_out):
                pass
            else:
                id_equiv_dict = cif2pdb(cif_in, pdb_out) # only generated if files do not exist
                cif2pdb_chain_dict[cif[:4]] = id_equiv_dict
                bio2asym_chain_dict[cif[:4]] = get_chain_dict(cif_in)

        print("Starting STAMP section!")
        
        domains_out = os.path.join(wd, "{}_{}_stamp.domains".format(prot, subdir))
        if os.path.isfile(domains_out):
                pass
        else:
            generate_domains(unsupp_pdbs_subdir, domains_out)

        prefix = "{}_{}_stamp".format(prot, subdir)
        n_strucs = len(os.listdir(unsupp_pdbs_subdir))
        matrix_file = prefix + "." + str(n_strucs-1)

        if os.path.isfile(os.path.join(stamp_out_dir, subdir, matrix_file)):
            pass
        else:
            ec = stamp(
                domains_out,
                prefix, os.path.join(wd, prefix + ".out")
            )
            if ec == 0:
                pass
            else:
                print("Something went wrong with STAMP :(")

        structure_files = os.listdir(unsupp_pdbs_subdir)
        c = 0
        
        for file in structure_files:
            if os.path.isfile(os.path.join(supp_pdbs_dir, subdir, file)): # only when they alaready have been transformed
                c += 1
        if c == len(structure_files):
            pass
        else:
            print("Proceeding to run TRANSFORM")
            if not os.path.isfile(matrix_file): # RUNNING TRANSFORM ONCE STAMP OUTPUT HAS BEEN MOVED TO STAMP_OUT_DIR
                matrix_file = os.path.join(stamp_out_dir, subdir, prefix + "." + str(n_strucs-1))
            ec = transform(matrix_file) #running transform with matrix on cwd
            if ec == 0:
                pass
            else:
                print("Something went wrong with TRANSFORM :(")
                
        move_supp_files(unsupp_pdbs_subdir, supp_pdbs_subdir)

        move_stamp_output(wd, prefix, stamp_out_subdir)
        
        ### PROCESSING OF STRUCTURES BEFORE PIPELINE ###

        lig_data_path = os.path.join(wd, "{}_{}_lig_data.csv".format(prot, subdir))
        if os.path.isfile(lig_data_path):
            ligs_df = pd.read_csv(lig_data_path)
        else:
            ligs_df = get_lig_data(supp_pdbs_subdir, lig_data_path)
        
        ### PDB - UNIPROT SEQUENCE MAPPING ###

        print("Starting UNIPROT-PDB mapping section!")
        
        pdb_mappings = []
        strucs = [file for file in os.listdir(unsupp_pdbs_subdir) if file.endswith("_bio.pdb")]
        for struc in strucs:
            pdb_id = struc[:4]
            dssp_csv = os.path.join(dssp_subdir, "dssp_" + struc.replace("pdb", "csv"))
            if os.path.isfile(dssp_csv):
                dssp_data = pd.read_csv(dssp_csv)
                pass
            else:
                dssp_data = run_dssp(struc, supp_pdbs_subdir, dssp_subdir)
            input_sifts = os.path.join(os.getcwd(), ".prointvar/sifts", "{}.xml".format(pdb_id))
            input_sifts_moved = os.path.join(sifts_subdir, "{}.xml".format(pdb_id))
            if os.path.isfile(input_sifts_moved):
                pass
            else:
                try:
                    download_sifts_from_ebi(pdb_id)
                except:
                    print("ERROR: SIFTS not found for {}!".format(pdb_id))
                    continue
            sifts_out1 = os.path.join(sifts_subdir, "sifts_" + pdb_id + ".csv")
            sifts_out2 = os.path.join(sifts_subdir, "sifts_mapping_" + pdb_id + ".csv")
            if os.path.isfile(sifts_out1) and os.path.isfile(sifts_out2):
                sifts_data = pd.read_csv(sifts_out1)
                sifts_mapping = pd.read_csv(sifts_out2)
            else:
                sifts_data, sifts_mapping = process_sifts_data(input_sifts, sifts_subdir, pdb_id)
                print("SIFTS output being processed for {}!".format(pdb_id))
            mapping = pd.merge(sifts_mapping, dssp_data, left_on = "PDB_ResNum", right_on = "PDB_ResNum")
            pdb_mappings.append(mapping)
    
        ### ARPEGGIO PART ###

        print("Starting ARPEGGIO section!")
        struc2ligs = {}
        for struc in strucs:
            struc2ligs[struc] = []
            struc_df = ligs_df[ligs_df.struc_name == struc]
            pdb_path = os.path.join(supp_pdbs_subdir, struc)
            clean_pdb_path = pdb_path.replace(".pdb", ".clean.pdb")
            if os.path.isfile(os.path.join(clean_pdbs_subdir, struc.replace(".pdb", ".clean.pdb"))) or os.path.isfile(os.path.join(supp_pdbs_subdir, struc.replace(".pdb", ".clean.pdb"))):
                pass
            else:
                ec = run_clean_pdb(pdb_path)
                if ec == 0:
                    pass
                else:
                    print("Something went wrong when cleaning {} :(".format(struc[:4]))
                    pass
            ligs = struc_df.label_comp_id.unique().tolist()

            for the_lig in ligs: # RUNs ARPEGGIO ONCE FOR EACH LIGAND
                struc2ligs[struc].append(the_lig)

                if not os.path.isfile(clean_pdb_path):
                    clean_pdb_path = os.path.join(clean_pdbs_subdir, struc.replace(".pdb", ".clean.pdb"))

                if os.path.isfile(os.path.join(arpeggio_subdir, struc[:-3] + "clean_{}.bs_contacts".format(the_lig))) or os.path.isfile(os.path.join(supp_pdbs_subdir, struc[:-3] + "clean_{}.bs_contacts".format(the_lig))):
                    continue

                ec = run_arpeggio(clean_pdb_path, the_lig)
                if ec == 0:
                    print("Arpeggio ran sucessfully for {} in {}!".format(the_lig, struc[:4]))
                    for arpeggio_suff in arpeggio_suffixes: # CHANGES ARPEGGIO OUTPUT FILENAMES SO THEY INCLUDE LIGAND NAME
                        arpeggio_file_old_name_supp = os.path.join(supp_pdbs_subdir, struc[:-3] + "clean" + "." + arpeggio_suff)
                        arpeggio_file_new_name_supp = os.path.join(supp_pdbs_subdir, struc[:-3] + "clean_{}".format(the_lig) + "." + arpeggio_suff)
                        arpeggio_file_old_name_clean = os.path.join(clean_pdbs_subdir, struc[:-3] + "clean" + "." + arpeggio_suff)
                        arpeggio_file_new_name_clean = os.path.join(clean_pdbs_subdir, struc[:-3] + "clean_{}".format(the_lig) + "." + arpeggio_suff)
                        if os.path.isfile(arpeggio_file_old_name_supp):
                            os.rename(arpeggio_file_old_name_supp, arpeggio_file_new_name_supp)
                        elif os.path.isfile(arpeggio_file_old_name_clean):
                            os.rename(arpeggio_file_old_name_clean, arpeggio_file_new_name_clean)
                else:
                    print("Something went wrong when running Arpeggio for {} :(".format(struc[:4]))
                    pass
        move_arpeggio_output(wd, subdir, strucs, supp_pdbs_subdir, clean_pdbs_subdir, struc2ligs)

        ligand_contact_list = []
        for struc in strucs:
            all_ligs = ligs_df[ligs_df.struc_name == struc].label_comp_id.unique().tolist()

            arpeggio_out1 = os.path.join(arpeggio_subdir, "arpeggio_all_cons_split_" + struc[:4] + ".csv") # output file 1
            arpeggio_out2 = os.path.join(arpeggio_subdir,  "arpeggio_lig_cons_" + struc[:4] + ".csv") # output file 2
            if os.path.isfile(arpeggio_out1) and os.path.isfile(arpeggio_out2):
                lig_cons_split = pd.read_csv(arpeggio_out1)
                arpeggio_lig_cons = pd.read_csv(arpeggio_out2)
            else:
                if len(all_ligs) == 0:
                    print("No LOIs in {}, so skipping!".format(struc[:4]))
                    continue
                else:
                    lig_cons_split, arpeggio_lig_cons = process_arpeggio(struc, all_ligs, clean_pdbs_subdir, arpeggio_subdir, sifts_subdir, bio2asym_chain_dict[struc[:4]], cif2pdb_chain_dict[struc[:4]]) ### NOW PROCESSES ALL LIGANDS ###
                    print("Arpeggio output being processed for {}!".format(struc[:4]))
            ligand_contact = arpeggio_lig_cons["PDB_ResNum"].astype(str)
            ligand_contact_list.append(ligand_contact)
    
        ### BINDING SITE DEFINITION AND CHIMERA SCRIPT WRITING ###

        print("Starting BS definition section!")
        
        pdb_paths = [os.path.join(clean_pdbs_subdir, file) for file in os.listdir(clean_pdbs_subdir)]

        ligs = ligs_df.label_comp_id.unique().tolist()
        string_name = "{}_{}_BS_def_OC_{}_{}_{}".format(prot, subdir, oc_method, oc_metric, oc_dist)
        bs_def_out = os.path.join(results_dir, "{}.csv".format(string_name))
        attr_out = os.path.join(results_dir, "{}.attr".format(string_name))
        chimera_script_out = os.path.join(results_dir, "{}.com".format(string_name))
        if os.path.isfile(bs_def_out) and os.path.isfile(attr_out) and os.path.isfile(chimera_script_out):
            pass
        else:
            def_bs_oc(results_dir, pdb_paths, prot, subdir, ligs, bs_def_out, attr_out, chimera_script_out, arpeggio_subdir, metric = oc_metric, dist = oc_dist, method = oc_method, alt_fmt = False)
            print("Binding sites were sucessfully defined!")
        
        bs_definition = pd.read_csv(bs_def_out)

        ### CONSERVATION AND VARIATION ANALYSIS

        example_struc = os.path.join(clean_pdbs_subdir, os.listdir(clean_pdbs_subdir)[0])
        fasta_path = os.path.join(varalign_subdir, "{}_{}.fa".format(prot, subdir))
        hits_aln = fasta_path.replace("fa", "sto")
        hits_aln_rf = fasta_path.replace(".fa", "_rf.sto")
        shenkin_out = os.path.join(varalign_subdir, "{}_{}_shenkin.csv".format(prot, subdir))
        shenkin_filt_out = os.path.join(varalign_subdir, "{}_{}_shenkin_filt.csv".format(prot, subdir))

        if os.path.isfile(hits_aln_rf):
            pass
        else:
            create_alignment_from_struc(example_struc, fasta_path)
            print("jackhmmer was generated correctly!")

        aln_obj = Bio.AlignIO.read(hits_aln_rf, "stockholm") #crashes if target protein is not human!
        aln_info_path = os.path.join(varalign_subdir, hits_aln_rf + "_info_table.p.gz")
        if os.path.isfile(aln_info_path):
            aln_info = pd.read_pickle(aln_info_path)
        else:
            aln_info = varalign.alignments.alignment_info_table(aln_obj)
            aln_info.to_pickle(aln_info_path)
            print("Aln info was correctly created and saved!")
        aln_info_human = aln_info[aln_info.species == "HUMAN"]
        
        if len(aln_info_human) == 0: #THERE ARE NOT ANY HUMAN HOMOLOGUES FOR THE TARGET PROTEIN
            print("There are {} sequences in the MSA for {} {}. Skipping to next protein sequence".format(len(aln_info_human), prot, subdir))
            continue
            
        human_hits_msa = os.path.join(hits_aln_rf[:-4] + "_human.sto")
        if os.path.isfile(human_hits_msa):
            pass
        else:
            get_human_subset_msa(hits_aln_rf, human_hits_msa)
            print("Human subset MSA generated correctly!")

        prot_cols = get_target_prot_cols(hits_aln)
        
        variant_table_path = os.path.join(varalign_subdir, human_hits_msa + "_human_variants.p.gz")
        if os.path.isfile(variant_table_path):
            variants_table = pd.read_pickle(variant_table_path)
        else:
            variants_table = varalign.align_variants.align_variants(aln_info_human, path_to_vcf = gnomad_vcf, include_other_info = False)
            variants_table.to_pickle(variant_table_path)
            print("Variant table was created and saved correctly!")

        indexed_mapping_path = os.path.join(varalign_subdir, hits_aln_rf + '_mappings.p.gz')
        if os.path.isfile(indexed_mapping_path):
            indexed_mapping_table = pd.read_pickle(indexed_mapping_path)
        else:
            indexed_mapping_table = varalign.align_variants._mapping_table(aln_info) # now contains all species
            indexed_mapping_table.to_pickle(indexed_mapping_path) # important for merging later on
            print("Mapping table was created and saved correctly!")

        if os.path.isfile(shenkin_out):
            shenkin = pd.read_csv(shenkin_out)
        else:
            shenkin = calculate_shenkin(hits_aln_rf, aln_fmt, shenkin_out)
            print("Shenkin dataframe was created and saved correctly!")

        if os.path.isfile(shenkin_filt_out):
            shenkin_filt = pd.read_csv(shenkin_filt_out)
        else:
            shenkin_filt = get_and_format_shenkin(shenkin, prot_cols, shenkin_filt_out)
            print("Shenkin dataframe was filtered and saved correctly!")

        human_miss_vars = format_variant_table(variants_table, prot_cols) # GET ONLY MISSENSE VARIANTS ROWS

        human_miss_vars_msa_out = os.path.join(varalign_subdir, hits_aln_rf[:-4] + "_human_missense_variants_seqs.sto")

        miss_df_out = os.path.join(varalign_subdir, "{}_{}_missense_df.csv".format(prot, subdir))
        if os.path.isfile(miss_df_out):
            missense_variants_df = pd.read_csv(miss_df_out)

        else:
            missense_variants_df = get_missense_df(
                hits_aln_rf, aln_fmt, human_miss_vars,
                shenkin_filt, prot_cols, human_miss_vars_msa_out
            )
            missense_variants_df = add_miss_class(
                missense_variants_df, miss_df_out,
                cons_col = "rel_norm_shenkin", thresholds = [25, 75]
            )
            print("Missense dataframe was created and saved correctly!")
        
        shenkin_filt["human_shenkin"] = missense_variants_df.shenkin
        shenkin_filt["human_occ"] = missense_variants_df.occ
        shenkin_filt["human_gaps"] = missense_variants_df.gaps
        shenkin_filt["human_occ_pct"] = missense_variants_df.occ_pct
        shenkin_filt["human_gaps_pct"] = missense_variants_df.gaps_pct
        shenkin_filt["variants"] = missense_variants_df.variants
        shenkin_filt["oddsratio"] = missense_variants_df.oddsratio
        shenkin_filt["log_oddsratio"] = missense_variants_df.log_oddsratio
        shenkin_filt["pvalue"] = missense_variants_df.pvalue
        shenkin_filt["ci_dist"] = missense_variants_df.ci_dist
        
        aln_ids = list(set([seqid[0] for seqid in indexed_mapping_table.index.tolist() if prot in seqid[0]])) # THIS IS EMPTY IF QUERY SEQUENCE IS NOT FOUND

        mapped_data = merge_shenkin_df_and_mapping(shenkin_filt, indexed_mapping_table, aln_ids) #does it need to be only human?
        
        contact_variation_list = []
        for pdb_mapping, ligand_contact in zip(pdb_mappings, ligand_contact_list): 
            mapped_data_pdb = mapped_data.merge(pdb_mapping, on = "UniProt_ResNum") # mapping conservation and variation data to pdb mapping for each structure
            contact_variation = mapped_data_pdb[mapped_data_pdb["PDB_ResNum"].isin(ligand_contact)] # subsetting table for ligand-interacting residues
            contact_variation = contact_variation.drop(axis = 1, labels = ["PDB_ChainID"])
            contact_variation = contact_variation.drop_duplicates("UniProt_ResNum") 
            contact_variation["UniProt_ResNum"].astype(int)
            contact_variation_list.append(contact_variation)

        for structure, contact_variation in zip(strucs, contact_variation_list): # concatenate all strucutre tables 
            contact_variation['structure'] = structure[:4] # adding column indicating structure id
        all_contact_variations = pd.concat(contact_variation_list) # concatenate all different structure tables ()
        bs_ids = sorted(bs_definition.binding_site.unique().tolist()) #if not sorted, mapping of residues and binding sites has labeling mixed up
        
        binding_site_res = get_bs_residues(bs_definition, sifts_subdir, arpeggio_subdir)
        bs_sig_cols = get_bs_sig_cols(strucs, bs_ids, binding_site_res, all_contact_variations)

        fragsys_df_path = os.path.join(results_dir, "{}_{}_fragsys_df.csv".format(prot, subdir))
        
        if os.path.isfile(fragsys_df_path):
            fragsys_df = pd.read_csv(fragsys_df_path)
        else:
            fragsys_df = add_bs_info2df(bs_sig_cols, all_contact_variations, fragsys_df_path)
            print("Fragsys results dataframe was created and saved successfully!")
    
        totals = get_totals(mapped_data, prot, sifts_subdir)

        ### BINDING SITE PLOTTING SECION ###
    
        mes_sgc_df_out = os.path.join(results_dir, "{}_{}_BS_df.csv".format(prot, subdir))
        if os.path.isfile(mes_sgc_df_out):
            mes_sgc_df = pd.read_csv(mes_sgc_df_out)
        else:
            mes_sgc_df = create_binding_site_df([fragsys_df, totals, binding_site_res])
            mes_sgc_df.to_csv(mes_sgc_df_out, index = False)
            print("Binding site dataframe was created and saved successfully for {}_{}!".format(prot, subdir))
        mes_sgc_df, bs_unique_res = add_bs_msa_coverage(mes_sgc_df, fragsys_df, binding_site_res)
        df_prot_miss = pd.read_csv(miss_df_out)
        df_prot_miss = df_prot_miss[df_prot_miss.occ > 0]
        mes_sgc_df_filt = mes_sgc_df[(mes_sgc_df.vars != 0) & (mes_sgc_df.occ != 0) & (mes_sgc_df.shenkin_ci != 0)&(mes_sgc_df.number_bs_res > 1)] # filter df to get rid of sites with 0 variants and other anomalies
        sample_colors = list(itertools.islice(rgbs(), 200)) # new_colours
        bs_color_dict = {}
        for i, color in enumerate(sample_colors):
            bs_color_dict["BS"+str(i)] = color

        mes_sgc_df_filt["bs_color"] = mes_sgc_df_filt.bs_id.map(bs_color_dict) # add binding site colour column

        for k, v in bs_unique_res.items():
            sorted_bs_res = sorted(fragsys_df[fragsys_df.UniProt_ResNum.isin(v)].drop_duplicates("UniProt_ResNum").alignment_column.tolist())
            plot_binding_site(
                df_prot_miss[df_prot_miss.col.isin(sorted_bs_res)],
                cons_col = "rel_norm_shenkin",
                error = False,
                thresholds = [25, 75],
                color = [bs_color_dict["BS"+str(k)]],
                ylims = [fragsys_df.log_oddsratio.min() - 0.1, fragsys_df.log_oddsratio.max() + 0.1], 
                pltitle = "BS{} of {}".format(str(k), prot),
                out = os.path.join(figs_subdir, "{}_{}_BS{}.png".format(prot, subdir,str(k))),
                show = False
        )
        
        bs_ids = mes_sgc_df_filt.bs_id.unique().tolist()
        
        out = os.path.join(figs_subdir, "{}_{}_bss.png".format(prot, subdir))

        plot_prot_bss(
            mes_sgc_df_filt, prot,
            bs_color_dict,
            out, show = False, override = False
        )

        get_overall_stats(wd, prot, subdir, lig_data_path, bs_ids, varalign_subdir, totals)
    
        
        print("Fragsys has finished running for group {} of {}!".format(subdir, prot))

    print("Fragsys has finished running for {}!".format(prot))

### THE END ###
