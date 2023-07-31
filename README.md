# FRAGSYS
This repository contains the fragment screeening analysis pipeline (**FRAGSYS**) used for the analysis of our manuscript [_Classification of likely functional state for ligand binding sites identified from fragment screening_](https://doi.org/10.21203/rs.3.rs-3185838/v1).

Our pipeline for the analysis of binding sites, **FRAGSYS**, can be executed from the jupyter notebook [`running_fragsys.ipynb`](running_fragsys.ipynb). The input for this pipeline is a list of protein [UniProt](https://www.uniprot.org/) accession identifiers.

## Pipeline methodology

This is how to run **FRAGSYS**: `main(main_dir, prot, panddas)`

Where `main_dir` is the directory where the output will be saved, `prot` is the query protein, and `panddas` is a pandas dataframe that has to contain at least two columns: `entry_uniprot_accession`, and `pdb_id`, for all protein structures in the data set.

For each structural segment of each protein in `panddas`, **FRAGSYS** will:
1. Download biological assemblies from [PDBe](https://www.ebi.ac.uk/pdbe/)
2. Structurally superimpose structures using [STAMP](http://www.compbio.dundee.ac.uk/downloads/stamp/)
3. Get accessibility and secondary structure elements from [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) via [ProIntVar](https://github.com/bartongroup/prointvar)<sup>*</sup>
4. Mapping PDB residues to UniProt using [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/)
5. Obtain protein-ligand interactions running [Arpeggio](https://github.com/harryjubb/arpeggio)
6. Cluster ligands into binding sites using [OC](http://www.compbio.dundee.ac.uk/downloads/oc/)
7. Generate visualisation scripts for [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/)
8. Generate multiple sequence alignment (MSA) with [jackhmmer](http://hmmer.org/)
9. Calculate Shenkin divergence score [[1](https://doi.org/10.1002/prot.340110408)]
10. Calculate missense enrichment scores with [VarAlign](https://github.com/bartongroup/SM_varalign)<sup>*</sup>

The final output of the pipeline consists of multiple tables for each structural segment collating the results from the different steps of the analysis for each residue, and for the defined ligand binding sites. These data include relative solvent accessibility (RSA), angles, secondary structure, PDB/UniProt residue number, alignment column, column occupancy, divergence score, missense enrichment score, p-value, etc.

These tables are concatenated into master tables, with data for all 37 structual segments, which form the input for the analyses carried out in the [`analysis`](analysis/) notebooks.

## Dependencies
The pipeline, as well as the whole of the analysis are run in an interactive manner in a series of jupyter notebooks, found in the [`analysis`](analysis/) folder.

Third party dependencies for these notebooks include:
- [Arpeggio](https://github.com/harryjubb/arpeggio)
- [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/)
- [Hmmer](http://hmmer.org/)
- [OC](http://www.compbio.dundee.ac.uk/downloads/oc/)
- [STAMP](http://www.compbio.dundee.ac.uk/downloads/stamp/)
- [ProIntVar](https://github.com/bartongroup/prointvar)<sup>*</sup>
- [ProteoFAV](https://github.com/bartongroup/ProteoFAV)<sup>*</sup>
- [VarAlign](https://github.com/bartongroup/SM_varalign)<sup>*</sup>

<sup>*</sup> This is an in-house python library that will be publicly available soon.

Other standard python libraries:
- [Biopython](https://biopython.org/)
- [Keras](https://keras.io/)
- [Matplotlib](https://matplotlib.org/)
- [Numpy](https://numpy.org/)
- [Pandas](https://pandas.pydata.org/)
- [Scipy](https://scipy.org/)
- [Seaborn](https://seaborn.pydata.org/)
- [Scikit-learn](https://scikit-learn.org/stable/)
- [Tensorflow](https://www.tensorflow.org/)

For more information, refere to the .yaml files in the [`envs`](envs/) directory.

## Files
There are 5 different files on this repository main directory. Three of these are python libraries, then a notebook, and a configuration file.
  +  [`fragsys_config.txt`](fragsys_config.txt) contains the default parameters to run **FRAGSYS** and it is read by [`fragsys.py`](fragsys.py).
  +  [`fragsys.py`](fragsys.py) contains all the function, lists and dictionaries needed to run the pipeline.
  +  [`fragsys_main.py`](fragsys_main.py) contains the main **FRAGSYS** function, where all functions in [`fragsys.py`](fragsys.py) are called. This script represents the pipeline itself.
  +  [`running_fragsys.ipynb`](running_fragsys.ipynb) is the notebook where the pipeline is executed in an interactive way.
## Directories
There are five directories in this repository.

### [`envs`](envs/)
The envs folder contains three .yml files describing the necessary packages and dependencies for the different parts of the pipeline and analysis.
  +  [varalign_env](envs/varalign_env.yml) is needed to run **FRAGSYS**.
  +  [myenv2](envs/myenv2.yml) supports all analysis notebooks, with the exception of number [11](analysis/11_fragsys_ML_create_models.ipynb), [12](analysis/12_fragsys_ML_test_models.ipynb), in which the machine learning models are executed.
  +  [deep_learning](envs/deep_learning.yml) contains the packages necessary to do the machine learning in notebooks [11](analysis/11_fragsys_ML_create_models.ipynb), and [12](analysis/12_fragsys_ML_test_models.ipynb).

### [`input`](input/)
The input folder contains the main input file which is used as input to run **FRAGSYS** on the [running_fragsys notebook](running_fragsys.ipynb).

### [`analysis`](analysis/)
The analysis folder contains all the notebooks used to carry out the analysis of the 37 fragment screening experiments. [myenv2](envs/myenv2.yml) is needed to run these notebooks.

### [`results`](results/)
The results folder contains all the results files generated by the notebooks in the analysis folder.

### [`figs`](figs/)
The figs folder contains the main figures generated and saved by the analysis notebooks.

## References
1. Shenkin PS, Erman B, Mastrandrea LD. Information-theoretical entropy as a measure of sequence variability.
Proteins. 1991; 11(4):297–313. Epub 1991/01/01. https://doi.org/10.1002/prot.340110408
PMID: 1758884.
