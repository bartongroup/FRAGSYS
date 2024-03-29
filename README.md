# FRAGSYS
This repository contains the fragment screeening analysis pipeline (**FRAGSYS**) used for the analysis of our manuscript [_Classification of likely functional class for ligand binding sites identified from fragment screening_](https://www.nature.com/articles/s42003-024-05970-8).

Our pipeline for the analysis of binding sites, **FRAGSYS**, can be executed from the jupyter notebook [`running_fragsys.ipynb`](running_fragsys.ipynb). The input for this pipeline is a table containing a series of [PDB](https://www.ebi.ac.uk/pdbe/) codes and their respective [UniProt](https://www.uniprot.org/) accession identifiers.

[![DOI](https://zenodo.org/badge/634598069.svg)](https://zenodo.org/badge/latestdoi/634598069)

## Installation
For complete installation instructions refer [here](INSTALL.md).

## Pipeline methodology

Refer to  run jupyter notebook [`running_fragsys.ipynb`](running_fragsys.ipynb) in order to run **FRAGSYS**. You can do so interactively in a notebook by running this command: `main(main_dir, prot, panddas)` using the appropriate environment: [varalign_env](envs/varalign_env.yml).

Where `main_dir` is the directory where the output will be saved, `prot` is the query protein, and `panddas` is a pandas dataframe that has to contain at least two columns: `entry_uniprot_accession`, and `pdb_id`, for all protein structures in the data set.

For another example, check this other [`notebook`](https://github.com/bartongroup/FRAGSYS/blob/main/running_fragsys_for_MPRO.ipynb) where we ran **FRAGSYS** for the main protease (MPro) of SARS-CoV-2 (P0DTD1).

For each structural segment of each protein in `panddas`, **FRAGSYS** will:
1. Download biological assemblies from [PDBe](https://www.ebi.ac.uk/pdbe/)
2. Structurally superimpose structures using [STAMP](http://www.compbio.dundee.ac.uk/downloads/stamp/)
3. Get accessibility and secondary structure elements from [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) via [ProIntVar](https://github.com/bartongroup/prointvar)
4. Mapping PDB residues to UniProt using [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/)
5. Obtain protein-ligand interactions running [Arpeggio](https://github.com/harryjubb/arpeggio)
6. Cluster ligands into binding sites using [OC](http://www.compbio.dundee.ac.uk/downloads/oc/)
7. Generate visualisation scripts for [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/)
8. Generate multiple sequence alignment (MSA) with [jackhmmer](http://hmmer.org/)
9. Calculate Shenkin divergence score [[1](https://doi.org/10.1002/prot.340110408)]
10. Calculate missense enrichment scores with [VarAlign](https://github.com/bartongroup/SM_varalign)

The final output of the pipeline consists of multiple tables for each structural segment collating the results from the different steps of the analysis for each residue, and for the defined ligand binding sites. These data include relative solvent accessibility (RSA), angles, secondary structure, PDB/UniProt residue number, alignment column, column occupancy, divergence score, missense enrichment score, p-value, etc.

These tables are concatenated into master tables, with data for all 37 structual segments, which form the input for the analyses carried out in the [`analysis`](analysis/) notebooks.

Refer to notebook [15](analysis/15_ML_predicting_rsa_labels.ipynb) to predict RSA cluster labels for your binding sites of interest.

## Dependencies
The pipeline, as well as the whole of the analysis are run in an interactive manner in a series of jupyter notebooks, found in the [`analysis`](analysis/) folder.

Third party dependencies for these notebooks include:
- [Arpeggio](https://github.com/harryjubb/arpeggio) [(GNU GPL v3.0 License)](https://github.com/harryjubb/arpeggio/blob/master/LICENSE)
- [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) [(Boost Software License)](https://swift.cmbi.umcn.nl/gv/dssp/)
- [Hmmer](http://hmmer.org/) [(BSD-3 Clause License)](http://eddylab.org/software/hmmer/Userguide.pdf)
- [OC](http://www.compbio.dundee.ac.uk/downloads/oc/)
- [STAMP](http://www.compbio.dundee.ac.uk/downloads/stamp/) [(GNU GPL v3.0 License)](http://www.compbio.dundee.ac.uk/manuals/stamp.4.4/stamp.html)
- [ProIntVar](https://github.com/bartongroup/prointvar) [(MIT License)](https://github.com/bartongroup/ProIntVar/blob/master/LICENSE.md)
- [ProteoFAV](https://github.com/bartongroup/ProteoFAV) [(MIT License)](https://github.com/bartongroup/ProteoFAV/blob/master/LICENSE.md)
- [VarAlign](https://github.com/bartongroup/SM_varalign) [(MIT License)](https://github.com/bartongroup/SM_VarAlign/blob/master/LICENSE)

Other standard python libraries:
- [Biopython](https://biopython.org/) [(BSD 3-Clause License)](https://github.com/biopython/biopython/blob/master/LICENSE.rst)
- [Keras](https://keras.io/) [(Apache v2.0 License)](https://github.com/keras-team/keras/blob/master/LICENSE)
- [Matplotlib](https://matplotlib.org/) [(PSF License)](https://github.com/matplotlib/matplotlib/blob/main/LICENSE/LICENSE)
- [Numpy](https://numpy.org/) [(BSD 3-Clause License)](https://github.com/numpy/numpy/blob/main/LICENSE.txt)
- [Pandas](https://pandas.pydata.org/) [(BSD 3-Clause License)](https://github.com/pandas-dev/pandas/blob/main/LICENSE)
- [Scipy](https://scipy.org/) [(BSD 3-Clause License)](https://github.com/scipy/scipy/blob/main/LICENSE.txt)
- [Seaborn](https://seaborn.pydata.org/) [(BSD 3-Clause License)](https://github.com/mwaskom/seaborn/blob/master/LICENSE.md)
- [Scikit-learn](https://scikit-learn.org/stable/) [(BSD 3-Clause License)](https://github.com/scikit-learn/scikit-learn/blob/main/COPYING)
- [Tensorflow](https://www.tensorflow.org/) [(Apache v2.0 License)](https://github.com/tensorflow/tensorflow/blob/master/LICENSE)

For more information on the dependencies, refere to the .yml files in the [`envs`](envs/) directory. To install all the dependencies, refer to the [installation manual](INSTALL.md).

## Files
Apart from the [INSTALL](INSTALL.md), [LICENSE](LICENSE) and [README](README.md) files, there are 5 other files on this repository main directory. Two of these are python libraries, a configuration file and two notebooks.
  +  [`fragsys_config.txt`](fragsys_config.txt) contains the default parameters to run **FRAGSYS** and it is read by [`fragsys.py`](fragsys.py).
  +  [`fragsys.py`](fragsys.py) contains all the function, lists and dictionaries needed to run the pipeline.
  +  [`fragsys_main.py`](fragsys_main.py) contains the main **FRAGSYS** function, where all functions in [`fragsys.py`](fragsys.py) are called. This script represents the pipeline itself.
  +  [`running_fragsys.ipynb`](running_fragsys.ipynb) is the notebook where the pipeline is executed in an interactive way.
  +  [`running_fragsys_for_MPRO.ipynb.ipynb`](running_fragsys_for_MPRO.ipynb.ipynb) is the notebook where the pipeline is executed in an interactive way for a case study of SARS-CoV-2 MPro.
## Directories
There are 6 directories in this repository.

### [`scripts`](scripts/clean_pdb.py)
This environment contains [clean_pdb.py](scriptscripts/clean_pdb.py), a python script grabbed from [here](https://github.com/harryjubb/pdbtools/blob/master/clean_pdb.py). This script will be used to pre-process the PDB files before running Arpeggio on them.

### [`envs`](envs/)
The envs folder contains three .yml files describing the necessary packages and dependencies for the different parts of the pipeline and analysis.
  +  [arpeggio_env](envs/arpeggio_env.yml) contains Arpeggio.
  +  [deep_learning_env](envs/deep_learning_env.yml) contains the packages necessary to do the machine learning in notebooks [11](analysis/11_fragsys_ML_create_models.ipynb), and [12](analysis/12_fragsys_ML_test_models.ipynb).
  +  [main_env](envs/main_env.yml) supports all analysis notebooks, with the exception of number [11](analysis/11_fragsys_ML_create_models.ipynb), [12](analysis/12_fragsys_ML_test_models.ipynb), in which the machine learning models are executed.
  +  [varalign_env](envs/varalign_env.yml) is needed to run **FRAGSYS**.

### [`input`](input/)
The input folder contains the main input file which is used as input to run **FRAGSYS** on the [running_fragsys notebook](running_fragsys.ipynb).

### [`analysis`](analysis/)
The analysis folder contains all the notebooks used to carry out the analysis of the 37 fragment screening experiments. [main_env](envs/main_env.yml) is needed to run these notebooks.

### [`results`](results/)
The results folder contains all the results files generated by the notebooks in the analysis folder.

### [`figs`](figs/)
The figs folder contains the main figures generated and saved by the analysis notebooks.

## Citation

If you use *FRAGSYS*, please cite:

Utgés, J.S. et al. Classification of likely functional class for ligand binding sites identified from fragment screening. Commun Biol 7, 320 (2024). https://doi.org/10.1038/s42003-024-05970-8

## References
1. Shenkin PS, Erman B, Mastrandrea LD. Information-theoretical entropy as a measure of sequence variability.
Proteins. 1991; 11(4):297–313. Epub 1991/01/01. https://doi.org/10.1002/prot.340110408
PMID: 1758884.
