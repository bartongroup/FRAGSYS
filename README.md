# FRAGSYS
This repository contains the fragment screeening analysis pipeline (**FRAGSYS**) used for the analysis of our manuscript [_Classification of likely functional state for ligand binding sites identified from fragment screening_](https://doi.org/10.21203/rs.3.rs-3185838/v1).

Our pipeline for the analysis of binding sites, **FRAGSYS**, can be executed from the jupyter notebook [`running_fragsys.ipynb`](running_fragsys.ipynb). The input for this pipeline is a table containing a series of [PDB](https://www.ebi.ac.uk/pdbe/) codes and their respective [UniProt](https://www.uniprot.org/) accession identifiers.

[![DOI](https://zenodo.org/badge/634598069.svg)](https://zenodo.org/badge/latestdoi/634598069)

## Installation

### Installation of OC

OC requires no installation, as the OC binary can be found on this repository, it is the filled called `oc`. Otherwise, the OC manual can be found [here](https://www.compbio.dundee.ac.uk/manuals/oc/oc_manual.txt), and OC can be downloaded [here](https://www.compbio.dundee.ac.uk/downloads/oc/).

### Installation of STAMP
The following instructions are to install STAMP. For more information refer to the [STAMP installation instructions](https://www.compbio.dundee.ac.uk/downloads/stamp/INSTALL).

```
# installing STAMP

## download STAMP
curl -O https://www.compbio.dundee.ac.uk/downloads/stamp/stamp.4.4.2.tar.gz

## decompress STAMP
tar -xzvf stamp.4.4.2.tar.gz
```
To install STAMP, run the BUILD script in this directory using:
```
## building STAMP
./BUILD <system-type>
```
where \<system-type\> is one of:

- linux
- osx 
- dec
- sgi
- sun

The executables will be installed in bin/\<system-type\>/.

For more information refer to the [STAMP manual](https://www.compbio.dundee.ac.uk/manuals/stamp.4.4/stamp.html)

### Installation of FRAGSYS

The first step to install **FRAGSYS** is to Git Clone the repository.

```
# git clone FRAGSYS from repository
git clone https://github.com/bartongroup/FRAGSYS.git
```

The next step is to install the three Conda environments needed to run the pipeline and analyse the results. This can be done with Conda using the .yml files in the [`envs`](envs/) directory.

```
## change directory to environments directory
cd FRAGSYS/envs

## install environments

### install main_env environment
conda env create -f main_env.yml 

### install deep_learning environment
conda env create -f deep_learning_env.yml

### install arpeggio environment
conda env create -f arpeggio_env.yml

#### activating arpeggio environment
conda activate arpeggio_env

#### change directory to main working directory
cd ../..

#### git clone Arpeggio from repository
git clone https://bitbucket.org/biomadeira/arpeggio

#### change directory to Arpeggio directory
cd arpeggio

#### test arpeggio with help function
python arpeggio.py -h
```

### Installation of VarAlign

The following instructions are to install VarAlign. Fore more information refer to the [VarAlign repository](https://github.com/bartongroup/SM_VarAlign/tree/JSU_branch).

```
#### change directory to FRAGSYS envs directory
cd ../FRAGSYS/envs/

### install varalign environment
conda env create -f varalign_env.yml

# VarAlign installation (Grabbed from URL)

#### change directory to main working directory
cd ../..

## git clone specific branch of VarAlign from repository
git clone -b JSU_branch https://github.com/bartongroup/SM_VarAlign.git

## change directory to VarAlign directory
cd SM_VarAlign

## activate varalign_env environment
conda activate varalign_env

## install VarAlign
pip install .
```

### Installation of ProIntVar

The following instructions are to install ProIntVar. Fore more information refer to the [ProIntVar repository](https://github.com/bartongroup/ProIntVar/tree/JSU_branch).

```
## ProIntVar installation (Grabbed from URL)

### change directory to main working directory
cd ..

### git clone specific branch of ProIntVar from repository
git clone -b JSU_branch https://github.com/bartongroup/ProIntVar.git

### change directory to ProIntVar directory
cd ProIntVar

### pip install ProIntVar dependencies
pip install -r requirements.txt

#then
python setup.py install

### set up ProIntVar configuration
ProIntVar-config-setup new_config.ini
```

## Pipeline methodology

Refer to  run jupyter notebook [`running_fragsys.ipynb`](running_fragsys.ipynb) in order to run **FRAGSYS**. You can do so interactively in a notebook by running this command: `main(main_dir, prot, panddas)` using the appropriate environment: [varalign_env](envs/varalign_env.yml).

Where `main_dir` is the directory where the output will be saved, `prot` is the query protein, and `panddas` is a pandas dataframe that has to contain at least two columns: `entry_uniprot_accession`, and `pdb_id`, for all protein structures in the data set.

For another example, check this other [`notebook`](https://github.com/bartongroup/FRAGSYS/blob/main/running_fragsys_for_MPRO.ipynb) where we ran **FRAGSYS** for the main protease (MPro) of SARS-CoV-2 (P0DTD1).

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
  +  [main_env](envs/main_env.yml) supports all analysis notebooks, with the exception of number [11](analysis/11_fragsys_ML_create_models.ipynb), [12](analysis/12_fragsys_ML_test_models.ipynb), in which the machine learning models are executed.
  +  [deep_learning_env](envs/deep_learning_env.yml) contains the packages necessary to do the machine learning in notebooks [11](analysis/11_fragsys_ML_create_models.ipynb), and [12](analysis/12_fragsys_ML_test_models.ipynb).

### [`input`](input/)
The input folder contains the main input file which is used as input to run **FRAGSYS** on the [running_fragsys notebook](running_fragsys.ipynb).

### [`analysis`](analysis/)
The analysis folder contains all the notebooks used to carry out the analysis of the 37 fragment screening experiments. [main_env](envs/main_env.yml) is needed to run these notebooks.

### [`results`](results/)
The results folder contains all the results files generated by the notebooks in the analysis folder.

### [`figs`](figs/)
The figs folder contains the main figures generated and saved by the analysis notebooks.

## References
1. Shenkin PS, Erman B, Mastrandrea LD. Information-theoretical entropy as a measure of sequence variability.
Proteins. 1991; 11(4):297–313. Epub 1991/01/01. https://doi.org/10.1002/prot.340110408
PMID: 1758884.
