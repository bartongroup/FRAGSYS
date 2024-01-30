# Installation

## Installation of DSSP

DSSP is incompatible with all other environments, and so must go on its environment of its own. You can install locally as well, all we will need is the location of the executable. The version of the libboost library must be specified to be this one, otherwise, dssp will not run.

```
conda create -n dssp_env salilab::dssp=3.0.0 libboost=1.73.0
```

## Installation of OC

The following instructions are to install OC. For more information refer to the OC manual [here](https://www.compbio.dundee.ac.uk/manuals/oc/oc_manual.txt), and OC can be downloaded [here](https://www.compbio.dundee.ac.uk/downloads/oc/).

```
# download OC
wget https://www.compbio.dundee.ac.uk/downloads/oc/oc-2.1a.tar.gz

# decompress OC
tar -xzvf oc-2.1a.tar.gz 

# change directory to OC main directory
cd oc-2.1a

# change getline function name
sed -i '1282s/getline/readline/; 1290s/getline/readline/; 1338s/getline/readline/' oc.c

# compile OC
gcc -o oc -O3  oc.c gjutil.c gjtimes.c -I./gjutil -I./ -lm

# check OC works
./oc
```

## Installation of STAMP
The following instructions are to install STAMP. For more information refer to the [STAMP installation instructions](https://www.compbio.dundee.ac.uk/downloads/stamp/INSTALL).

```
# download STAMP
curl -O https://www.compbio.dundee.ac.uk/downloads/stamp/stamp.4.4.2.tar.gz

# decompress STAMP
tar -xzvf stamp.4.4.2.tar.gz
```
To install STAMP, run the BUILD script in this directory using:
```
# building STAMP
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

## Installation of FRAGSYS

The first step to install **FRAGSYS** is to Git Clone the repository.

```
# git clone FRAGSYS from repository
git clone https://github.com/bartongroup/FRAGSYS.git
```

The next step is to install the three Conda environments needed to run the pipeline and analyse the results. This can be done with Conda using the .yml files in the [`envs`](envs/) directory.

```
# change directory to environments directory
cd FRAGSYS/envs

# install environments

# install main_env environment
conda env create -f main_env.yml 

# install deep_learning environment
conda env create -f deep_learning_env.yml

# install arpeggio environment
conda env create -f arpeggio_env.yml

# activating arpeggio environment
conda activate arpeggio_env

# change directory to main working directory
cd ../..

# git clone Arpeggio from repository
git clone https://bitbucket.org/biomadeira/arpeggio

# change directory to Arpeggio directory
cd arpeggio

# test arpeggio with help function
python arpeggio.py -h
```

## Installation of VarAlign

The following instructions are to install VarAlign. Fore more information refer to the [VarAlign repository](https://github.com/bartongroup/SM_VarAlign/tree/JSU_branch).

```
# change directory to FRAGSYS envs directory
cd ../FRAGSYS/envs/

# install varalign environment
conda env create -f varalign_env.yml

# VarAlign installation (Grabbed from URL)

# change directory to main working directory
cd ../..

# git clone specific branch of VarAlign from repository
git clone -b JSU_branch https://github.com/bartongroup/SM_VarAlign.git

# change directory to VarAlign directory
cd SM_VarAlign

# activate varalign_env environment
conda activate varalign_env

# install VarAlign
pip install .
```

## Installation of ProIntVar

The following instructions are to install ProIntVar. Fore more information refer to the [ProIntVar repository](https://github.com/bartongroup/ProIntVar/tree/JSU_branch).

```
# ProIntVar installation (Grabbed from URL)

# change directory to main working directory
cd ..

# git clone specific branch of ProIntVar from repository
git clone -b JSU_branch https://github.com/bartongroup/ProIntVar.git

### change directory to ProIntVar directory
cd ProIntVar

# pip install ProIntVar dependencies
pip install -r requirements.txt

#then
python setup.py install
```

## Configuration of ProIntVar

```
# change directory to VarAlign directory
cd ../SM_VarAlign

# set up ProIntVar configuration
ProIntVar-config-setup new_config.ini

# edit the following values in new_config.ini
# arpeggio_bin = /path/to/arpeggio/arpeggio.py
# python_exe = /path/to/anaconda/envs/arpeggio/bin/python
# python_path = /path/to/anaconda/envs/arpeggio/python/lib/site-packages/
# dssp_bin = /path/to/anaconda/envs/bin/mkdssp
```

## Downloading SwissProt

This is the database used for our analysis, but can be changed according to the user purposes, e.g. TrEMBL. What is important is to add the correct path in the [fragsys configuration file](fragsys_config.txt). To download SwissProt, follow the next steps.
```
# download SwissProt in fasta format
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# decrompress the file
gzip -d uniprot_sprot.fasta.gz
```

## Downloading gnomAD v2.1

This is the database used for our analysis, but can be changed according to the user purposes, e.g. v > 2.1. What is important is to add the correct path in the [fragsys configuration file](fragsys_config.txt). To download gnomAD v2.1, follow the next steps.
```
# download gnomAD Exomves vcf (large file 58GB)
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
```

For more information, refer to [gnomAD](https://gnomad.broadinstitute.org/).

After downloading gnomAD, it is required to run VEP on it, as VarAlign uses its annotations. Information on how to do so [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html).

## Configuration of FRAGSYS

Head to the [fragsys configuration file](fragsys_config.txt) and edit the following fields. Note: not all of them need changing, just showing the ones that do.
```
[binaries]
dssp_bin = /path/to/miniconda/envs/dssp_env/bin/mkdssp
stamp_bin = ./../stamp.4.4.2/bin/operative_system/stamp # change operative_system here
transform_bin = ./../stamp.4.4.2/bin/operative_system/transform # change operative_system here
clean_pdb_python_bin = /path/to/miniconda/envs/arpeggio_env/bin/python
arpeggio_python_bin = /path/to/miniconda/envs/arpeggio_env/bin/python
arpeggio_bin = ./../arpeggio/arpeggio.py

[dbs]
gnomad_vcf = /path/to/your/VEP/annotated/gnomAD/in/vcf/format
swissprot = ./../uniprot_sprot.fasta
```
