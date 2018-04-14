#!/bin/bash

# unload potentially conflicting anaconda instances
{ # try
    module unload python-anaconda2 &&
    module unload python-anaconda3
} || { # catch
    echo 'module unloading failed: maybe module does not exist'
}


# install miniconda for local independent package management
wget http://repo.continuum.io/miniconda/Miniconda2-4.3.21-Linux-x86_64.sh -O miniconda.sh
mkdir dependencies
chmod 775 miniconda.sh
chmod 775 dependencies
bash miniconda.sh -b -p ./dependencies/miniconda
rm miniconda.sh

# activate conda virtual environment
source ./dependencies/miniconda/bin/activate

pip install click

# add bioconda and r channel for easy dependency installations
conda install -c bioconda bbmap megahit prodigal prokka emirge

# install pyflow for automated task management
wget https://github.com/Illumina/pyflow/releases/download/v1.1.17/pyflow-1.1.17.tar.gz
pip install pyflow-1.1.17.tar.gz
rm pyflow-1.1.17.tar.gz

# Virsorter Install https://github.com/simroux/VirSorter
wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
tar -xvzf virsorter-data-v2.tar.gz -C dependencies/

wget http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz
tar -xvzf metagene/mga_x86_64.tar.gz -C dependencies/

git clone https://github.com/simroux/VirSorter.git
mv VirSorter dependencies/
cd dependencies/VirSorter/Scripts && make && cd ../../
ln -s dependencies/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl dependencies/miniconda/bin
ln -s dependencies/VirSorter/Scripts dependencies/miniconda/bin

# install VirHostMatcher
git clone https://github.com/jessieren/VirHostMatcher.git
mv VirHostMatcher dependencies/

# download adapters for trimming
wget https://raw.githubusercontent.com/BioInfoTools/BBMap/master/resources/adapters.fa -O dependencies/adapters.fa
