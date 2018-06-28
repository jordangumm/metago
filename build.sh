#!/bin/bash

# unload potentially conflicting anaconda instances
{ # try
    module unload python-anaconda2 &&
    module unload python-anaconda3
} || { # catch
    echo 'module unloading failed: maybe module does not exist'
}


# install miniconda for local independent package management
wget https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh -O miniconda.sh
mkdir dependencies
chmod 775 miniconda.sh
chmod 775 dependencies
bash miniconda.sh -b -p ./dependencies/miniconda
rm miniconda.sh

# activate conda virtual environment
source ./dependencies/miniconda/bin/activate

pip install click
pip install git+https://github.com/jordangumm/pyleup.git 

# add bioconda and r channel for easy dependency installations
conda install -c bioconda bbmap megahit prodigal prokka emirge pysam maxbin2 kallisto
conda install -c conda-forge conda-execute

# install pyflow for automated task management
pip install https://github.com/Illumina/pyflow/releases/download/v1.1.17/pyflow-1.1.17.tar.gz

# Virsorter Install https://github.com/simroux/VirSorter
wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
tar -xvzf virsorter-data-v2.tar.gz -C dependencies/
rm virsorter-data-v2.tar.gz

wget http://metagene.nig.ac.jp/metagene/mga_x86_64.tar.gz
tar -xvzf mga_x86_64.tar.gz -C dependencies/
rm mga_x86_64.tar.gz

git clone https://github.com/simroux/VirSorter.git
mv VirSorter dependencies/
cd dependencies/VirSorter/Scripts && make && cd ../../../
ln -s dependencies/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl dependencies/miniconda/bin
ln -s dependencies/VirSorter/Scripts dependencies/miniconda/bin

# install VirHostMatcher
git clone https://github.com/jessieren/VirHostMatcher.git
mv VirHostMatcher dependencies/

# download adapters for trimming
wget https://raw.githubusercontent.com/BioInfoTools/BBMap/master/resources/adapters.fa -O dependencies/adapters.fa

conda create -n py3 python=3.5
source activate py3

pip install -U pip
pip install -U Cython
pip install -U jupyter jupyter_client ipython pandas matplotlib scipy scikit-learn khmer

pip install -U https://github.com/dib-lab/sourmash/archive/master.zip
