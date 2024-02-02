#!/usr/bin/python3
import os

#reinstall pipeline
if os.path.isdir('/usr/local/bin/vorthology'):
    os.system('sudo rm -r /usr/local/bin/vorthology')
    try:
        os.system('sudo rm /usr/local/bin/viralorthology')
    except:
        pass

os.system('sudo apt-get install wget -y')
os.system('sudo apt-get install gzip -y')
os.system('sudo apt-get install unzip -y')

#pipeline
os.system('wget https://github.com/lfm-bio/viralorthology/archive/refs/heads/main.zip')
os.system('unzip main.zip')
os.chdir('viralorthology-main')
os.system('sudo mkdir /usr/local/bin/vorthology')
os.system('sudo cp -rp modules /usr/local/bin/vorthology/modules')
os.rename('viralorthology.py', 'viralorthology')
os.system('sudo cp viralorthology /usr/local/bin/vorthology')
os.rename('viralorthology', 'viralorthology.py')
os.system('sudo ln -s /usr/local/bin/vorthology/viralorthology /usr/local/bin/viralorthology') #symbolic link
os.system('sudo chmod a+rx /usr/local/bin/viralorthology')

#EntrezDirect
os.system('sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"')

#Blast
os.system('sudo apt-get install ncbi-blast+')

#ORFfinder [https://www.ncbi.nlm.nih.gov/orffinder/]
os.system('wget https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz')
os.system('gzip -d ORFfinder.gz')
os.system('sudo mv ORFfinder /usr/local/bin/ORFfinder')
os.system('sudo chmod a+rx /usr/local/bin/ORFfinder')

#proteinortho [https://gitlab.com/paulklemm_PHD/proteinortho]
os.system('sudo apt-get install proteinortho -y')

#clustal [http://www.clustal.org/]
os.system('sudo apt-get install clustalw -y')

#proteinortho_grab_proteins [https://gitlab.com/paulklemm_PHD/proteinortho]
os.system('sudo wget https://gitlab.com/paulklemm_PHD/proteinortho/-/raw/master/src/proteinortho_grab_proteins.pl -O /usr/local/bin/po_grab_proteins')
os.system('sudo chmod a+rx /usr/local/bin/po_grab_proteins')

#HMMER [http://hmmer.org/]
os.system('sudo apt-get install hmmer -y')

#Muscle [https://github.com/rcedgar/muscle]
os.system('wget https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.linux_intel64')
os.rename('muscle5.1.linux_intel64', 'muscle')
os.system('sudo mv muscle /usr/local/bin/muscle')
os.system('sudo chmod a+rx /usr/local/bin/muscle')

#hhsuite [https://github.com/soedinglab/hh-suite]
os.system('sudo apt-get install hhsuite -y')
os.system('sudo wget https://raw.githubusercontent.com/soedinglab/hh-suite/master/scripts/reformat.pl -O /usr/local/bin/reformat_msa')
os.system('sudo chmod a+rx /usr/local/bin/reformat_msa')

#python libraries
os.system('pip install biopython')
os.system('pip install tqdm')