#bash

#reinstall pipeline
if [ -d "/usr/local/bin/vorthology" ]; then
    sudo rm -r /usr/local/bin/vorthology
    sudo rm /usr/local/bin/viralorthology
fi

#basic linux tools
sudo apt-get install wget -y
sudo apt-get install gzip -y
sudo apt-get install unzip -y

#pipeline
wget https://github.com/lfm-bio/viralorthology/archive/refs/heads/main.zip
unzip main.zip
cd viralorthology-main
sudo mkdir /usr/local/bin/vorthology
sudo cp -rp modules /usr/local/bin/vorthology/modules
mv viralorthology.py viralorthology
sudo cp viralorthology /usr/local/bin/vorthology
mv viralorthology viralorthology.py
sudo ln -s /usr/local/bin/vorthology/viralorthology /usr/local/bin/viralorthology #symbolic link
sudo chmod a+rx /usr/local/bin/viralorthology

#EntrezDirect 
sudo apt install ncbi-entrez-direct
# wget https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh && bash install-edirect.sh

#Blast
sudo apt-get install ncbi-blast+ -y

#ORFfinder [https://www.ncbi.nlm.nih.gov/orffinder/]
wget https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz
gzip -d ORFfinder.gz
sudo mv ORFfinder /usr/local/bin/ORFfinder
sudo chmod a+rx /usr/local/bin/ORFfinder

#proteinortho [https://gitlab.com/paulklemm_PHD/proteinortho]
sudo apt-get install proteinortho -y

#clustal [http://www.clustal.org/]
sudo apt-get install clustalw -y

#proteinortho_grab_proteins [https://gitlab.com/paulklemm_PHD/proteinortho]
sudo wget https://gitlab.com/paulklemm_PHD/proteinortho/-/raw/master/src/proteinortho_grab_proteins.pl -O /usr/local/bin/po_grab_proteins
sudo chmod a+rx /usr/local/bin/po_grab_proteins

#HMMER [http://hmmer.org/]
sudo apt-get install hmmer -y

#Muscle [https://github.com/rcedgar/muscle]
wget https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.linux_intel64
mv muscle5.1.linux_intel64 muscle
sudo mv muscle /usr/local/bin/muscle
sudo chmod a+rx /usr/local/bin/muscle

#hhsuite [https://github.com/soedinglab/hh-suite]
sudo apt-get install hhsuite -y
sudo wget https://raw.githubusercontent.com/soedinglab/hh-suite/master/scripts/reformat.pl -O /usr/local/bin/reformat_msa
sudo chmod a+rx /usr/local/bin/reformat_msa

#python libraries
pip install biopython
pip install tqdm