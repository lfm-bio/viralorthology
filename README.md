# ViralOrthology

ViralOrthology es una herramienta para detectar ortólogos y parálogos con análisis de sintenia. Utiliza programas como proteinortho, blastp, HMMER, ORFfinder.
Input: genomas, orfeomas y proteomas. Output: Grupos de ortología, grupos de paralogía.

## Features

* descarga genomas, proteomas y orfeomas a partir de lista de numeros de acceso.
* actualiza la base de datos de genomas, orfeomas y proteomas. Descarga nuevas entradas de genbank, realiza el análisis de ortología de sus proteinas y las agrega a los grupos preexistentes.
* genera los archivos necesarios para hacer el análisis de kimura.

## Dependencies (TODAS SE INSTALAN AUTOMATICAMENTE (VER INSTALACIÓN))

* (EntrezDirect)[https://www.ncbi.nlm.nih.gov/books/NBK179288/]
* (BLAST)[https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html]
* (ORFfinder)[https://www.ncbi.nlm.nih.gov/orffinder/]
* (ProteinOrtho)[https://gitlab.com/paulklemm_PHD/proteinortho]
* (ProteinOrtho_grab_proteins)[https://gitlab.com/paulklemm_PHD/proteinortho]
* (HMMER)[http://hmmer.org/]
* (Muscle)[https://github.com/rcedgar/muscle]
* (hhsuite)[https://github.com/soedinglab/hh-suite]
* (BioPython)[https://biopython.org/]
* (tqdm)[https://tqdm.github.io/]

## Instalación

* Debian-based Linux distributions:

Correr las siguientes lineas en la consola:

```
wget https://raw.githubusercontent.com/lfm-bio/viralorthology/main/setup_pipeline.sh
./setup_pipeline.sh
```
## BEFORE USING

Run the following command:

```
viralorthology -check_dependencies
```

## Descripción

Todos los análisis se hacen utilizando las secuencias aminoacidicas.
(Grafico del pipeline)

## Parametros probados

* ORFfinder: -ml, -s, -n
* BASTP: -evalue, -qcov_hsp_perc
* Proteinortho: --e, --p

## Parametros que no se pueden tocar (son utilizados internamente por el pipeline)
* ORFfinder: -in, -out, -outfmt
* Proteinortho: --project
* SearchParalogs: -query, -db, -out 
* blastp: -query, -db, -out

## Ejemplo

```
viralorthology --orffinder -ml 90 -s 0 -n true --blastp -evalue 0.001 -qcov_hsp_perc 70 --proteinortho --p=blastp+
```

## Comandos

* -try_to_merge_groups fasta1 fasta2 fastaN: checkea si puede unir los grupos y los une [les pone el nombre del primero de la lista] 
* -prots_per_group: checkea en cada fasta si algun grupo (A, B, G, etc) está completo.
* -update_db (-query, -minlenght, -mindate)
* -download_seqs: descarga genomas, orfeomas y proteomas. (Lista de numers de acceso en ids.txt, uno por linea)
* -kimura: prepara el fasta con los orf para hacer análisis de k2p

## fastas iniciales:

* download_seqs

* Descarga a mano de genbank

genomes.fasta: send to > complete record > file > format: FASTA
orfeomes.fasta: send to > coding sequences > format: FASTA nucleotide
proteomes.fasta: send to > coding sequences > format: FASTA nucleotide
