# Antimicrobioresistance – Formation Africa CDC Juillet 2025

## Table des matières

- [Workflow](#workflow)
- [1) Mise à jour du système](#1-mise-à-jour-du-système)
- [2) Installer Miniconda](#2-installer-miniconda)
- [3) Créer un environnement conda](#3-créer-un-environnement-conda)
- [4) Installer les outils bioinfo](#4-installer-les-outils-bioinfo)
- [5) Contrôle qualité](#5-contrôle-qualité)
- [6) Correction des reads](#6-correction-des-reads)
- [7) Assemblage de novo](#7-assemblage-de-novo)
- [8) Identification de l'espèce](#8-identification-de-lespèce)
- [9) Évaluation de l'assemblage](#9-évaluation-de-lassemblage)
- [10) Détection de la contamination](#10-détection-de-la-contamination-avec-kraken2-et-bracken)
- [11) Évaluation de la complétude](#11-évaluation-de-la-complétude-du-génome)
- [12) Typage E.coli et Salmonella](#12-typage-ecoli-et-salmonella)
- [13) Installation de abricate](#13-installation-de-abricate-pour-recherche-de-gènes-de-résistance-et-facteurs-de-virulence)
- [14) Recherche avec abricate](#14-recherche-de-gènes-de-résistance-facteurs-de-virulence-et-plasmides-avec-abricate)
- [15) Recherche avec ResFinder](#15-recherche-de-gènes-de-résistance-et-mutations-avec-resfinder)
- [16) Analyse de la MLST](#16-analyse-de-la-mlst)
- [Exercice 2](#exercice-2--reprendre-le-workflow-bactérien-à-partir-des-reads-téléchargés)
- [18) Annotation du génome](#18-annotation-du-génome)
- [Commandes avec Docker](#commandes-avec-docker)

---

## Workflow

<img width="1542" height="731" alt="DIC1 drawio" src="https://github.com/user-attachments/assets/365e4b01-4ac4-4e0e-9239-6e72c6d3340f" />


---

## 1) Mise à jour du système

```bash
sudo apt update && sudo apt upgrade -y
```

**Explication des commandes :**
- **sudo** : permet de donner tous les droits
- **apt update** : Met à jour la liste des paquets disponibles
- **apt upgrade** : Met à jour les paquets installés s'il existe une nouvelle version

---

## 2) Installer Miniconda

**⚠️ UNIQUEMENT POUR PREMIÈRE UTILISATION DE LINUX**

```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

**Explication des commandes :**
- **wget** : Permet de télécharger un fichier grâce à un lien de téléchargement
- **bash** : Exécute le script d'installation

---

## 3) Créer un environnement conda

```bash
conda create -n amr -c bioconda -c conda-forge -y
conda create -n busco -c bioconda -c conda-forge -y
conda create -n abricate -c bioconda -c conda-forge -y
conda create -n mlst -c bioconda -c conda-forge -y
conda create -n seqsero2 -c bioconda -c conda-forge seqsero2=1.3.1 -y
```

**Explication des paramètres :**
- **conda create** : Permet de créer un environnement conda
- **-n** : donne le nom de l'environnement
- **-c** : ajoute un channel dans l'environnement
- **-y** : répondre "yes" à tout

---

## 4) Installer les outils bioinfo

```bash
# Active l'environnement amr
conda activate amr

# Installe les logiciels avec des versions spécifiques
conda install -c bioconda -c conda-forge fastp=1.0.1 mash=2.3 nanoplot=1.44.1 resfinder sra-tools=3.2.1 ectyper=2.0.0 -y

# Vérification des installations
fastp --version
mash --version
NanoPlot --version
ectyper --version
prefetch --version
python -m resfinder --version

# Désactive l'environnement amr
conda deactivate
```

> ** Note :** Il faut au préalable activer l'environnement conda dans lequel vous avez installé le logiciel

---

## 5) Contrôle qualité

### Pour lancer fastQC sur des fichiers spécifiques :

```bash
fastqc S1_R1.fastq.gz
fastqc S1_R2.fastq.gz
```

### Pour lancer fastQC sur tous les reads du dossier de travail :

```bash
fastqc *fastq.gz
```

** Ressources :** [Interprétation du rapport fastQC - Modules d'analyse de FastQC (site officiel)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)

---

## 6) Correction des reads

```bash
fastp -i S1_R1.fastq.gz -I S1_R2.fastq.gz -o trim_S1_R1.fastq.gz -O trim_S1_R2.fastq.gz -f 18 -F 18 -t 18 -T 18 -M 20 --detect_adapter_for_pe --dedup
```

**Explication des paramètres :**
- **-i** : input read 1
- **-I** : input read 2
- **-o** : nom de sortie du read1
- **-O** : nom de sortie du read2
- **-f** : nombre de bases à couper au début de chaque séquence read1
- **-F** : nombre de bases à couper au début de chaque séquence read2
- **-t** : nombre de bases à couper à la fin de chaque séquence read1
- **-T** : nombre de bases à couper à la fin de chaque séquence read2
- **-M** : the mean quality requirement
- **--detect_adapter_for_pe** : détection automatique des adaptateurs
- **--dedup** : supprimer les séquences dupliquées

---

## 7) Assemblage de novo

```bash
spades.py -1 trim_S1_R1.fastq.gz -2 trim_S1_R2.fastq.gz --only-assembler -o S1_assembled
```

**Explication des paramètres :**
- **-1** : Reads 1
- **-2** : Reads 2
- **--only-assembler** : Assemblage sans étape de correction des reads
- **-o** : Nom du dossier de sortie

---

## 8) Identification de l'espèce

### Création du sketch

```bash
conda activate amr
# Création du sketch
mash sketch S1_assembled/scaffolds.fasta -k 16

# Faire le mash sur plusieurs fichiers fasta
for i in *.fasta; do mash sketch $i; done
```

### Comparaison avec la base de données de référence

```bash
# Comparaison avec la base de référence
mash dist refseq.msh scaffolds.msh > mash_output.csv

# Comparaison pour partir plusieurs fichiers .msh
for i in *.fasta.msh; do mash dist refseq.msh $i > dist-$i.csv; done
```

### Identification de la référence la plus proche

```bash
# Afficher le numéro d'accession de la référence identifiée la plus proche
sort -k 3 -i mash_output.csv | head -n 1 | awk '{print $1}'
```

---

## 9) Évaluation de l'assemblage

```bash
quast S1_assembled/scaffolds.fasta -1 trim_S1_R1.fastq.gz -2 trim_S1_R2.fastq.gz -o quast_resultats -r S1_ref.fasta
```

**Explication des paramètres :**
- **-1** : Reads 1 (non obligatoire)
- **-2** : Reads 2 (non obligatoire)
- **-o** : Dossier des résultats
- **-r** : Référence au format fasta (non obligatoire)

---

## 10) Détection de la contamination avec Kraken2 et Bracken

```bash
# Classification taxonomique avec Kraken2
kraken2 --db /chemin/vers/kraken_db --paired \
        cleaned_R1.fastq.gz cleaned_R2.fastq.gz \
        --output kraken.out --report kraken.report

kraken2 --db /chemin/vers/kraken --output kraken.out --report kraken.report file.fasta

# Estimation des abondances avec Bracken (niveau espèce)
bracken -d /chemin/vers/kraken_db \
        -i kraken.report -o bracken.species.report \
        -l S -r 150
```

---

## 11) Évaluation de la complétude du génome

```bash
# Avec BUSCO
conda activate busco

# Installation de busco
conda install -c bioconda busco -y

# Lancer l'analyse busco
busco -i S1_assembled/scaffolds.fasta -o busco -l bacteria_odb10 -m genome
conda deactivate
```

**Explication des paramètres :**
- **-i** : Fichier fasta en input
- **-o** : Nom du dossier de sortie
- **-l** : Base de données à utiliser (pour avoir toutes les bases de données faire : `busco --list-datasets`)
- **-m** : Type de séquence (genome, proteins ou transcriptome)

---

## 12) Typage E.coli et Salmonella

```bash
# Activer l'environnement où se trouve ECtyper
conda activate amr

# Lancer l'analyse d'ectyper
ectyper -i S1_assembled/scaffolds.fasta -o ectyper
conda deactivate

# Lancer ECtyper sur plusieurs fichiers fasta
ectyper -i *.fasta -o ectyper

# Activer l'environnement seqsero2
conda activate seqsero2

# Lancer SeqSero sur des fasta
SeqSero2_package.py -i salmonella.fasta -m k -t 4
conda deactivate
```

**Pour seqsero2 (peut aussi marcher pour les reads) :**
- **-i** : fichier input fasta
- **-m** : recherche k-mer based
- **-t** : Type d'entrée attendu (ici génome donc =4)

---

## 13) Installation de abricate pour recherche de gènes de résistance et facteurs de virulence

```bash
conda activate abricate

# Installation d'abricate
conda install -c bioconda -c conda-forge abricate

# Mise à jour des bases de données d'abricate
abricate --setupdb

# Affiche les bases de données d'abricate
abricate --list
```

### Bases de données disponibles dans abricate

| DATABASE | SEQUENCES | DBTYPE |
|----------|-----------|--------|
| megares | 6635 | nucl |
| card | 2631 | nucl |
| argannot | 22231 | nucl |
| resfinder | 3077 | nucl |
| plasmidfinder | 460 | nucl |
| ncbi | 5386 | nucl |
| ecoli_vf | 2701 | nucl |
| vfdb | 2597 | nucl |
| ecoh | 597 | nucl |

---

## 14) Recherche de gènes de résistance facteurs de virulence et plasmides avec abricate

```bash
# Lance abricate avec la base de données par défaut (ncbi)
abricate S1_assembled/scaffolds.fasta > abricate.csv

# Lance abricate avec la base de données resfinder
abricate --db resfinder S1_assembled/scaffolds.fasta > abricate.csv

# Lancer abricate sur plusieurs fichiers fasta
abricate --db resfinder *.fasta > resfinder.csv

# Lançons abricate sur une souche multirésistante (MDR = Multi Drug Resistant)
abricate --db resfinder kpneumoniae_mdr.fasta

# Recherche la présence de plasmides
abricate --db plasmidfinder S1_assembled/scaffolds.fasta

# Lance abricate avec la base vfdb (virulence factors database)
abricate --db vfdb S1_assembled/scaffolds.fasta > abricate_virulence.csv
```

**Explication des paramètres :**
- **--db** : spécifie la base de données à utiliser
- **>** : permet de rediriger le résultat de sortie vers un fichier

---

## 15) Recherche de gènes de résistance et mutations avec ResFinder

```bash
# Télécharger la base de données de Resfinder
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder

# Télécharger la base de données de Pointfinder
git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git db_pointfinder

# Installer la base de données de Resfinder
cd db_resfinder
python3 INSTALL.py

# Installer la base de données de Pointfinder
cd ../db_pointfinder
python3 INSTALL.py

cd ../

# Analyse Resfinder
python -m resfinder -ifa S1_assembled/scaffolds.fasta -o resultats_resfinder --point --acquired --species "Klebsiella pneumoniae" -t 80 -l 0.5
```

---

## 16) Analyse de la MLST

```bash
conda activate mlst

# Installer la version 2.23.0 de mlst
conda install -c bioconda mlst=2.23.0 -y

# Lancer l'analyse de la mlst
mlst S1_assembled/scaffolds.fasta > mlst.csv

# Lancer l'analyse MLST sur plusieurs fichiers fasta
mlst *.fasta > mlst.csv
conda deactivate
```

---

## Exercice 1 : Reprendre le workflow bactérien à partir des reads téléchargés

```bash
conda activate amr

# Télécharger le fichier SRA avec suivi de progression
prefetch SRR34342573 --progress

# Extraire les fichiers FASTQ à partir du fichier .sra téléchargé
fasterq-dump SRR34342573/SRR34342573.sra --split-files --threads 4 --outdir ./

# Supprimer le dossier SRA après extraction
rm -r SRR34342573

# Compresser les fichiers FASTQ avec gzip
gzip SRR34342573_1.fastq
gzip SRR34342573_2.fastq
conda deactivate
```

---

## 18) Annotation du génome

**🔗 Outil d'annotation complète du génome :** [Bakta (site officiel)](https://bakta.computational.bio/submit)

---

## Commandes avec Docker

### Installation de Docker

```bash
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "${UBUNTU_CODENAME:-$VERSION_CODENAME}") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

# Télécharger la dernière version
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# Ajouter votre user au docker group
sudo usermod -aG docker $USER

# Changer de groupe actif
newgrp docker
```

### MLST avec Docker

```bash
# Lancer commande de docker interactif dans le dossier de travail
docker run -it --rm -v "$PWD":/data staphb/mlst /bin/bash
mlst fichier.fasta > mlst.csv
# Taper exit pour sortir du docker !
```

### ABRicate avec Docker

```bash
# Lancer commande de docker interactif d'abricate dans le dossier de travail
docker run -it --rm -v $PWD:/data staphb/abricate /bin/bash
abricate fichier.fasta --db resfinder
abricate fichier.fasta --db plasmidfinder
abricate fichier.fasta --db vfdb
# Taper exit pour sortir du docker !
```

### QUAST avec Docker

```bash
# Lancer commande de docker interactif de quast dans le dossier de travail
docker run -it --rm -v "$PWD":/data staphb/quast /bin/bash

# Lancer l'analyse de quast
quast S1_assembled/scaffolds.fasta -1 trim_S1_R1.fastq.gz -2 trim_S1_R2.fastq.gz -o quast_resultats -r S1_ref.fasta

# Taper exit pour sortir du docker !
```

### BUSCO avec Docker

```bash
# Lancer commande de docker interactif de busco dans le dossier de travail
docker run -it --rm -v $PWD:/data ezlabgva/busco:v5.8.2_cv1 /bin/bash

# Se déplacer dans le répertoire de travail (ce n'est pas automatique pour busco !)
cd ../data

# Lancer l'analyse de busco
busco -i fichier.fasta -o busco -l bacteria_odb10 -m genome
# Taper exit pour sortir du docker !
```

### SeqSero2 avec Docker

```bash
# Lancer commande de docker interactif de SeqSero2 dans le dossier de travail
docker run -it --rm -v "$PWD":/data staphb/seqsero2 /bin/bash

# Lancer l'analyse de SeqSero
SeqSero2_package.py -i salmonella.fasta -m k -t 4

# Taper exit pour sortir du docker !
```

---

*Formation Africa CDC - Module bactériologie niveau débutant - Juillet 2025*
