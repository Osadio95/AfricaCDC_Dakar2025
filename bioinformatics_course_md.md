# Formation Bioinformatique – Kraken & Bracken

## Table des matières

1. [Identification des espèces](#identification-des-espèces)
   - [1.1 Mash](#11-mash--estimation-rapide-des-similarités-génomiques)
   - [1.2 Kraken](#12-kraken2--classification-taxonomique)
   - [1.3 Bracken](#13-bracken--estimation-de-labondance)
   - [1.4 Busco](#14-busco--évaluation-de-la-complétude-dun-assemblage)
2. [Génotypage bacterien](#2-génotypage-des-souches-bactériennes)
   - [2.1 Serotypage](#21-sérotypage--classification-antigénique-des-souches)
   - [2.2 MLST](#22-mlst--multilocus-sequence-typing)
   - [2.3 CgMLST](#23-cgmlst--typage-par-génome-cœur)

---

## 1. Identification des espèces

L'identification des espèces est une étape fondamentale en génomique, notamment dans le cadre de la surveillance microbiologique, des investigations d'épidémies ou de la recherche clinique. En effet, déterminer à quelle espèce appartient un génome séquencé permet de :

- ✅ **Connaître le pathogène impliqué** : chaque espèce bactérienne possède des caractéristiques biologiques, écologiques et pathogéniques spécifiques. Savoir de quelle espèce il s'agit oriente directement les décisions cliniques et de santé publique.

- 🧬 **Mieux interpréter les résultats génomiques** : la recherche de gènes de résistance, de virulence ou de plasmides n'a de sens que si l'on connaît le contexte biologique de l'organisme étudié. Par exemple, un gène de résistance trouvé chez Salmonella n'a pas les mêmes implications que chez E. coli.

- 🧭 **Comparer les génomes à des références appropriées** : pour des analyses comme le typage (MLST, cgMLST), le calcul de l'ANI, ou la recherche de SNPs, il est indispensable d'avoir une référence de la même espèce.

- 🧪 **Éviter les erreurs d'interprétation dues à des contaminations** : l'identification rapide permet aussi de détecter la présence éventuelle de séquences exogènes, comme des contaminations croisées ou des co-infections.

- 🌍 **Contribuer à la veille sanitaire** : en identifiant les espèces en circulation, on peut suivre l'évolution des pathogènes dans une population ou un environnement donné (ex. : hôpital, élevage, alimentation, etc.).

En résumé, l'identification des espèces est un prérequis essentiel à toute analyse bioinformatique fiable. C'est une étape de tri et de validation qui permet d'assurer la pertinence biologique des résultats génomiques produits. Ici nous allons travailler avec l'outil **MASH**.

---

## 1.1 Mash : Estimation rapide des similarités génomiques

**Mash** (pour "MinHash") est un outil bioinformatique qui permet de comparer rapidement de grands ensembles de séquences, comme des génomes, sans avoir besoin d'un alignement classique. Il est particulièrement adapté à l'identification des espèces dans un contexte de séquençage à haut débit.

### 🔬 Principe de fonctionnement

- Mash réduit les séquences à un "sketch" composé d'un sous-ensemble représentatif de k-mers.
- Il mesure la similarité entre deux sketchs en calculant la distance de Jaccard approximée.
- Les distances sont utilisées pour estimer la proximité génomique entre un échantillon et une base de données de référence.

### ⚙️ Exemple de commande

```bash
# Créer un sketch à partir de vos contigs
mash sketch -o echantillon mash_contigs.fasta

# Comparer à une base de référence
mash dist refseq_sketch.msh echantillon.msh > distances.txt
```

### 🔗 Lien vers le projet GitHub

[https://github.com/marbl/Mash](https://github.com/marbl/Mash)

### ✅ Avantages

- Très rapide, même sur des milliers de génomes
- Indépendant du type de séquençage (long ou court)
- Peut être utilisé pour identifier, regrouper ou comparer des échantillons

---

## 1.2 Kraken2 : Classification taxonomique

**Kraken2** est un outil bioinformatique de classification taxonomique rapide et efficace, utilisé pour identifier les espèces présentes dans un échantillon métagénomique à partir de lectures (reads) de séquençage.

### 🔍 Principe de fonctionnement

- **Basé sur les k-mers** : Kraken découpe chaque séquence en petits segments de longueur fixe appelés *k-mers*.
- **Base de données pré-indexée** : Chaque k-mer est associé à un taxon (espèce, genre, etc.) dans une base de données construite à l'avance à partir de génomes de référence.
- **Attribution par LCA (Lowest Common Ancestor)** : Si un k-mer est partagé par plusieurs espèces, Kraken attribue le k-mer au plus petit ancêtre commun.
- **Vote majoritaire** : Chaque read est classé selon le taxon qui reçoit le plus de votes parmi les k-mers qui le composent.

### ⚙️ Exemple de commande Kraken2

```bash
kraken2 --db /chemin/vers/base_de_donnees \
  --threads 4 --use-names \
  --report rapport.txt --output resultats.txt \
  echantillon.fastq
```

### ✅ Avantages

- Très rapide, même sur de gros jeux de données
- Utilisation efficace de la mémoire
- Résultats hiérarchiques exploitables (espèce, genre, etc.)
- Peut être utilisé sur des lectures Illumina ou Nanopore et même sur des fichiers d'assemblage

---

## 1.3 Bracken : Estimation de l'abondance

**Bracken** (Bayesian Reestimation of Abundance after Classification with Kraken) est un outil complémentaire à Kraken2. Il permet de recalculer les abondances taxonomiques en corrigeant les biais introduits par les classifications initiales de Kraken.

### 🔍 Pourquoi utiliser Bracken ?

- Kraken2 classe les lectures individuellement, ce qui peut conduire à une **sous-estimation des espèces spécifiques**.
- Bracken utilise un modèle statistique bayésien pour **réévaluer l'abondance réelle des taxons**.

### ⚙️ Fonctionnement

- Utilise le fichier `--report` de Kraken2
- Redistribue les lectures classées à des taxons plus précis
- Génère un tableau d'abondance par taxon

### 🧪 Exemple de commande Bracken

```bash
bracken -d /chemin/vers/base_de_donnees \
  -i rapport.txt \
  -o abondance_species.txt \
  -r 150 -l S
```

**Paramètres clés :**
- `-r` : longueur des reads (ex: 150 pour Illumina)
- `-l` : niveau taxonomique (`S`=Species, `G`=Genus)

### ✅ Résultats

- Meilleure estimation des abondances
- Sortie tabulaire : nom taxon, rang, nombre de lectures, pourcentage

---

## 1.4 BUSCO : Évaluation de la complétude d'un assemblage

**BUSCO** (Benchmarking Universal Single-Copy Orthologs) est un outil bioinformatique qui permet d'évaluer la qualité et la complétude d'un génome ou d'un assemblage en recherchant des gènes orthologues universels et présents en copie unique dans un groupe taxonomique donné.

### 🔬 Principe de fonctionnement

- BUSCO utilise des bases de données de gènes orthologues spécifiques à un clade (bactéries, protistes, etc.).
- Il recherche la présence de ces gènes dans un génome ou un assemblage pour évaluer la complétude.
- Le résultat inclut le pourcentage de gènes : **complets**, **dupliqués**, **fragmentés** ou **absents**.

### ⚙️ Exemple de commande

```bash
busco -i mon_assemblage.fasta \
  -o resultat_busco \
  -l bacteria_odb10 \
  -m genome
```

### 🔗 Lien vers le projet GitLab

[https://gitlab.com/ezlab/busco](https://gitlab.com/ezlab/busco)

### ✅ Avantages

- Évaluation rapide et standardisée de la qualité d'un assemblage
- Permet de comparer différents assemblages entre eux
- Large choix de bases de données adaptées à différents groupes taxonomiques

---

## 2. Génotypage des souches bactériennes

Le **génotypage bactérien** regroupe l'ensemble des méthodes permettant de **caractériser finement une souche** à partir de son matériel génétique. Contrairement à l'identification de l'espèce, le génotypage vise à différencier les souches **au sein d'une même espèce**, ce qui est essentiel pour :

- 🔍 Suivre la propagation d'une souche lors d'une épidémie,
- 🧭 Établir des liens épidémiologiques entre différents cas,
- 🧪 Identifier des clones résistants ou virulents,
- 🧬 Mieux comprendre la diversité génétique des populations bactériennes.

### 🧰 Trois approches majeures

#### 1. Sérotypage in silico

Le sérotypage permet de prédire le **sérotype** d'une souche (par exemple pour *Salmonella* ou *E. coli*) à partir de la présence de gènes spécifiques (O, H, K). Cette approche est souvent utilisée pour remplacer les méthodes sérologiques traditionnelles.

#### 2. MLST (Multilocus Sequence Typing)

Le MLST repose sur l'analyse de **7 gènes de ménage (housekeeping genes)**. Chaque allèle est numéroté et une combinaison d'allèles définit un **ST (sequence type)**. C'est une méthode largement utilisée pour la surveillance internationale, car elle est standardisée.

#### 3. cgMLST (core genome MLST)

Le cgMLST est une extension du MLST classique, utilisant **des centaines voire des milliers de loci du génome core**. Il offre une **résolution beaucoup plus fine**, utile pour les investigations d'épidémies hospitalières ou alimentaires.

### 🎯 En résumé

Le génotypage génomique est un outil de **traçabilité**, de **surveillance** et de **compréhension évolutive** des pathogènes. Chaque méthode a son niveau de résolution et ses cas d'usage. Dans cette présentation, nous allons explorer successivement le sérotypage, le MLST et le cgMLST.

---

## 2.1 Sérotypage : Classification antigénique des souches

### 🧬 Définition et importance

Le **sérotypage** est une méthode de classification des bactéries basée sur les **antigènes de surface**, tels que :

- 🔸 **Antigènes O** : présents sur la paroi cellulaire (lipopolysaccharides), notamment chez les *Enterobacteriaceae*.
- 🔸 **Antigènes H** : flagellaires, liés à la mobilité.
- 🔸 **Antigènes K** : capsulaires, comme chez *E. coli* ou *Salmonella*.

### 🎯 Applications clés

- 🧪 **Épidémiologie** : Suivi des épidémies (ex. *Salmonella enterica* sérovar Typhi).
- 🏥 **Diagnostic médical** : Identification de pathogènes (ex. *E. coli* O157:H7).
- 💉 **Vaccinologie** : Conception de vaccins (ex. contre *Streptococcus pneumoniae*).

### 🧰 Outils bioinformatiques pour le sérotypage

Il existe un certain nombre d'outils bioinformatiques disponibles pour le sérotypage des souches bactériennes. La plupart sont spécifiques à une espèce :

- 🔬 **Pour Salmonella spp.** : [SeqSero](https://github.com/denglab/SeqSero2), [sistr_cmd](https://github.com/phac-nml/sistr_cmd), [SeroTools](https://github.com/SeroTools/SeroTools)
- 🔬 **Pour Neisseria meningitidis** : [meningotype](https://github.com/INNUENDOCON/meningotype)
- 🔬 **Pour Escherichia coli** : [ECTyper](https://github.com/phac-nml/ECTyper), [ecoli_serotyper](https://github.com/INNUENDOCON/ecoli_serotyper)
- 🔬 **Pour STEC** : [STECFinder](https://github.com/ssi-dk/STECFinder)
- 🔬 **Pour Vibrio parahaemolyticus** : [VPsero](https://github.com/lanl/VPsero)
- 🔬 **Pour Pseudomonas aeruginosa** : [PAst](https://github.com/Pseudomonas-Serotyping/PAst)
- 🔬 **Pour Listeria monocytogenes** : [LisSero](https://github.com/Listeria-genomics/LisSero)
- 🔬 **Pour Shigella** : [ShigEiFinder](https://github.com/phiweger/ShigEiFinder), [ShigaPass](https://github.com/phe-bioinformatics/ShigaPass)
- 🔬 **Pour Streptococcus pneumoniae** : [seroBA](https://github.com/sanger-pathogens/seroBA), [PneumoCaT](https://github.com/phe-bioinformatics/PneumoCaT), [SeroCall](https://github.com/CDCgov/SeroCall), [Serotyping](https://github.com/MDU-PHL/Serotyping), [seqSerotyper](https://github.com/shakya-lab/seqSerotyper)
- 🔬 **Pour Streptococcus suis** : [SsuisSerotyping_pipeline](https://github.com/cptobie/SsuisSerotyping_pipeline)

### 📘 Outils utilisés dans ce cours

Dans ce cours, nous allons travailler avec les outils suivants :

- 🔧 **SeqSero2** pour le sérotypage de *Salmonella*
- 🔧 **ECTyper** pour le sérotypage de *E. coli*

### ⚙️ Exemples de commande

**▶️ SeqSero2 avec un fichier d'assemblage FASTA :**

```bash
SeqSero2_package.py -m fasta -i echantillon.fasta -o resultat_seqsero
```

**▶️ ECTyper avec un fichier d'assemblage FASTA :**

```bash
ectyper -i echantillon.fasta -o resultat_ectyper
```

---

## 2.2 MLST : Multilocus Sequence Typing

### 🧬 Principe général

Le **typage par séquence multi-focus (MLST)** consiste à séquencer de manière systématique **plusieurs loci conservés** du génome bactérien (généralement sept gènes de ménage). Les séquences alléliques sont comparées à une base de données de référence pour attribuer un **type de séquence (ST)** à l'isolat.

Chaque combinaison unique d'allèles constitue un profil qui est ensuite utilisé pour :

- 🧭 Étudier les lignées bactériennes
- 🧬 Suivre l'évolution clonale
- 🌍 Comparer des souches à l'échelle mondiale

Les profils alléliques sont comparés à ceux référencés dans les bases de données hébergées sur le site [pubmlst.org](https://pubmlst.org).

### 📘 Spécificité par espèce

Le nombre et le type de gènes utilisés pour le MLST varient selon les espèces. Voici deux exemples :

**🔬 Staphylococcus aureus** (7 gènes de ménage) :

- `arcC` – carbamate kinase
- `aroE` – shikimate déshydrogénase
- `glpF` – glycérol kinase
- `gmk` – guanylate kinase
- `pta` – phosphate acétyltransférase
- `tpi` – triosephosphate isomérase
- `yqiL` – acétyl coenzyme A acétyltransférase

**🔬 Vibrio vulnificus** (10 gènes de ménage) :

- `glp` – glucose-6-phosphate isomérase
- `gyrB` – ADN gyrase, sous-unité B
- `mdh` – malate-lactate déshydrogénase
- `metG` – méthionyl-ARNt synthétase
- `purM` – phosphoribosylaminoimidazole synthétase
- `dtdS` – thréonine déshydrogénase
- `lysA` – diaminopimélate décarboxylase
- `pntA` – transhydrogénase alpha sous-unité
- `pyrC` – dihydroorotase
- `tnaA` – tryptophanase

### 🧰 Outil utilisé dans ce cours

Nous allons utiliser l'outil `mlst` développé par **Torsten Seemann**, disponible sur GitHub :

[https://github.com/tseemann/mlst](https://github.com/tseemann/mlst)

### ⚙️ Exemple de commande

Pour lancer une analyse MLST sur un fichier d'assemblage FASTA :

```bash
mlst echantillon.fasta
```

Le programme détecte automatiquement l'espèce et ensuite détermine la ST correspondante.

---

## 2.3 cgMLST : Typage par génome cœur

### 🧬 Qu'est-ce que le cgMLST ?

Le **cgMLST**, ou *core genome Multilocus Sequence Typing*, est une méthode de typage génomique qui étend le MLST classique. Alors que le MLST traditionnel analyse seulement 7 gènes dits "de ménage", le cgMLST prend en compte **des centaines à des milliers de gènes conservés** dans le génome (appelés *core genes*).

Chaque gène est comparé à une base de données pour détecter sa version (appelée **allèle**), et le profil allélique complet est ensuite utilisé pour comparer les souches entre elles.

👉 Cette approche permet une résolution beaucoup plus fine que le MLST traditionnel, ce qui la rend très utile pour :

- 🦠 Étudier des **clusters épidémiques** de manière très précise
- 🧬 **Comparer des souches clonales** dans un même hôpital ou une même région
- 🧭 Réaliser une **surveillance génomique à haute résolution**

### 🛠️ Outil : chewBBACA

**chewBBACA** (prononcé "Chewbacca") est un outil développé pour faire du cgMLST à grande échelle. Il permet :

- 📦 De créer ou utiliser un schéma cgMLST existant (par espèce)
- 🔍 D'assigner un allèle à chaque gène d'un assemblage
- 📊 De générer un fichier de profils alléliques prêt pour la visualisation

Chaque échantillon obtient un "code" basé sur des centaines de loci, ce qui permet une comparaison fine entre souches très proches.

### 🧾 Exemple de commande chewBBACA

```bash
chewBBACA.py AlleleCall \
  -i genomes_folder/ \
  -g schema_cgmlst/ \
  -o resultats_cgmlst/ \
  --ptf bacteria.ptf 
```

*NB: cette partie sera plus développée en exercices pratiques.*

### 📈 Visualisation avec GrapeTree

Les résultats du cgMLST peuvent être visualisés sous forme d'arbre avec **GrapeTree**, un outil interactif développé pour visualiser les distances génétiques entre souches.

Chaque nœud représente une souche, et les branches reflètent les différences d'allèles entre les profils. Cela permet de :

- 👀 Identifier visuellement les clusters de souches proches
- 🧪 Étudier la dynamique d'une épidémie
- 📌 Intégrer des métadonnées (origine, date, statut, etc.) pour l'analyse

Le fichier des profils alléliques produit par chewBBACA (ex. `allele_matrix.tsv`) peut être directement importé dans GrapeTree (en ligne ou via EnteroBase).

### 🔗 Liens utiles

- 🧰 [chewBBACA (GitHub)](https://github.com/B-UMMI/chewBBACA)
- 🌳 [GrapeTree Viewer](https://github.com/achtman-lab/GrapeTree)