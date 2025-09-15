## ExMODE -- Extremophile Multi-Omics DatabasE
Data analysis process：
### pipeline.sh 
#### 1. Quality control (Fastp)
#### 2. Taxonomic classification (Metaphlan4)
#### 3. Metagenomic assembly (MEGAHIT)
#### 4. Gene prediction (Metagenemark)
#### 5. Gene sets (MMseq2)
#### 6. Gene functional annotation (eggNOG-mapper)
#### 7. Structure sets (ESMfold)
#### 8. Metagenomic binning (MetaWRAP)
#### 9. Genomesets (dRep)
#### 10. Biosynthetic gene cluster analysis（Antismash BiG-SCAPE）
### metadata.xlsx
所有样本元信息表，包括BioSample、Taxonomy name、Isolation Source、Sample Title、Sample description、Project name、Project description等信息。
All sample metadata sheets, including BioSample, Taxonomy Name, Isolation Source, Sample Title, Sample Description, Project Name, Project Description, and other relevant information.
### sample.txt
All sample list.
### biome_list
List of samples from five habitats.
