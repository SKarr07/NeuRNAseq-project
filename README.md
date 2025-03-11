# Transcriptional heterogeneity of neutrophils in steady state: a meta-analysis (add Journal and authors)

## Description 
Works containing samples that followed all the inclusion criteria were selected to proceed with rawdata downloading, processing and analysis.

## Citing
Add APA format citing and DOI: URL.

## Pipeline

- **Step 1**. Download public data: - [`Download_Data.md`](https://github.com/SKarr07/NeuRNAseq-project/blob/main/Download_Data.md)

- **Step 2**. Alignment and gene counts generation. 

- **Step 3**. Import data from R (gene counts file) + Metadata. Create a count matrix with all the transcriptomic samples. - `Script load_data_inR.R`

- **Step 4**. Create dds file using DESeq2, run differential expression analysis, data normalization, batch effect correction and contrast conditions - `Script DEG_analysis.R`

- **Step 5**. Data visualization (PCA, heatmap)- `Script VisualizacionDatos.R`
  
    - Figures 1-6.
    - Supplementary Figures.
 
- **Step 6**.Análisis de Terminos funcionales (GO terms) - `Script GOterms_analysis.R`

## Contact

Corresponding mail.
Affiliation.
Social networks.


------------------

## RNA-seq analysis workflow

### 1.Realizar la descarga del rawdata (fastq).

a) Existe la opción de descarga de fastqs directamente desde el ENA browser usando wget.

b) O también más agilmente usando la opción del SRAtoolkit.

```
$ qlogin

$ module load sra/3.0.0
```

i. Se puede hacer de manera individual.

ii. O utilizando una lista con los SRR de las muestras.

```
## Para una lista con los SRR de interés.
$ prefetch --option-file ID_list_names.txt		
#y una vez terminado el prefetch,
$ fastq-dump --gzip --skip-technical --dumpbase --split-3 --clip --outdir /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/SRP *.sra

## Para una `lista` de SRR pertenecientes a distintos datasets, hicimos lo siguiente:
$nano list
#1.Preparar una lista con las carpetas SRP y sus muestras correspondientes SRR.sra (1 por línea).
#2.Un archivo .sge para automatizar el fastqdump y mediante qsub.
$ cat list | while read line; do c=$(echo $line|cut -d' ' -f1); m=$(echo $line|cut -d' ' -f2); fastq-dump --gzip --skip-technical --dumpbase --split-3 --clip --outdir /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/$c $m;done
Una vez descargados los SRR en sus carpetas correspondientes, proceder al FastQC.
```

### 2.FastQC-Trimming-FastQC.

Primeramente, organizar las muestras por tipo de librería: SE o PE para un análisis en conjunto.

#Hacer una lista con todas las muestras SE y otra con los PE.
```
$ls /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/SRP*/*fastq.gz >> newSEfastq.txt

$ls /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/SRP*/*fastq.gz >> newPEfastq.txt

$mkdir newSEapr2023 newPEapr2023

$cd newSEapr2023
$for i in $(cat /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/newSEfastq.txt); do echo $i; ln -s $i .; done

$cd newPEapr2023
$for i in $(cat /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/newPEfastq.txt); do echo $i; ln -s $i .; done
```
Una vez teniendo los symlinks podemos trabajar con ellos.

En este caso, hice un .sh con los 3 pasos de FastQC-Trimming-FastQc consecutivos.

#FastQCTrimmed_wflow.sh

#!/bin/bash

#USAGE= ./FastQCTrimmed_wflow.sh
sed -i 's/\r//' FastQCTrimmed_wflow.sh

#AUTHOR: Saul Karr (adaptado de Eve Coss)

###a partir de # /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/newSEapr2023

#PARTE 1.- FastQC y multiQC
```
fastqc ./*.fastq.gz -o ./newFQCout
multiqc ./newFQCout -o ./newFQCout
```

#PARTE 2.- Limpieza de adaptadores
#single-end
cd newSEapr2023
```
for i in *.fastq.gz;
do echo
trimmomatic SE -threads 8 -phred33 $i data_trimmed/"${i%.fastq}_trimmed.fastq.gz" ILLUMINACLIP:/mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:22
done
```
#PARTE 3.- FastQC y multiQC
cd newSEapr2023
```
fastqc ./data_trimmed/*.fastq.gz -o ./FastQC_trimmed
```
#Reporte en MultiQC
```
multiqc ./FastQC_trimmed -o ./FastQC_trimmed
```
#Finalmente crear un .sge para enviar por qsub.
SI se va a enviar por qsub, cuidar de pedir los nucleos adecuados:
```
#$ -pe openmp 8 #para trimmomatic SE -threads 8
```
#Repetir y adaptar para muestras PE.

### 3.Alignment with STAR y FastQC.

Una vez revisadas las muestras post-trimming con FastQC, proceder al alineamiento.

#para STAR, correr primero el índice, si no se ha hecho antes
```
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir  /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/STAR_index \
--genomeFastaFiles /mnt/Archives/genome/human/GRCh38/ensembl76/chromosomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /mnt/Archives/genome/human/GRCh38/ensembl76/GTF-file/Homo_sapiens.GRCh38.76.gtf \
--sjdbOverhang 99
```
#STAR para SE
```
$nano star.sh

index=/mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/STAR_index
FILES=/mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/newSEapr2023/data_trimmed/*.fastq.gz
for f in $FILES
do
    echo $f
    base=$(basename $f .fastq.gz)
    echo $base
    STAR --runThreadN 20 --genomeDir $index --readFilesIn $f --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts --readFilesCommand zcat --outFileNamePrefix /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/newSEapr2023/STAR_output/$base"_"
done
```
SI se va a enviar por qsub, hay que solicitar los nucleos adecuados:
#$ -pe openmp 20 #para STAR --runThreadN 20

#STAR para PE

Creé un .sh para el alineamiento con STAR y el posterior fastqc de los . bam generados.

#PARTE 4.- STAR paired end reads
```
index=/mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/STAR_index
FILES=/mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/newPEapr2023/data_trimmed/*1_trimmed.fastq.gz
for f in $FILES
do
    echo $f
    base=$(basename $f _1_trimmed.fastq.gz)
    echo $base
    STAR --runThreadN 20 --genomeDir $index --readFilesIn $f /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/newPEapr2023/data_trimmed/$base"_2_trimmed.fastq.gz" --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts --readFilesCommand zcat --outFileNamePrefix /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/newPEapr2023/STAR_output/$base"_"
done
```
#PARTE 5.- FastQC y multiQC
cd newPEapr2023
```
fastqc ./STAR_output/*.out.bam -o ./BamQC
```
#Reporte en MultiQC
```
multiqc ./BamQC -o ./BamQC
```
Se verificaron aquellos archivos en /BamQC correspondientes a muestras con conflictos post-trimming p.ej. contenido de adaptadores y secuencias sobrerrepresentadas.


### 4.Exportar data para R.

Junté todos los ReadsPerGene.out.tab de distintos destinos usando $ln -s a partir de una lista.

Después en R conectarse a la terminal.
```
setwd("/mnt/Citosina/amedina/skarr/neu/monorail/SRP114762/STAR_output")
files <- dir(pattern="ReadsPerGene.out.tab") # count files
counts <- c()
for( i in seq_along(files) ){
  x <- read.table(file=files[i], sep="\t", header=F, as.is=T)
  counts <- cbind(counts, x[,2])
}
head(counts)
tail(counts)
dim(counts)
class(counts)
class(x) # the original count tab file from star
counts <- as.data.frame(counts)

#Assigning the rows a name (the ID)
rownames(counts) <- x[,1]
colnames(counts)

#solo para eliminar alguna muestra, indicar el numero de columna al que pertenece
counts <- counts[,-c(3,5)]

#Assinging the columns a name
#set the column names
colnames(counts) <- sub("_ReadsPerGene.out.tab", "", files)
colnames(counts)

#o en su defecto
colnames(counts) <- sub(".fastq.gz_trimmed_ReadsPerGene.out.tab", "", files)
colnames (counts)
# saving the counts
save(counts, file = "/mnt/Citosina/amedina/skarr/neu/monorail/SRP*/SRP*-counts.RData")
```


```{r Load modules}
library(rlang)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(tidyverse)
library(sf)
library(AnnotationHub)
library(AnnotationDbi)
library(dbplyr)
library(vctrs)
library(VennDiagram)
library(pheatmap)
library(limma)
library(ggrepel)
library(EnhancedVolcano)
library(tidyr)
library(biomaRt)
library(org.Hs.eg.db)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(topGO)
library(dplyr)
library(ComplexHeatmap)
library(colorRamp2)
```

```{r 184 spl Matrix configuration}
setwd("C:/Users/SaulKarr/OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey/DBT_data/neubulkRNA")

#The original Metadata
OctMetadata184 <- read.csv("OctMetadata184.csv")

#nmfclustering <- read.delim("nmfconsensus.k.4.txt") #volver a correr en cybersort

#carga el archivo de conteos 184spl-counts.RData
load("C:/Users/SaulKarr/OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey/DBT_data/neubulkRNA/184spl-counts.RData")

Neu184 <-counts

# eliminar las primeras 4 filas
head(Neu184)
Neu184 <- Neu184[-c(1:4),]

#Agregar columna de clasificacion por edad
OctMetadata184$Age_group <- ifelse(OctMetadata184$Age <= 30, "Y", ifelse(OctMetadata184$Age > 50, "O", "A"))
head(OctMetadata184)

#Si hay NA dejarlo vacio
OctMetadata184$Age_group <- ifelse(is.na(OctMetadata184$Age_group), "", OctMetadata184$Age_group)

#Nombrar la muestra por sus variables
OctMetadata184$Names <- paste(OctMetadata184$P_Method, OctMetadata184$Sex, OctMetadata184$Age_group, 1:184, sep = "_")

#OctMetadata184 <- as.matrix(OctMetadata184)

rownames(OctMetadata184) <- OctMetadata184$Names

colnames(Neu184)<- OctMetadata184$Names

#Verificacion debe ser TRUE
all(rownames(OctMetadata184) == colnames(Neu184))

#Set as.factor
OctMetadata184$Beads      <- as.factor(OctMetadata184$Beads)
OctMetadata184$Density    <- as.factor(OctMetadata184$Density)
OctMetadata184$FACS       <- as.factor(OctMetadata184$FACS)
OctMetadata184$X2step     <- as.factor(OctMetadata184$X2step)
OctMetadata184$Sex        <- as.factor(OctMetadata184$Sex)
OctMetadata184$Type       <- as.factor(OctMetadata184$Type)
OctMetadata184$P_Method  <- as.factor(OctMetadata184$P_Method)
OctMetadata184$Batch      <- as.factor(OctMetadata184$Batch)
OctMetadata184$Age_group  <- as.factor(OctMetadata184$Age_group)

#Create a dds object
dds <- DESeqDataSetFromMatrix(countData = Neu184,
                              colData = OctMetadata184,
                              design = ~ Type + P_Method)

#test type(SE or PE) as batch in the formula
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,] #add the pval threshold somewhere
dds <- DESeq(dds)
dds

#guardar el dds original
save(dds, file = "dds.RData")
load("dds.RData")

# A) Plot PCA
vsd <- vst(dds, blind=FALSE) #normalizacion
write.table(assay(vsd),file ="Normalized_counts.tsv")

#Matriz de vsd
vsd_mad <- apply(assay(vsd), 1, mad)
```

```{r}
# BATCH - lo checamos y "no hay batch effect"
mat <- assay(vsd)
mm <- model.matrix(~P_Method + Type, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Type, design=mm)
assay(vsd) <- mat

pcaData <- plotPCA(vsd, intgroup=c("P_Method", "Type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=P_Method, shape=Type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
```

```{r}
#CIBERSORT: lets get the gene names on the rownames of normalized counts

normcounts <- assay(vsd)

listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

geneIDs <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'gene_biotype', 'entrezgene_id'),
              filters = "ensembl_gene_id",
              values = rownames(normcounts),
              mart = ensembl.con)

normcounts_df<- as.data.frame(normcounts)
normcounts_df<-rownames_to_column(normcounts_df, var="ensembl_gene_id")
newnormcounts <- left_join(normcounts_df,geneIDs, by="ensembl_gene_id")

save(newnormcounts, file = "newnormcountsGeneIDs.RData")
load("C:/Users/SaulKarr/OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey/DBT_data/neubulkRNA/newnormcountsGeneIDs.RData")

#drop na
lnc <- newnormcounts[(newnormcounts$gene_biotype == "lncRNA"), ]

NArnas <- newnormcounts[(newnormcounts$gene_biotype == "NA"), ]

nncounts <- newnormcounts[(newnormcounts$gene_biotype == "protein_coding"), ]
nncounts2 <- nncounts %>% drop_na(external_gene_name)
nncounts3 <- nncounts2[!duplicated(nncounts2$external_gene_name), ]
rownames(nncounts3) <- nncounts3$external_gene_name

#Remove the columns created to filter
neucounts <- nncounts3[,-c(1,186,187,188)]

dim(neucounts)
save(neucounts, file= "184neusplGeneIDs.RData")

load("184neusplGeneIDs.RData")

#import to the cibersort cell fraction module
write.csv(neucounts,file ="Matriz184splNcounts.csv")
```
```{r}

#install.packages('extrafont')
library(extrafont)
font_import()
fonts()
loadfonts(device="win")


cibersort184 <- read.csv("CIBERSORT_184results.csv")

Cell_matrix <- as.data.frame(cibersort184[,-(24:26)])

# --- Change barplot ------
# Convert to long data
cellfracs_tidy <- Cell_matrix %>% arrange(Neutrophils) %>% #order
  pivot_longer(!Mixture, names_to = "cell_types", values_to = "percent")
cellfracs_tidy$cell_types <- as.factor(cellfracs_tidy$cell_types)

# order Samples
names_order <- cellfracs_tidy %>% dplyr::filter(cell_types == "Neutrophils") %>%
  arrange(percent)
names_order <- names_order$Mixture
cellfracs_tidy$Mixture <- factor(cellfracs_tidy$Mixture,
                              levels = names_order)
# order cells
cellfracs_tidy$name <- factor(cellfracs_tidy$Mixture,
                              levels = c("B.cells.naive", "B.cells.memory", "Plasma.cells",                
"T.cells.CD8",                  "T.cells.CD4.naive",            "T.cells.CD4.memory.resting",  
"T.cells.CD4.memory.activated", "T.cells.follicular.helper",    "T.cells.regulatory..Tregs.",  
"T.cells.gamma.delta",          "NK.cells.resting"        ,     "NK.cells.activated",          
"Monocytes"   ,                 "Macrophages.M0"           ,    "Macrophages.M1",              
"Macrophages.M2",               "Dendritic.cells.resting"   ,   "Dendritic.cells.activated",   
"Mast.cells.resting",           "Mast.cells.activated"       ,  "Eosinophils",                 
"Neutrophils"))

# Add label
cellfracs_tidy <- cellfracs_tidy %>%
  mutate(Category =
  case_when(cell_types == "Neutrophils" & percent*100 >=75 ~ "More than 75 %",
    .default = "Less than 75 %",
  ))
cellfracs_tidy$Category <- as.factor(cellfracs_tidy$Category)
# Last sample
last_sample <- cellfracs_tidy %>% dplyr::filter(cell_types == "Neutrophils" & Category == "Less than 75 %") %>% tail(1)
last_sample <- as.character(last_sample$Mixture)

library(RColorBrewer)
display.brewer.all()

# plot

png(filename = "deconv184neus.png", 
    res = 300, width = 50, height = 25, units = "cm")

quanti <- ggplot(cellfracs_tidy, aes(x = Mixture, y = percent, fill = cell_types)) +
  geom_bar(position = "fill", stat = "identity") +
  #geom_abline(xintercept = "D_F__83") + 
scale_fill_manual(values = c('gold1', 'yellow', 'gold4', 'sandybrown', 'tomato3', 'wheat', 'sienna', 'palegreen', 'olivedrab3', 'olivedrab', 'dodgerblue', 'darkblue', 'purple1', 'purple4', 'deeppink4', 'firebrick1', 'hotpink', 'deeppink1', 'azure2', 'azure3', 'azure4', 'black')) + # manual
  xlab("Samples") +
  ylab("Cell Fractions") +
  theme_minimal() + 
  theme(axis.text=element_text(size=9, angle = 90, hjust = 1),
        text=element_text(family = "Arial"),
        axis.title.x = element_text(size = 14, face="bold"),
        axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(),
      axis.title.y = element_text(size = 14, face="bold"),  # Change y-axis title size
  axis.text.y = element_text(size = 12),     # Change y-axis label size
  legend.position="top",
  legend.title = element_text(color = "black", size = 10),
          legend.text = element_text(color = "black", size = 15))
  
quanti

dev.off()




```

```{r}
#Import cibersort cell fraction results

#152 samples filtering after running CIBERSORTX

#We made a subset of samples with >75% Neutrophil Identity
cibersortxResults <- read.csv("CIBERSORTx_Job11.csv")

Metadata152spl <- left_join(OctMetadata184,cibersortxResults, by="Names")

rownames(Metadata152spl) <- Metadata152spl$Names

#Delete Identity NAs
Metadata152 <- Metadata152spl %>% drop_na(Identity)

#select columns for the samples present in Metadata152
Neu152 <- Neu184[,rownames(Metadata152)]
dim(Neu152)

colnames(Neu152)
rownames(Metadata152)

#Verify this is TRUE
all(rownames(Metadata152) == colnames(Neu152))

#orden
Metadata152 <- Metadata152[,c(14,1:13,15)]

write.csv(Metadata152, file= "Metadata152.csv",quote = F, col.names = T, row.names = F)

save(Neu152, Metadata152, file= "Neutrophils152.RData")

```

```{r}

load("Neutrophils152.RData")

str(Metadata152)

#cambiar contenido de celdas
Metadata152$P_Method  <- as.character(Metadata152$P_Method)

Metadata152$P_Method <- ifelse(Metadata152$P_Method == 'DF', '2steps', Metadata152$P_Method)
Metadata152$P_Method <- ifelse(Metadata152$P_Method == 'DBNF', '2steps', Metadata152$P_Method)
Metadata152$P_Method <- ifelse(Metadata152$P_Method == 'DBP', '2steps', Metadata152$P_Method)
Metadata152$P_Method <- ifelse(Metadata152$P_Method == 'DBN', '2steps', Metadata152$P_Method)

Metadata152$P_Method  <- as.factor(Metadata152$P_Method)

dds152 <- DESeqDataSetFromMatrix(countData = Neu152,
                              colData = Metadata152,
                              design = ~ Type + P_Method)


#test type(SE or PE) as batch in the formula
keep152 <- rowSums(counts(dds152)) >= 10
dds152 <- dds152[keep152,] #add the pval threshold somewhere

counts <- assays(dds152)$counts
head(counts(dds152))

# ---Compute z-score transformation---
var_non_zero <- apply(counts, 1, var) !=0 #  filter out genes that have low variance across samples
filtered_counts <- counts[var_non_zero, ]
zscores <- t(scale(t(filtered_counts)))
dim(zscores) #[1] 18782   165
zscore_mat <- as.matrix(zscores)

# Save melted_norm_counts with z-scores
save(zscore_mat, file = "zscore_mat.RData")
load("zscore_mat.RData")

head(zscore_mat, 3)


dds152 <- DESeq(dds152)
dds152

#guardar el dds original
save(dds152, file = "dds152.RData")
load("dds152.RData")

# A) Plot PCA
vsd152 <- vst(dds152, blind=FALSE) #normalizacion

#library(limma)
#plotDensities(assay(vsd), legend=F)
vsd_mad152 <- apply(assay(vsd152), 1, mad)

write.table(assay(vsd152),file ="Normalized152_counts.tsv")

# BATCH - lo checamos y "no hay batch effect"
mat152 <- assay(vsd152)
mm152 <- model.matrix(~P_Method + Type, colData(vsd152))
mat152 <- limma::removeBatchEffect(mat152, batch=vsd152$Type, design=mm152)
assay(vsd152) <- mat152

#FigS1. PCA

png(filename = "FigS1_PCA.png", res = 300, height = 20, width = 20, units = "cm")

pcaData <- plotPCA(vsd152, intgroup=c("P_Method", "Type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=P_Method, shape=Type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

dev.off()


pdf("FigS1_PCA.pdf")

pcaData <- plotPCA(vsd152, intgroup=c("P_Method", "Type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=P_Method, shape=Type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

dev.off()



```


```{r}

#lets get the gene names on the rownames of normalized counts

normcounts152 <- assay(vsd152)

listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

geneIDs <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'gene_biotype', 'entrezgene_id'),
              filters = "ensembl_gene_id",
              values = rownames(normcounts152),
              mart = ensembl.con)

#salvar con mi vida!
#write.csv(geneIDs, file= "normcounts152_geneIDs.csv",quote = F, col.names = T, row.names = F)
geneIDs <-read.csv("normcounts152_geneIDs.csv")

normcounts152_df<- as.data.frame(normcounts152)
normcounts152_df<-rownames_to_column(normcounts152_df, var="ensembl_gene_id")
newnormcounts152 <- left_join(normcounts152_df,geneIDs, by="ensembl_gene_id")

save(newnormcounts152, file = "newnormcounts152GeneIDs.RData")
load("C:/Users/SaulKarr/OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey/DBT_data/neubulkRNA/newnormcounts152GeneIDs.RData")

zscore_mat

zscore_mat_df<- as.data.frame(zscore_mat)
zscore_mat_df<-rownames_to_column(zscore_mat_df, var="ensembl_gene_id")
newzscore_mat <- left_join(zscore_mat_df,geneIDs, by="ensembl_gene_id")

save(newzscore_mat, file = "newzscore152GeneIDs.RData")
load("C:/Users/SaulKarr/OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey/DBT_data/neubulkRNA/newzscore152GeneIDs.RData")


#look for lncs 20/06/24
lnccounts152 <- newnormcounts152[(newnormcounts152$gene_biotype == "lncRNA"), ]
lnccounts152_2 <- lnccounts152 %>% drop_na(external_gene_name)
lnccounts152_3 <- lnccounts152_2[!duplicated(lnccounts152_2$external_gene_name), ]
rownames(lnccounts152_3) <- lnccounts152_3$external_gene_name #5147 rows

#drop na
nncounts152 <- newnormcounts152[(newnormcounts152$gene_biotype == "protein_coding"), ]
nncounts152_2 <- nncounts152 %>% drop_na(external_gene_name)
nncounts152_3 <- nncounts152_2[!duplicated(nncounts152_2$external_gene_name), ]
rownames(nncounts152_3) <- nncounts152_3$external_gene_name #18908 rows

neucounts152 <- nncounts152_3[,-c(1,154,155,156)]

dim(neucounts152)
save(neucounts152, file= "152neusplGeneIDs.RData") #para correr el NMF
load("152neusplGeneIDs.RData")

write.csv(neucounts152,file ="Matriz152splNcounts.csv") #para CIBERSORTX
#write.csv(neucounts,file ="Matriz152splNcounts.csv") 

zcounts152 <- newzscore_mat[(newzscore_mat$gene_biotype == "protein_coding"), ]
zcounts152_2 <- zcounts152 %>% drop_na(external_gene_name)
zcounts152_3 <- zcounts152_2[!duplicated(zcounts152_2$external_gene_name), ]
rownames(zcounts152_3) <- zcounts152_3$external_gene_name #18908 rows

zneucounts152 <- zcounts152_3[,-c(1,154,155,156)]

dim(zneucounts152)
save(zneucounts152, file= "z152neusplGeneIDs.RData") #para correr el NMF
load("z152neusplGeneIDs.RData")

write.csv(zneucounts152,file ="zMatriz152splNcounts.csv") #para CIBERSORTX
#write.csv(neucounts,file ="Matriz152splNcounts.csv") 


#load dds y neucounts

# TOP 20 genes
select <- order(rowMeans(neucounts152),
                decreasing=TRUE)[1:20]

#selectgenes <- (counts(dds,normalized=TRUE))[c('ENSG00000126012', 'ENSG00000147050', 'ENSG00000229807','ENSG00000067048','ENSG00000183878'),] 
absolutedf <- as.data.frame(colData(dds152)[,c("Age_group","Sex","P_Method")])

matneucounts152 <- as.matrix(neucounts152)

col_fun = colorRamp2(c(5, 15, 25), c("blue", "white", "red"))

pheatmap(matneucounts152[select,], 
         cluster_rows=T, 
         show_rownames=T, 
         show_colnames = F,
         cluster_cols=T,
         col = col_fun,
         name = "Norm counts",
         #annotation_row= ifelse(rownames(assay(vsd)) %in% BeadsonlyDEGs, BeadsonlyDEGs$external_gene_name, "" ),
         annotation_col = absolutedf,
         row_names_side = "left",
         show_row_dend = F,
         show_column_dend = T)

pdf("hmap152.pdf")

hmap152 <- pheatmap(matneucounts152[select,], 
         cluster_rows=T, 
         show_rownames=T, 
         show_colnames = F,
         cluster_cols=T,
         col = col_fun,
         name = "Norm counts",
         #annotation_row= ifelse(rownames(assay(vsd)) %in% BeadsonlyDEGs, BeadsonlyDEGs$external_gene_name, "" ),
         annotation_col = absolutedf,
         row_names_side = "left",
         show_row_dend = F,
         show_column_dend = T)

hmap152

dev.off()


# TOP 20 genes z-scores
selectZ <- order(rowMeans(zneucounts152),
                decreasing=TRUE)[1:20]

#selectgenes <- (counts(dds,normalized=TRUE))[c('ENSG00000126012', 'ENSG00000147050', 'ENSG00000229807','ENSG00000067048','ENSG00000183878'),] 
absolutedfZ <- as.data.frame(colData(dds152)[,c("Age_group","Sex","P_Method")])

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))


png(filename = "zTop20genes.png", 
    res = 300, width = 50, height = 25, units = "cm")

zTop20genes <- pheatmap(zneucounts152[selectZ,], 
         cluster_rows=T, 
         show_rownames=T, 
         show_colnames = F,
         cluster_cols=T,
         col = col_fun,
         name = "Norm counts",
         #annotation_row= ifelse(rownames(assay(vsd)) %in% BeadsonlyDEGs, BeadsonlyDEGs$external_gene_name, "" ),
         annotation_col = absolutedfZ,
         row_names_side = "left",
         show_row_dend = F,
         show_column_dend = T)

zTop20genes

dev.off()

```
```{r DiffExp and EVolcano for Density only vs Beads only}
ddsDvB <- DESeqDataSetFromMatrix(countData = Neu152,
                              colData = Metadata152,
                              design = ~ P_Method)

ddsDvB

keepDvB <- rowSums(counts(ddsDvB)) >= 10

ddsDvB <- ddsDvB[keepDvB,]

ddsDvB <- DESeq(ddsDvB)
ddsDvB

save(ddsDvB, file= "ddsDvB.RData")

#Expresion diferencial

resultsNames(ddsDvB)
DvB_res <-results(ddsDvB, contrast=c("P_Method","D","BN"))
DvB_res

save(DvB_res, file = "DensityvsBeadsRes.RData")

DvB_res_df <- as.data.frame(DvB_res)

#en values usar la lista de Ensembl IDs que se quiere convertir
DvB_resGenes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'gene_biotype'),
              filters = "ensembl_gene_id",
              values = rownames(DvB_res_df),
              mart = ensembl.con)

DvB_res_df<-rownames_to_column(DvB_res_df, var="ensembl_gene_id")
DvB_res <- left_join(DvB_res_df,DvB_resGenes, by="ensembl_gene_id")

dim(DvB_res) #52595 genes

DvB_res$gene_biotype <- as.factor(DvB_res$gene_biotype)
levels(DvB_res$gene_biotype) #37 diferentes gene_biotypes

DvB_res2 <- DvB_res[(DvB_res$gene_biotype == "protein_coding"), ]
dim(DvB_res2) #21312 filas o genes codificantes
DvB_res3 <- DvB_res2 %>% drop_na(c(padj, gene_biotype))
dim(DvB_res3) #18298 genes con padj

#remover duplicates y guardar uniques
DvB_uniqgenes <- DvB_res3[!duplicated(DvB_res3$external_gene_name), ] #18212 unique genes
rownames(DvB_uniqgenes) <- DvB_uniqgenes$external_gene_name

save(DvB_uniqgenes, file = "DBvB_uniqprotcodgenes.RData")
load("DBvB_uniqprotcodgenes.RData")

#clasificar genes significativos por padj
DvB_uniqgenes$diffexpressed <- ""
DvB_uniqgenes$diffexpressed[DvB_uniqgenes$log2FoldChange >= 1 & DvB_uniqgenes$padj < 0.05] <- "UP"
DvB_uniqgenes$diffexpressed[DvB_uniqgenes$log2FoldChange <= -1 & DvB_uniqgenes$padj < 0.05] <- "DOWN"

#make a df with only top and down genes using drop na then filter by head and tails
DvB_onlyDEGs <- DvB_uniqgenes[!(DvB_uniqgenes$diffexpressed == ""), ] #4358 DEGs
DvBUpgenes <- DvB_onlyDEGs[(DvB_onlyDEGs$diffexpressed == "UP"), ] #3158
DvBDowngenes <- DvB_onlyDEGs[(DvB_onlyDEGs$diffexpressed == "DOWN"), ] #1200

# Identify the top 20 or desired hits, arrange() apparently goes from AtoZ
DvB10Upgenes <- tail(arrange(DvB_onlyDEGs,log2FoldChange),10)
DvB10Downgenes <- head(arrange(DvB_onlyDEGs,log2FoldChange),10)

# Add column label, containing the gene name for the top hits or nothing for all others
DvB_uniqgenes$genetags <- ifelse(DvB_uniqgenes$external_gene_name %in% DvB10Upgenes$external_gene_name,
                                     DvB_uniqgenes$external_gene_name, 
                                     ifelse(DvB_uniqgenes$external_gene_name %in% DvB10Downgenes$external_gene_name,
                                            DvB_uniqgenes$external_gene_name, ""))

## Enhanced volcano

keyvalsDvB <- ifelse(
    DvB_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(DvB_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))
  keyvalsDvB[is.na(keyvalsDvB)] <- 'gray'
  names(keyvalsDvB)[keyvalsDvB == 'red'] <- 'Highest'
  names(keyvalsDvB)[keyvalsDvB == 'gray'] <- 'Mid'
  names(keyvalsDvB)[keyvalsDvB == 'blue'] <- 'Lowest'

  keyvals.shapeDvB <- ifelse(
    DvB_uniqgenes$external_gene_name %in% DvB10Downgenes$external_gene_name, 19,
      ifelse(DvB_uniqgenes$external_gene_name %in% DvB10Upgenes$external_gene_name, 19,
        20))
  
  keyvals.shapeDvB[is.na(keyvals.shapeDvB)] <- 20
  names(keyvals.shapeDvB)[keyvals.shapeDvB == 20] <- 'Mid'
  names(keyvals.shapeDvB)[keyvals.shapeDvB == 19] <- 'Lowest'
  names(keyvals.shapeDvB)[keyvals.shapeDvB == 19] <- 'Highest'

keyvals.colourDvB <- ifelse(
    DvB_uniqgenes$pval < 0.05 & DvB_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(DvB_uniqgenes$pval < 0.05 & DvB_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))

keyvals.colourDvB[is.na(keyvals.colourDvB)] <- 'gray'
names(keyvals.colourDvB)[keyvals.colourDvB == 'blue'] <- 'Negative selection by beads'
names(keyvals.colourDvB)[keyvals.colourDvB == 'gray'] <- 'Mid'
names(keyvals.colourDvB)[keyvals.colourDvB == 'red'] <- 'Density separation by centrifugation'

#library(ggtext)

DvB1 <- EnhancedVolcano(DvB_uniqgenes,
    lab = DvB_uniqgenes$genetags,
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = DvB_uniqgenes$external_gene_name[which(names(keyvalsDvB) %in% c('Highest', 'Lowest'))],
    title = "Differential gene expression by purification method",
    subtitle = NULL,
    caption = NULL,
  pCutoff = 0.05,
    FCcutoff = 1,
  pointSize = c(ifelse((DvB_uniqgenes$padj<0.05 & 
                          (DvB_uniqgenes$log2FoldChange >= 1 | DvB_uniqgenes$log2FoldChange <= -1)),3, 0.5)),
    labSize = 5,
  #col=c('gray2', 'green2', 'blue2', 'red'),
  shapeCustom = keyvals.shapeDvB,
  colCustom = keyvals.colourDvB,
  #shape = c(1, 4, 23, 25),
    colAlpha = 0.5,
    shadeAlpha = 1/2,
    cutoffLineType = 'twodash',
    cutoffLineCol = 'black',
    cutoffLineWidth = 0.8,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
  legendLabels=c('Not sig.','Log (base 2) FC','p-value',
      'p-value & Log (base 2) FC'),
   legendPosition = 'top',
    legendLabSize = 14,
    legendIconSize = 2,
  drawConnectors = TRUE,
  arrowheads = FALSE,
    endsConnectors = "first",
   widthConnectors = 0.25,
      labFace = 'bold',
    boxedLabels = F,
 max.overlaps = 10,
 maxoverlapsConnectors = Inf
  )

DvB1 +
     ggplot2::coord_cartesian(xlim=c(-7.5,7.5)) +
    ggplot2::scale_x_continuous(
   breaks=seq(-5,5, 5))+
    theme(plot.subtitle = element_text(hjust = 0.5)+
            theme(legend.position="top", legend.box = "horizontal"))+
           guides(shape=FALSE)+
       scale_colour_manual(values = keyvals.colourDvB, breaks=c('Negative selection by beads', '', 'Density separation by centrifugation'))

pdf("EV_Pm86.pdf",width = 8, height=6)

DvB1 +
     ggplot2::coord_cartesian(xlim=c(-7.5,7.5)) +
    ggplot2::scale_x_continuous(
   breaks=seq(-5,5, 5))+
    theme(plot.subtitle = element_text(hjust = 0.5)+
            theme(legend.position="top", legend.box = "horizontal"))+
           guides(shape=FALSE)+
       scale_colour_manual(values = keyvals.colourDvB, breaks=c('Negative selection by beads', '', 'Density separation by centrifugation'))

dev.off()

```

```{r DiffExp and EVolcano for Sex}
ddsSex <- DESeqDataSetFromMatrix(countData = Neu152,
                              colData = Metadata152,
                              design = ~ Sex)

ddsSex

keepSex <- rowSums(counts(ddsSex)) >= 10

ddsSex <- ddsSex[keepSex,]

ddsSex <- DESeq(ddsSex)
ddsSex

save(ddsSex, file = "ddsSex.RData")

#Expresion diferencial

resultsNames(ddsSex)
SexFvM_res <-results(ddsSex, contrast=c("Sex","F","M"))
SexFvM_res

save(SexFvM_res, file = "SexFvMRes.RData")
load("SexFvMRes.RData")

SexFvM_res_df <- as.data.frame(SexFvM_res)

#en values usar la lista de Ensembl IDs que se quiere convertir
SexFvM_resGenes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'gene_biotype'),
              filters = "ensembl_gene_id",
              values = rownames(SexFvM_res_df),
              mart = ensembl.con)

SexFvM_res_df<-rownames_to_column(SexFvM_res_df, var="ensembl_gene_id")
SexFvM_res <- left_join(SexFvM_res_df,SexFvM_resGenes, by="ensembl_gene_id")

dim(SexFvM_res) #52595 genes

SexFvM_res$gene_biotype <- as.factor(SexFvM_res$gene_biotype)
levels(SexFvM_res$gene_biotype) #37 diferentes gene_biotypes

SexFvM_res2 <- SexFvM_res[(SexFvM_res$gene_biotype == "protein_coding"), ]
dim(SexFvM_res2) #21312 filas o genes codificantes
SexFvM_res3 <- SexFvM_res2 %>% drop_na(c(padj, gene_biotype))
dim(SexFvM_res3) #15802 genes con padj

#remover duplicates y guardar uniques
SexFvM_uniqgenes <- SexFvM_res3[!duplicated(SexFvM_res3$external_gene_name), ] #15757 unique genes
rownames(SexFvM_uniqgenes) <- SexFvM_uniqgenes$external_gene_name

save(SexFvM_uniqgenes, file = "SexFvM_uniqprotcodgenes.RData")
load("SexFvM_uniqprotcodgenes.RData")

#clasificar genes significativos por padj
SexFvM_uniqgenes$diffexpressed <- ""
SexFvM_uniqgenes$diffexpressed[SexFvM_uniqgenes$log2FoldChange >= 1 & SexFvM_uniqgenes$padj < 0.05] <- "UP"
SexFvM_uniqgenes$diffexpressed[SexFvM_uniqgenes$log2FoldChange <= -1 & SexFvM_uniqgenes$padj < 0.05] <- "DOWN"

#make a df with only top and down genes using drop na then filter by head and tails
SexFvM_onlyDEGs <- SexFvM_uniqgenes[!(SexFvM_uniqgenes$diffexpressed == ""), ] #5 DEGs

# Identify the top 20 or desired hits, arrange() apparently goes from AtoZ
SexFvM10Upgenes <- tail(arrange(SexFvM_onlyDEGs,log2FoldChange),2)
SexFvM10Downgenes <- head(arrange(SexFvM_onlyDEGs,log2FoldChange),3)

# Add column label, containing the gene name for the top hits or nothing for all others
SexFvM_uniqgenes$genetags <- ifelse(SexFvM_uniqgenes$external_gene_name %in% SexFvM10Upgenes$external_gene_name,
                                     SexFvM_uniqgenes$external_gene_name, 
                                     ifelse(SexFvM_uniqgenes$external_gene_name %in% SexFvM10Downgenes$external_gene_name,
                                            SexFvM_uniqgenes$external_gene_name, ""))

keyvalsS <- ifelse(
    SexFvM_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(SexFvM_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))
  keyvalsS[is.na(keyvalsS)] <- 'gray'
  names(keyvalsS)[keyvalsS == 'red'] <- 'Highest'
  names(keyvalsS)[keyvalsS == 'gray'] <- 'Mid'
  names(keyvalsS)[keyvalsS == 'blue'] <- 'Lowest'

keyvals.shapeS <- ifelse(
    SexFvM_uniqgenes$external_gene_name %in% SexFvM10Downgenes$external_gene_name, 19,
      ifelse(SexFvM_uniqgenes$external_gene_name %in% SexFvM10Upgenes$external_gene_name, 19,
        20))
  
  keyvals.shapeS[is.na(keyvals.shapeS)] <- 20
  names(keyvals.shapeS)[keyvals.shapeS == 20] <- 'Mid'
  names(keyvals.shapeS)[keyvals.shapeS == 19] <- 'Lowest'
  names(keyvals.shapeS)[keyvals.shapeS == 19] <- 'Highest'

keyvals.colourS <- ifelse(
    SexFvM_uniqgenes$pval < 0.05 & SexFvM_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(SexFvM_uniqgenes$pval < 0.05 & SexFvM_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))
keyvals.colourS[is.na(keyvals.colourS)] <- 'gray'
names(keyvals.colourS)[keyvals.colourS == 'red'] <- 'Female'
names(keyvals.colourS)[keyvals.colourS == 'gray'] <- 'Mid'
names(keyvals.colourS)[keyvals.colourS == 'blue'] <- 'Male'

SexFvM <- EnhancedVolcano(SexFvM_uniqgenes,
    lab = SexFvM_uniqgenes$genetags,
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = SexFvM_uniqgenes$external_gene_name[which(names(keyvalsS) %in% c('Highest', 'Lowest'))],
    title = "Differential expression by sex",
  subtitle = NULL,
  caption=NULL,
  pCutoff = 0.05,
    FCcutoff = 1,
  pointSize = c(ifelse((SexFvM_uniqgenes$padj<0.05 & 
                          (SexFvM_uniqgenes$log2FoldChange >= 1 | SexFvM_uniqgenes$log2FoldChange <= -1)),3, 0.5)),
    labSize = 5,
  #col=c('gray2', 'green2', 'blue2', 'red'),
  shapeCustom = keyvals.shapeS,
  colCustom = keyvals.colourS,
  #shape = c(1, 4, 23, 25),
    colAlpha = 0.5,
    shadeAlpha = 1/2,
    cutoffLineType = 'twodash',
    cutoffLineCol = 'black',
    cutoffLineWidth = 0.8,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
  legendLabels=c('Not sig.','Log (base 2) FC','p-value',
      'p-value & Log (base 2) FC'),
    legendPosition = 'top',
    legendLabSize = 14,
    legendIconSize = 2,
  drawConnectors = TRUE,
  arrowheads = FALSE,
    endsConnectors = "first",
   widthConnectors = 0.25,
  labFace = 'bold',
      boxedLabels = F,
  max.overlaps = 10,
  maxoverlapsConnectors = Inf
  )

SexFvM +
    ggplot2::coord_cartesian(xlim=c(-4, 4), ylim=c(0,8)) +
    ggplot2::scale_x_continuous(
      breaks=seq(-4,4, 2))+
    theme(plot.subtitle = element_text(hjust = 0.5)+
            theme(legend.position="top", legend.box = "horizontal"))+
           guides(shape=FALSE)+
       scale_colour_manual(values = keyvals.colourS, breaks=c('Male', '', 'Female'))


pdf("EV_Sex85.pdf", width = 8, height = 5)

SexFvM +
    ggplot2::coord_cartesian(xlim=c(-4, 4), ylim=c(0,8)) +
    ggplot2::scale_x_continuous(
      breaks=seq(-4,4, 2))+
    theme(plot.subtitle = element_text(hjust = 0.5)+
            theme(legend.position="top", legend.box = "horizontal"))+
           guides(shape=FALSE)+
       scale_colour_manual(values = keyvals.colourS, breaks=c('Male', '', 'Female'))


dev.off()


```

```{r DiffExp and EVolcano for Age}

ddsAge <- DESeqDataSetFromMatrix(countData = Neu152,
                              colData = Metadata152,
                              design = ~ Age_group)

ddsAge

keepAge <- rowSums(counts(ddsAge)) >= 10

ddsAge <- ddsAge[keepAge,]

ddsAge <- DESeq(ddsAge)
ddsAge

save(ddsAge, file = "ddsAge.RData")

#Expresion diferencial

resultsNames(ddsAge)

AgeOvY_res <-results(ddsAge, contrast=c("Age_group","O","Y"))
AgeOvY_res

AgeYvO_res <-results(ddsAge, contrast=c("Age_group","Y","O"))
AgeYvO_res

save(AgeOvY_res, file = "AgeOvYRes.RData")
load("AgeOvYRes.RData")

save(AgeYvO_res, file = "AgeYvORes.RData")
load("AgeYvORes.RData")

AgeOvY_res_df <- as.data.frame(AgeOvY_res)
AgeYvO_res_df <- as.data.frame(AgeYvO_res)

#en values usar la lista de Ensembl IDs que se quiere convertir
AgeOvY_resGenes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'gene_biotype'),
              filters = "ensembl_gene_id",
              values = rownames(AgeOvY_res_df),
              mart = ensembl.con)

AgeOvY_res_df<-rownames_to_column(AgeOvY_res_df, var="ensembl_gene_id")
AgeOvY_res <- left_join(AgeOvY_res_df,AgeOvY_resGenes, by="ensembl_gene_id")

dim(AgeOvY_res) #52595 genes

AgeOvY_res$gene_biotype <- as.factor(AgeOvY_res$gene_biotype)
levels(AgeOvY_res$gene_biotype) #37 diferentes gene_biotypes

AgeOvY_res2 <- AgeOvY_res[(AgeOvY_res$gene_biotype == "protein_coding"), ]
dim(AgeOvY_res2) #21312 filas o genes codificantes
AgeOvY_res3 <- AgeOvY_res2 %>% drop_na(c(padj, gene_biotype))
dim(AgeOvY_res3) #17380 genes con padj

#remover duplicates y guardar uniques
AgeOvY_uniqgenes <- AgeOvY_res3[!duplicated(AgeOvY_res3$external_gene_name), ] #17313 unique genes
rownames(AgeOvY_uniqgenes) <- AgeOvY_uniqgenes$external_gene_name

save(AgeOvY_uniqgenes, file = "AgeOvY_uniqprotcodgenes.RData")
load("AgeOvY_uniqprotcodgenes.RData")

#clasificar genes significativos por padj
AgeOvY_uniqgenes$diffexpressed <- ""
AgeOvY_uniqgenes$diffexpressed[AgeOvY_uniqgenes$log2FoldChange >= 1 & AgeOvY_uniqgenes$padj < 0.05] <- "UP"
AgeOvY_uniqgenes$diffexpressed[AgeOvY_uniqgenes$log2FoldChange <= -1 & AgeOvY_uniqgenes$padj < 0.05] <- "DOWN"

#make a df with only top and down genes using drop na then filter by head and tails
AgeOvY_onlyDEGs <- AgeOvY_uniqgenes[!(AgeOvY_uniqgenes$diffexpressed == ""), ] #3 DEGs

# Identify the top 20 or desired hits, arrange() apparently goes from AtoZ
AgeOvY10Upgenes <- tail(arrange(AgeOvY_onlyDEGs,log2FoldChange),2)
AgeOvY10Downgenes <- head(arrange(AgeOvY_onlyDEGs,log2FoldChange),1)

# Add column label, containing the gene name for the top hits or nothing for all others
AgeOvY_uniqgenes$genetags <- ifelse(AgeOvY_uniqgenes$external_gene_name %in% AgeOvY10Upgenes$external_gene_name,
                                     AgeOvY_uniqgenes$external_gene_name, 
                                     ifelse(AgeOvY_uniqgenes$external_gene_name %in% AgeOvY10Downgenes$external_gene_name,
                                            AgeOvY_uniqgenes$external_gene_name, ""))


keyvalsA <- ifelse(
    AgeOvY_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(AgeOvY_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))
  keyvalsA[is.na(keyvalsA)] <- 'gray'
  names(keyvalsA)[keyvalsA == 'red'] <- 'Highest'
  names(keyvalsA)[keyvalsA == 'gray'] <- 'Mid'
  names(keyvalsA)[keyvalsA == 'blue'] <- 'Lowest'

keyvals.shapeA <- ifelse(
    AgeOvY_uniqgenes$external_gene_name %in% AgeOvY10Downgenes$external_gene_name, 19,
      ifelse(AgeOvY_uniqgenes$external_gene_name %in% AgeOvY10Upgenes$external_gene_name, 19,
        20))
  
  keyvals.shapeA[is.na(keyvals.shapeA)] <- 20
  names(keyvals.shapeA)[keyvals.shapeA == 20] <- 'Mid'
  names(keyvals.shapeA)[keyvals.shapeA == 19] <- 'Lowest'
  names(keyvals.shapeA)[keyvals.shapeA == 19] <- 'Highest'

keyvals.colourA <- ifelse(
    AgeOvY_uniqgenes$pval < 0.05 & AgeOvY_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(AgeOvY_uniqgenes$pval < 0.05 & AgeOvY_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))
keyvals.colourA[is.na(keyvals.colourA)] <- 'gray'
names(keyvals.colourA)[keyvals.colourA == 'red'] <- 'Older than 50 years old'
names(keyvals.colourA)[keyvals.colourA == 'gray'] <- 'Mid'
names(keyvals.colourA)[keyvals.colourA == 'blue'] <- '30 years old and younger'

AgeOvY <- EnhancedVolcano(AgeOvY_uniqgenes,
    lab = AgeOvY_uniqgenes$genetags,
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = AgeOvY_uniqgenes$external_gene_name[which(names(keyvalsA) %in% c('Highest', 'Lowest'))],
    title = "Differential gene expression by age groups",
  subtitle = NULL,
  caption=NULL,
  pCutoff = 0.05,
    FCcutoff = 1,
  pointSize = c(ifelse((AgeOvY_uniqgenes$padj<0.05 & 
                          (AgeOvY_uniqgenes$log2FoldChange >= 1 | AgeOvY_uniqgenes$log2FoldChange <= -1)),3, 0.5)),
    labSize = 5,
  #col=c('gray2', 'green2', 'blue2', 'red'),
  shapeCustom = keyvals.shapeA,
  colCustom = keyvals.colourA,
  #shape = c(1, 4, 23, 25),
    colAlpha = 0.5,
    shadeAlpha = 1/2,
    cutoffLineType = 'twodash',
    cutoffLineCol = 'black',
    cutoffLineWidth = 0.8,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
  legendLabels=c('Not sig.','Log (base 2) FC','p-value',
      'p-value & Log (base 2) FC'),
    legendPosition = 'top',
    legendLabSize = 14,
    legendIconSize = 2.0,
  drawConnectors = TRUE,
  arrowheads = FALSE,
    endsConnectors = "first",
   widthConnectors = 0.25,
  labFace = 'bold',
      boxedLabels = F,
  max.overlaps = 10,
  maxoverlapsConnectors = Inf
  )

AgeOvY +
    ggplot2::coord_cartesian(xlim=c(-3.5, 3.5), ylim=c(0,7)) +
    ggplot2::scale_x_continuous(
      breaks=seq(-3,3, 1))+
      theme(plot.subtitle = element_text(hjust = 0.5)+
            theme(legend.position="top", legend.box = "horizontal"))+
           guides(shape=FALSE)+
       scale_colour_manual(values = keyvals.colourA, breaks=c('30 years old and younger', '', 'Older than 50 years old'))

pdf("EV_Ages877.pdf", width = 8, height = 7)

AgeOvY +
    ggplot2::coord_cartesian(xlim=c(-3.5, 3.5), ylim=c(0,7)) +
    ggplot2::scale_x_continuous(
      breaks=seq(-3,3, 1))+
      theme(plot.subtitle = element_text(hjust = 0.5)+
            theme(legend.position="top", legend.box = "horizontal"))+
           guides(shape=FALSE)+
       scale_colour_manual(values = keyvals.colourA, breaks=c('30 years old and younger', '', 'Older than 50 years old'))


dev.off()




```

```{r DE and EVolcano for clusters}

clustering3 <- read.delim("spl152_nmfconsensus.k.3.txt")

#cambiar nombre de columna
clustering3 <- clustering3 %>% rename("Names" = "Name")
clustering3 <- clustering3 %>% rename("Cluster" = "membership.ordered")

Neu3Clusters <- left_join(Metadata152, clustering3, by="Names")

Neu3Clusters$Cluster <- as.factor(Neu3Clusters$Cluster)

Neu3Clusters$Cluster1 <- ""
Neu3Clusters$Cluster1[Neu3Clusters$Cluster == 1] <- "1"
Neu3Clusters$Cluster1[!Neu3Clusters$Cluster == 1] <- "0"

Neu3Clusters$Cluster2 <- ""
Neu3Clusters$Cluster2[Neu3Clusters$Cluster == 2] <- "1"
Neu3Clusters$Cluster2[!Neu3Clusters$Cluster == 2] <- "0"

Neu3Clusters$Cluster3 <- ""
Neu3Clusters$Cluster3[Neu3Clusters$Cluster == 3] <- "1"
Neu3Clusters$Cluster3[!Neu3Clusters$Cluster == 3] <- "0"

Neu3Clusters$Cluster1 <- as.factor(Neu3Clusters$Cluster1)
Neu3Clusters$Cluster2 <- as.factor(Neu3Clusters$Cluster2)
Neu3Clusters$Cluster3 <- as.factor(Neu3Clusters$Cluster3)

ddsCluster1 <- DESeqDataSetFromMatrix(countData = Neu152,
                              colData = Neu3Clusters,
                              design = ~ Cluster1)
ddsCluster1
keepCluster1 <- rowSums(counts(ddsCluster1)) >= 10
ddsCluster1 <- ddsCluster1[keepCluster1,]
ddsCluster1 <- DESeq(ddsCluster1)

ddsCluster2 <- DESeqDataSetFromMatrix(countData = Neu152,
                              colData = Neu3Clusters,
                              design = ~ Cluster2)
ddsCluster2
keepCluster2 <- rowSums(counts(ddsCluster2)) >= 10
ddsCluster2 <- ddsCluster2[keepCluster2,]
ddsCluster2 <- DESeq(ddsCluster2)

ddsCluster3 <- DESeqDataSetFromMatrix(countData = Neu152,
                              colData = Neu3Clusters,
                              design = ~ Cluster3)
ddsCluster3
keepCluster3 <- rowSums(counts(ddsCluster3)) >= 10
ddsCluster3 <- ddsCluster3[keepCluster3,]
ddsCluster3 <- DESeq(ddsCluster3)


#Expresion diferencial

resultsNames(ddsCluster1)
resultsNames(ddsCluster2)
resultsNames(ddsCluster3)

Cluster1_res <-results(ddsCluster1, contrast=c("Cluster1","1","0"))
Cluster2_res <-results(ddsCluster2, contrast=c("Cluster2","1","0"))
Cluster3_res <-results(ddsCluster3, contrast=c("Cluster3","1","0"))

save(Cluster1_res, file = "Cluster1res.RData")
save(Cluster2_res, file = "Cluster2res.RData")
save(Cluster3_res, file = "Cluster3res.RData")
load("Cluster1res.RData")
load("Cluster2res.RData")
load("Cluster3res.RData")


Cluster1_res_df <- as.data.frame(Cluster1_res)
Cluster2_res_df <- as.data.frame(Cluster2_res)
Cluster3_res_df <- as.data.frame(Cluster3_res)

#en values usar la lista de Ensembl IDs que se quiere convertir
Cluster1_resGenes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'gene_biotype'),
              filters = "ensembl_gene_id",
              values = rownames(Cluster1_res_df),
              mart = ensembl.con)

Cluster1_res_df<-rownames_to_column(Cluster1_res_df, var="ensembl_gene_id")
Cluster1_res <- left_join(Cluster1_res_df,Cluster1_resGenes, by="ensembl_gene_id")

Cluster2_resGenes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'gene_biotype'),
              filters = "ensembl_gene_id",
              values = rownames(Cluster2_res_df),
              mart = ensembl.con)
Cluster2_res_df<-rownames_to_column(Cluster2_res_df, var="ensembl_gene_id")
Cluster2_res <- left_join(Cluster2_res_df,Cluster2_resGenes, by="ensembl_gene_id")

Cluster3_resGenes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'gene_biotype'),
              filters = "ensembl_gene_id",
              values = rownames(Cluster3_res_df),
              mart = ensembl.con)
Cluster3_res_df<-rownames_to_column(Cluster3_res_df, var="ensembl_gene_id")
Cluster3_res <- left_join(Cluster3_res_df,Cluster3_resGenes, by="ensembl_gene_id")

dim(Cluster1_res) #52998 genes
dim(Cluster2_res) #52998 genes
dim(Cluster3_res) #52998 genes

Cluster1_res$gene_biotype <- as.factor(Cluster1_res$gene_biotype)
Cluster1_res2 <- Cluster1_res[(Cluster1_res$gene_biotype == "protein_coding"), ]
dim(Cluster1_res2) #21411 filas o genes codificantes
Cluster1_res3 <- Cluster1_res2 %>% drop_na(c(padj, gene_biotype))
dim(Cluster1_res3) #18530 genes con padj
#remover duplicates y guardar uniques
Cluster1_uniqgenes <- Cluster1_res3[!duplicated(Cluster1_res3$external_gene_name), ] #18424 unique genes
rownames(Cluster1_uniqgenes) <- Cluster1_uniqgenes$external_gene_name

Cluster2_res$gene_biotype <- as.factor(Cluster2_res$gene_biotype)
Cluster2_res2 <- Cluster2_res[(Cluster2_res$gene_biotype == "protein_coding"), ]
dim(Cluster2_res2) #21411 filas o genes codificantes
Cluster2_res3 <- Cluster2_res2 %>% drop_na(c(padj, gene_biotype))
dim(Cluster2_res3) #18530 genes con padj
#remover duplicates y guardar uniques
Cluster2_uniqgenes <- Cluster2_res3[!duplicated(Cluster2_res3$external_gene_name), ] #18652 unique genes
rownames(Cluster2_uniqgenes) <- Cluster2_uniqgenes$external_gene_name

Cluster3_res$gene_biotype <- as.factor(Cluster3_res$gene_biotype)
Cluster3_res2 <- Cluster3_res[(Cluster3_res$gene_biotype == "protein_coding"), ]
dim(Cluster3_res2) #21411 filas o genes codificantes
Cluster3_res3 <- Cluster3_res2 %>% drop_na(c(padj, gene_biotype))
dim(Cluster3_res3) #18530 genes con padj
#remover duplicates y guardar uniques
Cluster3_uniqgenes <- Cluster3_res3[!duplicated(Cluster3_res3$external_gene_name), ] #18645 unique genes
rownames(Cluster3_uniqgenes) <- Cluster3_uniqgenes$external_gene_name

save(Cluster1_uniqgenes, file = "Cluster1_uniqprotcodgenes.RData")
save(Cluster2_uniqgenes, file = "Cluster2_uniqprotcodgenes.RData")
save(Cluster3_uniqgenes, file = "Cluster3_uniqprotcodgenes.RData")
load("Cluster1_uniqprotcodgenes.RData")
load("Cluster2_uniqprotcodgenes.RData")
load("Cluster3_uniqprotcodgenes.RData")

#clasificar genes significativos por padj
Cluster1_uniqgenes$diffexpressed <- ""
Cluster1_uniqgenes$diffexpressed[Cluster1_uniqgenes$log2FoldChange >= 1 & Cluster1_uniqgenes$padj < 0.05] <- "UP"
Cluster1_uniqgenes$diffexpressed[Cluster1_uniqgenes$log2FoldChange <= -1 & Cluster1_uniqgenes$padj < 0.05] <- "DOWN"
Cluster1_onlyDEGs <- Cluster1_uniqgenes[!(Cluster1_uniqgenes$diffexpressed == ""), ] #6225 DEGs

Cluster2_uniqgenes$diffexpressed <- ""
Cluster2_uniqgenes$diffexpressed[Cluster2_uniqgenes$log2FoldChange >= 1 & Cluster2_uniqgenes$padj < 0.05] <- "UP"
Cluster2_uniqgenes$diffexpressed[Cluster2_uniqgenes$log2FoldChange <= -1 & Cluster2_uniqgenes$padj < 0.05] <- "DOWN"
Cluster2_onlyDEGs <- Cluster2_uniqgenes[!(Cluster2_uniqgenes$diffexpressed == ""), ] #3656 DEGs

Cluster3_uniqgenes$diffexpressed <- ""
Cluster3_uniqgenes$diffexpressed[Cluster3_uniqgenes$log2FoldChange >= 1 & Cluster3_uniqgenes$padj < 0.05] <- "UP"
Cluster3_uniqgenes$diffexpressed[Cluster3_uniqgenes$log2FoldChange <= -1 & Cluster3_uniqgenes$padj < 0.05] <- "DOWN"
Cluster3_onlyDEGs <- Cluster3_uniqgenes[!(Cluster3_uniqgenes$diffexpressed == ""), ] #5465 DEGs

# UP-regulated genes
Cluster1_upGenes <- subset(Cluster1_onlyDEGs, log2FoldChange >= 1)
rownames(Cluster1_upGenes) <- Cluster1_upGenes$ensembl_gene_id
Cluster1_upGenes_names <- rownames(as.data.frame(Cluster1_upGenes))
length(Cluster1_upGenes_names) #877
write.table(Cluster1_upGenes,file ="./Cluster1_upGenes.tsv", quote=FALSE, sep="\t")

Cluster2_upGenes <- subset(Cluster2_onlyDEGs, log2FoldChange >= 1)
rownames(Cluster2_upGenes) <- Cluster2_upGenes$ensembl_gene_id
Cluster2_upGenes_names <- rownames(as.data.frame(Cluster2_upGenes))
length(Cluster2_upGenes_names) #1209
write.table(Cluster2_upGenes,file ="./Cluster2_upGenes.tsv", quote=FALSE, sep="\t")

Cluster3_upGenes <- subset(Cluster3_onlyDEGs, log2FoldChange >= 1)
rownames(Cluster3_upGenes) <- Cluster3_upGenes$ensembl_gene_id
Cluster3_upGenes_names <- rownames(as.data.frame(Cluster3_upGenes))
length(Cluster3_upGenes_names) #751
write.table(Cluster3_upGenes,file ="./Cluster3_upGenes.tsv", quote=FALSE, sep="\t")


Clusterbar <- as.data.frame(Neu3Clusters$Cluster)
colorCluster=c("blue3","blue4","magenta4")

Clusterbar %>% 
  ggplot(aes(x = `Neu3Clusters$Cluster`)) +
geom_bar(color="white", fill=colorCluster) +
  theme(legend.position="none") +
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_fill(vjust=2),
             colour="black")+
  labs(y = " ", x= "Clusters")

```

```{r EVolcano for clusters}

load("Cluster1_uniqprotcodgenes.RData")
load("Cluster2_uniqprotcodgenes.RData")
load("Cluster3_uniqprotcodgenes.RData")

#clasificar genes significativos por padj
Cluster1_uniqgenes$diffexpressed <- ""
Cluster1_uniqgenes$diffexpressed[Cluster1_uniqgenes$log2FoldChange >= 1 & Cluster1_uniqgenes$padj < 0.05] <- "UP"
Cluster1_uniqgenes$diffexpressed[Cluster1_uniqgenes$log2FoldChange <= -1 & Cluster1_uniqgenes$padj < 0.05] <- "DOWN"
Cluster1_onlyDEGs <- Cluster1_uniqgenes[!(Cluster1_uniqgenes$diffexpressed == ""), ] #6225 DEGs

Cluster2_uniqgenes$diffexpressed <- ""
Cluster2_uniqgenes$diffexpressed[Cluster2_uniqgenes$log2FoldChange >= 1 & Cluster2_uniqgenes$padj < 0.05] <- "UP"
Cluster2_uniqgenes$diffexpressed[Cluster2_uniqgenes$log2FoldChange <= -1 & Cluster2_uniqgenes$padj < 0.05] <- "DOWN"
Cluster2_onlyDEGs <- Cluster2_uniqgenes[!(Cluster2_uniqgenes$diffexpressed == ""), ] #3656 DEGs

Cluster3_uniqgenes$diffexpressed <- ""
Cluster3_uniqgenes$diffexpressed[Cluster3_uniqgenes$log2FoldChange >= 1 & Cluster3_uniqgenes$padj < 0.05] <- "UP"
Cluster3_uniqgenes$diffexpressed[Cluster3_uniqgenes$log2FoldChange <= -1 & Cluster3_uniqgenes$padj < 0.05] <- "DOWN"
Cluster3_onlyDEGs <- Cluster3_uniqgenes[!(Cluster3_uniqgenes$diffexpressed == ""), ] #5465 DEGs


# Identify the top 20 or desired hits, arrange() apparently goes from AtoZ
Cluster1_10Upgenes <- tail(arrange(Cluster1_onlyDEGs,log2FoldChange),10)
Cluster1_10Downgenes <- head(arrange(Cluster1_onlyDEGs,log2FoldChange),10)

# Add column label, containing the gene name for the top hits or nothing for all others
Cluster1_uniqgenes$genetags <- ifelse(Cluster1_uniqgenes$external_gene_name %in% Cluster1_10Upgenes$external_gene_name,
                                     Cluster1_uniqgenes$external_gene_name, 
                                     ifelse(Cluster1_uniqgenes$external_gene_name %in% Cluster1_10Downgenes$external_gene_name,
                                            Cluster1_uniqgenes$external_gene_name, ""))

keyvalsC1 <- ifelse(
    Cluster1_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(Cluster1_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))
  keyvalsC1[is.na(keyvalsC1)] <- 'gray'
  names(keyvalsC1)[keyvalsC1 == 'red'] <- 'Highest'
  names(keyvalsC1)[keyvalsC1 == 'gray'] <- 'Mid'
  names(keyvalsC1)[keyvalsC1 == 'blue'] <- 'Lowest'

keyvals.shapeC1 <- ifelse(
    Cluster1_uniqgenes$external_gene_name %in% Cluster1_10Downgenes$external_gene_name, 19,
      ifelse(Cluster1_uniqgenes$external_gene_name %in% Cluster1_10Upgenes$external_gene_name, 19,
        20))
  
  keyvals.shapeC1[is.na(keyvals.shapeC1)] <- 20
  names(keyvals.shapeC1)[keyvals.shapeC1 == 20] <- 'Mid'
  names(keyvals.shapeC1)[keyvals.shapeC1 == 19] <- 'Lowest'
  names(keyvals.shapeC1)[keyvals.shapeC1 == 19] <- 'Highest'

  
keyvals.colourC1 <- ifelse(
    Cluster1_uniqgenes$pval < 0.05 & Cluster1_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(Cluster1_uniqgenes$pval < 0.05 & Cluster1_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))

keyvals.colourC1[is.na(keyvals.colourC1)] <- 'gray'
names(keyvals.colourC1)[keyvals.colourC1 == 'red'] <- 'Upregulated'
names(keyvals.colourC1)[keyvals.colourC1 == 'gray'] <- 'Mid'
names(keyvals.colourC1)[keyvals.colourC1 == 'blue'] <- 'Downregulated'


C1EV <- EnhancedVolcano(Cluster1_uniqgenes,
    lab = Cluster1_uniqgenes$genetags,
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = Cluster1_uniqgenes$external_gene_name[which(names(keyvalsC1) %in% c('Highest', 'Lowest'))],
    title = "Differentially expressed genes in Cluster 1",
  subtitle = NULL,
  caption = NULL,
  pCutoff = 0.05,
    FCcutoff = 1,
  pointSize = c(ifelse((Cluster1_uniqgenes$padj<0.05 & 
                          (Cluster1_uniqgenes$log2FoldChange >= 1 | Cluster1_uniqgenes$log2FoldChange <= -1)),3, 0.5)),
    labSize = 5,
  #col=c('gray2', 'green2', 'blue2', 'red'),
  shapeCustom = keyvals.shapeC1,
  colCustom = keyvals.colourC1,
  #shape = c(1, 4, 23, 25),
    colAlpha = 0.5,
    shadeAlpha = 1/2,
    cutoffLineType = 'twodash',
    cutoffLineCol = 'black',
    cutoffLineWidth = 0.8,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
  legendLabels=c('Not sig.','Log (base 2) FC','p-value',
      'p-value & Log (base 2) FC'),
    legendPosition = 'top',
    legendLabSize = 14,
    legendIconSize = 2,
  drawConnectors = TRUE,
  arrowheads = FALSE,
    endsConnectors = "first",
   widthConnectors = 0.25,
  labFace = 'bold',
      boxedLabels = F,
  max.overlaps = 10,
  maxoverlapsConnectors = Inf
  )


pdf("EV_Cluster1.pdf", width = 8, height = 5)

C1EV +
    ggplot2::coord_cartesian(xlim=c(-10, 10), ylim=c(0,80)) +
    ggplot2::scale_x_continuous(
      breaks=seq(-10,10, 5))+
    theme(plot.subtitle = element_text(hjust = 0.5)+
            theme(legend.position="top", legend.box = "horizontal"))+
           guides(shape=FALSE)+
       scale_colour_manual(values = keyvals.colourC1, breaks=c('Downregulated', '', 'Upregulated'))

dev.off()



# Identify the top 20 or desired hits, arrange() apparently goes from AtoZ
Cluster2_10Upgenes <- tail(arrange(Cluster2_onlyDEGs,log2FoldChange),10)
Cluster2_10Downgenes <- head(arrange(Cluster2_onlyDEGs,log2FoldChange),10)

# Add column label, containing the gene name for the top hits or nothing for all others
Cluster2_uniqgenes$genetags <- ifelse(Cluster2_uniqgenes$external_gene_name %in% Cluster2_10Upgenes$external_gene_name,
                                     Cluster2_uniqgenes$external_gene_name, 
                                     ifelse(Cluster2_uniqgenes$external_gene_name %in% Cluster2_10Downgenes$external_gene_name,
                                            Cluster2_uniqgenes$external_gene_name, ""))

keyvalsC2 <- ifelse(
    Cluster2_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(Cluster2_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))
  keyvalsC2[is.na(keyvalsC2)] <- 'gray'
  names(keyvalsC2)[keyvalsC2 == 'red'] <- 'Highest'
  names(keyvalsC2)[keyvalsC2 == 'gray'] <- 'Mid'
  names(keyvalsC2)[keyvalsC2 == 'blue'] <- 'Lowest'

keyvals.shapeC2 <- ifelse(
    Cluster2_uniqgenes$external_gene_name %in% Cluster2_10Downgenes$external_gene_name, 19,
      ifelse(Cluster2_uniqgenes$external_gene_name %in% Cluster2_10Upgenes$external_gene_name, 19,
        20))
  
  keyvals.shapeC2[is.na(keyvals.shapeC2)] <- 20
  names(keyvals.shapeC2)[keyvals.shapeC2 == 20] <- 'Mid'
  names(keyvals.shapeC2)[keyvals.shapeC2 == 19] <- 'Lowest'
  names(keyvals.shapeC2)[keyvals.shapeC2 == 19] <- 'Highest'

  
keyvals.colourC2 <- ifelse(
    Cluster2_uniqgenes$pval < 0.05 & Cluster2_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(Cluster2_uniqgenes$pval < 0.05 & Cluster2_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))

keyvals.colourC2[is.na(keyvals.colourC2)] <- 'gray'
names(keyvals.colourC2)[keyvals.colourC2 == 'red'] <- 'Upregulated'
names(keyvals.colourC2)[keyvals.colourC2 == 'gray'] <- 'Mid'
names(keyvals.colourC2)[keyvals.colourC2 == 'blue'] <- 'Downregulated'


C2EV <- EnhancedVolcano(Cluster2_uniqgenes,
    lab = Cluster2_uniqgenes$genetags,
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = Cluster2_uniqgenes$external_gene_name[which(names(keyvalsC2) %in% c('Highest', 'Lowest'))],
    title = "Differentially expressed genes in Cluster 2",
  subtitle = NULL,
  caption = NULL,
  pCutoff = 0.05,
    FCcutoff = 1,
  pointSize = c(ifelse((Cluster2_uniqgenes$padj<0.05 & 
                          (Cluster2_uniqgenes$log2FoldChange >= 1 | Cluster2_uniqgenes$log2FoldChange <= -1)),3, 0.5)),
    labSize = 5,
  #col=c('gray2', 'green2', 'blue2', 'red'),
  shapeCustom = keyvals.shapeC2,
  colCustom = keyvals.colourC2,
  #shape = c(1, 4, 23, 25),
    colAlpha = 0.5,
    shadeAlpha = 1/2,
    cutoffLineType = 'twodash',
    cutoffLineCol = 'black',
    cutoffLineWidth = 0.8,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
  legendLabels=c('Not sig.','Log (base 2) FC','p-value',
      'p-value & Log (base 2) FC'),
    legendPosition = 'top',
    legendLabSize = 14,
    legendIconSize = 2,
  drawConnectors = TRUE,
  arrowheads = FALSE,
    endsConnectors = "first",
   widthConnectors = 0.25,
  labFace = 'bold',
      boxedLabels = F,
  max.overlaps = 10,
  maxoverlapsConnectors = Inf
  )


pdf("EV_Cluster2.pdf", width = 8, height = 5)

C2EV +
    ggplot2::coord_cartesian(xlim=c(-8,8), ylim=c(0,150)) +
    ggplot2::scale_x_continuous(
      breaks=seq(-8,8,4))+
    theme(plot.subtitle = element_text(hjust = 0.5)+
            theme(legend.position="top", legend.box = "horizontal"))+
           guides(shape=FALSE)+
       scale_colour_manual(values = keyvals.colourC2, breaks=c('Downregulated', '', 'Upregulated'))

dev.off()




# Identify the top 20 or desired hits, arrange() apparently goes from AtoZ
Cluster3_10Upgenes <- tail(arrange(Cluster3_onlyDEGs,log2FoldChange),10)
Cluster3_10Downgenes <- head(arrange(Cluster3_onlyDEGs,log2FoldChange),10)

# Add column label, containing the gene name for the top hits or nothing for all others
Cluster3_uniqgenes$genetags <- ifelse(Cluster3_uniqgenes$external_gene_name %in% Cluster3_10Upgenes$external_gene_name,
                                     Cluster3_uniqgenes$external_gene_name, 
                                     ifelse(Cluster3_uniqgenes$external_gene_name %in% Cluster3_10Downgenes$external_gene_name,
                                            Cluster3_uniqgenes$external_gene_name, ""))

keyvalsC3 <- ifelse(
    Cluster3_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(Cluster3_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))
  keyvalsC3[is.na(keyvalsC3)] <- 'gray'
  names(keyvalsC3)[keyvalsC3 == 'red'] <- 'Highest'
  names(keyvalsC3)[keyvalsC3 == 'gray'] <- 'Mid'
  names(keyvalsC3)[keyvalsC3 == 'blue'] <- 'Lowest'

keyvals.shapeC3 <- ifelse(
    Cluster3_uniqgenes$external_gene_name %in% Cluster3_10Downgenes$external_gene_name, 19,
      ifelse(Cluster3_uniqgenes$external_gene_name %in% Cluster3_10Upgenes$external_gene_name, 19,
        20))
  
  keyvals.shapeC3[is.na(keyvals.shapeC3)] <- 20
  names(keyvals.shapeC3)[keyvals.shapeC3 == 20] <- 'Mid'
  names(keyvals.shapeC3)[keyvals.shapeC3 == 19] <- 'Lowest'
  names(keyvals.shapeC3)[keyvals.shapeC3 == 19] <- 'Highest'

  
keyvals.colourC3 <- ifelse(
    Cluster3_uniqgenes$pval < 0.05 & Cluster3_uniqgenes$log2FoldChange <= -1, 'blue',
      ifelse(Cluster3_uniqgenes$pval < 0.05 & Cluster3_uniqgenes$log2FoldChange >= 1, 'red',
        'gray'))

keyvals.colourC3[is.na(keyvals.colourC3)] <- 'gray'
names(keyvals.colourC3)[keyvals.colourC3 == 'red'] <- 'Upregulated'
names(keyvals.colourC3)[keyvals.colourC3 == 'gray'] <- 'Mid'
names(keyvals.colourC3)[keyvals.colourC3 == 'blue'] <- 'Downregulated'


C3EV <- EnhancedVolcano(Cluster3_uniqgenes,
    lab = Cluster3_uniqgenes$genetags,
    x = 'log2FoldChange',
    y = 'pvalue',
    selectLab = Cluster3_uniqgenes$external_gene_name[which(names(keyvalsC3) %in% c('Highest', 'Lowest'))],
    title = "Differentially expressed genes in Cluster 3",
  subtitle = NULL,
  caption = NULL,
  pCutoff = 0.05,
    FCcutoff = 1,
  pointSize = c(ifelse((Cluster3_uniqgenes$padj<0.05 & 
                          (Cluster3_uniqgenes$log2FoldChange >= 1 | Cluster3_uniqgenes$log2FoldChange <= -1)),3, 0.5)),
    labSize = 5,
  #col=c('gray2', 'green2', 'blue2', 'red'),
  shapeCustom = keyvals.shapeC3,
  colCustom = keyvals.colourC3,
  #shape = c(1, 4, 23, 25),
    colAlpha = 0.5,
    shadeAlpha = 1/2,
    cutoffLineType = 'twodash',
    cutoffLineCol = 'black',
    cutoffLineWidth = 0.8,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
  legendLabels=c('Not sig.','Log (base 2) FC','p-value',
      'p-value & Log (base 2) FC'),
    legendPosition = 'top',
    legendLabSize = 14,
    legendIconSize = 2,
  drawConnectors = TRUE,
  arrowheads = FALSE,
    endsConnectors = "first",
   widthConnectors = 0.25,
  labFace = 'bold',
      boxedLabels = F,
  max.overlaps = 10,
  maxoverlapsConnectors = Inf
  )


pdf("EV_Cluster3.pdf", width = 8, height = 5)

C3EV +
    ggplot2::coord_cartesian(xlim=c(-8,8), ylim=c(0,80)) +
    ggplot2::scale_x_continuous(
      breaks=seq(-8,8,4))+
    theme(plot.subtitle = element_text(hjust = 0.5)+
            theme(legend.position="top", legend.box = "horizontal"))+
           guides(shape=FALSE)+
       scale_colour_manual(values = keyvals.colourC3, breaks=c('Downregulated', '', 'Upregulated'))

dev.off()





```

```{r}
#save(DvB_uniqgenes, file = "DvB_uniqgenes.RData")

setwd("C:/Users/SaulKarr/OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey/DBT_data/neubulkRNA")

#https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

#Usaremos el df con los geneIDs de ENSEMBLE en rownames

library(clusterProfiler)
library(gprofiler2)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(cowplot)
library(fgsea)
library(patchwork)

# SET THE DESIRED ORGANISM
#https://bioconductor.org/packages/release/BiocViews.html#___OrgDb

organism = "org.Hs.eg.db"

#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#Prepare Input

# data from deseq2
load("DBvB_uniqprotcodgenes.RData")
load("SexFvM_uniqprotcodgenes.RData")
load("AgeOvY_uniqprotcodgenes.RData")
load("Cluster1_uniqprotcodgenes.RData")
load("Cluster2_uniqprotcodgenes.RData")
load("Cluster3_uniqprotcodgenes.RData")

DvB_uniqgenes <- DBvB_uniqgenes


# we want the log2 fold change 
original_gene_listDvB    <- DvB_uniqgenes$log2FoldChange
original_gene_listSexFvM <- SexFvM_uniqgenes$log2FoldChange
original_gene_listAgeOvY <- AgeOvY_uniqgenes$log2FoldChange
original_gene_listC1 <- Cluster1_uniqgenes$log2FoldChange
original_gene_listC2 <- Cluster2_uniqgenes$log2FoldChange
original_gene_listC3 <- Cluster3_uniqgenes$log2FoldChange

# name the vector
names(original_gene_listDvB) <- DvB_uniqgenes$external_gene_name
names(original_gene_listSexFvM) <- SexFvM_uniqgenes$external_gene_name
names(original_gene_listAgeOvY) <- AgeOvY_uniqgenes$external_gene_name
names(original_gene_listC1) <- Cluster1_uniqgenes$external_gene_name
names(original_gene_listC2) <- Cluster2_uniqgenes$external_gene_name
names(original_gene_listC3) <- Cluster3_uniqgenes$external_gene_name

# omit any NA values 
gene_listDvB<-na.omit(original_gene_listDvB)
gene_listSexFvM<-na.omit(original_gene_listSexFvM)
gene_listAgeOvY<-na.omit(original_gene_listAgeOvY)
gene_listC1<-na.omit(original_gene_listC1)
gene_listC2<-na.omit(original_gene_listC2)
gene_listC3<-na.omit(original_gene_listC3)


# sort the list in decreasing order (required for clusterProfiler)
gene_listDvB = sort(gene_listDvB, decreasing = TRUE)
gene_listSexFvM = sort(gene_listSexFvM, decreasing = TRUE)
gene_listAgeOvY = sort(gene_listAgeOvY, decreasing = TRUE)
gene_listC1 = sort(gene_listC1, decreasing = TRUE)
gene_listC2 = sort(gene_listC2, decreasing = TRUE)
gene_listC3 = sort(gene_listC3, decreasing = TRUE)

## Gene Set Enrichment
Params:  
  
**keyType** This is the source of the annotation (gene ids). The options vary for each annotation. In the example of *org.Dm.eg.db*, the options are:   
  
"ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"      
"ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "FLYBASE"      "FLYBASECG"    "FLYBASEPROT"   
"GENENAME"     "GO"           "GOALL"        "MAP"          "ONTOLOGY"     "ONTOLOGYALL"   
"PATH"         "PMID"         "REFSEQ"       "SYMBOL"       "UNIGENE"      "UNIPROT"  
  
Check which options are available with the `keytypes` command, for example `keytypes(org.Hs.eg.db)`. 

keytypes(org.Hs.eg.db)
  
**ont** one of "BP", "MF", "CC" or "ALL"  
**nPerm** permutation numbers, the higher the number of permutations you set, the more accurate your results is, but it will also cost longer time for running permutation.  
**minGSSize** minimal size of each geneSet for analyzing.   
**maxGSSize** maximal size of genes annotated for testing.   
**pvalueCutoff** pvalue Cutoff.   
**pAdjustMethod** one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" 

gseDvB <-    gseGO(gene_listDvB, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db")
gseSexFvM <- gseGO(gene_listSexFvM, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db")
gseAgeOvY <- gseGO(gene_listAgeOvY, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db")
gseC1 <-     gseGO(gene_listC1, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db")
gseC2 <-     gseGO(gene_listC2, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db")
gseC3 <-     gseGO(gene_listC3, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db")

gseDvBfdr <-    gseGO(gene_listDvB, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", pAdjustMethod = "fdr")
gseSexFvMfdr <- gseGO(gene_listSexFvM, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", pAdjustMethod = "fdr")
gseAgeOvYfdr <- gseGO(gene_listAgeOvY, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", pAdjustMethod = "fdr")
gseC1fdr <-     gseGO(gene_listC1, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", pAdjustMethod = "fdr")
gseC2fdr <-     gseGO(gene_listC2, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", pAdjustMethod = "fdr")
gseC3fdr <-     gseGO(gene_listC3, ont = "ALL", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", pAdjustMethod = "fdr")


save(gseDvB, file =   "gseDvB.RData")
save(gseSexFvM, file ="gseSexFvM.RData")
save(gseAgeOvY, file ="gseAgeOvY.RData")
save(gseC1, file =    "gseC1.RData")
save(gseC2, file =    "gseC2.RData")
save(gseC3, file =    "gseC3.RData")

load(file =   "gseDvB.RData")
load(file ="gseSexFvM.RData")
load(file ="gseAgeOvY.RData")
load(file =    "gseC1.RData")
load(file =    "gseC2.RData")
load(file =    "gseC3.RData")


```

##Dotplot

```{r echo=TRUE, fig.width=15, fig.height=8}
require(DOSE)

png("dptest1.png", 
    res = 300, width = 30, height = 35, units = "cm")

dp1 <- dotplot(gseDvB, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
#+ scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100))

dp1

dev.off()




pdf("dp1test100.pdf", width = 10, height = 10)          #w and h in inches

# Creamos un gráfico
dp1 <- dotplot(gseDvB, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
#+ scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100))

dp1

# Cerramos el dispositivo gráfico
dev.off() 


#trying to convert to ggplot format, didn't work


png("dp2.png", 
    res = 300, width = 30, height = 35, units = "cm")

dp2 <- dotplot(gseSexFvM, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
#+ scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100))

dp2

dev.off()


pdf("dp2test100.pdf", width = 10, height = 10)

dp2 <- dotplot(gseSexFvM, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
dp2

dev.off()


png("dp3.png", 
    res = 300, width = 30, height = 40, units = "cm")

dp3 <- dotplot(gseAgeOvY, showCategory=10, split=".sign", font.size = 16, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
#+ scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100))

dp3

dev.off()


pdf("dp3test100.pdf", width = 10, height = 10)

dp3 <- dotplot(gseAgeOvY, showCategory=10, split=".sign", font.size = 16, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
#+ scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100))

dp3

dev.off()

png("dp4.png", 
    res = 300, width = 30, height = 35, units = "cm")

dp4 <- dotplot(gseC1, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
dp4

dev.off()

pdf("dp4test100.pdf", width = 10, height = 10)

dp4 <- dotplot(gseC1, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
dp4

dev.off()



png("dp5.png", 
    res = 300, width = 30, height = 40, units = "cm")

dp5 <- dotplot(gseC2, showCategory=10, split=".sign", font.size = 16, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
#+ scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100))

dp5

dev.off()

pdf("dp5test100.pdf", width = 10, height = 10)

dp5 <- dotplot(gseC2, showCategory=10, split=".sign", font.size = 16, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
dp5

dev.off()



png("dp6.png", 
    res = 300, width = 30, height = 40, units = "cm")

dp6 <- dotplot(gseC3, showCategory=10, split=".sign", font.size = 16, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
#+ scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100))

dp6

dev.off()


pdf("dp6test100.pdf", width = 10, height = 10)

dp6 <- dotplot(gseC3, showCategory=10, split=".sign", font.size = 16, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
dp6

dev.off()



pdf("dp6testfdr.pdf", width = 10, height = 10)

dp6 <- dotplot(gseC3fdr, showCategory=10, split=".sign", font.size = 16, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "GO terms", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
dp6

dev.off()

```

##Encrichment plot map: Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. In this way, mutually overlapping gene sets are tend to cluster together, making it easy to identify functional modules.

```{r echo=TRUE}
#library(GOSemSim)
emapplot(gse, showCategory = 1)

#pairwise_termsim(x, method = "JC", semData = NULL, showCategory = 200)

#gene <- names(gene_list)[abs(gene_list) > 2]
#ego <- enrichGO(gene = gene,
#universe = names(gene_list),
#OrgDb = org.Hs.eg.db,
#ont = "BP",
#pAdjustMethod = "BH",
#pvalueCutoff = 0.01,
#qvalueCutoff = 0.05,
#readable = TRUE)
#d <- godata('org.Hs.eg.db', ont="BP")
#ego2 <- pairwise_termsim(ego, method="Wang", semData = d)
#emapplot(ego2)
#emapplot_cluster(ego2)

```

##Category Netplot The cnetplot depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network (helpful to see which genes are involved in enriched pathways and genes that may belong to multiple annotation categories).

```{r fig.width=18}

#install.packages("ggnewscale",dependencies = TRUE)
# categorySize can be either 'pvalue' or 'geneNum'

png("cnet1.png", 
    res = 300, width = 26, height = 28, units = "cm")

cnet1 <- cnetplot(gseDvB, categorySize="pvalue", foldChange=gene_listDvB, 
         showCategory = 3,
        cex_category = 1,
  cex_gene = 1,
  node_label_size = NULL,
  cex_label_category = 1.5,
  cex_label_gene = 1.8,
         circular = F, colorEdge = TRUE,
          )+
theme(legend.position = "bottom", legend.box="vertical",
      legend.title = element_text(size = 15, face="bold"),
          legend.text = element_text(size=12, face="bold"))+
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 3),
         category = guide_legend(order = 2))+
  guides(size = FALSE)

cnet1

dev.off()

pdf("cnet1test100.pdf", width = 10, height = 10)

cnet1 <- cnetplot(gseDvB, categorySize="pvalue", foldChange=gene_listDvB, 
         showCategory = 3,
        cex_category = 1,
  cex_gene = 1,
  node_label_size = NULL,
  cex_label_category = 1.5,
  cex_label_gene = 1.8,
         circular = F, colorEdge = TRUE,
          )+
theme(legend.position = "bottom", legend.box="vertical",
      legend.title = element_text(size = 15, face="bold"),
          legend.text = element_text(size=12, face="bold"))+
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 3),
         category = guide_legend(order = 2))+
  guides(size = FALSE)

cnet1

dev.off()


heatplot1 <- heatplot(gseDvB, foldChange=gene_listDvB, showCategory=5)
heatplot1


png("cnet2.png", 
    res = 300, width = 26, height = 28, units = "cm")

cnet2 <- cnetplot(gseSexFvM, categorySize="pvalue", foldChange=gene_listSexFvM, 
         showCategory = 3,
        cex_category = 1,
  cex_gene = 1,
  node_label_size = NULL,
  cex_label_category = 1.2,
  cex_label_gene = 1,
         circular = F, colorEdge = TRUE,
          )+
theme(legend.position = "bottom", legend.box="vertical",
      legend.title = element_text(size = 15, face="bold"),
          legend.text = element_text(size=12, face="bold"))+
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 3),
         category = guide_legend(order = 2))+
  guides(size = FALSE)

cnet2

dev.off()

pdf("cnet2test100.pdf", width = 10, height = 10)

cnet2 <- cnetplot(gseSexFvM, categorySize="pvalue", foldChange=gene_listSexFvM, 
         showCategory = 3,
        cex_category = 1,
  cex_gene = 1,
  node_label_size = NULL,
  cex_label_category = 1.2,
  cex_label_gene = 1,
         circular = F, colorEdge = TRUE,
          )+
theme(legend.position = "bottom", legend.box="vertical",
      legend.title = element_text(size = 15, face="bold"),
          legend.text = element_text(size=12, face="bold"))+
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 3),
         category = guide_legend(order = 2))+
  guides(size = FALSE)

cnet2

dev.off()



png("cnet3.png", 
    res = 300, width = 28, height = 30, units = "cm")

cnet3 <- cnetplot(gseAgeOvY, categorySize="pvalue", foldChange=gene_listAgeOvY, 
         showCategory = 3,
        cex_category = 1,
  cex_gene = 1,
  node_label_size = NULL,
  cex_label_category = 1.2,
  cex_label_gene = 1.5,
         circular = F, colorEdge = TRUE,
          )+
theme(legend.position = "bottom", legend.box="vertical",
      legend.title = element_text(size = 12, face="bold"),
          legend.text = element_text(size=10, face="bold"))+
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 3),
         category = guide_legend(order = 2))+
  guides(size = FALSE)

cnet3

dev.off()

pdf("cnet3test100.pdf", width = 10, height = 10)

cnet3 <- cnetplot(gseAgeOvY, categorySize="pvalue", foldChange=gene_listAgeOvY, 
         showCategory = 3,
        cex_category = 1,
  cex_gene = 1,
  node_label_size = NULL,
  cex_label_category = 1.2,
  cex_label_gene = 1.5,
         circular = F, colorEdge = TRUE,
          )+
theme(legend.position = "bottom", legend.box="vertical",
      legend.title = element_text(size = 12, face="bold"),
          legend.text = element_text(size=10, face="bold"))+
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 3),
         category = guide_legend(order = 2))+
  guides(size = FALSE)

cnet3

dev.off()

```

## Ridgeplot

Helpful to interpret up/down-regulated pathways.

```{r fig.width=18, fig.height=12}
install.packages("ggridges", dependencies = T)

ridgeplot(gse) + labs(x = "enrichment distribution")
ridgeplot(gse) + labs(x = "enrichment distribution") + theme(text = element_text(size=1)) + theme_minimal()

```

## GSEA Plot

Traditional method for visualizing GSEA result.

Params:\
**Gene Set** Integer. Corresponds to gene set in the gse object. The first gene set is 1, second gene set is 2, etc.

```{r fig.height=6}
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
```

## PubMed trend of enriched terms

Plots the number/proportion of publications trend based on the query result from PubMed Central.

```{r fig.width=10}
install.packages("europepmc", dependencies = T)
terms <- gse$Description[1:3]
pmcplot(terms, 2010:2018, proportion=FALSE)
```

```{r}
library(DOSE)
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
idsDvB<-bitr(names(original_gene_listDvB), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

idsSexFvM<-bitr(names(original_gene_listSexFvM), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
idsAgeOvY<-bitr(names(original_gene_listAgeOvY), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
idsC1<-bitr(names(original_gene_listC1), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
idsC2<-bitr(names(original_gene_listC2), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
idsC3<-bitr(names(original_gene_listC3), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_idsDvB = idsDvB[!duplicated(idsDvB[c("ENTREZID")]),]

dedup_idsSexFvM = idsSexFvM[!duplicated(idsSexFvM[c("ENTREZID")]),]
dedup_idsAgeOvY = idsAgeOvY[!duplicated(idsAgeOvY[c("ENTREZID")]),]
dedup_idsC1 = idsC1[!duplicated(idsC1[c("ENTREZID")]),]
dedup_idsC2 = idsC2[!duplicated(idsC2[c("ENTREZID")]),]
dedup_idsC3 = idsC3[!duplicated(idsC3[c("ENTREZID")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
DvB_uniqgenes2 = DvB_uniqgenes[DvB_uniqgenes$external_gene_name %in% dedup_idsDvB$SYMBOL,]
DvB_uniqgenes2$SYMBOL = DvB_uniqgenes2$external_gene_name
DvB_genesEntrez <- left_join(DvB_uniqgenes2,dedup_idsDvB, by="SYMBOL")

SexFvM_uniqgenes2 = SexFvM_uniqgenes[SexFvM_uniqgenes$external_gene_name %in% dedup_idsSexFvM$SYMBOL,]
SexFvM_uniqgenes2$SYMBOL = SexFvM_uniqgenes2$external_gene_name
SexFvM_genesEntrez <- left_join(SexFvM_uniqgenes2,dedup_idsSexFvM, by="SYMBOL")

AgeOvY_uniqgenes2 = AgeOvY_uniqgenes[AgeOvY_uniqgenes$external_gene_name %in% dedup_idsAgeOvY$SYMBOL,]
AgeOvY_uniqgenes2$SYMBOL = AgeOvY_uniqgenes2$external_gene_name
AgeOvY_genesEntrez <- left_join(AgeOvY_uniqgenes2,dedup_idsAgeOvY, by="SYMBOL")

C1_uniqgenes2 = Cluster1_uniqgenes[Cluster1_uniqgenes$external_gene_name %in% dedup_idsC1$SYMBOL,]
C1_uniqgenes2$SYMBOL = C1_uniqgenes2$external_gene_name
C1_genesEntrez <- left_join(C1_uniqgenes2,dedup_idsC1, by="SYMBOL")

C2_uniqgenes2 = Cluster2_uniqgenes[Cluster2_uniqgenes$external_gene_name %in% dedup_idsC2$SYMBOL,]
C2_uniqgenes2$SYMBOL = C2_uniqgenes2$external_gene_name
C2_genesEntrez <- left_join(C2_uniqgenes2,dedup_idsC2, by="SYMBOL")

C3_uniqgenes2 = Cluster3_uniqgenes[Cluster3_uniqgenes$external_gene_name %in% dedup_idsC3$SYMBOL,]
C3_uniqgenes2$SYMBOL = C3_uniqgenes2$external_gene_name
C3_genesEntrez <- left_join(C3_uniqgenes2,dedup_idsC3, by="SYMBOL")

# Create a new column in df2 with the corresponding ENTREZ IDs
#DvB_uniqgenes2$Y = dedup_idsDvB$ENTREZID

# Create a vector of the gene universe
kegg_gene_listDvB <- DvB_genesEntrez$log2FoldChange

kegg_gene_listSexFvM <- SexFvM_genesEntrez$log2FoldChange
kegg_gene_listAgeOvY <- AgeOvY_genesEntrez$log2FoldChange
kegg_gene_listC1 <- C1_genesEntrez$log2FoldChange
kegg_gene_listC2 <- C2_genesEntrez$log2FoldChange
kegg_gene_listC3 <- C3_genesEntrez$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_listDvB) <- DvB_genesEntrez$ENTREZID

names(kegg_gene_listSexFvM) <- SexFvM_genesEntrez$ENTREZID
names(kegg_gene_listAgeOvY) <- AgeOvY_genesEntrez$ENTREZID
names(kegg_gene_listC1) <- C1_genesEntrez$ENTREZID
names(kegg_gene_listC2) <- C2_genesEntrez$ENTREZID
names(kegg_gene_listC3) <- C3_genesEntrez$ENTREZID

# omit any NA values 
kegg_gene_listDvB<-na.omit(kegg_gene_listDvB)

kegg_gene_listSexFvM<-na.omit(kegg_gene_listSexFvM)
kegg_gene_listAgeOvY<-na.omit(kegg_gene_listAgeOvY)
kegg_gene_listC1<-na.omit(kegg_gene_listC1)
kegg_gene_listC2<-na.omit(kegg_gene_listC2)
kegg_gene_listC3<-na.omit(kegg_gene_listC3)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_listDvB = sort(kegg_gene_listDvB, decreasing = TRUE)

kegg_gene_listSexFvM = sort(kegg_gene_listSexFvM, decreasing = TRUE)
kegg_gene_listAgeOvY = sort(kegg_gene_listAgeOvY, decreasing = TRUE)
kegg_gene_listC1 = sort(kegg_gene_listC1, decreasing = TRUE)
kegg_gene_listC2 = sort(kegg_gene_listC2, decreasing = TRUE)
kegg_gene_listC3 = sort(kegg_gene_listC3, decreasing = TRUE)

## Create gseKEGG object
 
**organism** KEGG Organism Code: The full list is here: https://www.genome.jp/kegg/catalog/org_list.html (need the 3 letter code). I define this as `kegg_organism` first, because it is used again below when making the pathview plots.  
**nPerm** permutation numbers, the higher the number of permutations you set, the more accurate your results is, but it will also cost longer time for running permutation.  
**minGSSize** minimal size of each geneSet for analyzing.   
**maxGSSize** maximal size of genes annotated for testing.   
**pvalueCutoff** pvalue Cutoff.   
**pAdjustMethod** one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".  
**keyType** one of 'kegg', 'ncbi-geneid', 'ncib-proteinid' or 'uniprot'.  

kegg_organism = "hsa"
kk2DvB <- gseKEGG(geneList     = kegg_gene_listDvB,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

kk2SexFvM <- gseKEGG(geneList     = kegg_gene_listSexFvM,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
kk2AgeOvY <- gseKEGG(geneList     = kegg_gene_listAgeOvY,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
kk2C1 <- gseKEGG(geneList     = kegg_gene_listC1,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
kk2C2 <- gseKEGG(geneList     = kegg_gene_listC2,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
kk2C3 <- gseKEGG(geneList     = kegg_gene_listC3,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

save(kk2DvB, file =   "kk2DvB.RData")
save(kk2SexFvM, file ="kk2SexFvM.RData")
save(kk2AgeOvY, file ="kk2AgeOvY.RData")
save(kk2C1, file =    "kk2C1.RData")
save(kk2C2, file =    "kk2C2.RData")
save(kk2C3, file =    "kk2C3.RData")

load("kk2DvB.RData")
load("kk2SexFvM.RData")
load("kk2AgeOvY.RData")
load("kk2C1.RData")
load("kk2C2.RData")
load("kk2C3.RData")


png("k1.png", 
    res = 300, width = 30, height = 35, units = "cm")

k1 <- dotplot(kk2DvB, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "Enriched Pathways", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
k1

dev.off()



png("k2.png", 
    res = 300, width = 30, height = 35, units = "cm")

k2 <- dotplot(kk2SexFvM, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "Enriched Pathways", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
k2

dev.off()



png("k3.png", 
    res = 300, width = 30, height = 35, units = "cm")

k3 <- dotplot(kk2AgeOvY, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "Enriched Pathways", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
k3

dev.off()



png("k4.png", 
    res = 300, width = 30, height = 35, units = "cm")

k4 <- dotplot(kk2C1, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "Enriched Pathways", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
k4

dev.off()





png("k5.png", 
    res = 300, width = 30, height = 35, units = "cm")

k5 <- dotplot(kk2C2, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "Enriched Pathways", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
k5

dev.off()



png("k6.png", 
    res = 300, width = 30, height = 40, units = "cm")

k6 <- dotplot(kk2C3, showCategory=10, split=".sign", font.size = 20, title = "") + facet_grid(.~.sign)+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20, face="bold.italic"))+ 
  theme(legend.position="right", 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) + 
  labs(title = "Enriched Pathways", y = "", x="Gene ratio") + 
  theme(plot.title = element_text(size=20, hjust=0.5), axis.title.x = element_text(size=20, hjust=0.5))
  
k6

dev.off()


```











