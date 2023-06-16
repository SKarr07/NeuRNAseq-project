# NeuRNAseq-project

Habiendo seleccionado las muestras bajo los impuestos criterios de inclusión, proceder al análisis de los datos de secuenciación.

## RNA-seq analysis workflow

### 1.Realizar la descarga del rawdata (fastq).

a)Existe la opción de descarga de fastqs directamente desde el ENA browser usando wget.

b)O también más agilmente usando la opción del SRAtoolkit.

```
$ qlogin

$ module load sra/3.0.0
```

i.Se puede hacer de manera individual.

ii.O utilizando una lista con los SRR de las muestras.

$ prefetch --option-file ID_list_names.txt		## Para una lista con los SRR de interés.

#una vez terminado el prefetch,

$ fastq-dump --gzip --skip-technical --dumpbase --split-3 --clip --outdir /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/SRP *.sra

#para una lista de SRR pertenecientes a distintos datasets, hicimos lo siguiente:

$nano list #1.Preparar una lista con las carpetas SRP y sus muestras correspondientes SRR.sra (1 por línea).

#2.Un archivo .sge para automatizar el fastqdump y mediante qsub.

$cat list | while read line; do c=$(echo $line|cut -d' ' -f1); m=$(echo $line|cut -d' ' -f2); fastq-dump --gzip --skip-technical --dumpbase --split-3 --clip --outdir /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/$c $m;done

Una vez descargados los SRR en sus carpetas correspondientes, proceder al FastQC.

### 2.FastQC-Trimming-FastQC.

Primeramente, organizar las muestras por tipo de librería: SE o PE para un análisis en conjunto.

#Hacer una lista con todas las muestras SE y otra con los PE.

$ls /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/SRP*/*fastq.gz >> newSEfastq.txt

$ls /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/SRP*/*fastq.gz >> newPEfastq.txt

$mkdir newSEapr2023 newPEapr2023

$cd newSEapr2023
$for i in $(cat /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/newSEfastq.txt); do echo $i; ln -s $i .; done

$cd newPEapr2023
$for i in $(cat /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/newPEfastq.txt); do echo $i; ln -s $i .; done

Una vez teniendo los symlinks podemos trabajar con ellos.

En este caso, hice un .sh con los 3 pasos de FastQC-Trimming-FastQc consecutivos.

#FastQCTrimmed_wflow.sh

#!/bin/bash

#USAGE= ./FastQCTrimmed_wflow.sh
sed -i 's/\r//' FastQCTrimmed_wflow.sh

#AUTHOR: Saul Karr (adaptado de Eve Coss)

###a partir de # /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/newSEapr2023

#PARTE 1.- FastQC y multiQC
fastqc ./*.fastq.gz -o ./newFQCout
multiqc ./newFQCout -o ./newFQCout

#PARTE 2.- Limpieza de adaptadores
#single-end
cd newSEapr2023
for i in *.fastq.gz;
do echo
trimmomatic SE -threads 8 -phred33 $i data_trimmed/"${i%.fastq}_trimmed.fastq.gz" ILLUMINACLIP:/mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:22
done

#PARTE 3.- FastQC y multiQC
cd newSEapr2023
fastqc ./data_trimmed/*.fastq.gz -o ./FastQC_trimmed

#Reporte en MultiQC
multiqc ./FastQC_trimmed -o ./FastQC_trimmed

#Finalmente crear un .sge para enviar por qsub.
SI se va a enviar por qsub, cuidar de pedir los nucleos adecuados:
#$ -pe openmp 8 #para trimmomatic SE -threads 8

#Repetir y adaptar para muestras PE.

### 3.Alignment with STAR y FastQC.

Una vez revisadas las muestras post-trimming con FastQC, proceder al alineamiento.

#para STAR, correr primero el índice, si no se ha hecho antes

STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir  /mnt/Citosina/amedina/skarr/neu/GEO_bulkRNAseq/STAR_index \
--genomeFastaFiles /mnt/Archives/genome/human/GRCh38/ensembl76/chromosomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /mnt/Archives/genome/human/GRCh38/ensembl76/GTF-file/Homo_sapiens.GRCh38.76.gtf \
--sjdbOverhang 99

#STAR para SE

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

SI se va a enviar por qsub, hay que solicitar los nucleos adecuados:
#$ -pe openmp 20 #para STAR --runThreadN 20

#STAR para PE

Creé un .sh para el alineamiento con STAR y el posterior fastqc de los . bam generados.

#PARTE 4.- STAR paired end reads

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

#PARTE 5.- FastQC y multiQC
cd newPEapr2023
fastqc ./STAR_output/*.out.bam -o ./BamQC

#Reporte en MultiQC
multiqc ./BamQC -o ./BamQC

Se verificaron aquellos archivos en /BamQC correspondientes a muestras con conflictos post-trimming p.ej. contenido de adaptadores y secuencias sobrerrepresentadas.


### 4.Exportar data para R.

Junté todos los ReadsPerGene.out.tab de distintos destinos usando $ln -s a partir de una lista.

Después en R conectarse a la terminal.

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
