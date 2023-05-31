# NeuRNAseq-project

Habiendo seleccionado las muestras bajo los impuestos criterios de inclusión, proceder al análisis de los datos de secuenciación.

## RNA-seq analysis workflow

### 1.Realizar la descarga del rawdata (fastq).

a)Existe la opción de descarga de fastqs directamente desde el ENA browser usando wget.

b)O también más agilmente usando la opción del SRAtoolkit.

$ qlogin

$ module load sra/3.0.0

i.Se puede hacer de manera individual.

ii.O utilizando una lista con los SRR de las muestras.

$ prefetch --option-file ID_list_names.txt		## Para una lista con los SRR de interés.

#una vez terminado el prefetch,

$ fastq-dump --gzip --skip-technical --dumpbase --split-3 --clip --outdir /mnt/Citosina/amedina/skarr/neu/monorail/SRP *.sra

#para una lista de SRR pertenecientes a distintos datasets, hicimos lo siguiente:

$nano list #1.Preparar una lista con las carpetas SRP y sus muestras correspondientes SRR.sra (1 por línea).

#2.Un archivo .sge para automatizar el fastqdump y mediante qsub.

$cat list | while read line; do c=$(echo $line|cut -d' ' -f1); m=$(echo $line|cut -d' ' -f2); fastq-dump --gzip --skip-technical --dumpbase --split-3 --clip --outdir /mnt/Citosina/amedina/skarr/neu/monorail/$c $m;done

Una vez descargados los SRR en sus carpetas correspondientes, proceder al FastQC.

### 2.FastQC-Trimming-FastQC.

Primeramente, organizar las muestras por tipo de librería: SE o PE para un análisis en conjunto.

#Hacer una lista con todas las muestras SE y otra con los PE.

$ls /mnt/Citosina/amedina/skarr/neu/monorail/SRP*/*fastq.gz >> newSEfastq.txt
$ls /mnt/Citosina/amedina/skarr/neu/monorail/SRP*/*fastq.gz >> newPEfastq.txt

$mkdir newSEapr2023 newPEapr2023

$cd newSEapr2023
$for i in $(cat /mnt/Citosina/amedina/skarr/neu/monorail/newSEfastq.txt); do echo $i; ln -s $i .; done

$cd newPEapr2023
$for i in $(cat /mnt/Citosina/amedina/skarr/neu/monorail/newPEfastq.txt); do echo $i; ln -s $i .; done

Una vez teniendo los symlinks podemos trabajar con ellos.

En este caso, hice un .sh con los 3 pasos de FastQC-Trimming-FastQc consecutivos.

#FastQCTrimmed_wflow.sh

#!/bin/bash

#USAGE= ./FastQCTrimmed_wflow.sh
sed -i 's/\r//' FastQCTrimmed_wflow.sh

#AUTHOR: Saul Karr (adaptado de Eve Coss)

###a partir de # /mnt/Citosina/amedina/skarr/neu/monorail/newSEapr2023

#PARTE 1.- FastQC y multiQC
fastqc ./*.fastq.gz -o ./newFQCout
multiqc ./newFQCout -o ./newFQCout

#PARTE 2.- Limpieza de adaptadores
#single-end
cd newSEapr2023
for i in *.fastq.gz;
do echo
trimmomatic SE -threads 8 -phred33 $i data_trimmed/"${i%.fastq}_trimmed.fastq.gz" ILLUMINACLIP:/mnt/Citosina/amedina/skarr/neu/monorail/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:22
done

#PARTE 3.- FastQC y multiQC
cd newSEapr2023
fastqc ./data_trimmed/*.fastq.gz -o ./FastQC_trimmed

#Reporte en MultiQC
multiqc ./FastQC_trimmed -o ./FastQC_trimmed

#Finalmente crear un .sge para enviar por qsub.

#Repetir y adaptar para muestras PE.

### 3.Alignment with STAR y FastQC.
### 4.Exportar data para R.
