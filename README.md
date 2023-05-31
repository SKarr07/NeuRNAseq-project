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

Una vez descargados los SRR, proceder al FastQC.

### 2.FastQC-Trimming-FastQC.

Primeramente organizar las muestras por tipo de librería: SE o PE para un análisis en conjunto.

Para las muestras PE.



### 3.Alignment with STAR y FastQC.
### 4.Exportar data para R.
