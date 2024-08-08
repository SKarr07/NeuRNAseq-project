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
