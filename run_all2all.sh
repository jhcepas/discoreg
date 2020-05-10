#!/usr/bin/bash 

if [[ $# -eq 0 ]] ; then
    echo 'provide: FASTAFILE NCPUS'
    exit 0
fi

TARGET_FASTA=`realpath $1`
CORES=$2
#../diamond makedb --in cdhit_cluster_reps.faa -d cdhit_cluster_reps.dmnd
#../diamond blastp -q cdhit_cluster_reps.faa -d cdhit_cluster_reps.dmnd -o cdhit_cluster_hits.dmnd.m8 --sallseqid -k 0 -f 6 qseqid sseqid evalue bitscore length qstart qend sstart send -e 0.1 --more-sensitive --comp-based-stats 0

#~/mmseqs createdb  cdhit_cluster_reps.faa cdhit_cluster_reps.db
#~/mmseqs search cdhit_cluster_reps.db cdhit_cluster_reps.db cdhit_cluster_hits.mmseqs.db tmp/ --start-sens 2 -s 7 --sens-steps 3
#~/mmseqs convertalis  cdhit_cluster_reps.db cdhit_cluster_reps.db cdhit_cluster_hits.mmseqs.db cdhit_cluster_hits.mmseqs.m8  --format-output query,target,evalue,bits,alnlen,qstart,qend,tstart,tend

makeblastdb -in $TARGET_FASTA -dbtype prot
blastp -num_threads $CORES -query $TARGET_FASTA -db $TARGET_FASTA -outfmt "6 qseqid sseqid evalue bitscore length pident qstart qend sstart send"  > $TARGET_FASTA.blast.m8
