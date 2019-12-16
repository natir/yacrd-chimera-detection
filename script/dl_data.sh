#!/bin/bash

# references

mkdir -p references

## E. coli CFT073
curl "https://www.ebi.ac.uk/ena/data/view/AE014075&display=fasta" | seqtk seq -A - > references/ref_e_coli_cft073.fasta

## H. sapiens chr1
curl "ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz" | seqtk seq -A - > references/h_sapiens_chr1_ref.fasta

## D. melanogaster
curl "https://www.ebi.ac.uk/ena/data/view/AE014298.5,AE014134.6,AE013599.5,AE014296.5,AE014297.3,AE014135.4,CP007106.1,KJ947872.2&display=fasta" | seqtk seq -A - > references/d_melanogaster_ref.fasta

## C. elegans
curl "https://www.ebi.ac.uk/ena/data/view/BX284601.5,BX284602.5,BX284603.4,BX284604.4,BX284605.5,BX284606.5&display=fasta" | seqtk seq -A - > references/c_elegans_ref.fasta

# reads

mkdir -p reads

## E. coli ONT
curl "ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/000/SRR8494940/SRR8494940_1.fastq.gz" | seqtk seq -A - > reads/SRR8494940_ont.fasta

## E. coli RSII
curl "ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/001/SRR8494911/SRR8494911_subreads.fastq.gz" | seqtk seq -A - > reads/SRR8494911_pb.fasta

## H. sapiens ONT
curl "http://s3.amazonaws.com/nanopore-human-wgs/chr1.sorted.bam" | samtools bam2fq - | seqtk seq -A - > reads/h_sapiens_chr1_ont.fasta

## D. melanogaster ONT
curl "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR670/003/SRR6702603/SRR6702603_1.fastq.gz" | seqtk seq -A - > reads/d_melanogaster_reads_ont.fasta

## C. elegans RSII
rm reads/c_elegans_pb.fasta
for i in $(curl http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/wget.html | grep "fasta" | cut -d\" -f2 | grep "40X")
do
    curl "http://datasets.pacb.com.s3.amazonaws.com${i}" >> reads/c_elegans_pb.fasta
done
    


