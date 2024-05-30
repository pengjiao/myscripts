#!/bin/bash

# Check the number of command-line arguments
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <raw_dir> <project_dir> <nproc>"
    exit 1
fi

# Get variables from command-line arguments
raw_dir=$1
project_dir=$2
nproc=$3

# Quality control and directory setup
cd ${raw_dir} || { echo "Cannot enter raw data directory"; exit 1; }
fastqc *fq.gz
multiqc *zip
gunzip *fq.gz

cd ${project_dir} || { echo "Cannot enter project directory"; exit 1; }
mkdir -p preprocessed/{trimmed_first,trimmed_primer,annotated_reads,consensus_seq,syn_reads,assembled_reads,quality_reads} logs

# Filtering and preprocessing

## Remove low quality reads: reads with mean Phred quality scores less than 10 (-q 10) are removed
for file in ${raw_dir}*fq; do
    file_name=$(basename $file .fq)
    echo "Processing ${file_name}"
    FilterSeq.py quality -s ${file} -q 10 --nproc ${nproc} --outname ${file_name}_filtered --outdir preprocessed --log logs/${file_name}.log
done

## Identify Primers and UMI
### Trim first few bases before AAGCAGTGGTATCAA
for file in ${raw_dir}*L1_1.fq; do
    file_name=$(basename $file .fq)
    echo $file_name
    awk '{
        if (NR % 4 == 2) {
            idx = index($0, "AAGCAGTGGTATCAA");
            if (idx) {
                seq = substr($0, idx);
                print seq;
                getline;
                print;
                getline;
                print substr($0, idx);
            } else {
                print;
                getline;
                print;
                getline;
                print;
            }
        } else {
            print;
        }
    }' ${project_dir}/preprocessed/${file_name}_filtered_quality-pass.fastq > ${project_dir}/preprocessed/trimmed_first/${file_name}_filtered_quality-pass_trimmed.fastq
done

### Identify primers and UMI in R1
for file in ${raw_dir}*L1_1.fq; do
    file_name=$(basename $file .fq)
    echo $file_name
    MaskPrimers.py extract -s ${project_dir}/preprocessed/trimmed_first/${file_name}_filtered_quality-pass_trimmed.fastq \
        --start 0 --len 25 --mode cut \
        --log ${project_dir}/logs/${file_name}_primers.log \
        --outname ${file_name} --outdir ${project_dir}/preprocessed/trimmed_primer
    
    MaskPrimers.py extract -s ${project_dir}/preprocessed/trimmed_primer/${file_name}_primers-pass.fastq \
        --start 10 --len 5 --mode cut \
        --bf UMI --pf TS --barcode \
        --log ${project_dir}/logs/${file_name}_primers_UMI.log \
        --outname ${file_name}_UMI --outdir ${project_dir}/preprocessed/trimmed_primer
done

## Annotate R2 with Internal C-Region
for file in ${raw_dir}*L1_2.fq; do
    file_name=$(basename $file .fq)
    echo $file_name
    MaskPrimers.py align -s ${project_dir}/preprocessed/${file_name}_filtered_quality-pass.fastq \
        -p ${raw_dir}/Mouse_IG_CRegion_RC.fasta \
        --mode cut --pf C_CALL --skiprc --failed \
        --log ${project_dir}/logs/${file_name}_primers_UMI.log \
        --outname ${file_name}_CRegion --outdir ${project_dir}/preprocessed/trimmed_primer
done

## Copy Annotations Between Reads
for file in ${raw_dir}*L1_1.fq; do
    file_name=$(basename $file _L1_1.fq)
    echo $file_name
    PairSeq.py -1 ${project_dir}/preprocessed/trimmed_primer/${file_name}_L1_1_UMI_primers-pass.fastq \
        -2 ${project_dir}/preprocessed/trimmed_primer/${file_name}_L1_2_CRegion_primers-pass.fastq \
        --1f UMI --2f C_CALL --coord illumina \
        --outdir ${project_dir}/preprocessed/annotated_reads
done

## Generation of UMI Consensus Sequences
for file in ${project_dir}/preprocessed/annotated_reads/*fastq; do
    file_name=$(basename $file _primers-pass_pair-pass.fastq)
    echo $file_name
    BuildConsensus.py -s ${project_dir}/preprocessed/annotated_reads/${file_name}_primers-pass_pair-pass.fastq \
        --bf UMI --pf C_CALL --prcons 0.6 \
        -n 1 -q 0 --maxerror 0.1 --maxgap 0.5 \
        --nproc ${nproc} \
        --log ${project_dir}/logs/${file_name}_consensus.log \
        --outdir ${project_dir}/preprocessed/consensus_seq
done

## Synchronize Reads
for file in ${project_dir}/preprocessed/consensus_seq/*_L1_1*fastq; do
    file_name=$(basename $file _L1_1_UMI_primers-pass_pair-pass_consensus-pass.fastq)
    echo $file_name
    PairSeq.py -1 ${project_dir}/preprocessed/consensus_seq/${file_name}_L1_1_UMI_primers-pass_pair-pass_consensus-pass.fastq \
        -2 ${project_dir}/preprocessed/consensus_seq/${file_name}_L1_2_CRegion_primers-pass_pair-pass_consensus-pass.fastq \
        --coord presto --outdir ${project_dir}/preprocessed/syn_reads --outname $file_name
done

## Assemble Pairs
for file in ${raw_dir}*L1_1.fq; do
    file_name=$(basename $file _L1_1.fq)
    echo $file_name
    AssemblePairs.py sequential -1 ${project_dir}/preprocessed/syn_reads/${file_name}-1_pair-pass.fastq \
        -2 ${project_dir}/preprocessed/syn_reads/${file_name}-2_pair-pass.fastq \
        -r ~/share/igblast/fasta/imgt_mouse_ig_v.fasta \
        --coord presto --rc tail --1f CONSCOUNT --2f PRCONS CONSCOUNT \
        --minlen 8 --maxerror 0.3 --alpha 1e-5 --scanrev \
        --minident 0.5 --evalue 1e-5 --maxhits 100 --aligner blastn \
        --log ${project_dir}/logs/${file_name}_assemble.log \
        --outname $file_name --outdir ${project_dir}/preprocessed/assembled_reads
done

## Mask Low-Quality positions
for file in ${project_dir}/preprocessed/assembled_reads/*fastq; do
    file_name=$(basename $file _assemble-pass.fastq)
    echo $file_name
    FilterSeq.py maskqual -s ${project_dir}/preprocessed/assembled_reads/${file_name}_assemble-pass.fastq -q 30 \
        --outname $file_name --outdir ${project_dir}/preprocessed/quality_reads/ \
        --log ${project_dir}/logs/${file_name}_maskqual.log
done

## Track the Number of Sequences that Contributed to the Consensus
for file in ${project_dir}/preprocessed/quality_reads/*_maskqual-pass.fastq; do
    file_name=$(basename $file _maskqual-pass.fastq)
    echo $file_name
    ParseHeaders.py collapse -s ${project_dir}/preprocessed/quality_reads/${file_name}_maskqual-pass.fastq \
        -f CONSCOUNT --act min \
        --outname ${file_name} --outdir ${project_dir}/preprocessed/quality_reads/
done
  
## Collapse Duplicates
for file in ${project_dir}/preprocessed/quality_reads/*_reheader.fastq; do
    file_name=$(basename $file _reheader.fastq)
    echo $file_name
    CollapseSeq.py -s ${project_dir}/preprocessed/quality_reads/${file_name}_reheader.fastq -n 0 \
        --uf PRCONS --cf CONSCOUNT --act sum --inner \
        --keepmiss --outname ${file_name} --outdir ${project_dir}/preprocessed/quality_reads/ \
        --log ${project_dir}/logs/${file_name}_colduplicate.log
done

# Move final files to the final directory
mkdir -p preprocessed/final
for file in ${project_dir}/preprocessed/quality_reads/*_collapse-unique.fastq; do
    file_name=$(basename $file _collapse-unique.fastq)
    echo $file_name
    mv ${file} preprocessed/final/${file_name}_final.fastq
done


mkdir ${project_dir}/annotation_changeo

###################
## Convert FASTQ to FASTA
for file in ${project_dir}/preprocessed/assembled_reads/*fastq
do
  file_name=$(basename $file _assemble-pass.fastq)
  echo $file_name
  ParseHeaders.py rename -s ${project_dir}/preprocessed/final/${file_name}_final.fastq --fasta -f PRCONS -k C_CALL
done

### The wrapper tool AssignGenes.py, from Change-O, uses IgBLAST, and a reference database created with germlines from IMGT, to make V(D)J allele calls
for file in ${project_dir}/preprocessed/assembled_reads/*fastq
do
  file_name=$(basename $file _assemble-pass.fastq)
  echo $file_name
  AssignGenes.py igblast -s ${project_dir}/preprocessed/final/${file_name}_final_reheader.fasta \
 --organism mouse --loci ig --nproc 60\
 -b ~/share/igblast --format airr \
 --outdir  ${project_dir}/annotation_changeo --outname ${file_name}_airr
done