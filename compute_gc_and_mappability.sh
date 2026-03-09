#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --job-name=gc_map
***REMOVED***
***REMOVED***
#SBATCH --mail-type=BEGIN,END,FAIL

#this script will take a reference genome (fasta file) and will compute both the gc content and the mappability of specified bins
#for this to work we need art_illumina, bwa, samtools and bedtools. 
#art ilumina is in the scDNAseq container.

module --force purge
module load calcua/2024a            # base toolchain (GCC 13.3.0 family)

# Load matching builds (all with GCC-13.3.0)
module load BWA/0.7.18-GCCcore-13.3.0
module load BEDTools/2.31.1-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0

#load the container containing art_ilumina
unset PYTHONPATH
export PATH="/scratch/antwerpen/205/vsc20542/containers/scDNAseq/bin:$PATH"

set -e
# ======== Parameters ===========
REF_PATH="/scratch/antwerpen/205/vsc20542/reference_genomes/Leishmania/L_donovani/LdCL/LdCL_original.fasta" # reference genome file path
OUTPUT_DIR="/scratch/antwerpen/205/vsc20542/projects/BPK081_subclones/"
BIN_SIZE=1000 # bin size in bp
SIM_LIBRARY_SIZE=400 # simulated library size (in bp)
SIM_LIBRARY_SIZE_SE=60 # simulated library size standard error
SIM_READ_LENGTH=150 # simulated read length (bp)
SIM_COVERAGE=1 # simulated coverage (per base pair)
NTHREADS=8 # number of threads to use

# =========== Create artificial Library ====================
#now we create an artificial fastq file which is used te calculate bin mappability
wd="$OUTPUT_DIR"
wd="$wd/gc_and_mappability/$BIN_SIZE/"
mkdir -p "$wd"

#copy the genome to local dir (so I can have write rights on it)
cp "$REF_PATH" "$wd/ref_genome.fa"
REF_PATH="$wd/ref_genome.fa"

art_illumina -ss HS25 -i "$REF_PATH" -p -l $SIM_READ_LENGTH -f $SIM_COVERAGE -m $SIM_LIBRARY_SIZE -s $SIM_LIBRARY_SIZE_SE -o "$wd/simulated_reads"
rm "$wd"/*.aln


# ========== Map Artificial library to reference genome =======================
#now we map the artificial reads to the reference using bwa
# Check if BWA index files exist
if [[ ! -f "${REF_PATH}.amb" || ! -f "${REF_PATH}.ann" || ! -f "${REF_PATH}.bwt" || \
      ! -f "${REF_PATH}.pac" || ! -f "${REF_PATH}.sa" ]]; then  
    echo "BWA index not found. Running bwa index..."
    bwa index "$REF_PATH"
else
    echo "BWA index already exists. Skipping indexing."
fi

#map the artificial reads to the reference genome
bwa mem -t $NTHREADS "$REF_PATH" "$wd/simulated_reads1.fq" "$wd/simulated_reads2.fq" | \
samtools fixmate -u -m - - | \
samtools sort -u -@ $NTHREADS - | \
samtools markdup -@ $NTHREADS - "$wd/mapped_simulated_reads_sorted.bam"

#keep only reads with mapping quality higher than 30
samtools view -b -q 30 "$wd/mapped_simulated_reads_sorted.bam" > "$wd/mapped_simulated_reads_sorted_filtered.bam"


# =============== Calculate mappability and GC content ========================

#to get coverage per bin we first need to create a genome size file
samtools faidx "$REF_PATH" > "$REF_PATH.fai"
cut -f1,2 "$REF_PATH.fai" > "$wd/genome_sizes.txt"

#then we crete a bed file to store the information
bedtools makewindows -g "$wd/genome_sizes.txt" -w $BIN_SIZE > "$wd/genome_bins.bed"

#get total bin coverage
bedtools coverage -a "$wd/genome_bins.bed" -b "$wd/mapped_simulated_reads_sorted.bam" > "$wd/temp.txt"
#get bin coverage with unique maps only
bedtools coverage -a "$wd/genome_bins.bed" -b "$wd/mapped_simulated_reads_sorted_filtered.bam" > "$wd/temp2.txt"

#combine both results into a single file
paste <(cut -f1-4 "$wd/temp.txt") <(cut -f4-7 "$wd/temp2.txt") > "$wd/temp3.txt"

#now calculate the ratio between total reads mapping and uniquely mapped reads (the mappability). Note: This part was written by ChatGPT. 
awk 'BEGIN {FS="\t"; OFS="\t"}
{
  if ($4 > 0) {
    m = $5 / $4
  } else {
    m = 0
  }
  print $0, m
}' "$wd/temp3.txt" > "$wd/temp3_with_mappability.txt" 

echo -e "chrom\tstart\tend\tsim_total_read_count\tsim_unique_read_count\tsim_bases_covered\tBIN_SIZE\tsim_fraction_covered\tmappability" > "$wd/header.txt"
cat "$wd/header.txt" "$wd/temp3_with_mappability.txt" > "$wd/reads_per_bin.txt"
#get bin gc content
bedtools nuc -fi "$REF_PATH" -bed "$wd/genome_bins.bed" > "$wd/gc_content_per_bin.txt" 
#extract only the column with the gc content
tail -n +1 "$wd/gc_content_per_bin.txt" | cut -f5 > "$wd/temp.txt" 

# =============== Save final results ===================
#append the gc content to the final file
paste "$wd/reads_per_bin.txt" "$wd/temp.txt" > "$OUTPUT_DIR/gc_and_mappability/bin_gc_and_mappability_$BIN_SIZE.tsv"

#remove unecessary files
rm -rf $wd