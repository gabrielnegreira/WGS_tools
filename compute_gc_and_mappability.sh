#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --job-name=gc_map
#SBATCH --mail-type=BEGIN,END,FAIL

#this script will take a reference genome (fasta file) and will compute both the gc content and the mappability of specified bins
#for this to work we need art_illumina, bwa, samtools and bedtools. 
#art ilumina is in the scDNAseq container.
set -euo pipefail
IFS=$'\n\t'

# --- preflight --- 
##set usage helper
set -euo pipefail

usage() {
  cat <<EOF
Usage: sbatch $(basename "$0") [options]

Required:
  --ref <path>              Reference genome FASTA
  --output-dir <path>       Output directory
  --bin-size <int>          Bin size in bp

Optional:
  --lib-size <int>          Library size (default: 400)
  --lib-size-se <int>       Library size SE (default: 60)
  --read-length <int>       Read length (default: 150)
  --coverage <float>        Coverage (default: 1)
  -h, --help                Show this help
EOF
  exit 1
}

# --- arguments ---

# Defaults
REF_PATH=""
OUTPUT_DIR=""
BIN_SIZE=""
SIM_LIBRARY_SIZE=400
SIM_LIBRARY_SIZE_SE=60
SIM_READ_LENGTH=150
SIM_COVERAGE=1

# Parse
while [[ $# -gt 0 ]]; do
  case $1 in
    --ref)
      REF_PATH="$2"
      shift 2
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --bin-size)
      BIN_SIZE="$2"
      shift 2
      ;;
    --lib-size)
      SIM_LIBRARY_SIZE="$2"
      shift 2
      ;;
    --lib-size-se)
      SIM_LIBRARY_SIZE_SE="$2"
      shift 2
      ;;
    --read-length)
      SIM_READ_LENGTH="$2"
      shift 2
      ;;
    --coverage)
      SIM_COVERAGE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "ERROR: Unknown option '$1'" >&2
      echo "Use --help for usage information" >&2
      exit 1
      ;;
  esac
done

# Validate
if [[ -z "$REF_PATH" || -z "$OUTPUT_DIR" || -z "$BIN_SIZE" ]]; then
  echo "ERROR: --ref, --output-dir, and --bin-size are required" >&2
  usage
fi

echo ">>> Parameters:"
echo "    Reference:    $REF_PATH"
echo "    Output dir:   $OUTPUT_DIR"
echo "    Bin size:     $BIN_SIZE bp"
echo "    Lib size:     $SIM_LIBRARY_SIZE ± $SIM_LIBRARY_SIZE_SE bp"
echo "    Read length:  $SIM_READ_LENGTH bp"
echo "    Coverage:     ${SIM_COVERAGE}x"
echo

# --- get packages ---
#test if it is running as a slurm array
if [[ ${SLURM_JOB_ID+x} ]]; then
  is_slurm=1
  threads=${SLURM_CPUS_PER_TASK:-1}
  echo "Script is running as a slurm job"
    if [[ ${SLURM_ARRAY_TASK_ID}+x ]]; then
      is_slurm_array=1
      echo "Script is running as a part of a slurm array"
      else
        echo "Script is not running as a part of a slurm array"
        is_slurm_array=0
        SLURM_ARRAY_TASK_ID=1 #to make the rest of the code run
    fi
  else
  is_slurm=0
  SLURM_ARRAY_TASK_ID=1 #to make the rest of the code run
  threads=$(nproc)
  threads=$(( threads < 16 ? threads : 16 )) #cap to a maximum of 8 CPU cores
fi

echo This script will use $threads CPU cores to run.

#if [[ is_slurm == 1 ]]; then

  module --force purge
  module load calcua/2024a            # base toolchain (GCC 13.3.0 family)

  # Load matching builds (all with GCC-13.3.0)
  module load BWA/0.7.18-GCCcore-13.3.0
  module load BEDTools/2.31.1-GCC-13.3.0
  module load SAMtools/1.21-GCC-13.3.0

  #load the container containing art_ilumina
  unset PYTHONPATH
  export PATH="${VSC_SCRATCH}/containers/scDNAseq/bin:$PATH"
#fi

# --- actual run ---
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
bwa mem -t $threads "$REF_PATH" "$wd/simulated_reads1.fq" "$wd/simulated_reads2.fq" | \
samtools fixmate -u -m - - | \
samtools sort -u -@ $threads - | \
samtools markdup -@ $threads - "$wd/mapped_simulated_reads_sorted.bam"

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