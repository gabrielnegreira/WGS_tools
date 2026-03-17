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
REF_GENOME=""
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
      REF_GENOME="$2"
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
if [[ -z "$REF_GENOME" || -z "$OUTPUT_DIR" || -z "$BIN_SIZE" ]]; then
  echo "ERROR: --ref, --output-dir, and --bin-size are required" >&2
  usage
fi

echo ">>> Parameters:"
echo "    Reference:    $REF_GENOME"
echo "    Output dir:   $OUTPUT_DIR"
echo "    Bin size:     $BIN_SIZE bp"
echo "    Lib size:     $SIM_LIBRARY_SIZE ± $SIM_LIBRARY_SIZE_SE bp"
echo "    Read length:  $SIM_READ_LENGTH bp"
echo "    Coverage:     ${SIM_COVERAGE}x"
echo

# --- get packages ---
#test if it is running as a slurm array
if [[ -z "SLURM_JOB_ID" ]]; then
  is_slurm=1
  threads=${SLURM_CPUS_PER_TASK:-1}
  echo "Script is running as a slurm job"
    if [[ -z "$SLURM_ARRAY_TASK_ID" ]]; then
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
mkdir -p "$OUTPUT_DIR"

art_illumina -ss HS25 -i "$REF_GENOME" -p -l $SIM_READ_LENGTH -f $SIM_COVERAGE -m $SIM_LIBRARY_SIZE -s $SIM_LIBRARY_SIZE_SE -o "$OUTPUT_DIR/simulated_reads"
rm "$OUTPUT_DIR"/*.aln


# ========== Map Artificial library to reference genome =======================
#now we map the artificial reads to the reference using bwa
# Check if BWA index files exist
if [[ ! -f "${REF_GENOME}.amb" || ! -f "${REF_GENOME}.ann" || ! -f "${REF_GENOME}.bwt" || \
      ! -f "${REF_GENOME}.pac" || ! -f "${REF_GENOME}.sa" ]]; then  
    echo "BWA index not found. Running bwa index..."
    bwa index "$REF_GENOME"
else
    echo "BWA index already exists. Skipping indexing."
fi

#map the artificial reads to the reference genome
bwa mem -t $threads "$REF_GENOME" "$OUTPUT_DIR/simulated_reads1.fq" "$OUTPUT_DIR/simulated_reads2.fq" | \
samtools fixmate -u -m - - | \
samtools sort -u -@ $threads - | \
samtools markdup -@ $threads - "$OUTPUT_DIR/mapped_simulated_reads_sorted.temp"

#keep only reads with mapping quality higher than 30
samtools view -b -q 30 "$OUTPUT_DIR/mapped_simulated_reads_sorted.temp" > "$OUTPUT_DIR/mapped_simulated_reads_sorted_filtered.temp"


# =============== Calculate mappability and GC content ========================

#to get coverage per bin we first need to create a genome size file
samtools faidx "$REF_GENOME" > "$REF_GENOME.fai"
cut -f1,2 "$REF_GENOME.fai" > "$OUTPUT_DIR/genome_sizes.temp"

#then we crete a bed file to store the information
bedtools makewindows -g "$OUTPUT_DIR/genome_sizes.temp" -w $BIN_SIZE > "$OUTPUT_DIR/genome_bins.temp"

#get total bin coverage
bedtools coverage -a "$OUTPUT_DIR/genome_bins.temp" -b "$OUTPUT_DIR/mapped_simulated_reads_sorted.temp" > "$OUTPUT_DIR/total_cov.temp"
#get bin coverage with unique maps only
bedtools coverage -a "$OUTPUT_DIR/genome_bins.temp" -b "$OUTPUT_DIR/mapped_simulated_reads_sorted_filtered.temp" > "$OUTPUT_DIR/unique_cov.temp"

#combine both results into a single file
paste <(cut -f1-4 "$OUTPUT_DIR/total_cov.temp") <(cut -f4-7 "$OUTPUT_DIR/unique_cov.temp") > "$OUTPUT_DIR/total_vs_unique.temp"

#now calculate the ratio between total reads mapping and uniquely mapped reads (the mappability). Note: This part was written by ChatGPT. 
awk 'BEGIN {FS="\t"; OFS="\t"}
{
  if ($4 > 0) {
    m = $5 / $4
  } else {
    m = 0
  }
  print $0, m
}' "$OUTPUT_DIR/total_vs_unique.temp" > "$OUTPUT_DIR/total_vs_unique_with_mappability.temp" 

echo -e "chrom\tstart\tend\tsim_total_read_count\tsim_unique_read_count\tsim_bases_covered\tBIN_SIZE\tsim_fraction_covered\tmappability" > "$OUTPUT_DIR/header.temp"
cat "$OUTPUT_DIR/header.temp" "$OUTPUT_DIR/total_vs_unique_with_mappability.temp" > "$OUTPUT_DIR/reads_per_bin.temp"
#get bin gc content
bedtools nuc -fi "$REF_GENOME" -bed "$OUTPUT_DIR/genome_bins.temp" > "$OUTPUT_DIR/gc_content_per_bin.temp" 
#extract only the column with the gc content
tail -n +1 "$OUTPUT_DIR/gc_content_per_bin.temp" | cut -f5 > "$OUTPUT_DIR/total_cov.temp" 

# =============== Save final results ===================
#append the gc content to the final file
ref_name=$( basename "$REF_GENOME" )
ref_name="${ref_name%.*}"
paste "$OUTPUT_DIR/reads_per_bin.temp" "$OUTPUT_DIR/total_cov.temp" > "$OUTPUT_DIR/${ref_name}_${BIN_SIZE}bp_mappability.tsv"

#remove unecessary files
rm -rf *.temp
rm -rf simulated_reads1.fq simulated_reads2.fq