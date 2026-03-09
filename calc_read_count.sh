#!/bin/bash

#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --job-name=count_reads
***REMOVED***
***REMOVED***
#SBATCH --mail-type=BEGIN,END,FAIL

#this is a simple bash script that will compute binned read counts from a given bam file.
# --- internal variables and bash settings --- 

set -euo pipefail
IFS=$'\n\t'

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

# --- parameters --- 
##set usage helper
usage() {
  cat <<EOF
Usage: sbatch calc_read_count.sh [options]

Options:
  -b <bam file>           the path to the bam file (required)
  -s <bin size>           the size of the genomic bins in basepairs
  -r <reference genome>   path to reference genome file in fasta format (required)
  -o <output_dir>         Output directory (required)
  -l <file list>          a tsv file containing a list of files, one per line. Meant to be used with slurm arrays.
  -h                      Show this help
EOF
  exit 1
}

# parse options (note new -l)
while getopts ":b:s:r:o:l:h" opt; do
  case "$opt" in
    b) BAM_FILE="$OPTARG" ;;
    s) BIN_SIZE="$OPTARG" ;;
    r) REF_GENOME="$OPTARG" ;;
    o) OUTPUTS_DIR="$OPTARG" ;;
    l) LISTFILE="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done


# --- tests ---
# If -l provided, ignore -b (explicit)
if [[ ${LISTFILE+x} ]]; then
    if [[ ! $is_slurm_array ]]; then
        echo "[ERROR] Option \`-l\` is meant to be used with slurm arrays only!" >&2
        usage
    else
        unset BAM_FILE >/dev/null 2>&1 || true
        BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LISTFILE") #tells wich line of the tsv file to read based on the SLURM_ARRAY_TASK_ID
    fi
fi

# Basic validation
if [[ ! ${BAM_FILE+x} || ! ${OUTPUTS_DIR+x} || ! ${REF_GENOME+x} || ! ${BIN_SIZE+x} ]]; then
  echo "ERROR: -b <bam file> -r <REF_GENOME> -s <bin size> and -o <output_dir> are required." >&2
  usage
fi

echo ">>> Parameters:"
echo "    bam file:  $BAM_FILE"
echo "    Rererence genome: $REF_GENOME"
echo "    Bin Size: ${BIN_SIZE} bp"
echo "    Output directory: $OUTPUTS_DIR"
echo
#make sure we have the full path of all needed files
REF_GENOME=$(realpath "$REF_GENOME")
BAM_FILE=$(realpath "$BAM_FILE") #get's the full path to the file
#get sample name (basically the name of the bam file without the file extension)
sample=$(basename "$BAM_FILE") #get's the file name
sample="${sample%%.*}" #remove any extension.

# --- modules ---
module --force purge
module load calcua/2024a calcua/all
module load BEDTools/2.31.1-GCC-13.3.0 SAMtools/1.21-GCC-13.3.0 parallel/20240722-GCCcore-13.3.0

# --- run ---
##store the time the run started
start_time=$(date +%s)

##define output dir (same as input)
mkdir -p "$OUTPUTS_DIR"

##test if the reference genome was indexed by samtools. If not, task 1 will index it while the other tasks wait
###to make sure indexing was finished, task 1 will create a `.index_samtools_done` dummy file to signal the conclusion of indexing. 
index_done="${REF_GENOME}.index_samtools_done"

if [[ ! -f "$index_done" ]]; then
  if (( SLURM_ARRAY_TASK_ID == 1 )); then
    echo "[INFO] Index not found. Task 1 (this task) will build it..."
    samtools faidx "$REF_GENOME"
    touch "$index_done"                # marker created only when bwa index finishes
    echo "[INFO] Index built."
  else
    echo "[INFO] Waiting for index to finish..."
    while [[ ! -f "$index_done" ]]; do
      sleep 20
    done
    echo "[INFO] Index detected."
  fi
fi

##create the bed file for the genomic bins
temp_dir="${OUTPUTS_DIR}/tmp_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$temp_dir"
echo "$temp_dir"
bedtools makewindows -g "${REF_GENOME}.fai" -w "$BIN_SIZE" > "${temp_dir}/genome_${BIN_SIZE}bp_bins.bed"
bins_bed="${temp_dir}/genome_${BIN_SIZE}bp_bins.bed"

##count reads
printf "Counting reads for sample %s with bin size %s bp based on file:\n\n%s\n" "$sample" "$BIN_SIZE" "$BAM_FILE"

    samtools sort -n -@ "${threads}" -O BAM -T "${temp_dir}.${sample}.nsort" "${BAM_FILE}" \
  | samtools fixmate -m -@ "${threads}" - - \
  | samtools sort    -@ "${threads}" -O BAM -T "${temp_dir}.${sample}.csort" - \
  | samtools markdup -r -@ "${threads}" - - \
  | samtools view -b -F 4 -q 30 -@ "${threads}" - \
  | bedtools coverage -a "${bins_bed}" -b stdin -counts -sorted -g "${REF_GENOME}.fai" \
    > "${OUTPUTS_DIR}/${sample}.counts.bed"

##remove temporary files
rm -rf ${temp_dir}/

#display how long the run took
end_time=$(date +%s)
echo "Job completed in $(( (end_time - start_time)/60 )) minutes."