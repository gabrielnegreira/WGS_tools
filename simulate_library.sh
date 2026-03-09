#!/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --job-name=gc_map
#SBATCH --mail-type=BEGIN,END,FAIL

#this script will take as input a mapped bam file and will simulate and artificial library based on the main metrics of that file. 
#for this to work we need art_illumina, bwa, samtools, picard, and bedtools. 
#art ilumina is in the scDNAseq container.

#=========== Prepare the Environment ====================
module --force purge
module load calcua/2024a

# Load matching builds (all with GCC-13.3.0)
module load BWA/0.7.18-GCCcore-13.3.0
module load BEDTools/2.31.1-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load R/4.4.2-gfbf-2024a #required for picard to run the estimated insert distribution
module load picard/3.3.0-Java-17

#load the container containing art_ilumina
unset PYTHONPATH
export PATH="/scratch/antwerpen/205/vsc20542/containers/scDNAseq/bin:$PATH"

#create a function to facilitate running picard
picard() {
  java -jar "$EBROOTPICARD/picard.jar" "$@"
}

set -e

# ======== User-defined Parameters ===========
REF_PATH="/scratch/antwerpen/205/vsc20542/projects/BPK081_subclones/inputs/reference_genome/LdCL_original.fasta" # reference genome file path
BAM_FILE="/scratch/antwerpen/205/vsc20542/projects/BPK081_subclones/inputs/WGS/subclones/bam_files/H3GWTDRXY-104487-001-049.bam" #path to bam file
OUTPUT_DIR="/scratch/antwerpen/205/vsc20542/projects/BPK081_subclones/inputs/WGS/subclones/bam_files/" #
NTHREADS=8

# =========== Get the metrics from the bam file ============
## get the file name
bam_name="$(basename "$BAM_FILE")"
bam_name="${bam_name%.*}"

## get the mean and standard deviation for insert sizes
printf "Using picard to determine mean and standard deviation of the insert sizes...\n"
picard CollectInsertSizeMetrics \
  -I "$BAM_FILE" \
  -O "bam_metrics.temp" \
  -H "histogram.temp"

picard CollectAlignmentSummaryMetrics \
  -I "$BAM_FILE" \
  -R "$REF_PATH" \
  -O "bam_alignment_metrics.temp"

### assign the mean and standard deviation insert sizes to variables
#### first we need to get the index of the desired columns
col_names=$(head bam_metrics.temp | awk '/^MEDIAN_INSERT_SIZE/ {print}') #get all column names in the summary table created by picard (MEDIAN_INSERT_SIZE is the first column)
read -a col_names <<< "$col_names" #convert to array
#### convert to key-index pair
declare -A col_index
for i in "${!col_names[@]}"; do
    col_index["${col_names[$i]}"]=$i
done
#### map the values to the variables
mean_col=$(( col_index["MEAN_INSERT_SIZE"] + 1 )) #bash index is 0-based, awk is 1-based
sd_col=$(( col_index["STANDARD_DEVIATION"] + 1 ))
mean_size=$(awk -v col="$mean_col" '/^MEDIAN_INSERT_SIZE/{getline; print $col}' bam_metrics.temp) #`getline` gets the value in the next line following the line that contains the `MEDIAN_INSERT_SIZE` header
sd_size=$(awk -v col="$sd_col" '/^MEDIAN_INSERT_SIZE/{getline; print $col}' bam_metrics.temp)

## assign the mean coverage to a variable
coverage_depth=$(samtools depth -aa -q 20 -Q 20 $BAM_FILE | awk '{sum+=$3; n++} END{print sum/n}')

## assign the mode read length
read_length=$(
  samtools view -F 0x904 "$BAM_FILE" \
  | awk '{print length($10)}' \
  | sort -n | uniq -c | sort -nr \
  | awk 'NR==1{print $2}'
)

echo "Mean insert size:" $mean_size
echo "Insert size standard deviation:" $sd_size
echo "Read length:" $read_length
echo "Mean coverage depth:" $coverage_depth

# =========== Create artificial Library ====================
#now we create an artificial fastq file which is used te calculate bin mappability
mkdir -p "$OUTPUT_DIR"

art_illumina \
  -i "$REF_PATH" \
  -p \
  -l $read_length \
  -f $coverage_depth \
  -m $mean_size \
  -s $sd_size \
  -o "$OUTPUT_DIR/${bam_name}_simulated_reads"

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
bwa mem \
  -t $NTHREADS \
  "$REF_PATH" \
  "$OUTPUT_DIR/${bam_name}_simulated_reads1.fq" \
  "$OUTPUT_DIR/${bam_name}_simulated_reads2.fq" | \
samtools fixmate -u -m - - | \
samtools sort -u -@ $NTHREADS - | \
samtools markdup -@ $NTHREADS - "$OUTPUT_DIR/${bam_name}_simulated_reads_sorted.bam"

#keep only reads with mapping quality higher than 30
samtools view -b -q 30 "$OUTPUT_DIR/${bam_name}_simulated_reads_sorted.bam" > "$OUTPUT_DIR/${bam_name}_simulated_reads_sorted_mapq30.bam"