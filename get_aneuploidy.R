#This is an R script that takes as input a bedgraph file containing binned counts and outputs aneuploidy estimations
#set libraries (this vector is then used by the `check_packages()` function)
packages <- c("tidyverse", "mgcv", "ggalign")

#set parameters
plot_scale <- 1.8

#custom functions####
##function to check and automatically install packages if needed####
# Function to install packages from CRAN, Bioconductor, or fail
check_packages <- function(package_list) {
  
  # Check if BiocManager is installed (needed for Bioconductor)
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager for Bioconductor package management...")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  
  # Track failed packages
  failed_packages <- character(0)
  
  # Iterate through each package
  for (pkg in package_list) {
    
    # Check if package is already installed
    if (requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("✓ Package '%s' is already installed", pkg))
      next
    }
    
    message(sprintf("Package '%s' not found. Attempting to install...", pkg))
    
    # Try CRAN first
    message(sprintf("  Trying CRAN for '%s'...", pkg))
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      
      # Verify installation
      if (requireNamespace(pkg, quietly = TRUE)) {
        message(sprintf("✓ Successfully installed '%s' from CRAN", pkg))
        next
      }
    }, error = function(e) {
      message(sprintf("  Failed to install '%s' from CRAN", pkg))
    })
    
    # Try Bioconductor if CRAN failed
    message(sprintf("  Trying Bioconductor for '%s'...", pkg))
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
      
      # Verify installation
      if (requireNamespace(pkg, quietly = TRUE)) {
        message(sprintf("✓ Successfully installed '%s' from Bioconductor", pkg))
        next
      }
    }, error = function(e) {
      message(sprintf("  Failed to install '%s' from Bioconductor", pkg))
    })
    
    # If we reach here, both methods failed
    failed_packages <- c(failed_packages, pkg)
  }
  
  # Stop with error if any packages failed
  if (length(failed_packages) > 0) {
    error_msg <- sprintf(
      "Failed to install the following package(s): %s\nPlease install them manually.",
      paste(failed_packages, collapse = ", ")
    )
    stop(error_msg)
  }
  
  message("\n✓ All packages are installed, loading them...")
  for(pkg in package_list){
    library(pkg, character.only = TRUE)
    message(sprintf("✓ Loaded '%s'", pkg))
  }
}
##function to print a usage helper####
helper <- function(){
  cat("[INFO] Ussage: Rscript get_aneuploidy.R [options]\n")
  cat("\n")
  cat("Options:\n")
  cat("   --inputs_dir      (required) path to the directory containing the bedgraph files.\n")
  cat("   --outputs_dir     (required) path to the directory where outputs should be stored.\n")
  cat("   --ploidy          (optional) The expected ploidy of the samples. Defaults to 2.\n")
  if(interactive()){
    stop()
  }else{
    quit(save = "no", status = 0)    
  }
}
##function to validate the CLI arguments####
validate_args <- function(x){
  expected_args <- c("inputs_dir", "outputs_dir", "ploidy")
  
  for(arg in args){
    parts <- strsplit(arg, "=")[[1]]
    key <- gsub("--", "", parts[1])
    if(!key %in% expected_args){
      cat("[ERROR]: \`", key, "\` is not a valid argument.\n", sep = "")
      helper()
    }else{
      value <- parts[2]
      cat(key, ": ", value, "\n", sep = "")
      assign(key, value, envir = .GlobalEnv) 
    }
  }
}

##function to validate bed files####
validate_file <- function(x){
  error_message <- function(){
    cat("[ERROR] The input bed graph must be a table containing 6 columns (no column names!): \n")
    cat("\n")
    cat("Column 1: contig name (character)\n")
    cat("Column 2: start (integer)\n")
    cat("Column 3: end (integer)\n")
    cat("Column 4: read count (integer)\n")
    cat("Column 5: GC content (numeric between 0 and 1)\n")
    cat("Column 6: mappability (numeric between 0 and 1)\n")
    if(interactive()){
      stop()
    }else{
      quit(save = "no", status = 0)    
    }
  }
  
  #test number of columns
  if(ncol(x) != 6){
    error_message()
  }
  
  if(any(sapply(x, class) != c("character", "numeric", "numeric", "numeric", "numeric", "numeric"))){
    error_message()
  }
  
  if(any(x[,2:4] != sapply(x[,2:4], as.integer))){
    error_message()
  }
  
  if(min(x[,5:6]) < 0 | max(x[,5:6]) > 1){
    error_message()
  }
  colnames(x) <- c("chromosome", "start", "end", 'read_count', "gc_content", "mappability")
  return(x)
}

##function to compensate gc_bias####
correct_gc_bias <- function(x){
  #determine baseline chromosomes
  x <- x %>%
    group_by(chromosome) %>%
    mutate(is_empty = read_count < 3) %>%
    mutate(is_outlier = read_count < quantile(read_count, 0.05) | read_count > quantile(read_count, 0.95)) %>%
    rowwise() %>%
    mutate(to_remove = any(is_outlier, is_empty, mappability < 1)) %>%
    ungroup()

  #to calculate gc bias we need to exclude aneuploid chromosomes. So we first calculate roughly which are the disomic chromosomes
  base_chromo <- x %>% 
    filter(!to_remove) %>%
    mutate(norm_read_count = read_count/median(read_count)) %>%
    group_by(chromosome) %>%
    summarise(mean_norm_read_count = median(norm_read_count)) %>%
    mutate(is_aneuploid = mean_norm_read_count < 0.9 | mean_norm_read_count > 1.1) %>%
    filter(!is_aneuploid) %>%
    pull(chromosome)
  
  #now we use gam modeling to determine the relationship between GC content and read count
  to_model <- x %>%
    filter(!to_remove) %>%
    filter(chromosome %in% base_chromo) %>%
    filter(gc_content < quantile(gc_content, 0.99) & gc_content > quantile(gc_content, 0.01)) %>%
    select(chromosome, read_count, gc_content)
  
  model <- gam(read_count ~ s(gc_content), data = to_model)
  
  #finally, we apply the model to correct the counts
  x <- x %>%
    mutate(gc_predicted = as.numeric(predict(model, newdata = data.frame(gc_content = gc_content)))) %>%
    mutate(gc_correction_factor = median(gc_predicted) / gc_predicted) %>%
    mutate(read_count_corrected = read_count * gc_correction_factor) %>%
    mutate(read_count_corrected = ifelse(read_count_corrected < 0, 0, read_count_corrected)) %>%
    mutate(norm_read_count = ifelse(to_remove, NA, read_count_corrected)) %>%
    mutate(norm_read_count = norm_read_count/median(norm_read_count, na.rm = TRUE))
  
  return(x)
}

#actual code####
##load libraries####
check_packages(packages)
##set defaults####
ploidy <- 2

##get arguments####
##args <- c("--inputs_dir=inputs/WGS/E4_G3_pTB007/count_files/", "--outputs_dir=outputs/WGS/E4_G3_pTB007") #meant for testing
args <- commandArgs(trailingOnly = TRUE)
validate_args(args)

##load input files#####
file_names <- list.files(inputs_dir, pattern = "*.bed")
files <- paste0(inputs_dir, "/", file_names)
cat("Loading a total of", length(files), "files.\n")

files <- files %>%
  lapply(read_delim, col_names = FALSE) %>%
  lapply(validate_file)

sample_names <- gsub(".counts|.bed", "", file_names)
names(files) <- sample_names

##now compensate gc_bias####
files <- files %>%
  lapply(correct_gc_bias)

##calculate somies####
combined_data <- files %>%
  bind_rows(.id = "sample") 

somies <- combined_data %>%
  filter(!to_remove) %>%
  group_by(sample, chromosome) %>%
  summarise(somy = median(norm_read_count) * ploidy)

#create output dir
dir.create(paste0(outputs_dir, "/", "individual_plots"), recursive = TRUE)

##save the outputs####
###save the list of bed files as an RDS file###
write_rds(files, paste0(outputs_dir, "/bedfile_list.rds"))

###save the somies as a matrix####
somies %>%
  pivot_wider(names_from = chromosome, values_from = somy) %>%
  column_to_rownames("sample") %>%
  write.table(paste0(outputs_dir, "/somy_matrix.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

#plot individual plots####
for(sample in names(files)){
  to_plot <- files[[sample]]
  plot_list <- list()
  
  ###plot gc bias####
  #plot bias before correction
  plot_list[[1]] <- to_plot %>%
    mutate(norm_read_count = read_count/median(read_count)) %>% #to include noisy bins
    ggplot(aes(x = gc_content, y = norm_read_count, color = mappability == 1))+
    labs(title = "before correction\n(red = low-mappability bin)", x = "GC %", y = "Normalized read count")
  
  #plot it after correction
  plot_list[[2]] <- to_plot %>%
    ggplot(aes(x = gc_content, y = norm_read_count))+
    labs(title = "after correction", x = "GC %", y = "Normalized read count (corrected)")
  
  #add shared settings to both plots
  plot_list  <- plot_list %>%
    lapply(function(x){
      x <- x+
        geom_point(size = 0.1)+
        geom_smooth(data = ~ filter(.x, mappability == 1), color = "blue")+
        scale_x_continuous(limits = c(0.4, 0.7))+
        scale_y_continuous(limits = c(0, 3))+
        scale_color_manual(values = c("darkred", "black"))+
        guides(color = "none")
      return(x)
    })
  
  #combine the plots
  gc_bias_plot <- ggalign::align_plots(!!!plot_list, ncol = 1)
  
  ###plot genome coverage####
  ####plot outlier bins
  to_plot <- to_plot %>%
    group_by(chromosome) %>%
    mutate(median_count = median(read_count_corrected, na.rm = TRUE))
  
  plot_list[[1]] <- to_plot %>%
    ggplot(aes(x = start, y = read_count_corrected, color = is_empty))+
    scale_color_manual(values = c("lightgrey", "darkred"))+
    ggtitle("\nEmpty bins = red")
  
  plot_list[[2]] <- to_plot %>%
    ggplot(aes(x = start, y = read_count_corrected, color = is_outlier))+
    scale_color_manual(values = c("lightgrey", "darkred"))+
    ggtitle("\nOutlier bins = red")
  
  ####plot low_mappability bins
  plot_list[[3]] <- to_plot %>%
    ggplot(aes(x = start, y = read_count_corrected, color = mappability < 1))+
    scale_color_manual(values = c("lightgrey", "darkred"))+
    ggtitle("\nLow mappability bins = red")
  
  ####plot the kept bins only
  plot_list[[4]] <- to_plot %>%
    filter(!to_remove) %>%
    ggplot(aes(x = start, y = read_count_corrected))+
    scale_color_manual(values = "black")+
    ggtitle("\nfiltered bins")
  
  #set y-axis scale limits
  y_limits <- 2 * ceiling(c(0, quantile(to_plot$read_count_corrected,  0.99)))
  
  ####add common settings
  plot_list <- plot_list %>%
    lapply(function(x){
      x <- x+
        geom_point(size = 0.1)+
        facet_grid(cols = vars(chromosome), scale = "free_x", switch = "x")+
        scale_x_continuous(breaks = NULL, labels = NULL, name = NULL)+
        scale_y_continuous(limits = y_limits, name = "corrected read count")+
        guides(col = "none")+
        theme_void()+
        theme(
          strip.clip = "off", 
          strip.text.x = element_text(angle = 45, hjust = 1),
          strip.text.y = element_text(angle = 270),
          axis.text.y = element_text(),
          axis.title.y = element_text(angle = 90),
        )
      return(x)
    })
  coverage_plot <- ggalign::align_plots(!!!plot_list, ncol = 1)
  
  ###plot somies####
  somy_plot <- to_plot %>%
    ungroup() %>%
    filter(!to_remove) %>%
    mutate(norm_read_count = norm_read_count/median(norm_read_count)) %>%
    ggplot(aes(x = chromosome, y = norm_read_count * ploidy))+
    geom_boxplot()+
    scale_y_continuous(name = "Somy")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  top <- ggalign::align_plots(gc_bias_plot, somy_plot, nrow = 1, widths = c(0.25, 0.75))
  final_plot <- ggalign::align_plots(top, coverage_plot, ncol = 1, heights = c(0.3, 0.7))
  ggsave(paste0(outputs_dir, "/individual_plots/", sample, "_coverage_plot.png"), final_plot, height = 11.6*plot_scale, width = 8.27*plot_scale, dpi = 600)
}

##plot grouped plots####
##plot somies
heat_col <- c("#001221", "#002342", "#002342", "#014175", "#035ba3", "#00c3ff", "#00ffee", "#33ff00", "#ccff00", "#fffa00","#ffa600", "#D73027", "#A50026", "#541b1b", "#4d0600")

plot_list <- list()
plot_list[["somies"]] <- somies %>%
  ggplot(aes(x = sample, y = fct_rev(chromosome), fill = somy, label = round(somy, 1)))+
  geom_tile()+
  geom_text(color = "black", size = 3.2)+
  geom_text(color = "white", size = 3)+
  scale_fill_gradientn(colors = heat_col, limits = c(0, 8), breaks = c(0:8))+
  labs(y = "Chromosome", title = "bulk somies")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#visualize coverage
plot_list[["coverage"]] <- combined_data %>%
  filter(!to_remove) %>%
  group_by(sample, chromosome) %>%
  mutate(norm_read_count = norm_read_count * ploidy) %>%
  mutate(median_count = median(norm_read_count)) %>%
  ggplot(aes(x = start, y = norm_read_count))+
  geom_point(size = 0.1)+
  geom_line(aes(y = median_count), color = "darkred")+
  facet_grid(rows = vars(sample), cols = vars(chromosome), scale = "free_x", switch = "x")+
  scale_x_continuous(breaks = NULL, labels = NULL, name = NULL)+
  scale_y_continuous(breaks = seq(0, 10, length = 6))+
  labs(y = "Normalized read count * ploidy", title = "coverage per 1kb bin")+
  theme_void()+
  theme(
    strip.clip = "off", 
    strip.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 270),
    axis.text.y = element_text(),
    axis.title.y = element_text(angle = 90),
  )

final_plot <- ggalign::align_plots(!!!plot_list, ncol = 1, heights = c(0.3, 0.7))
ggsave(paste0(outputs_dir, "/", "all_samples_coverage_and_somy.png"), 
       final_plot, 
       height = 11.6*plot_scale, 
       width = 8.27*plot_scale, 
       dpi = 600)

