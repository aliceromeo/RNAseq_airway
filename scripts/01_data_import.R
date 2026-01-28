# ==== Loading libraries ====
# Suppress logs
pkgs <- c("tximport", "readr", "txdbmaker", "AnnotationDbi")
invisible(lapply(pkgs, function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
}))

# ==== Load input files via Snakemake ====
# Define Salmon input and metadata files
input_files <- snakemake@input[["salmon_files"]]
samples_files <- snakemake@input[["samples_meta"]]
samples <- read.csv(samples_files, header = TRUE)
# Map filenames to SRA ID
names(input_files) <- samples$Run..SRA.ID

if (!all(file.exists(input_files))) {
  stop(paste("ERROR! Salmon quant.sf files missing!"))
}

# ==== Annotations and Mapping ====
gtf_path <- snakemake@input[["gtf"]]
# Create path for database SQLite in the same directory as GTF
sqlite_path <- gsub("\\.gtf\\.gz$|\\.gtf$", ".sqlite", gtf_path)
# If sqlite exists, it is loaded, otherwise it is created
if (!file.exists(sqlite_path)) {
  message("Generating SQLite database from GTF")
  txdb <- txdbmaker::makeTxDbFromGFF(gtf_path)
  saveDb(txdb, file = sqlite_path)
} else {
  message("Loading existing SQLite database...")
  txdb <- loadDb(sqlite_path)
}

# Mapping transcript -> gene
k <- keys(txdb, keytype = "TXNAME")
# Conversion table
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
# Clean IDs
tx2gene$TXNAME <- gsub("\\..*", "", tx2gene$TXNAME)
tx2gene$GENEID <- gsub("\\..*", "", tx2gene$GENEID)

# ==== Salmon data import ====
txi <- tximport(input_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# ==== Save output object ====
output_data <- list(
  txi = txi,
  metadata = samples
)

saveRDS(output_data, file = snakemake@output[["rds"]])
message("Success: Data imported and saved to ", snakemake@output[["rds"]])
